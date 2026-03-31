 @testset "Chain solve" begin
    function build_chain_test_config()
        config = terra.nitrogen_10ev_config(; isothermal = false)
        config = terra.with_case_path(config, mktempdir())
        config = terra.with_time(config;
                                 dt = 5e-12,
                                 dt_output = 5e-7,
                                 duration = 5e-7,
                                 nstep = 200000,
                                 method = 2)
        config = terra.with_runtime(config;
                                    validate_species_against_terra = false,
                                    print_source_terms = false)
        return config
    end

    function profile_inlet(config;
                           source_compact_index::Int = 1,
                           mole_fractions = config.reactor.composition.mole_fractions,
                           total_number_density_m3 = nothing,
                           thermal = config.reactor.thermal)
        n_total_m3 = config.runtime.unit_system == :SI ?
                     config.reactor.composition.total_number_density :
                     terra.convert_number_density_cgs_to_si(config.reactor.composition.total_number_density)
        n_total_val = total_number_density_m3 === nothing ? n_total_m3 :
                      Float64(total_number_density_m3)
        composition = terra.ChainProfileInletComposition(;
                                                         species = config.reactor.composition.species,
                                                         mole_fractions = mole_fractions,
                                                         total_number_density_m3 = n_total_val,)
        return terra.ChainProfileInlet(;
                                       composition = composition,
                                       thermal = thermal,
                                       source_compact_index = source_compact_index,)
    end

    function cleanup_terra!()
        try
            terra.finalize_api_wrapper()
        catch
            # ignore cleanup errors in tests
        end
        try
            terra.close_terra_library()
        catch
            # ignore cleanup errors in tests
        end
    end

    function species_velocity_profile(;
                                      n_segments::Int,
                                      neutral_base::Float64,
                                      ion_base::Float64,)
        return Dict("N" => fill(neutral_base + 20.0, n_segments),
                    "N2" => fill(neutral_base, n_segments),
                    "N+" => collect(range(ion_base, step = 500.0, length = n_segments)),
                    "N2+" => collect(range(ion_base + 1000.0, step = 500.0,
                                           length = n_segments)))
    end

    function wall_profile(; n_segments::Int)
        return terra.ChainWallProfile(;
                                      a_wall_over_v_m_inv = fill(2.0 / 0.0155, n_segments),
                                      channel_gap_m = fill(0.0155, n_segments),)
    end

    @testset "Single segment run" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))
        marching = terra.AxialMarchingConfig()

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.failed_cell === nothing
        @test length(chain.cells) == 1
        @test chain.cells[1].reactor.success == true
        @test chain.cells[1].species_u_m_s["N+"] == 18000.0

        segment_result = chain.cells[1].reactor
        @test !isempty(segment_result.t)
        @test segment_result.frames[1].temperatures.tt ≈ config.reactor.thermal.Tt
        @test segment_result.frames[1].temperatures.tv ≈ config.reactor.thermal.Tv
        @test all(abs(frame.temperatures.te - profile.te_K[1]) <= 1e-6
                  for
                  frame in segment_result.frames)
        segment1_config = terra._chain_segment_config(config,
                                                      profile,
                                                      1,
                                                      terra._profile_inlet_reactor(profile,
                                                                                   config.runtime.unit_system),
                                                      marching)
        @test segment1_config.models.physics.is_isothermal_teex == true
        chain_log_path = joinpath(config.runtime.case_path, "output", "logs", "chain.log")
        @test isfile(chain_log_path)
        chain_log = read(chain_log_path, String)
        @test occursin("chain:", chain_log)
        @test occursin("segments: 1", chain_log)
        @test occursin("segment 1/1:", chain_log)
        @test occursin("species_u_m_s:", chain_log)
        @test occursin("number_density_m3:", chain_log)
        @test occursin("result:", chain_log)
        run_log_path = joinpath(config.runtime.case_path, "output", "logs", "run.log")
        @test isfile(run_log_path)
        run_log = read(run_log_path, String)
        @test occursin("========= TERRA 1D Chain Simulation =========", run_log)
        @test occursin("segments: 1", run_log)
        @test occursin("===== Segment 1/1 0D Simulation =====", run_log)
        @test occursin("chain success!", run_log)
        @test occursin("============================================", run_log)
        @test !occursin("Initializing TERRA", run_log)

        cleanup_terra!()
    end

     @testset "Chain segment runtime logging stays file-only under segment case path" begin
        config = terra.nitrogen_10ev_config()
        config = terra.with_case_path(config, mktempdir())
        config = terra.with_logging(config;
                                    native_stream_mode = :both,
                                    integration_detail_mode = :both,
                                    chain_detail_mode = :both,
                                    log_dir = "custom_logs")
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = Dict("N" => [150.0],
                                                               "N2" => [120.0],
                                                               "N+" => [18000.0],
                                                               "N2+" => [19000.0]),
                                          inlet = profile_inlet(config))
        marching = terra.AxialMarchingConfig()
        inlet_reactor = terra._profile_inlet_reactor(profile, config.runtime.unit_system)
        segment_config = terra._chain_segment_config(config, profile, 1, inlet_reactor,
                                                     marching)

        expected_segment_log_dir = normpath(joinpath(terra.log_dir(config.runtime),
                                                     "chain_segments",
                                                     "segment_0001"))
        @test segment_config.runtime.logging.native_stream_mode == :file
        @test segment_config.runtime.logging.integration_detail_mode == :file
        @test segment_config.runtime.logging.chain_detail_mode == :off
        @test segment_config.runtime.logging.log_dir == expected_segment_log_dir
        @test terra.log_dir(segment_config.runtime) == expected_segment_log_dir
    end

     @testset "Chain detail logging can be disabled" begin
        config = build_chain_test_config()
        config = terra.with_logging(config;
                                    console_mode = :quiet,
                                    chain_detail_mode = :off)
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))

        chain = terra.solve_terra_chain_steady(config, profile)
        @test chain.success == true
        @test !isfile(joinpath(config.runtime.case_path, "output", "logs", "chain.log"))
        run_log_path = joinpath(config.runtime.case_path, "output", "logs", "run.log")
        @test isfile(run_log_path)
        run_log = read(run_log_path, String)
        @test occursin("========= TERRA 1D Chain Simulation =========", run_log)
        @test occursin("===== Segment 1/1 0D Simulation =====", run_log)
        @test occursin("chain success!", run_log)

        cleanup_terra!()
    end

     @testset "Chain segment presentation uses reactor-owned 0D seam" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))
        marching = terra.AxialMarchingConfig()
        inlet_reactor = terra._profile_inlet_reactor(profile,
                                                     config.runtime.unit_system)
        segment_config = terra._chain_segment_config(config, profile, 1,
                                                     inlet_reactor, marching)

        result_ref = Ref{Any}(nothing)
        cache_ref = Ref{Any}(nothing)
        capture_pipe = Pipe()
        redirect_stdout(capture_pipe) do
            @test_nowarn terra.initialize_terra(segment_config,
                                                segment_config.runtime.case_path;
                                                lifecycle_console = :never,
                                                preserve_active_runtime = false)
            result_ref[], cache_ref[] = terra._solve_terra_0d_internal(segment_config;
                                                                       presentation = terra.CHAIN_SEGMENT_PRESENTATION)
            nothing
        end
        close(Base.pipe_writer(capture_pipe))

        console_text = read(capture_pipe, String)
        result = result_ref[]::terra.ReactorResult
        cache = cache_ref[]::terra.ReactorStateCache

        @test result.success == true
        @test !isempty(result.frames)
        @test cache.species == segment_config.reactor.composition.species
        @test occursin("starting ODE integration...", console_text)
        @test !occursin("TERRA 0D Simulation", console_text)

        cleanup_terra!()
    end

     @testset "Chain detail logging can mirror to console and file" begin
        config = build_chain_test_config()
        config = terra.with_logging(config;
                                    console_mode = :quiet,
                                    chain_detail_mode = :both)
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))

        chain_ref = Ref{Any}(nothing)
        capture_pipe = Pipe()
        redirect_stdout(capture_pipe) do
            chain_ref[] = terra.solve_terra_chain_steady(config, profile)
            nothing
        end
        close(Base.pipe_writer(capture_pipe))

        console_text = read(capture_pipe, String)
        chain = chain_ref[]::terra.ChainSimulationResult
        @test chain.success == true
        @test occursin("chain:", console_text)
        @test occursin("segment 1/1:", console_text)
        @test occursin("result:", console_text)
        @test isfile(joinpath(config.runtime.case_path, "output", "logs", "chain.log"))

        cleanup_terra!()
    end

     @testset "Chain solve does not emit standalone 0D banner" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))

        chain_ref = Ref{Any}(nothing)
        capture_pipe = Pipe()
        redirect_stdout(capture_pipe) do
            chain_ref[] = terra.solve_terra_chain_steady(config, profile)
            nothing
        end
        close(Base.pipe_writer(capture_pipe))

        console_text = read(capture_pipe, String)
        chain = chain_ref[]::terra.ChainSimulationResult
        @test chain.success == true
        @test occursin("========= TERRA 1D Chain Simulation =========", console_text)
        @test occursin("===== Segment 1/1 0D Simulation =====", console_text)
        @test occursin("starting ODE integration...", console_text)
        @test occursin("chain success!", console_text)
        @test occursin("============================================", console_text)
        @test !occursin("TERRA 0D Simulation", console_text)

        cleanup_terra!()
    end

     @testset "Chain solve restores the top-level active runtime" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))

        chain = terra.solve_terra_chain_steady(config, profile)

        @test chain.success == true
        active_runtime = terra._active_runtime_for_logging()
        @test active_runtime !== nothing
        @test active_runtime.case_path == config.runtime.case_path

        terra.finalize_terra()

        run_log_path = joinpath(config.runtime.case_path, "output", "logs", "run.log")
        run_log = read(run_log_path, String)
        @test occursin("Finalizing TERRA", run_log)
        @test occursin("TERRA finalized successfully", run_log)
    end

     @testset "Two-segment handoff and Te profile enforcement" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0, 0.01],
                                          dx_m = [0.01, 0.01],
                                          te_K = [115000.0, 90000.0],
                                          species_u_m_s = Dict("N" => [220.0, 225.0],
                                                               "N2" => [200.0, 205.0],
                                                               "N+" => [18000.0, 22000.0],
                                                               "N2+" => [19000.0, 23000.0]),
                                          inlet = profile_inlet(config))
        marching = terra.AxialMarchingConfig()

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test all(cell.success for cell in chain.cells)

        first_result = chain.cells[1].reactor
        second_result = chain.cells[2].reactor
        first_endpoint = chain.cells[1].endpoint_reactor
        @test first_result.success == true
        @test second_result.success == true
        @test all(abs(frame.temperatures.te - profile.te_K[1]) <= 1e-6
                  for
                  frame in first_result.frames)
        @test all(abs(frame.temperatures.te - profile.te_K[2]) <= 1e-6
                  for
                  frame in second_result.frames)

        @test first_endpoint !== nothing
        @test first_endpoint.thermal.Tee ≈ first_endpoint.thermal.Te
        @test second_result.frames[1].temperatures.tt ≈ first_endpoint.thermal.Tt
        @test second_result.frames[1].temperatures.tv ≈ first_endpoint.thermal.Tv
        @test second_result.frames[1].temperatures.te ≈ profile.te_K[2]
        segment2_config = terra._chain_segment_config(config, profile, 2,
                                                      first_endpoint, marching)
        @test segment2_config.models.physics.is_isothermal_teex == true
        @test segment2_config.reactor.thermal.Tee ≈ profile.te_K[2]
        @test chain.cells[2].species_u_m_s["N2+"] == 23000.0

        cleanup_terra!()
    end

     @testset "Two-segment full-state rho_ex handoff" begin
        config = build_chain_test_config()
        config = terra.with_time(config;
                                 dt = 5e-12,
                                 dt_output = 5e-8,
                                 duration = 5e-8,
                                 nstep = 20000,
                                 method = 2)
        profile = terra.AxialChainProfile(z_m = [0.0, 0.01],
                                          dx_m = [0.01, 0.01],
                                          te_K = [115000.0, 90000.0],
                                          species_u_m_s = Dict("N" => [220.0, 225.0],
                                                               "N2" => [200.0, 205.0],
                                                               "N+" => [18000.0, 22000.0],
                                                               "N2+" => [19000.0, 23000.0]),
                                          inlet = profile_inlet(config;
                                                                thermal = terra.ReactorThermalState(;
                                                                                                    Tt = 500.0,
                                                                                                    Tv = 650.0,
                                                                                                    Tee = 2000.0,
                                                                                                    Te = 115000.0)))
        marching = terra.AxialMarchingConfig(;
                                             handoff_policy = terra.FullStateHandoff(),
                                             override_tt_K = 800.0,
                                             override_tv_K = 900.0,)

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.metadata.diagnostics["segment_rho_ex_handoff_applied"] == [false, true]
        @test chain.metadata.diagnostics["full_state_rho_ex_handoff_supported"] == true
        @test !haskey(chain.metadata.diagnostics, "tee_policy")
        @test chain.cells[1].reactor.frames[1].temperatures.tt ≈ 800.0
        @test chain.cells[1].reactor.frames[1].temperatures.tv ≈ 900.0
        @test chain.cells[1].reactor.frames[1].temperatures.te ≈ profile.te_K[1]

        inlet_reactor = terra._profile_inlet_reactor(profile,
                                                     config.runtime.unit_system)
        segment1_config = terra._chain_segment_config(config, profile, 1,
                                                      inlet_reactor, marching)
        @test segment1_config.models.physics.is_isothermal_teex == true
        @test segment1_config.sources.residence_time !== nothing
        @test_nowarn terra.initialize_terra(segment1_config,
                                            segment1_config.runtime.case_path)
        segment1_result, segment1_cache = terra._solve_terra_0d_internal(segment1_config)

        segment1_endpoint = terra._segment_endpoint_reactor(config, segment1_result)
        segment2_config = terra._chain_segment_config(config, profile, 2,
                                                      segment1_endpoint, marching)
        state_from_cache = terra.config_to_initial_state(segment2_config;
                                                         state_cache = segment1_cache)
        state_from_boltz = terra.config_to_initial_state(segment2_config)
        boltz_rho_ex_active = state_from_boltz.rho_ex[:,
                                                      1:length(profile.inlet.composition.species)]
        @test any(abs.(state_from_cache.rho_ex .- boltz_rho_ex_active) .> 0.0)
        @test segment2_config.models.physics.is_isothermal_teex == true

        chain_segment2 = chain.cells[2].reactor
        @test chain_segment2.frames[1].temperatures.tt ≈ 800.0
        @test chain_segment2.frames[1].temperatures.tv ≈ 900.0
        @test chain_segment2.frames[1].temperatures.te ≈ profile.te_K[2]

        cleanup_terra!()
    end

     @testset "Chain marching can opt out of isothermal Teex enforcement" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0, 0.01],
                                          dx_m = [0.01, 0.01],
                                          te_K = [115000.0, 90000.0],
                                          species_u_m_s = Dict("N" => [220.0, 225.0],
                                                               "N2" => [200.0, 205.0],
                                                               "N+" => [18000.0, 22000.0],
                                                               "N2+" => [19000.0, 23000.0]),
                                          inlet = profile_inlet(config))
        marching = terra.AxialMarchingConfig(; is_isothermal_teex = false)

        inlet_reactor = terra._profile_inlet_reactor(profile,
                                                     config.runtime.unit_system)
        segment1_config = terra._chain_segment_config(config, profile, 1,
                                                      inlet_reactor, marching)
        segment2_config = terra._chain_segment_config(config, profile, 2,
                                                      inlet_reactor, marching)

        @test segment1_config.models.physics.is_isothermal_teex == false
        @test segment2_config.models.physics.is_isothermal_teex == false

        cleanup_terra!()
    end

     @testset "Metadata propagation and compact/source mapping" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0, 0.01],
                                          dx_m = [0.01, 0.01],
                                          te_K = [115000.0, 90000.0],
                                          species_u_m_s = Dict("N" => [220.0, 225.0],
                                                               "N2" => [200.0, 205.0],
                                                               "N+" => [18000.0, 22000.0],
                                                               "N2+" => [19000.0, 23000.0]),
                                          inlet = profile_inlet(config),
                                          generator = Dict("tool" => "unit-test",
                                                           "tool_version" => "1.0.0",
                                                           "created_utc" => "2026-03-05T00:00:00Z"),
                                          selection = Dict("average_start_time_s" => 5e-4,
                                                           "exported_species" => [
                                                               "N",
                                                               "N2",
                                                               "N+",
                                                               "N2+",
                                                           ],
                                                           "ion_velocity_policy" => "trim_to_first_positive",
                                                           "u_ion_floor" => 0.0,
                                                           "min_consecutive_positive" => 3,
                                                           "trim_start_index" => 13,
                                                           "trim_start_z_m" => 0.0025,
                                                           "trimmed_point_count" => 12,
                                                           "original_point_count" => 102),
                                          schema_version = "terra_chain_profile_v4",
                                          source_snapshot = Dict("enabled" => true,
                                                                 "source_type" => "unit-test"))

        chain = terra.solve_terra_chain_steady(config, profile)

        @test chain.success == true
        @test chain.metadata.schema_version == "terra_chain_profile_v4"
        @test chain.metadata.generator["tool"] == "unit-test"
        @test chain.metadata.selection["trim_start_index"] == 13
        @test chain.metadata.source_snapshot !== nothing
        @test chain.metadata.source_snapshot["source_type"] == "unit-test"
        @test chain.metadata.compact_to_source_index == [13, 14]
        @test chain.metadata.original_point_count == 102
        @test chain.metadata.retained_point_count == 2
        @test chain.cells[1].compact_cell_index == 1
        @test chain.cells[2].compact_cell_index == 2
        @test chain.cells[1].source_cell_index == 13
        @test chain.cells[2].source_cell_index == 14
        @test !isempty(chain.cells[1].reactor.t)
        @test !isempty(chain.cells[2].reactor.t)

        cleanup_terra!()
    end

     @testset "Wall losses change the chain solution and expose diagnostics" begin
        base_config = build_chain_test_config()
        wall_cfg = terra.WallLossConfig(;
                                        species_models = Dict("N+" => terra.IonNeutralizationWallModel(;
                                                                                                        bohm_scale = 1.25,
                                                                                                        products = Dict("N" => 1.0)),
                                                              "N2+" => terra.IonNeutralizationWallModel(;
                                                                                                         products = Dict("N2" => 1.0)),
                                                              "N" => terra.ConstantNeutralRecombinationWallModel(;
                                                                                                                 k_wall_1_s = 4.0e5,
                                                                                                                 products = Dict("N2" => 0.5))),)
        wall_config = terra.Config(;
                                   reactor = base_config.reactor,
                                   models = base_config.models,
                                   sources = terra.SourceTermsConfig(;
                                                                     wall_losses = wall_cfg),
                                   numerics = base_config.numerics,
                                   runtime = base_config.runtime,)

        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [base_config.reactor.thermal.Te],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          wall_profile = wall_profile(; n_segments = 1),
                                          inlet = profile_inlet(base_config))

        chain_base = terra.solve_terra_chain_steady(base_config, profile)
        chain_wall = terra.solve_terra_chain_steady(wall_config, profile)

        @test chain_wall.success == true
        @test chain_wall.metadata.diagnostics["wall_profile"]["channel_gap_m"] == [0.0155]
        @test chain_wall.metadata.diagnostics["wall_profile"]["a_wall_over_v_m_inv"] ≈
              [2.0 / 0.0155]

        wall_segment_inputs = terra._segment_wall_inputs(profile, 1, wall_cfg)
        @test wall_segment_inputs !== nothing
        @test wall_segment_inputs.channel_gap_m ≈ 0.0155

        base_frame = chain_base.cells[1].reactor.frames[end]
        wall_frame = chain_wall.cells[1].reactor.frames[end]
        @test any(abs.(wall_frame.species_densities .- base_frame.species_densities) .> 0.0)
        @test wall_frame.temperatures.te≈base_frame.temperatures.te atol=0.0 rtol=1e-12
        @test wall_frame.source_terms !== nothing
        @test haskey(wall_frame.source_terms.wall_losses.species_mass_density_rates, "N2+")
        @test wall_frame.source_terms.wall_losses.species_mass_density_rates["N2+"] < 0.0
        @test wall_frame.source_terms.wall_losses.species_mass_density_rates["N2"] > 0.0
        @test haskey(chain_wall.cells[1].reactor.metadata, "wall_losses")
        species_model = chain_wall.cells[1].reactor.metadata["wall_losses"]["species_models"]["N+"]
        @test species_model["model_type"] == "IonNeutralizationWallModel"
        @test sort!(collect(keys(chain_wall.cells[1].reactor.metadata["wall_losses"]))) ==
              ["segment_inputs", "species_models", "species_order"]
        @test sort!(collect(keys(species_model))) ==
              ["charge_state", "model_type", "parameters", "product_indices", "products",
               "reactant_ground_index", "reactant_indices"]

        missing_wall_profile = terra.AxialChainProfile(z_m = [0.0],
                                                       dx_m = [0.01],
                                                       te_K = [base_config.reactor.thermal.Te],
                                                       species_u_m_s = species_velocity_profile(;
                                                                                                n_segments = 1,
                                                                                                neutral_base = 180.0,
                                                                                                ion_base = 18000.0),
                                                       inlet = profile_inlet(base_config))
        @test_throws ArgumentError terra._segment_wall_inputs(missing_wall_profile, 1,
                                                              wall_cfg)

        cleanup_terra!()
    end

     @testset "Failed segment returns partial chain result" begin
        base_config = build_chain_test_config()
        wall_cfg = terra.WallLossConfig(;
                                        species_models = Dict("N+" => terra.IonNeutralizationWallModel(;
                                                                                                        products = Dict("N" => 1.0))),)
        wall_config = terra.Config(;
                                   reactor = base_config.reactor,
                                   models = base_config.models,
                                   sources = terra.SourceTermsConfig(; wall_losses = wall_cfg),
                                   numerics = base_config.numerics,
                                   runtime = base_config.runtime,)
        profile = terra.AxialChainProfile(z_m = [0.0, 0.01],
                                          dx_m = [0.01, 0.01],
                                          te_K = [115000.0, 90000.0],
                                          species_u_m_s = Dict("N" => [220.0, 225.0],
                                                               "N2" => [200.0, 205.0],
                                                               "N+" => [18000.0, 22000.0],
                                                               "N2+" => [19000.0, 23000.0]),
                                          inlet = profile_inlet(base_config))

        chain = terra.solve_terra_chain_steady(wall_config, profile)

        @test chain.success == false
        @test chain.failed_cell == 1
        @test chain.message == "Chain failed at segment 1 during setup/integration."
        @test occursin("wall_profile is missing", chain.cells[1].message)
        @test chain.cells[1].reactor.success == false
        @test chain.cells[2].message == "Not executed."
        @test chain.metadata.compact_to_source_index == [1, 2]

        cleanup_terra!()
    end

     @testset "Rejects Missing or Extra Profile Species" begin
        config = build_chain_test_config()

        missing_profile = terra.AxialChainProfile(z_m = [0.0],
                                                  dx_m = [0.01],
                                                  te_K = [115000.0],
                                                  species_u_m_s = Dict("N" => [220.0],
                                                                       "N2" => [200.0],
                                                                       "N+" => [18000.0]),
                                                  inlet = profile_inlet(config))
        @test_throws ArgumentError terra.solve_terra_chain_steady(config, missing_profile)

        extra_profile = terra.AxialChainProfile(z_m = [0.0],
                                                dx_m = [0.01],
                                                te_K = [115000.0],
                                                species_u_m_s = Dict("N" => [220.0],
                                                                     "N2" => [200.0],
                                                                     "N+" => [18000.0],
                                                                     "N2+" => [19000.0],
                                                                     "Ar" => [1000.0]),
                                                inlet = profile_inlet(config))
        @test_throws ArgumentError terra.solve_terra_chain_steady(config, extra_profile)

        cleanup_terra!()
    end

     @testset "Unsupported termination policy" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(z_m = [0.0],
                                          dx_m = [0.01],
                                          te_K = [115000.0],
                                          species_u_m_s = species_velocity_profile(;
                                                                                   n_segments = 1,
                                                                                   neutral_base = 180.0,
                                                                                   ion_base = 18000.0),
                                          inlet = profile_inlet(config))

        @test_throws ArgumentError terra.solve_terra_chain_steady(config,
                                                                  profile;
                                                                  marching = terra.AxialMarchingConfig(;
                                                                                                       termination_policy = terra.SteadyStateTermination()),)

        cleanup_terra!()
    end
end
