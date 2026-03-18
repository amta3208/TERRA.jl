@testset "Chain CSTR solver" begin
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
            print_source_terms = false,
            write_native_outputs = false,
            print_integration_output = false)
        return config
    end

    function profile_inlet(config;
            source_compact_index::Int = 1,
            mole_fractions = config.reactor.composition.mole_fractions,
            total_number_density_m3 = nothing,
            thermal = config.reactor.thermal)
        n_total_m3 = config.runtime.unit_system == :SI ?
                     config.reactor.composition.total_number_density :
                     terra.convert_number_density_cgs_to_si(
            config.reactor.composition.total_number_density)
        n_total_val = total_number_density_m3 === nothing ? n_total_m3 : Float64(total_number_density_m3)
        composition = terra.ChainProfileInletComposition(;
            species = config.reactor.composition.species,
            mole_fractions = mole_fractions,
            total_number_density_m3 = n_total_val,
        )
        return terra.ChainProfileInlet(;
            composition = composition,
            thermal = thermal,
            source_compact_index = source_compact_index,
        )
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

    function species_velocity_profile(; n_segments::Int, neutral_base::Float64, ion_base::Float64)
        return Dict(
            "N" => fill(neutral_base + 20.0, n_segments),
            "N2" => fill(neutral_base, n_segments),
            "N+" => collect(range(ion_base, step = 500.0, length = n_segments)),
            "N2+" => collect(range(ion_base + 1000.0, step = 500.0, length = n_segments)),
        )
    end

    @testset "Single segment run" begin
        config = build_chain_test_config()
        inlet_thermal = terra.ReactorThermalState(; Tt = 500.0, Tv = 500.0,
            Tee = config.reactor.thermal.Te, Te = config.reactor.thermal.Te)
        profile = terra.AxialChainProfile(
            z_m = [0.0],
            dx_m = [0.01],
            te_K = [config.reactor.thermal.Te],
            species_u_m_s = species_velocity_profile(; n_segments = 1, neutral_base = 180.0,
                ion_base = 18000.0),
            inlet = profile_inlet(config;
                mole_fractions = [0.2, 0.65, 0.05, 0.05, 0.05],
                thermal = inlet_thermal),
        )
        marching = terra.AxialMarchingConfig()

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.failed_cell === nothing
        @test length(chain.cells) == 1
        @test chain.cells[1].reactor.success == true
        @test chain.cells[1].species_u_m_s["N+"] == 18000.0

        segment_result = chain.cells[1].reactor
        @test !isempty(segment_result.t)
        @test segment_result.frames[1].temperatures.tt ≈ inlet_thermal.Tt
        @test segment_result.frames[1].temperatures.tv ≈ inlet_thermal.Tv
        @test all(abs(frame.temperatures.te - profile.te_K[1]) <= 1e-6 for frame in segment_result.frames)
        segment1_config = terra._build_chain_segment_config(
            config,
            profile,
            1,
            terra._build_profile_inlet_reactor(profile, config.runtime.unit_system),
            marching,
        )
        @test segment1_config.models.physics.is_isothermal_teex == true

        cleanup_terra!()
    end

    @testset "Two-segment handoff and Te profile enforcement" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0, 0.01],
            dx_m = [0.01, 0.01],
            te_K = [115000.0, 90000.0],
            species_u_m_s = Dict(
                "N" => [220.0, 225.0],
                "N2" => [200.0, 205.0],
                "N+" => [18000.0, 22000.0],
                "N2+" => [19000.0, 23000.0],
            ),
            inlet = profile_inlet(config),
        )
        marching = terra.AxialMarchingConfig()

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test all(cell.success for cell in chain.cells)

        first_result = chain.cells[1].reactor
        second_result = chain.cells[2].reactor
        first_endpoint = chain.cells[1].endpoint_reactor
        @test first_result.success == true
        @test second_result.success == true
        @test all(abs(frame.temperatures.te - profile.te_K[1]) <= 1e-6 for frame in first_result.frames)
        @test all(abs(frame.temperatures.te - profile.te_K[2]) <= 1e-6 for frame in second_result.frames)

        @test first_endpoint !== nothing
        @test first_endpoint.thermal.Tee ≈ first_endpoint.thermal.Te
        @test second_result.frames[1].temperatures.tt ≈ first_endpoint.thermal.Tt
        @test second_result.frames[1].temperatures.tv ≈ first_endpoint.thermal.Tv
        @test second_result.frames[1].temperatures.te ≈ profile.te_K[2]
        segment2_config = terra._build_chain_segment_config(
            config, profile, 2, first_endpoint, marching)
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
        profile = terra.AxialChainProfile(
            z_m = [0.0, 0.01],
            dx_m = [0.01, 0.01],
            te_K = [115000.0, 90000.0],
            species_u_m_s = Dict(
                "N" => [220.0, 225.0],
                "N2" => [200.0, 205.0],
                "N+" => [18000.0, 22000.0],
                "N2+" => [19000.0, 23000.0],
            ),
            inlet = profile_inlet(config;
                thermal = terra.ReactorThermalState(;
                    Tt = 500.0, Tv = 650.0, Tee = 2000.0, Te = 115000.0)),
        )
        marching = terra.AxialMarchingConfig(;
            handoff_mode = :full_state,
            override_tt_K = 800.0,
            override_tv_K = 900.0,
        )

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.metadata.diagnostics["segment_rho_ex_handoff_applied"] == [false, true]
        @test chain.metadata.diagnostics["full_state_rho_ex_handoff_supported"] == true
        @test !haskey(chain.metadata.diagnostics, "tee_policy")
        @test chain.cells[1].reactor.frames[1].temperatures.tt ≈ 800.0
        @test chain.cells[1].reactor.frames[1].temperatures.tv ≈ 900.0
        @test chain.cells[1].reactor.frames[1].temperatures.te ≈ profile.te_K[1]

        inlet_reactor = terra._build_profile_inlet_reactor(profile, config.runtime.unit_system)
        segment1_config = terra._build_chain_segment_config(
            config, profile, 1, inlet_reactor, marching)
        @test segment1_config.models.physics.is_isothermal_teex == true
        @test segment1_config.sources.residence_time !== nothing
        @test_nowarn terra.initialize_terra(segment1_config, segment1_config.runtime.case_path)
        segment1_result, segment1_cache = terra._solve_terra_0d_internal(segment1_config)

        segment1_endpoint = terra._extract_segment_endpoint_reactor(config, segment1_result)
        segment2_config = terra._build_chain_segment_config(
            config, profile, 2, segment1_endpoint, marching)
        state_from_cache = terra.config_to_initial_state(segment2_config;
            state_cache = segment1_cache)
        state_from_boltz = terra.config_to_initial_state(segment2_config)
        boltz_rho_ex_active = state_from_boltz.rho_ex[:, 1:length(profile.inlet.composition.species)]
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
        profile = terra.AxialChainProfile(
            z_m = [0.0, 0.01],
            dx_m = [0.01, 0.01],
            te_K = [115000.0, 90000.0],
            species_u_m_s = Dict(
                "N" => [220.0, 225.0],
                "N2" => [200.0, 205.0],
                "N+" => [18000.0, 22000.0],
                "N2+" => [19000.0, 23000.0],
            ),
            inlet = profile_inlet(config),
        )
        marching = terra.AxialMarchingConfig(; is_isothermal_teex = false)

        inlet_reactor = terra._build_profile_inlet_reactor(profile, config.runtime.unit_system)
        segment1_config = terra._build_chain_segment_config(
            config, profile, 1, inlet_reactor, marching)
        segment2_config = terra._build_chain_segment_config(
            config, profile, 2, inlet_reactor, marching)

        @test segment1_config.models.physics.is_isothermal_teex == false
        @test segment2_config.models.physics.is_isothermal_teex == false

        cleanup_terra!()
    end

    @testset "Metadata propagation and compact/source mapping" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0, 0.01],
            dx_m = [0.01, 0.01],
            te_K = [115000.0, 90000.0],
            species_u_m_s = Dict(
                "N" => [220.0, 225.0],
                "N2" => [200.0, 205.0],
                "N+" => [18000.0, 22000.0],
                "N2+" => [19000.0, 23000.0],
            ),
            inlet = profile_inlet(config),
            generator = Dict(
                "tool" => "unit-test",
                "tool_version" => "1.0.0",
                "created_utc" => "2026-03-05T00:00:00Z",
            ),
            selection = Dict(
                "average_start_time_s" => 5e-4,
                "exported_species" => ["N", "N2", "N+", "N2+"],
                "ion_velocity_policy" => "trim_to_first_positive",
                "u_ion_floor" => 0.0,
                "min_consecutive_positive" => 3,
                "trim_start_index" => 13,
                "trim_start_z_m" => 0.0025,
                "trimmed_point_count" => 12,
                "original_point_count" => 102,
            ),
            schema_version = "terra_chain_profile_v3",
            source_snapshot = Dict(
                "enabled" => true,
                "source_type" => "unit-test",
            ),
        )

        chain = terra.solve_terra_chain_steady(config, profile)

        @test chain.success == true
        @test chain.metadata.schema_version == "terra_chain_profile_v3"
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

    @testset "Rejects Missing or Extra Profile Species" begin
        config = build_chain_test_config()

        missing_profile = terra.AxialChainProfile(
            z_m = [0.0],
            dx_m = [0.01],
            te_K = [115000.0],
            species_u_m_s = Dict(
                "N" => [220.0],
                "N2" => [200.0],
                "N+" => [18000.0],
            ),
            inlet = profile_inlet(config),
        )
        @test_throws ArgumentError terra.solve_terra_chain_steady(config, missing_profile)

        extra_profile = terra.AxialChainProfile(
            z_m = [0.0],
            dx_m = [0.01],
            te_K = [115000.0],
            species_u_m_s = Dict(
                "N" => [220.0],
                "N2" => [200.0],
                "N+" => [18000.0],
                "N2+" => [19000.0],
                "Ar" => [1000.0],
            ),
            inlet = profile_inlet(config),
        )
        @test_throws ArgumentError terra.solve_terra_chain_steady(config, extra_profile)

        cleanup_terra!()
    end

    @testset "Unsupported termination mode" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0],
            dx_m = [0.01],
            te_K = [115000.0],
            species_u_m_s = species_velocity_profile(; n_segments = 1, neutral_base = 180.0,
                ion_base = 18000.0),
            inlet = profile_inlet(config),
        )

        @test_throws ArgumentError terra.solve_terra_chain_steady(
            config,
            profile;
            marching = terra.AxialMarchingConfig(; termination_mode = :steady_state),
        )

        cleanup_terra!()
    end
end
