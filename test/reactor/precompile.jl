@testset "Precompile helpers" begin
    @testset "_native_precompile_ready" begin
        withenv(terra.TERRA_ENV_VAR_NAME => nothing) do
            @test terra._native_precompile_ready() == false
        end
    end

    @testset "_precompile_0d_config" begin
        for isothermal in (false, true)
            config = terra._precompile_0d_config(isothermal)
            try
                @test config.models.physics.is_isothermal_teex == isothermal
                @test config.numerics.time.dt == 5e-12
                @test config.numerics.time.dt_output == 1e-8
                @test config.numerics.time.duration == 1e-8
                @test config.numerics.time.nstep == 1000
                @test config.numerics.solver.saveat_count == 8
                @test config.runtime.validate_species_against_terra == false
                @test config.runtime.print_source_terms == false
                @test config.runtime.write_native_state_files == false
                @test config.runtime.logging.console_mode == :quiet
                @test config.runtime.logging.progress_mode == :off
                @test config.runtime.logging.native_stream_mode == :off
                @test config.runtime.logging.integration_detail_mode == :off
                @test config.runtime.logging.chain_detail_mode == :off
                @test isdir(config.runtime.case_path)
            finally
                rm(config.runtime.case_path; recursive = true, force = true)
            end
        end
    end

    @testset "_precompile_chain_profile" begin
        config = terra._precompile_0d_config(false)
        try
            profile = terra._precompile_chain_profile(config)
            marching = terra.AxialMarchingConfig()

            @test length(profile.z_m) == 1
            @test profile.z_m == [0.0]
            @test profile.dx_m == [0.01]
            @test profile.te_K == [config.reactor.thermal.Te]
            @test profile.species_u_m_s == Dict("N" => [200.0],
                                                "N2" => [180.0],
                                                "N+" => [18000.0],
                                                "N2+" => [19000.0])
            @test profile.inlet.source_compact_index == 1
            @test profile.inlet.composition.species == config.reactor.composition.species
            @test profile.inlet.composition.total_number_density_m3 ≈
                  terra.convert_number_density_cgs_to_si(config.reactor.composition.total_number_density)
            @test profile.inlet.thermal.Tt == 500.0
            @test profile.inlet.thermal.Tv == 500.0
            @test profile.inlet.thermal.Tee == config.reactor.thermal.Te
            @test profile.inlet.thermal.Te == config.reactor.thermal.Te
            @test_nowarn terra._validate_chain_inputs(config, profile, marching)
        finally
            rm(config.runtime.case_path; recursive = true, force = true)
        end
    end
end
