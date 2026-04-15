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
                @test isdir(config.runtime.case_path)
            finally
                rm(config.runtime.case_path; recursive = true, force = true)
            end
        end
    end
end
