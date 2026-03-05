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

    @testset "Single segment run" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0],
            dx_m = [0.01],
            te_K = [config.reactor.thermal.Te],
            u_neutral_m_s = [200.0],
            u_ion_m_s = [18000.0],
        )
        marching = terra.AxialMarchingConfig(; store_segment_histories = true)

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.failed_segment === nothing
        @test length(chain.segment_end_reactors) == 1
        @test chain.segment_success == [true]
        @test chain.segment_results !== nothing
        @test chain.segment_results[1].success == true

        segment_result = chain.segment_results[1]
        @test !isempty(segment_result.time)
        @test maximum(abs.(segment_result.temperatures.te .- profile.te_K[1])) <= 1e-6

        cleanup_terra!()
    end

    @testset "Two-segment handoff and Te profile enforcement" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0, 0.01],
            dx_m = [0.01, 0.01],
            te_K = [115000.0, 90000.0],
            u_neutral_m_s = [200.0, 200.0],
            u_ion_m_s = [18000.0, 22000.0],
        )
        marching = terra.AxialMarchingConfig(; store_segment_histories = true)

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.segment_success == [true, true]
        @test chain.segment_results !== nothing

        first_result = chain.segment_results[1]
        second_result = chain.segment_results[2]
        @test first_result.success == true
        @test second_result.success == true
        @test maximum(abs.(first_result.temperatures.te .- profile.te_K[1])) <= 1e-6
        @test maximum(abs.(second_result.temperatures.te .- profile.te_K[2])) <= 1e-6

        # In reinitialize mode without thermal overrides, segment 2 starts from
        # segment 1 endpoint Tt/Tv while Te is overwritten by the profile.
        @test second_result.temperatures.tt[1] ≈ chain.segment_end_reactors[1].thermal.Tt
        @test second_result.temperatures.tv[1] ≈ chain.segment_end_reactors[1].thermal.Tv
        @test second_result.temperatures.te[1] ≈ profile.te_K[2]

        cleanup_terra!()
    end

    @testset "Unsupported handoff/termination modes" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0],
            dx_m = [0.01],
            te_K = [115000.0],
            u_neutral_m_s = [200.0],
            u_ion_m_s = [18000.0],
        )

        @test_throws ArgumentError terra.solve_terra_chain_steady(
            config,
            profile;
            marching = terra.AxialMarchingConfig(; handoff_mode = :full_state),
        )

        @test_throws ArgumentError terra.solve_terra_chain_steady(
            config,
            profile;
            marching = terra.AxialMarchingConfig(; termination_mode = :steady_state),
        )

        cleanup_terra!()
    end
end
