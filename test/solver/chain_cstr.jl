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
        marching = terra.AxialMarchingConfig()

        chain = terra.solve_terra_chain_steady(config, profile; marching = marching)

        @test chain.success == true
        @test chain.failed_cell === nothing
        @test length(chain.cells) == 1
        @test chain.cells[1].reactor.success == true

        segment_result = chain.cells[1].reactor
        @test !isempty(segment_result.t)
        @test all(abs(frame.temperatures.te - profile.te_K[1]) <= 1e-6 for frame in segment_result.frames)

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

        # In reinitialize mode without thermal overrides, segment 2 starts from
        # segment 1 endpoint Tt/Tv while Te is overwritten by the profile.
        @test first_endpoint !== nothing
        @test second_result.frames[1].temperatures.tt ≈ first_endpoint.thermal.Tt
        @test second_result.frames[1].temperatures.tv ≈ first_endpoint.thermal.Tv
        @test second_result.frames[1].temperatures.te ≈ profile.te_K[2]

        cleanup_terra!()
    end

    @testset "Metadata propagation and compact/source mapping" begin
        config = build_chain_test_config()
        profile = terra.AxialChainProfile(
            z_m = [0.0, 0.01],
            dx_m = [0.01, 0.01],
            te_K = [115000.0, 90000.0],
            u_neutral_m_s = [200.0, 200.0],
            u_ion_m_s = [18000.0, 22000.0],
            generator = Dict(
                "tool" => "unit-test",
                "tool_version" => "1.0.0",
                "created_utc" => "2026-03-05T00:00:00Z",
            ),
            selection = Dict(
                "neutral_species" => "N2",
                "ion_species" => "N2",
                "ion_charge_state" => 1,
                "average_start_time_s" => 5e-4,
                "trim_start_index" => 13,
                "trimmed_point_count" => 12,
                "original_point_count" => 102,
            ),
            schema_version = "terra_chain_profile_v1",
            source_snapshot = Dict(
                "enabled" => true,
                "source_type" => "unit-test",
            ),
        )

        chain = terra.solve_terra_chain_steady(config, profile)

        @test chain.success == true
        @test chain.metadata.schema_version == "terra_chain_profile_v1"
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
