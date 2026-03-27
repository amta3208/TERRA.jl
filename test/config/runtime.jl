@testset "LoggingConfig" begin
    logging = terra.LoggingConfig()
    @test logging.console_mode == :minimal
    @test logging.progress_mode == :auto
    @test logging.native_stream_mode == :file
    @test logging.integration_detail_mode == :file
    @test logging.chain_detail_mode == :file
    @test logging.log_dir === nothing

    @test_throws ArgumentError terra.LoggingConfig(; console_mode = :loud)
    @test_throws ArgumentError terra.LoggingConfig(; native_stream_mode = :mirror)
    @test_throws ArgumentError terra.LoggingConfig(; log_dir = "   ")
end

@testset "RuntimeConfig rejects legacy aliases" begin
    @test_throws MethodError terra.RuntimeConfig(; case_path = pwd(),
                                                 write_native_outputs = true)
    @test_throws MethodError terra.RuntimeConfig(; case_path = pwd(),
                                                 print_integration_output = true)

    runtime = terra.RuntimeConfig(; case_path = pwd(),
                                  write_native_state_files = true,
                                  logging = terra.LoggingConfig(;
                                                                integration_detail_mode = :console))
    @test runtime.write_native_state_files == true
    @test runtime.logging.integration_detail_mode == :console
    @test_throws FieldError runtime.write_native_outputs
    @test_throws FieldError runtime.print_integration_output
    @test_throws MethodError terra.with_runtime(runtime; write_native_outputs = true)
    @test_throws MethodError terra.with_runtime(runtime; print_integration_output = true)
end

@testset "Chain Segment Runtime Logging" begin
    config = terra.nitrogen_10ev_config()
    config = terra.with_case_path(config, mktempdir())
    config = terra.with_logging(config;
                                native_stream_mode = :both,
                                integration_detail_mode = :both,
                                chain_detail_mode = :both,
                                log_dir = "custom_logs")

    n_total_m3 = config.runtime.unit_system == :SI ?
                 config.reactor.composition.total_number_density :
                 terra.convert_number_density_cgs_to_si(config.reactor.composition.total_number_density)
    inlet = terra.ChainProfileInlet(;
                                    composition = terra.ChainProfileInletComposition(;
                                                                                     species = config.reactor.composition.species,
                                                                                     mole_fractions = config.reactor.composition.mole_fractions,
                                                                                     total_number_density_m3 = n_total_m3),
                                    thermal = config.reactor.thermal,
                                    source_compact_index = 1)
    profile = terra.AxialChainProfile(z_m = [0.0],
                                      dx_m = [0.01],
                                      te_K = [config.reactor.thermal.Te],
                                      species_u_m_s = Dict("N" => [150.0],
                                                           "N2" => [120.0],
                                                           "N+" => [18000.0],
                                                           "N2+" => [19000.0]),
                                      inlet = inlet)
    marching = terra.AxialMarchingConfig()
    inlet_reactor = terra._build_profile_inlet_reactor(profile, config.runtime.unit_system)
    segment_config = terra._build_chain_segment_config(config, profile, 1, inlet_reactor,
                                                       marching)

    expected_segment_log_dir = normpath(joinpath(terra.log_dir(config.runtime),
                                                 "chain_segments", "segment_0001"))
    @test segment_config.runtime.logging.native_stream_mode == :file
    @test segment_config.runtime.logging.integration_detail_mode == :file
    @test segment_config.runtime.logging.chain_detail_mode == :off
    @test segment_config.runtime.logging.log_dir == expected_segment_log_dir
    @test terra.log_dir(segment_config.runtime) == expected_segment_log_dir
end
