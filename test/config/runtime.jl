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
