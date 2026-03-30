 @progress_testset "Active Runtime Context" begin
    runtime_a = terra.RuntimeConfig(; case_path = mktempdir())
    runtime_b = terra.RuntimeConfig(; case_path = mktempdir())

    try
        @test terra._active_runtime_for_logging() === nothing

        terra.with_active_runtime_for_logging(runtime_a; restore = nothing) do
            active = terra._active_runtime_for_logging()
            @test active !== nothing
            @test active.case_path == runtime_a.case_path

            terra.with_active_runtime_for_logging(runtime_b; restore = runtime_a) do
                inner = terra._active_runtime_for_logging()
                @test inner !== nothing
                @test inner.case_path == runtime_b.case_path
            end

            restored = terra._active_runtime_for_logging()
            @test restored !== nothing
            @test restored.case_path == runtime_a.case_path
        end

        @test terra._active_runtime_for_logging() === nothing
    finally
        rm(runtime_a.case_path; recursive = true, force = true)
        rm(runtime_b.case_path; recursive = true, force = true)
    end
end

 @progress_testset "Native Logging Preparation" begin
    temp_dir = mktempdir()
    custom_log_dir = mktempdir()
    try
        default_runtime = terra.RuntimeConfig(; case_path = temp_dir)
        default_logging = terra.prepare!(terra.NATIVE_LOG, default_runtime)
        @test default_logging.console_level == terra.API_NATIVE_LOG_OFF
        @test default_logging.file_level == terra.API_NATIVE_LOG_VERBOSE
        @test default_logging.log_path == joinpath(temp_dir, "output", "logs", "native.log")
        @test isfile(default_logging.log_path)

        custom_runtime = terra.RuntimeConfig(;
                                             case_path = temp_dir,
                                             logging = terra.LoggingConfig(;
                                                                           native_stream_mode = :off,
                                                                           log_dir = custom_log_dir))
        custom_logging = terra.prepare!(terra.NATIVE_LOG, custom_runtime)
        @test custom_logging.file_level == terra.API_NATIVE_LOG_OFF
        @test custom_logging.log_path === nothing
        @test isdir(custom_log_dir)
    finally
        rm(temp_dir; recursive = true, force = true)
        rm(custom_log_dir; recursive = true, force = true)
    end
end
