 @progress_testset "CaseLayout" begin
    temp_dir = mktempdir()
    try
        runtime = terra.RuntimeConfig(; case_path = temp_dir)
        layout = terra.case_layout(runtime)

        @test layout.case_path == normpath(temp_dir)
        @test layout.input_dir == joinpath(temp_dir, "input")
        @test layout.output_dir == joinpath(temp_dir, "output")
        @test layout.sources_dir == joinpath(temp_dir, "output", "sources")
        @test layout.states_dir == joinpath(temp_dir, "output", "states")
        @test layout.log_dir == joinpath(temp_dir, "output", "logs")
        @test terra.log_dir(runtime) == layout.log_dir
        @test layout.run_log_path == joinpath(temp_dir, "output", "logs", "run.log")
        @test layout.native_log_path == joinpath(temp_dir, "output", "logs", "native.log")
        @test layout.chain_log_path == joinpath(temp_dir, "output", "logs", "chain.log")

        terra.ensure_case_layout!(runtime)
        @test isdir(layout.input_dir)
        @test isdir(layout.output_dir)
        @test isdir(layout.sources_dir)
        @test isdir(layout.states_dir)
        @test isdir(layout.log_dir)
    finally
        rm(temp_dir; recursive = true, force = true)
    end
end

 @progress_testset "CaseLayout Custom Log Dir" begin
    temp_dir = mktempdir()
    custom_log_dir = mktempdir()
    try
        runtime_rel = terra.RuntimeConfig(;
                                          case_path = temp_dir,
                                          logging = terra.LoggingConfig(; log_dir = "custom_logs"))
        relative_layout = terra.case_layout(runtime_rel)
        @test relative_layout.log_dir == joinpath(temp_dir, "custom_logs")
        @test terra.log_dir(runtime_rel) == relative_layout.log_dir

        runtime_abs = terra.RuntimeConfig(;
                                          case_path = temp_dir,
                                          logging = terra.LoggingConfig(; log_dir = custom_log_dir))
        absolute_layout = terra.case_layout(runtime_abs)
        @test absolute_layout.log_dir == normpath(custom_log_dir)

        terra.ensure_case_layout!(runtime_abs)
        @test isdir(custom_log_dir)
    finally
        rm(temp_dir; recursive = true, force = true)
        rm(custom_log_dir; recursive = true, force = true)
    end
end
