function _logging_test_config(case_path::AbstractString;
        console_mode::Symbol = :minimal,
        native_stream_mode::Symbol = :file,
        integration_detail_mode::Symbol = :file,
        chain_detail_mode::Symbol = :file)
    config = terra.nitrogen_10ev_config()
    config = terra.with_case_path(config, case_path)
    return terra.with_logging(config;
                              console_mode = console_mode,
                              native_stream_mode = native_stream_mode,
                              integration_detail_mode = integration_detail_mode,
                              chain_detail_mode = chain_detail_mode)
end

function _logging_test_profile(config::terra.Config)
    total_number_density_m3 = config.runtime.unit_system == :SI ?
                              config.reactor.composition.total_number_density :
                              terra.convert_number_density_cgs_to_si(config.reactor.composition.total_number_density)
    inlet = terra.ChainProfileInlet(;
                                    composition = terra.ChainProfileInletComposition(;
                                                                                     species = config.reactor.composition.species,
                                                                                     mole_fractions = config.reactor.composition.mole_fractions,
                                                                                     total_number_density_m3 = total_number_density_m3),
                                    thermal = config.reactor.thermal,
                                    source_compact_index = 1)
    return terra.AxialChainProfile(; z_m = [0.0],
                                   dx_m = [0.01],
                                   te_K = [config.reactor.thermal.Te],
                                   species_u_m_s = Dict("N" => [150.0],
                                                        "N2" => [120.0],
                                                        "N+" => [18000.0],
                                                        "N2+" => [19000.0]),
                                   inlet = inlet)
end

function _capture_stdout(f::Function)
    capture_pipe = Pipe()
    redirect_stdout(capture_pipe) do
        f()
        nothing
    end
    close(Base.pipe_writer(capture_pipe))
    return read(capture_pipe, String)
end

 @testset "Logging" begin
     @testset "Run log routes to file and console" begin
        temp_dir = mktempdir()
        try
            runtime = terra.RuntimeConfig(; case_path = temp_dir,
                                          logging = terra.LoggingConfig(; console_mode = :minimal))
            entry = terra.EventEntry(:info, "run event"; console = :minimal, :step => 1)
            console_text = _capture_stdout() do
                terra.emit!(terra.RUN_LOG, runtime, entry)
            end

            @test occursin("run event", console_text)
            @test isfile(terra.run_log_path(runtime))

            run_log = read(terra.run_log_path(runtime), String)
            @test occursin("INFO run event", run_log)
            @test occursin("step=1", run_log)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

     @testset "Exception entries keep exception detail in file output" begin
        temp_dir = mktempdir()
        try
            runtime = terra.RuntimeConfig(; case_path = temp_dir,
                                          logging = terra.LoggingConfig(; console_mode = :quiet))
            terra.emit!(terra.RUN_LOG, runtime,
                        terra.ExceptionEntry(:error, "boom", ArgumentError("bad input")))

            run_log = read(terra.run_log_path(runtime), String)
            @test occursin("ERROR boom", run_log)
            @test occursin("bad input", run_log)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

     @testset "Chain detail routing honors file, console, and off modes" begin
        temp_file = mktempdir()
        temp_console = mktempdir()
        temp_off = mktempdir()
        try
            config_file = _logging_test_config(temp_file; console_mode = :quiet,
                                               chain_detail_mode = :file)
            profile_file = _logging_test_profile(config_file)
            header_file = terra.ChainHeaderEntry(config_file, profile_file,
                                                 terra.AxialMarchingConfig(), [1])

            terra.prepare!(terra.CHAIN_LOG, config_file.runtime)
            terra.emit!(terra.CHAIN_LOG, config_file.runtime, header_file)
            @test isfile(terra.chain_log_path(config_file.runtime))
            @test occursin("chain:", read(terra.chain_log_path(config_file.runtime), String))

            config_console = _logging_test_config(temp_console; console_mode = :quiet,
                                                  chain_detail_mode = :console)
            profile_console = _logging_test_profile(config_console)
            header_console = terra.ChainHeaderEntry(config_console, profile_console,
                                                    terra.AxialMarchingConfig(), [1])

            console_text = _capture_stdout() do
                terra.emit!(terra.CHAIN_LOG, config_console.runtime, header_console)
            end
            @test occursin("chain:", console_text)
            @test !isfile(terra.chain_log_path(config_console.runtime))

            config_off = _logging_test_config(temp_off; console_mode = :quiet,
                                              chain_detail_mode = :off)
            profile_off = _logging_test_profile(config_off)
            header_off = terra.ChainHeaderEntry(config_off, profile_off,
                                                terra.AxialMarchingConfig(), [1])

            console_off = _capture_stdout() do
                terra.emit!(terra.CHAIN_LOG, config_off.runtime, header_off)
            end
            @test isempty(console_off)
            @test !isfile(terra.chain_log_path(config_off.runtime))
        finally
            rm(temp_file; recursive = true, force = true)
            rm(temp_console; recursive = true, force = true)
            rm(temp_off; recursive = true, force = true)
        end
    end

     @testset "Progress gating follows logging policy" begin
        summary_logging = terra.LoggingConfig(; progress_mode = :summary)
        quiet_auto_logging = terra.LoggingConfig(; console_mode = :quiet, progress_mode = :auto)
        off_logging = terra.LoggingConfig(; progress_mode = :off)

        @test terra._progress_mode_enabled(summary_logging) == true
        @test terra._progress_mode_enabled(quiet_auto_logging) == false
        @test terra._progress_mode_enabled(off_logging) == false

        summary_runtime = terra.RuntimeConfig(; case_path = mktempdir(),
                                              logging = summary_logging)
        off_runtime = terra.RuntimeConfig(; case_path = mktempdir(),
                                          logging = off_logging)
        try
            @test terra._progress_reporter(summary_runtime, 1.0) isa terra.FractionProgressReporter
            @test terra._progress_reporter(off_runtime, 1.0) === nothing
        finally
            rm(summary_runtime.case_path; recursive = true, force = true)
            rm(off_runtime.case_path; recursive = true, force = true)
        end
    end

     @testset "Chain scalar rendering covers common value families" begin
        @test terra._chain_scalar(nothing) == "null"
        @test terra._chain_scalar(true) == "true"
        @test terra._chain_scalar(:full_state) == "full_state"
        @test terra._chain_scalar(7) == "7"
        @test terra._chain_scalar(0.0) == "0.0"
        @test terra._chain_scalar(1.5e5) == "1.500000e+05"
        @test terra._chain_scalar("chain") == "chain"
        @test terra._chain_scalar([1, 2]) == repr([1, 2])
        @test terra._chain_scalar(Dict("N" => 1.0)) == repr(Dict("N" => 1.0))
        @test terra._chain_scalar((handoff = "full_state",)) ==
              repr((handoff = "full_state",))
    end
end
