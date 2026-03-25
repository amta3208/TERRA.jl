@testset "nitrogen_10ev_config" begin
    @testset "Default Configuration" begin
        config = terra.nitrogen_10ev_config()

        @test config.reactor.composition.species == ["N", "N2", "N+", "N2+", "E-"]
        @test config.reactor.composition.mole_fractions ==
              [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
        @test config.reactor.composition.total_number_density == 1.0e13
        @test config.reactor.thermal.Tt == 750.0
        @test config.reactor.thermal.Tv == 750.0
        @test config.reactor.thermal.Te == 115000.0
        @test config.reactor.thermal.Tee == 750.0
        # Stored in seconds (prob_setup writes microseconds)
        @test config.numerics.time.dt ≈ 5e-12
        @test config.numerics.time.dt_output ≈ 5e-6
        @test config.numerics.time.duration ≈ 1e-3
        @test config.numerics.time.nstep == 500000
        @test config.numerics.time.method == 2
        @test config.sources.residence_time === nothing
        @test config.sources.wall_losses === nothing
    end
end

@testset "Nested Config solver controls" begin
    base = terra.nitrogen_10ev_config(; isothermal = false)
    case_path = mktempdir()

    config = terra.Config(;
                          reactor = base.reactor,
                          models = base.models,
                          numerics = terra.NumericsConfig(;
                                                          time = terra.TimeConfig(;
                                                                                  dt = 5e-12,
                                                                                  dt_output = 1e-6,
                                                                                  duration = 5e-7,
                                                                                  nstep = 1000,
                                                                                  method = 2),
                                                          solver = terra.ODESolverConfig(;
                                                                                         saveat_count = 7,
                                                                                         reltol = 1e-8,
                                                                                         abstol_density = 1e-10,
                                                                                         ramp_understep_ratio = inv(64),
                                                                                         ramp_history_steps = 4),
                                                          space = base.numerics.space),
                          runtime = terra.RuntimeConfig(;
                                                        database_path = base.runtime.database_path,
                                                        case_path = case_path,
                                                        unit_system = base.runtime.unit_system,
                                                        validate_species_against_terra = false,
                                                        print_source_terms = false,
                                                        write_native_outputs = false,
                                                        print_integration_output = false))

    @test_nowarn reset_and_init!(case_path; config = config)
    initial_state = terra.config_to_initial_state(config)
    results = terra.integrate_0d_system(config, initial_state)

    @test results.success == true
    @test length(results.t) == config.numerics.solver.saveat_count

    run_log_path = joinpath(case_path, "output", "logs", "run.log")
    @test isfile(run_log_path)
    run_log = read(run_log_path, String)
end

@testset "Direct 0D solve rejects wall-loss configs without profile inputs" begin
    base = terra.nitrogen_10ev_config(; isothermal = false)
    wall_cfg = terra.WallLossConfig(;
                                    species_models = Dict("N+" => terra.SpeciesWallModel(;
                                                                                         class = :ion_neutralization,
                                                                                         rate_model = :bohm_gap,
                                                                                         products = Dict("N" => 1.0),)),)
    config = terra.Config(;
                          reactor = base.reactor,
                          models = base.models,
                          sources = terra.SourceTermsConfig(; wall_losses = wall_cfg),
                          numerics = base.numerics,
                          runtime = base.runtime,)

    @test_throws ArgumentError terra.solve_terra_0d(config; sources = config.sources)
end

@testset "Native Output Generation" begin
    base_config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir(cleanup = false)
    config = terra.with_case_path(base_config, temp_case_path)
    config = terra.with_time(config;
                             dt = 5e-12, dt_output = 1e-6, duration = 5e-7, nstep = 1000,
                             method = 2)
    config = terra.with_runtime(config;
                                validate_species_against_terra = false,
                                print_source_terms = false,
                                write_native_outputs = true,
                                print_integration_output = false)

    # Fresh initialization that preserves the generated case directory
    try
        terra.finalize_api_wrapper()
    catch
        # ignore
    end
    try
        terra.close_terra_library()
    catch
        # ignore
    end

    terra.load_terra_library!()
    try
        terra.set_api_finalize_mpi_wrapper(false)
    catch
        # ignore if unavailable
    end

    terra.generate_input_files(config, temp_case_path)
    terra.initialize_api_wrapper(case_path = temp_case_path)

    case_path_used = terra.TERRA_CASE_PATH[]
    @test case_path_used == temp_case_path
    @test isdir(case_path_used)

    initial_state = terra.config_to_initial_state(config)
    results = terra.integrate_0d_system(config, initial_state)
    @test results.t[end] >= results.t[1]

    output_dir = joinpath(case_path_used, "output")
    @test isdir(output_dir)

    for fname in ("result-time.dat", "result-flow.dat", "result-temp.dat")
        fpath = joinpath(output_dir, fname)
        @test isfile(fpath)
        open(fpath, "r") do io
            header = readline(io)
            @test startswith(header, "TITLE=")
        end
    end

    states_dir = joinpath(output_dir, "states")
    @test isdir(states_dir)
    enex_files = filter(f -> occursin("enex", f), readdir(states_dir))
    @test !isempty(enex_files)

    terra.finalize_terra()
end

@testset "Run Log Detail Output" begin
    base_config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    config = terra.with_case_path(base_config, temp_case_path)
    config = terra.with_time(config;
                             dt = 5e-12, dt_output = 1e-6, duration = 5e-7, nstep = 1000,
                             method = 2)
    config = terra.with_runtime(config;
                                validate_species_against_terra = false,
                                print_source_terms = false,
                                write_native_outputs = false)
    config = terra.with_logging(config;
                                console_mode = :quiet,
                                progress_mode = :summary,
                                integration_detail_mode = :file)

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    quiet_results = Ref{Any}(nothing)
    quiet_pipe = Pipe()
    redirect_stdout(quiet_pipe) do
        quiet_results[] = terra.integrate_0d_system(config, initial_state)
        nothing
    end
    close(Base.pipe_writer(quiet_pipe))
    quiet_console_text = read(quiet_pipe, String)
    results = quiet_results[]::terra.ReactorResult
    @test results.success == true
    @test !occursin("Ytot,err", quiet_console_text)

    run_log_path = joinpath(temp_case_path, "output", "logs", "run.log")
    @test isfile(run_log_path)
    run_log = read(run_log_path, String)
    @test occursin("=========== TERRA 0D Simulation ===========", run_log)
    @test occursin("success!", run_log)
    @test occursin("===========================================", run_log)
end

@testset "0D Console Progress Routing" begin
    base_config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    config = terra.with_case_path(base_config, temp_case_path)
    config = terra.with_time(config;
                             dt = 5e-12, dt_output = 1e-6, duration = 5e-7, nstep = 1000,
                             method = 2)
    config = terra.with_runtime(config;
                                validate_species_against_terra = false,
                                print_source_terms = false,
                                write_native_outputs = false)
    config = terra.with_logging(config;
                                console_mode = :minimal,
                                progress_mode = :summary,
                                integration_detail_mode = :file)

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    routed_results = Ref{Any}(nothing)
    routed_pipe = Pipe()
    redirect_stdout(routed_pipe) do
        routed_results[] = terra.integrate_0d_system(config, initial_state)
        nothing
    end
    close(Base.pipe_writer(routed_pipe))
    console_text = read(routed_pipe, String)
    results = routed_results[]::terra.ReactorResult
    @test results.success == true
    @test occursin("=========== TERRA 0D Simulation ===========", console_text)
    @test occursin("starting ODE integration...", console_text)
    @test occursin("success!", console_text)
    @test occursin("===========================================", console_text)
end
