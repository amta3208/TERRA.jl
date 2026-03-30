 @progress_testset "Config (Nested)" begin
    reactor = terra.ReactorConfig(;
                                  composition = terra.ReactorComposition(;
                                                                         species = ["N",
                                                                             "N2",
                                                                             "E-"],
                                                                         mole_fractions = [0.1,
                                                                             0.8,
                                                                             0.1],
                                                                         total_number_density = 1e13),
                                  thermal = terra.ReactorThermalState(; Tt = 300.0,
                                                                      Tv = 350.0,
                                                                      Tee = 360.0,
                                                                      Te = 10000.0))
    models = terra.ModelConfig()
    numerics = terra.NumericsConfig(;
                                    time = terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4,
                                                            duration = 1e-3))
    sources = terra.SourceTermsConfig()
    runtime = terra.RuntimeConfig(; case_path = pwd(), unit_system = :CGS)

    config = terra.Config(;
                          reactor = reactor,
                          models = models,
                          sources = sources,
                          numerics = numerics,
                          runtime = runtime)

    @test config.reactor == reactor
    @test config.models == models
    @test config.sources.residence_time === sources.residence_time
    @test config.sources.wall_losses === nothing
    @test config.numerics == numerics
    @test config.runtime == runtime
end

 @progress_testset "Config Modifiers" begin
    wall_cfg = terra.WallLossConfig(;
                                    species_models = Dict("N+" => terra.IonNeutralizationWallModel(;
                                                                                                    products = Dict("N" => 1.0))),)
    config = terra.Config(;
                          reactor = terra.ReactorConfig(;
                                                        composition = terra.ReactorComposition(;
                                                                                               species = ["N",
                                                                                                   "N2",
                                                                                                   "E-"],
                                                                                               mole_fractions = [0.1,
                                                                                                   0.8,
                                                                                                   0.1],
                                                                                               total_number_density = 1e13),
                                                        thermal = terra.ReactorThermalState(;
                                                                                            Tt = 300.0,
                                                                                            Tv = 350.0,
                                                                                            Tee = 360.0,
                                                                                            Te = 10000.0)),
                          numerics = terra.NumericsConfig(;
                                                          time = terra.TimeConfig(;
                                                                                  dt = 1e-6,
                                                                                  dt_output = 1e-4,
                                                                                  duration = 1e-3,
                                                                                  nstep = 1234,
                                                                                  method = 2)),
                          sources = terra.SourceTermsConfig(; wall_losses = wall_cfg),
                          runtime = terra.RuntimeConfig(;
                                                        database_path = ".",
                                                        case_path = pwd(),
                                                        unit_system = :CGS,
                                                        validate_species_against_terra = false,
                                                        print_source_terms = false,
                                                        write_native_state_files = false,
                                                        logging = terra.LoggingConfig(;
                                                                                      integration_detail_mode = :off)))

    temp_case = mktempdir()
    config_case = terra.with_case_path(config, temp_case)
    @test config_case.runtime.case_path == temp_case
    @test config.runtime.case_path == pwd()

    config_time = terra.with_time(config; dt = 2e-6, duration = 2e-3, method = 1)
    @test config_time.numerics.time.dt == 2e-6
    @test config_time.numerics.time.duration == 2e-3
    @test config_time.numerics.time.dt_output == config.numerics.time.dt_output
    @test config_time.numerics.time.method == 1
    @test config.numerics.time.dt == 1e-6
    @test config_time.sources.residence_time === config.sources.residence_time
    @test config_time.sources.wall_losses === config.sources.wall_losses

    config_runtime = terra.with_runtime(config;
                                        unit_system = :SI,
                                        print_source_terms = true,
                                        write_native_state_files = true)
    @test config_runtime.runtime.unit_system == :SI
    @test config_runtime.runtime.print_source_terms == true
    @test config_runtime.runtime.write_native_state_files == true
    @test config_runtime.sources.residence_time === config.sources.residence_time
    @test config_runtime.sources.wall_losses === config.sources.wall_losses

    config_logging = terra.with_logging(config;
                                        console_mode = :verbose,
                                        progress_mode = :summary,
                                        integration_detail_mode = :both,
                                        chain_detail_mode = :both,
                                        log_dir = "logs/custom")
    @test config_logging.runtime.logging.console_mode == :verbose
    @test config_logging.runtime.logging.progress_mode == :summary
    @test config_logging.runtime.logging.integration_detail_mode == :both
    @test config_logging.runtime.logging.chain_detail_mode == :both
    @test config_logging.runtime.logging.log_dir == normpath("logs/custom")
    @test config.runtime.logging.console_mode == :minimal

    @test_throws ArgumentError terra.with_case_path(config, joinpath(temp_case, "missing"))
    @test_throws ArgumentError terra.with_time(config; method = 9)
    @test_throws MethodError terra.NumericsConfig(;
                                                  time = terra.TimeConfig(; dt = 1e-6,
                                                                          dt_output = 1e-4,
                                                                          duration = 1e-3),
                                                  residence_time = nothing)
end
