 @testset "Integrate 0D (adiabatic)" begin
    # Initialize using the config-driven input to ensure the selected
    # database and options are honored (rather than a stale case file).

    config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    # Shorten integration to keep tests fast
    config = terra.with_case_path(config, temp_case_path)
    config = terra.with_time(config;
        dt = 5e-12, dt_output = 1e-6, duration = 1e-6, nstep = 500000, method = 2)
    config = terra.with_runtime(config;
        validate_species_against_terra = false,
        print_source_terms = false)

    # Initialize the Fortran API using a temporary case generated from this config
    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    results = @time terra.integrate_0d_system(config, initial_state)
    @test results isa terra.ReactorResult
    densities = terra.species_density_matrix(results)
    temperatures = terra.temperature_history(results)
    @test results.t[end] > results.t[1]
    @test size(densities, 1) == length(config.reactor.composition.species)
    @test all(isfinite, temperatures.tt)
    @test all(isfinite, temperatures.te)
    @test all(isfinite, temperatures.tv)
end

 @testset "Integrate 0D (isothermal)" begin
    config = terra.nitrogen_10ev_config(; isothermal = true)
    temp_case_path = mktempdir()
    config = terra.with_case_path(config, temp_case_path)
    config = terra.with_time(config;
        dt = 5e-12, dt_output = 1e-6, duration = 1e-6, nstep = 500000, method = 2)
    config = terra.with_runtime(config;
        validate_species_against_terra = false,
        print_source_terms = false)

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    results = @time terra.integrate_0d_system(config, initial_state)
    @test results isa terra.ReactorResult

    densities = terra.species_density_matrix(results)
    temperatures = terra.temperature_history(results)
    @test results.t[end] > results.t[1]
    @test size(densities, 1) == length(config.reactor.composition.species)
    @test all(isfinite, temperatures.tt)
    @test all(isfinite, temperatures.te)
    @test maximum(abs.(temperatures.te .- config.reactor.thermal.Te)) <= 1e-6
    @test all(isfinite, temperatures.tv)
end
