@testset "Integrate 0D (adiabatic)" begin
    # Initialize using the config-driven input to ensure the selected
    # database and options are honored (rather than a stale case file).

    config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    # Shorten integration to keep tests fast
    config = terra.TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 1e-6, 500000, 2),
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = temp_case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = false,
        print_source_terms = false,
        print_integration_output = false
    )

    # Initialize the Fortran API using a temporary case generated from this config
    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    results = @time terra.integrate_0d_system(config, initial_state)
    @test results.time[end] > results.time[1]
    @test size(results.species_densities, 1) == length(config.species)
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test all(isfinite, results.temperatures.tv)
end

@testset "Integrate 0D (isothermal)" begin
    config = terra.nitrogen_10ev_config(; isothermal = true)
    temp_case_path = mktempdir()
    config = terra.TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 1e-6, 500000, 2),
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = temp_case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = false,
        print_source_terms = false,
        print_integration_output = false
    )

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    results = @time terra.integrate_0d_system(config, initial_state)

    @test results.time[end] > results.time[1]
    @test size(results.species_densities, 1) == length(config.species)
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test maximum(abs.(results.temperatures.te .- config.temperatures.Te)) <= 1e-6
    @test all(isfinite, results.temperatures.tv)
end

@testset "Benchmark with Fortran Solver - [0D Adiabatic Nitrogen 10eV]" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    config = terra.TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 1e-3, 500000, 2),
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = temp_case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = false,
        print_source_terms = false,
        print_integration_output = false
    )

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    results = @time terra.integrate_0d_system(config, initial_state)

    @test results.success == true
    @test length(results.time) >= 2
    @test size(results.species_densities, 1) == 5
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test all(isfinite, results.temperatures.tv)
    @test terra.validate_results(results)

    # Approximate final temperature values (update as needed)
    @test results.temperatures.tt[end]≈754.6 rtol=0.03
    @test results.temperatures.tv[end]≈759.2 rtol=0.03
    @test results.temperatures.te[end]≈2474.0 rtol=0.03

    # Approximate final species densities in CGS (update as needed)
    @test results.species_densities[1, end]≈4.182e-14 rtol=0.05 # N
    @test results.species_densities[2, end]≈4.650e-10 rtol=0.03 # N₂
    @test results.species_densities[3, end]≈5.513e-19 rtol=0.10 # N⁺
    @test results.species_densities[4, end]≈4.943e-14 rtol=0.05 # N₂⁺
    @test results.species_densities[5, end]≈9.680e-19 rtol=0.05 # E⁻
end

@testset "Benchmark with Fortran Solver - [0D Isothermal Nitrogen 10eV for 50us]" begin
    config = terra.nitrogen_10ev_config(; isothermal = true)
    temp_case_path = mktempdir()
    config = terra.TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 5e-5, 500000, 2),
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = temp_case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = false,
        print_source_terms = false,
        print_integration_output = false
    )

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    initial_state = terra.config_to_initial_state(config)
    results = @time terra.integrate_0d_system(config, initial_state)

    @test results.time[end] > results.time[1]
    @test size(results.species_densities, 1) == length(config.species)
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test maximum(abs.(results.temperatures.te .- config.temperatures.Te)) <= 1e-6
    @test all(isfinite, results.temperatures.tv)
    @test all(isfinite, results.total_energy)

    # Approximate final temperature values (update as needed)
    @test results.temperatures.tt[end]≈1782.0 rtol=0.05
    @test results.temperatures.tv[end]≈751.6 rtol=0.01
    @test results.temperatures.te[end]≈115000.0 rtol=0.01

    # Approximate final species densities in CGS (update as needed)
    @test results.species_densities[1, end]≈4.461e-11 rtol=0.03 # N
    @test results.species_densities[2, end]≈4.100e-10 rtol=0.03 # N₂
    @test results.species_densities[3, end]≈9.969e-13 rtol=0.10 # N⁺
    @test results.species_densities[4, end]≈9.487e-12 rtol=0.03 # N₂⁺
    @test results.species_densities[5, end]≈2.248e-16 rtol=0.05 # E⁻
end
