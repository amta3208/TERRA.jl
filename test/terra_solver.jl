@testset "State Vector Layout" begin
    # Ensure Fortran-side library is loaded and initialized for dimension queries
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    # Use canonical nitrogen config for consistent dimensions
    config = terra.nitrogen_10ev_config()
    dims = terra.get_state_dimensions(config)

    # Basic dimension sanity
    @test dims.n_species == length(config.species)
    @test prod(dims.rho_ex_size) >= 0
    @test prod(dims.rho_vx_size) >= 0

    # Build consistent arrays and verify pack/unpack round-trip
    rho_sp = collect(1.0:(1.0 + dims.n_species - 1))
    rho_etot = 5.0
    rho_ex = fill(0.1, dims.rho_ex_size)
    rho_vx = fill(0.01, dims.rho_vx_size)

    u = terra.pack_state_vector(rho_sp, rho_etot, dims;
        rho_ex = rho_ex, rho_vx = rho_vx,
        rho_u = 0.0, rho_v = 0.0, rho_w = 0.0,
        rho_erot = 0.0, rho_eeex = 0.0, rho_evib = 0.0)

    expected_len = dims.n_species + 1 + prod(dims.rho_ex_size) + prod(dims.rho_vx_size) + 6
    @test length(u) == expected_len

    state = terra.unpack_state_vector(u, dims)
    @test state.rho_sp ≈ rho_sp
    @test state.rho_etot ≈ rho_etot
    @test state.rho_ex ≈ rho_ex
    @test state.rho_vx ≈ rho_vx
    @test state.rho_u == 0.0 && state.rho_v == 0.0 && state.rho_w == 0.0
    @test state.rho_erot == 0.0 && state.rho_eeex == 0.0 && state.rho_evib == 0.0

    # Derivative packing with and without optional arrays
    du = zeros(length(u))
    derivs = (
        drho_sp = ones(dims.n_species) * 0.1,
        drho_etot = 0.2,
        drho_ex = nothing,
        drho_vx = nothing,
        drho_erot = nothing,
        drho_eeex = nothing,
        drho_evib = nothing
    )
    @test_nowarn terra.pack_derivative_vector!(du, derivs, dims)
    @test du[1:(dims.n_species)] ≈ derivs.drho_sp
    @test du[dims.n_species + 1] ≈ derivs.drho_etot

    derivs_full = (
        drho_sp = ones(dims.n_species) * 0.1,
        drho_etot = 0.2,
        drho_ex = fill(0.01, dims.rho_ex_size),
        drho_vx = fill(0.001, dims.rho_vx_size),
        drho_erot = 0.05,
        drho_eeex = 0.02,
        drho_evib = 0.03
    )
    @test_nowarn terra.pack_derivative_vector!(du, derivs_full, dims)
    @test du[end - 2] ≈ 0.05
    @test du[end - 1] ≈ 0.02
    @test du[end] ≈ 0.03
end

@testset "RHS and Temperatures (initialized Fortran)" begin
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    config = terra.nitrogen_10ev_config()
    state = terra.config_to_initial_state(config)
    dims = terra.get_state_dimensions(config)
    species = config.species
    molecular_weights = terra.get_molecular_weights(species)
    gas_constants = state.gas_constants
    teex_const_vec = fill(config.temperatures.Te, dims.n_species)
    electron_index = findfirst(==("E-"), species)

    @test haskey(state, :rho_rem)
    electron_enthalpy = electron_index === nothing ? 0.0 :
                        state.rho_sp[electron_index] * gas_constants[electron_index] *
                        config.temperatures.Te
    @test isapprox(state.rho_rem, state.rho_etot - state.rho_eeex - electron_enthalpy;
        rtol = 1e-12, atol = 0.0)

    # Pack current state
    u = terra.pack_state_vector(state.rho_sp, state.rho_etot, dims;
        rho_ex = state.rho_ex,
        rho_eeex = state.rho_eeex,
        rho_evib = state.rho_evib)
    du = zeros(length(u))
    p = (
        dimensions = dims,
        config = config,
        rho_etot0 = state.rho_etot,
        molecular_weights = molecular_weights,
        species = species,
        gas_constants = gas_constants,
        teex_const = state.teex_const,
        teex_const_vec = teex_const_vec,
        electron_index = electron_index,
        energy_cache = Ref(state.rho_energy),
        pressure_cache = Ref(state.pressure)
    )

    # Compute RHS; expect finite derivatives
    @test_nowarn terra.terra_ode_system!(du, u, p, 0.0)
    @test all(isfinite, du)
    @test length(du) == length(u)

    # Temperature calculation at current state (charge-balanced electrons optional)
    temps = terra.calculate_temperatures_wrapper(state.rho_sp, state.rho_etot;
        rho_ex = state.rho_ex, rho_eeex = state.rho_eeex, rho_evib = state.rho_evib)
    @test temps.tt > 0 && temps.teex > 0 && temps.tvib > 0
end

@testset "Isothermal TeeX Handling" begin
    temp_case_path = mktempdir()
    config = terra.nitrogen_10ev_config(; isothermal = true)
    @test_nowarn reset_and_init!(temp_case_path; config = config)

    state = terra.config_to_initial_state(config)
    dims = terra.get_state_dimensions(config)
    teex_const_vec_local = fill(config.temperatures.Te, dims.n_species)
    species_iso = config.species
    molecular_weights_iso = terra.get_molecular_weights(species_iso)
    gas_constants_iso = state.gas_constants
    electron_index_iso = findfirst(==("E-"), species_iso)

    u = terra.pack_state_vector(state.rho_sp, state.rho_rem, dims;
        rho_ex = state.rho_ex,
        rho_eeex = state.rho_eeex,
        rho_evib = state.rho_evib)
    du = zeros(length(u))
    p = (
        dimensions = dims,
        config = config,
        rho_etot0 = state.rho_etot,
        molecular_weights = molecular_weights_iso,
        species = species_iso,
        gas_constants = gas_constants_iso,
        teex_const = state.teex_const,
        teex_const_vec = teex_const_vec_local,
        electron_index = electron_index_iso,
        energy_cache = Ref(state.rho_energy),
        pressure_cache = Ref(state.pressure)
    )

    @test_nowarn terra.terra_ode_system!(du, u, p, 0.0)
    @test all(isfinite, du)

    state_unpacked = terra.unpack_state_vector(u, dims)
    energy = terra.reconstruct_energy_components(state_unpacked, config;
        teex_const = config.temperatures.Te,
        teex_vec = teex_const_vec_local,
        molecular_weights = molecular_weights_iso,
        species_names = species_iso,
        gas_constants = gas_constants_iso,
        has_electronic_sts = dims.has_electronic_sts,
        has_vibrational_sts = dims.has_vibrational_sts)
    @test isapprox(energy.rho_rem, state.rho_rem; rtol = 1e-12)
    @test isapprox(energy.rho_eeex, state.rho_eeex; rtol = 1e-12)
    @test isapprox(energy.rho_etot, state.rho_energy; rtol = 1e-6)

    temps = terra.calculate_temperatures_wrapper(state_unpacked.rho_sp, energy.rho_etot;
        rho_ex = state_unpacked.rho_ex,
        rho_eeex = energy.rho_eeex,
        rho_evib = state_unpacked.rho_evib)
    @test isapprox(temps.teex, config.temperatures.Te; rtol = 1e-6)
end

@testset "Native Output Generation" begin
    base_config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir(cleanup = false)
    println(temp_case_path)
    config = terra.TERRAConfig(
        species = base_config.species,
        mole_fractions = base_config.mole_fractions,
        total_number_density = base_config.total_number_density,
        temperatures = base_config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 5e-7, 1000, 2),
        physics = base_config.physics,
        processes = base_config.processes,
        database_path = base_config.database_path,
        case_path = temp_case_path,
        unit_system = base_config.unit_system,
        validate_species_against_terra = false,
        print_source_terms = false,
        write_native_outputs = true
    )

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
    @test results.time[end] >= results.time[1]

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
        print_source_terms = false
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

@testset "End-to-end Example (0D Adiabatic Nitrogen 10eV)" begin
    # Run the high-level example wrapper and verify structure and success
    results = @time terra.nitrogen_10ev_example()

    @test results.success == true
    @test length(results.time) >= 2
    @test size(results.species_densities, 1) == 5
    @test all(isfinite, results.temperatures.tt)
    @test all(isfinite, results.temperatures.te)
    @test all(isfinite, results.temperatures.tv)
    @test terra.validate_results(results)

    # Approximate final temperature values (update as needed)
    @test results.temperatures.tt[end]≈750.13 rtol=0.03
    @test results.temperatures.tv[end]≈758.49 rtol=0.03
    @test results.temperatures.te[end]≈2300.30 rtol=0.03

    # Approximate final species densities in CGS (update as needed)
    @test results.species_densities[1, end]≈4.341e-14 rtol=0.05 # N
    @test results.species_densities[2, end]≈4.650e-10 rtol=0.03 # N₂
    @test results.species_densities[3, end]≈4.889e-19 rtol=0.10 # N⁺
    @test results.species_densities[4, end]≈4.760e-14 rtol=0.05 # N₂⁺
    @test results.species_densities[5, end]≈9.322e-19 rtol=0.05 # E⁻
end

@testset "End-to-end Example (0D Isothermal Nitrogen 10eV)" begin
    config = terra.nitrogen_10ev_config(; isothermal = true)
    temp_case_path = mktempdir()
    config = terra.TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 1e-4, 500000, 2),
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = temp_case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = false,
        print_source_terms = false
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
    @test results.temperatures.tt[end]≈1252.0 rtol=0.03
    @test results.temperatures.tv[end]≈751.1 rtol=0.03
    @test results.temperatures.te[end]≈115000.0 rtol=0.01

    # Approximate final species densities in CGS (update as needed)
    @test results.species_densities[1, end]≈4.988e-11 rtol=0.03 # N
    @test results.species_densities[2, end]≈4.107e-10 rtol=0.03 # N₂
    @test results.species_densities[3, end]≈8.219e-13 rtol=0.10 # N⁺
    @test results.species_densities[4, end]≈3.754e-12 rtol=0.05 # N₂⁺
    @test results.species_densities[5, end]≈1.057e-16 rtol=0.05 # E⁻
end
