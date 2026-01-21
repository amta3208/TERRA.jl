@testset "State Vector Layout (rhs_api)" begin
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    layout = terra.get_api_layout()
    @test layout.neq > 0
    @test layout.neq ==
          layout.n_eq_vib + layout.n_eq_elec + layout.n_eq_sp + layout.n_eq_mom +
          layout.n_eq_energy

    nsp = layout.nsp
    rho_sp = collect(range(1.0, length = nsp))

    rho_ex = nothing
    if layout.is_elec_sts
        rho_ex = zeros(Float64, layout.mnex, nsp)
        @inbounds for isp in 1:nsp
            if layout.ies[isp] == 0
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 1
                    continue
                end
                rho_ex[iex, isp] = 1.0e-3 * (100 * isp + iex)
            end
        end
    end

    rho_vx = nothing
    if layout.n_eq_vib > 0
        rho_vx = zeros(Float64, layout.mnv + 1, layout.mmnex, nsp)
        @inbounds for isp in 1:nsp
            if layout.ies[isp] == 0 || layout.ih[isp] != 2
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 0
                    continue
                end
                mv_isp_iex = layout.mv[isp, iex]
                for ivx in 0:mv_isp_iex
                    rho_vx[ivx + 1, iex, isp] = 1.0e-6 * (10000 * isp + 100 * iex + ivx)
                end
            end
        end
    end

    rho_etot = 123.0
    rho_eeex = layout.eex_noneq ? 7.0 : nothing
    rho_erot = layout.rot_noneq ? 5.0 : nothing
    rho_evib = layout.vib_noneq ? 9.0 : nothing

    y = terra.pack_state_vector(layout, rho_sp, rho_etot;
        rho_ex = rho_ex,
        rho_vx = rho_vx,
        rho_u = layout.nd >= 1 ? 0.0 : nothing,
        rho_v = layout.nd >= 2 ? 0.0 : nothing,
        rho_w = layout.nd >= 3 ? 0.0 : nothing,
        rho_eeex = rho_eeex,
        rho_erot = rho_erot,
        rho_evib = rho_evib)
    @test length(y) == layout.neq

    st = terra.unpack_state_vector(y, layout)
    y2 = terra.pack_state_vector(layout, st.rho_sp, st.rho_etot;
        rho_ex = st.rho_ex,
        rho_vx = st.rho_vx,
        rho_u = layout.nd >= 1 ? st.rho_u : nothing,
        rho_v = layout.nd >= 2 ? st.rho_v : nothing,
        rho_w = layout.nd >= 3 ? st.rho_w : nothing,
        rho_eeex = layout.eex_noneq ? st.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? st.rho_erot : nothing,
        rho_evib = layout.vib_noneq ? st.rho_evib : nothing)
    @test y2 ≈ y
end

@testset "RHS (rhs_api)" begin
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    config = terra.nitrogen_10ev_config(; isothermal = false)
    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()

    y0 = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    du = zeros(length(y0))
    p = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.temperatures.Te, layout.nsp),
        work_y = similar(y0),
        work_dy = similar(y0),
        work_rho_sp = zeros(Float64, layout.nsp),
        work_rho_ex = layout.is_elec_sts ? zeros(Float64, layout.mnex, layout.nsp) :
                      zeros(Float64, 0, 0)
    )

    @test_nowarn terra.terra_ode_system!(du, y0, p, 0.0)
    @test all(isfinite, du)
    @test length(du) == length(y0)
end

@testset "Native Output Generation" begin
    base_config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir(cleanup = false)
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
        write_native_outputs = true,
        print_integration_output = false
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

@testset "Integrate 0D (adiabatic)" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
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

@testset "Benchmark with Fortran Solver - [0D Isothermal Nitrogen 10eV for 10us]" begin
    config = terra.nitrogen_10ev_config(; isothermal = true)
    temp_case_path = mktempdir()
    config = terra.TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = terra.TimeIntegrationConfig(5e-12, 1e-6, 1e-5, 500000, 2),
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
    @test results.temperatures.tt[end]≈756.0 rtol=0.03
    @test results.temperatures.tv[end]≈750.0 rtol=0.01
    @test results.temperatures.te[end]≈115000.0 rtol=0.01

    # Approximate final species densities in CGS (update as needed)
    @test results.species_densities[1, end]≈4.236e-13 rtol=0.03 # N
    @test results.species_densities[2, end]≈4.646e-10 rtol=0.03 # N₂
    @test results.species_densities[3, end]≈4.811e-17 rtol=0.10 # N⁺
    @test results.species_densities[4, end]≈1.083e-13 rtol=0.05 # N₂⁺
    @test results.species_densities[5, end]≈2.123e-18 rtol=0.05 # E⁻
end
