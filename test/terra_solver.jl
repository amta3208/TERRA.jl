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

    # If the solver probes intermediate states with small negative components,
    # the wrapper should still compute a meaningful RHS by clamping the state
    # passed to Fortran and applying a Shampine-style positivity rule.
    density_indices = Int[]
    append!(density_indices, collect(layout.vib_range))
    append!(density_indices, collect(layout.elec_range))
    append!(density_indices, collect(layout.sp_range))
    unique!(density_indices)
    @test !isempty(density_indices)
    idx = first(density_indices)

    y_neg = copy(y0)
    y_neg[idx] = -abs(y0[idx]) - 1.0e-30

    density_stop = layout.mom_range.start - 1
    @test density_stop >= 1
    rho_target = sum(@view y_neg[1:density_stop])

    y_proj = copy(y_neg)
    @test y_proj[idx] < 0.0
    y_proj[idx] = 0.0
    rho_clamped = sum(@view y_proj[1:density_stop])
    @test rho_target > 0.0
    @test rho_clamped > 0.0
    scale = rho_target / rho_clamped
    @test scale > 0.0
    @inbounds y_proj[1:density_stop] .*= scale

    du_proj = zeros(length(y0))
    @test_nowarn terra.terra_ode_system!(du_proj, y_proj, p, 0.0)
    @test all(isfinite, du_proj)

    du_neg = zeros(length(y0))
    @test_nowarn terra.terra_ode_system!(du_neg, y_neg, p, 0.0)
    @test all(isfinite, du_neg)
    @test du_neg[idx] >= 0.0

    du_expected = copy(du_proj)
    du_expected[idx] = max(du_expected[idx], 0.0)
    @test du_neg == du_expected
end

@testset "RHS (rhs_api) with residence time terms" begin
    test_case_path = joinpath(@__DIR__, "test_case")
    @test_nowarn reset_and_init!(test_case_path)

    config = terra.nitrogen_10ev_config(; isothermal = false)
    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()

    u_in = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    rt_cfg = terra.ResidenceTimeConfig(; L = 1.0, U_neutral = 1.0, U_ion = 2.0)
    rt = terra._prepare_residence_time_data(layout, config, u_in, rt_cfg)

    @test !isempty(rt.neutral_indices)
    @test !isempty(rt.ion_indices)
    @test !isempty(rt.energy_indices)

    rho_n = sum(rt.u_in[rt.neutral_indices])
    rho_i = sum(rt.u_in[rt.ion_indices])
    expected_inv_tau_energy = (rho_n * rt_cfg.U_neutral + rho_i * rt_cfg.U_ion) / (rho_n + rho_i)
    @test rt.inv_tau_energy ≈ expected_inv_tau_energy

    u = copy(u_in)
    i_neu = rt.neutral_indices[argmax(u_in[rt.neutral_indices])]
    i_ion = rt.ion_indices[argmax(u_in[rt.ion_indices])]
    i_E = rt.energy_indices[1]

    u[i_neu] *= 0.5
    u[i_ion] *= 2.0
    u[i_E] += 1.0

    du_base = zeros(length(u))
    du_rt = zeros(length(u))

    p_base = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.temperatures.Te, layout.nsp)
    )
    p_rt = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.temperatures.Te, layout.nsp),
        residence_time = rt
    )

    @test_nowarn terra.terra_ode_system!(du_base, u, p_base, 0.0)
    @test_nowarn terra.terra_ode_system!(du_rt, u, p_rt, 0.0)

    diff = du_rt .- du_base

    @test diff[i_neu] ≈ (rt.u_in[i_neu] - u[i_neu]) * rt.inv_tau_neutral
    @test diff[i_ion] ≈ (rt.u_in[i_ion] - u[i_ion]) * rt.inv_tau_ion
    @test diff[i_E] ≈ (rt.u_in[i_E] - u[i_E]) * rt.inv_tau_energy
end

@testset "Isothermal Teex: residence time skips rho_eeex" begin
    config = terra.nitrogen_10ev_config(; isothermal = true)
    @test_nowarn reset_and_init!(tempname(); config = config)

    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()

    @test layout.eex_noneq == true
    @test layout.idx_eeex != 0

    u_in = terra.pack_state_vector(layout, state.rho_sp, state.rho_rem;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    rt_cfg = terra.ResidenceTimeConfig(; L = 1.0, U_neutral = 1.0, U_ion = 2.0)
    rt = terra._prepare_residence_time_data(layout, config, u_in, rt_cfg)

    @test !(layout.idx_eeex in rt.energy_indices)

    u = copy(u_in)
    u[layout.idx_eeex] += 123.0

    du_base = zeros(length(u))
    du_rt = zeros(length(u))

    p_base = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.temperatures.Te, layout.nsp)
    )
    p_rt = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.temperatures.Te, layout.nsp),
        residence_time = rt
    )

    @test_nowarn terra.terra_ode_system!(du_base, u, p_base, 0.0)
    @test_nowarn terra.terra_ode_system!(du_rt, u, p_rt, 0.0)

    diff = du_rt .- du_base
    @test diff[layout.idx_eeex] == 0.0
    @test du_rt[layout.idx_eeex] == 0.0
end

@testset "Residence-time override toggle" begin
    rt_cfg = terra.ResidenceTimeConfig(; enabled = true, L = 1.0, U_neutral = 1.0, U_ion = 1.0)
    @test terra._resolve_residence_time(rt_cfg, nothing) === rt_cfg
    @test terra._resolve_residence_time(rt_cfg, false) === nothing

    forced_rt = terra._resolve_residence_time(nothing, true)
    @test forced_rt !== nothing
    @test forced_rt.enabled == true
end

@testset "Residence-time inlet_reactor support" begin
    config = terra.to_config(terra.nitrogen_10ev_config(; isothermal = false))
    @test_nowarn reset_and_init!(tempname(); config = config)

    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()
    u_base = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    inlet_reactor = terra.ReactorConfig(;
        composition = terra.ReactorComposition(;
            species = config.reactor.composition.species,
            mole_fractions = [0.2, 0.65, 0.05, 0.05, 0.05],
            total_number_density = config.reactor.composition.total_number_density),
        thermal = config.reactor.thermal)
    rt_cfg = terra.ResidenceTimeConfig(; L = 1.0, U_neutral = 1.0, U_ion = 2.0,
        inlet_reactor = inlet_reactor)
    rt = terra._prepare_residence_time_data(layout, config, u_base, rt_cfg)

    @test any(abs.(rt.u_in .- u_base) .> 0.0)
end

@testset "Nested Config solver controls" begin
    base = terra.to_config(terra.nitrogen_10ev_config(; isothermal = false))
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
            space = base.numerics.space,
            residence_time = nothing),
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
    @test length(results.time) == config.numerics.solver.saveat_count
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
