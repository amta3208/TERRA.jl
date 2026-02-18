@testset "RHS (rhs_api) with residence time terms" begin
    test_case_path = TEST_CASE_PATH
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
