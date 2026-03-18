@testset "RHS (rhs_api) with species-resolved residence time terms" begin
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

    u_species = Dict("N" => 1.0, "N2" => 1.5, "N+" => 2.0, "N2+" => 2.5)
    rt_cfg = terra.ResidenceTimeConfig(; L = 1.0, U_species = u_species)
    rt = terra._prepare_residence_time_data(layout, config, u_in, rt_cfg)

    @test sort(rt.species_order) == ["N", "N+", "N2", "N2+"]
    @test !isempty(rt.species_indices["N"])
    @test !isempty(rt.species_indices["N2"])
    @test !isempty(rt.species_indices["N+"])
    @test !isempty(rt.species_indices["N2+"])
    @test !isempty(rt.energy_indices)

    rho_weighted_u = 0.0
    rho_total = 0.0
    for (i, species_name) in pairs(config.reactor.composition.species)
        species_name == "E-" && continue
        rho_weighted_u += state.rho_sp[i] * u_species[species_name]
        rho_total += state.rho_sp[i]
    end
    expected_inv_tau_energy = rho_weighted_u / rho_total
    @test rt.inv_tau_energy ≈ expected_inv_tau_energy

    u = copy(u_in)
    i_n = first(rt.species_indices["N"])
    i_n2 = first(rt.species_indices["N2"])
    i_n_plus = first(rt.species_indices["N+"])
    i_n2_plus = first(rt.species_indices["N2+"])
    i_E = rt.energy_indices[1]

    u[i_n] *= 0.5
    u[i_n2] *= 0.5
    u[i_n_plus] *= 2.0
    u[i_n2_plus] *= 2.0
    u[i_E] += 1.0

    du_base = zeros(length(u))
    du_rt = zeros(length(u))

    p_base = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp)
    )
    p_rt = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp),
        sources = (residence_time = rt,)
    )

    @test_nowarn terra.terra_ode_system!(du_base, u, p_base, 0.0)
    @test_nowarn terra.terra_ode_system!(du_rt, u, p_rt, 0.0)

    diff = du_rt .- du_base

    @test diff[i_n] ≈ (rt.u_in[i_n] - u[i_n]) * rt.inv_tau_species["N"] atol=1e-20
    @test diff[i_n2] ≈ (rt.u_in[i_n2] - u[i_n2]) * rt.inv_tau_species["N2"] atol=1e-20
    @test diff[i_n_plus] ≈ (rt.u_in[i_n_plus] - u[i_n_plus]) * rt.inv_tau_species["N+"] atol=1e-20
    @test diff[i_n2_plus] ≈ (rt.u_in[i_n2_plus] - u[i_n2_plus]) * rt.inv_tau_species["N2+"] atol=1e-20
    @test diff[i_E] ≈ (rt.u_in[i_E] - u[i_E]) * rt.inv_tau_energy

    for species_name in rt.species_order
        for i in rt.species_indices[species_name]
            expected = (rt.u_in[i] - u[i]) * rt.inv_tau_species[species_name]
            @test diff[i] ≈ expected atol=1e-20
        end
    end
end

@testset "Residence time rejects missing or extra species keys" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
    @test_nowarn reset_and_init!(tempname(); config = config)

    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()
    u_in = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    missing_rt = terra.ResidenceTimeConfig(;
        L = 1.0,
        U_species = Dict("N" => 1.0, "N2" => 1.0, "N+" => 2.0))
    @test_throws ArgumentError terra._prepare_residence_time_data(layout, config, u_in, missing_rt)

    extra_rt = terra.ResidenceTimeConfig(;
        L = 1.0,
        U_species = Dict(
            "N" => 1.0,
            "N2" => 1.0,
            "N+" => 2.0,
            "N2+" => 2.0,
            "Ar" => 3.0,
        ))
    @test_throws ArgumentError terra._prepare_residence_time_data(layout, config, u_in, extra_rt)
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

    rt_cfg = terra.ResidenceTimeConfig(;
        L = 1.0,
        U_species = Dict("N" => 1.0, "N2" => 1.5, "N+" => 2.0, "N2+" => 2.5))
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
        teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp)
    )
    p_rt = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp),
        sources = (residence_time = rt,)
    )

    @test_nowarn terra.terra_ode_system!(du_base, u, p_base, 0.0)
    @test_nowarn terra.terra_ode_system!(du_rt, u, p_rt, 0.0)

    diff = du_rt .- du_base
    @test diff[layout.idx_eeex] == 0.0
    @test du_rt[layout.idx_eeex] == 0.0
end

@testset "Source-term preparation respects empty and populated sources" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
    @test_nowarn reset_and_init!(tempname(); config = config)

    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()
    u_base = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    rt_cfg = terra.ResidenceTimeConfig(;
        enabled = true,
        L = 1.0,
        U_species = Dict("N" => 1.0, "N2" => 1.0, "N+" => 1.0, "N2+" => 1.0))

    prepared = terra._prepare_source_terms_data(
        layout, config, u_base, terra.SourceTermsConfig(; residence_time = rt_cfg))
    @test prepared.residence_time !== nothing

    disabled = terra._prepare_source_terms_data(
        layout, config, u_base, terra.SourceTermsConfig())
    @test disabled.residence_time === nothing
end

@testset "Residence-time inlet_reactor support" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
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
    rt_cfg = terra.ResidenceTimeConfig(;
        L = 1.0,
        U_species = Dict("N" => 1.0, "N2" => 1.5, "N+" => 2.0, "N2+" => 2.5),
        inlet_reactor = inlet_reactor)
    rt = terra._prepare_residence_time_data(layout, config, u_base, rt_cfg)

    @test any(abs.(rt.u_in .- u_base) .> 0.0)
end

@testset "Residence-time inlet state cache preserves rho_ex handoff" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    config = terra.with_case_path(config, temp_case_path)
    config = terra.with_time(config;
        dt = 5e-12, dt_output = 5e-7, duration = 5e-7, nstep = 200000, method = 2)
    config = terra.with_runtime(config;
        validate_species_against_terra = false,
        print_source_terms = false,
        write_native_outputs = false,
        print_integration_output = false)

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()
    u_base = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    _, final_state_cache = terra._solve_terra_0d_internal(config)
    @test final_state_cache !== nothing
    @test final_state_cache.rho_ex_cgs !== nothing

    rt_cfg = terra.ResidenceTimeConfig(;
        L = 1.0,
        U_species = Dict("N" => 1.0, "N2" => 1.5, "N+" => 2.0, "N2+" => 2.5),
        inlet_reactor = config.reactor)

    rt_boltz = terra._prepare_residence_time_data(layout, config, u_base, rt_cfg)
    rt_cached = terra._prepare_residence_time_data(
        layout, config, u_base, rt_cfg;
        inlet_state_cache = final_state_cache)

    boltz_inlet_state = terra.unpack_state_vector(rt_boltz.u_in, layout)
    cached_inlet_state = terra.unpack_state_vector(rt_cached.u_in, layout)

    @test any(abs.(rt_cached.u_in .- rt_boltz.u_in) .> 0.0)
    @test cached_inlet_state.rho_ex ≈ final_state_cache.rho_ex_cgs
    @test any(abs.(cached_inlet_state.rho_ex .- boltz_inlet_state.rho_ex) .> 0.0)
end
