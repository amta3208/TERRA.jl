function _with_reactor_thermal(config::terra.Config;
        Tt = config.reactor.thermal.Tt,
        Tv = config.reactor.thermal.Tv,
        Tee = config.reactor.thermal.Tee,
        Te = config.reactor.thermal.Te)
    thermal = terra.ReactorThermalState(; Tt = Tt, Tv = Tv, Tee = Tee, Te = Te)
    reactor = terra.ReactorConfig(;
        composition = config.reactor.composition,
        thermal = thermal,
    )
    return terra.Config(;
        reactor = reactor,
        models = config.models,
        numerics = config.numerics,
        runtime = config.runtime,
    )
end

 @progress_testset "Initial-state rho_ex handoff" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
    temp_case_path = mktempdir()
    config = terra.with_case_path(config, temp_case_path)
    config = terra.with_runtime(config;
        validate_species_against_terra = false,
        print_source_terms = false)

    @test_nowarn reset_and_init!(temp_case_path; config = config)

    base_state = terra.config_to_initial_state(config)
    state_cache = terra.ReactorStateCache(;
        species = config.reactor.composition.species,
        rho_sp_cgs = base_state.rho_sp,
        rho_ex_cgs = base_state.rho_ex,
    )

    downstream_config_a = _with_reactor_thermal(config;
        Tt = 950.0,
        Tv = 1250.0,
        Tee = 3500.0,
        Te = 82000.0,
    )
    downstream_config_b = _with_reactor_thermal(downstream_config_a; Tee = 90000.0)

    cached_state_a = terra.config_to_initial_state(
        downstream_config_a; state_cache = state_cache)
    cached_state_b = terra.config_to_initial_state(
        downstream_config_b; state_cache = state_cache)
    boltz_state = terra.config_to_initial_state(downstream_config_a)
    boltz_rho_ex_active = boltz_state.rho_ex[:, 1:length(config.reactor.composition.species)]

    @test cached_state_a.state_cache_used == true
    @test cached_state_b.state_cache_used == true
    @test cached_state_a.rho_sp ≈ state_cache.rho_sp_cgs
    @test cached_state_a.rho_ex ≈ state_cache.rho_ex_cgs
    @test cached_state_b.rho_ex ≈ state_cache.rho_ex_cgs
    @test cached_state_a.rho_ex ≈ cached_state_b.rho_ex
    @test cached_state_a.rho_eeex ≈ cached_state_b.rho_eeex
    @test cached_state_a.rho_evib ≈ cached_state_b.rho_evib
    @test !isapprox(cached_state_a.rho_eeex, base_state.rho_eeex)
    @test !isapprox(cached_state_a.rho_evib, base_state.rho_evib)
    @test any(abs.(boltz_rho_ex_active .- cached_state_a.rho_ex) .> 0.0)
end
