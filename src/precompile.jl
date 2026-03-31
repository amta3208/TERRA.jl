const _PRECOMPILE_0D_DT = 5e-12
const _PRECOMPILE_0D_DT_OUTPUT = 1e-8
const _PRECOMPILE_0D_DURATION = 1e-8
const _PRECOMPILE_0D_NSTEP = 1000
const _PRECOMPILE_0D_SAVEAT_COUNT = 8
const _PRECOMPILE_CHAIN_TT = 500.0
const _PRECOMPILE_CHAIN_TV = 500.0
const _PRECOMPILE_CHAIN_DX_M = 0.01

function _precompile_chain_profile(config::Config)
    total_number_density_m3 = config.runtime.unit_system == :SI ?
                              config.reactor.composition.total_number_density :
                              convert_number_density_cgs_to_si(config.reactor.composition.total_number_density)
    inlet_composition = ChainProfileInletComposition(;
                                                     species = config.reactor.composition.species,
                                                     mole_fractions = [0.2, 0.65, 0.05, 0.05, 0.05],
                                                     total_number_density_m3 = total_number_density_m3)
    inlet_thermal = ReactorThermalState(; Tt = _PRECOMPILE_CHAIN_TT,
                                        Tv = _PRECOMPILE_CHAIN_TV,
                                        Tee = config.reactor.thermal.Te,
                                        Te = config.reactor.thermal.Te)
    inlet = ChainProfileInlet(; composition = inlet_composition,
                              thermal = inlet_thermal,
                              source_compact_index = 1)
    return AxialChainProfile(; z_m = [0.0],
                             dx_m = [_PRECOMPILE_CHAIN_DX_M],
                             te_K = [config.reactor.thermal.Te],
                             species_u_m_s = Dict("N" => [200.0],
                                                  "N2" => [180.0],
                                                  "N+" => [18000.0],
                                                  "N2+" => [19000.0]),
                             inlet = inlet)
end

function _native_precompile_ready()
    database_path = _nitrogen_10ev_database_path()
    chemistry_path = joinpath(database_path, "chemistry.dat")

    if !isdir(database_path) || !isfile(chemistry_path)
        return false
    end

    try
        return isfile(resolve_terra_library_path())
    catch
        return false
    end
end

function _precompile_0d_config(isothermal::Bool)
    config = nitrogen_10ev_config(; isothermal = isothermal)
    config = with_case_path(config, mktempdir())
    config = with_time(config;
                       dt = _PRECOMPILE_0D_DT,
                       dt_output = _PRECOMPILE_0D_DT_OUTPUT,
                       duration = _PRECOMPILE_0D_DURATION,
                       nstep = _PRECOMPILE_0D_NSTEP,
                       method = config.numerics.time.method)
    config = with_runtime(config;
                          validate_species_against_terra = false,
                          print_source_terms = false,
                          write_native_state_files = false)
    config = with_logging(config;
                          console_mode = :quiet,
                          progress_mode = :off,
                          native_stream_mode = :off,
                          integration_detail_mode = :off,
                          chain_detail_mode = :off)

    solver = ODESolverConfig(; reltol = config.numerics.solver.reltol,
                             abstol_density = config.numerics.solver.abstol_density,
                             saveat_count = _PRECOMPILE_0D_SAVEAT_COUNT,
                             ramp_understep_ratio = config.numerics.solver.ramp_understep_ratio,
                             ramp_history_steps = config.numerics.solver.ramp_history_steps)

    return Config(; reactor = config.reactor,
                  models = config.models,
                  sources = config.sources,
                  numerics = NumericsConfig(; time = config.numerics.time,
                                            solver = solver,
                                            space = config.numerics.space),
                  runtime = config.runtime)
end

function _run_native_0d_precompile_case(isothermal::Bool)
    config = _precompile_0d_config(isothermal)
    try
        if !is_terra_loaded()
            load_terra_library!()
        end
        try
            set_api_finalize_mpi_wrapper(false)
        catch
            # Older libraries may not expose this control.
        end

        initialize_terra(config)
        results = solve_terra_0d(config)
        results.success || error("Native precompile warmup failed: $(results.message)")
        return nothing
    finally
        try
            finalize_terra()
        catch
            try
                finalize_api_wrapper()
            catch
            end
            try
                close_terra_library()
            catch
            end
        end

        if isdir(config.runtime.case_path)
            rm(config.runtime.case_path; recursive = true, force = true)
        end
    end
end

function _run_native_chain_precompile_case()
    config = _precompile_0d_config(false)
    profile = _precompile_chain_profile(config)
    marching = AxialMarchingConfig()
    try
        results = solve_terra_chain_steady(config, profile; marching = marching)
        results.success || error("Native chain precompile warmup failed: $(results.message)")
        return nothing
    finally
        try
            finalize_terra()
        catch
            try
                finalize_api_wrapper()
            catch
            end
            try
                close_terra_library()
            catch
            end
        end

        if isdir(config.runtime.case_path)
            rm(config.runtime.case_path; recursive = true, force = true)
        end
    end
end

@setup_workload begin
    native_precompile_ready = _native_precompile_ready()

    @compile_workload begin
        adiabatic_config = _precompile_0d_config(false)
        isothermal_config = _precompile_0d_config(true)
        chain_config = _precompile_0d_config(false)
        chain_profile = _precompile_chain_profile(chain_config)
        chain_marching = AxialMarchingConfig()

        if native_precompile_ready
            try
                _run_native_0d_precompile_case(false)
            catch
            end

            try
                _run_native_0d_precompile_case(true)
            catch
            end

            try
                _run_native_chain_precompile_case()
            catch
            end
        end

        for config in (adiabatic_config, isothermal_config, chain_config)
            if isdir(config.runtime.case_path)
                rm(config.runtime.case_path; recursive = true, force = true)
            end
        end
    end
end
