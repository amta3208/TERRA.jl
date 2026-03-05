const CHAIN_SEGMENT_DIRNAME = "chain_segments"

function _segment_case_path(base_case_path::AbstractString, segment_index::Integer)
    return normpath(joinpath(base_case_path, CHAIN_SEGMENT_DIRNAME,
        @sprintf("segment_%04d", segment_index)))
end

function _failed_simulation_result(message::AbstractString)
    return SimulationResult(
        Float64[],
        zeros(0, 0),
        (tt = Float64[], te = Float64[], tv = Float64[]),
        Float64[],
        nothing,
        false,
        String(message)
    )
end

function _build_chain_models(base_models::ModelConfig)
    physics = base_models.physics
    chain_physics = PhysicsConfig(
        bbh_model = physics.bbh_model,
        esc_model = physics.esc_model,
        ar_et_model = physics.ar_et_model,
        eex_noneq = physics.eex_noneq,
        ev_relax_set = physics.ev_relax_set,
        et_relax_set = physics.et_relax_set,
        radiation_length = physics.radiation_length,
        get_electron_density_by_charge_balance = physics.get_electron_density_by_charge_balance,
        min_sts_frac = physics.min_sts_frac,
        is_isothermal_teex = true,
        energy_loss_per_eii = physics.energy_loss_per_eii,
    )
    return ModelConfig(; physics = chain_physics, processes = base_models.processes)
end

@inline function _resolve_segment_tee(inlet_reactor::ReactorConfig,
        te_profile::Float64,
        marching::AxialMarchingConfig)
    if marching.tee_policy == :match_te
        return te_profile
    end
    return inlet_reactor.thermal.Tee
end

function _build_segment_reactor(inlet_reactor::ReactorConfig,
        te_profile::Float64,
        marching::AxialMarchingConfig)
    tt = marching.override_tt_K === nothing ? inlet_reactor.thermal.Tt : marching.override_tt_K
    tv = marching.override_tv_K === nothing ? inlet_reactor.thermal.Tv : marching.override_tv_K
    tee = _resolve_segment_tee(inlet_reactor, te_profile, marching)

    thermal = ReactorThermalState(; Tt = tt, Tv = tv, Tee = tee, Te = te_profile)
    return ReactorConfig(; composition = inlet_reactor.composition, thermal = thermal)
end

function _build_segment_residence_time(base_config::Config,
        profile::AxialChainProfile,
        segment_index::Integer,
        inlet_reactor::ReactorConfig)
    base_rt = base_config.numerics.residence_time
    u_energy = base_rt === nothing ? nothing : base_rt.U_energy

    return ResidenceTimeConfig(
        enabled = true,
        L = profile.dx_m[segment_index],
        U_neutral = profile.u_neutral_m_s[segment_index],
        U_ion = profile.u_ion_m_s[segment_index],
        U_energy = u_energy,
        inlet_reactor = inlet_reactor,
    )
end

function _build_segment_runtime(base_runtime::RuntimeConfig,
        segment_case_path::AbstractString)
    return RuntimeConfig(
        database_path = base_runtime.database_path,
        case_path = String(segment_case_path),
        unit_system = base_runtime.unit_system,
        validate_species_against_terra = base_runtime.validate_species_against_terra,
        print_source_terms = base_runtime.print_source_terms,
        write_native_outputs = base_runtime.write_native_outputs,
        print_integration_output = base_runtime.print_integration_output,
    )
end

function _build_chain_segment_config(base_config::Config,
        profile::AxialChainProfile,
        segment_index::Integer,
        inlet_reactor::ReactorConfig,
        marching::AxialMarchingConfig)
    segment_case_path = _segment_case_path(base_config.runtime.case_path, segment_index)
    mkpath(segment_case_path)

    reactor = _build_segment_reactor(inlet_reactor, profile.te_K[segment_index], marching)
    models = _build_chain_models(base_config.models)
    rt = _build_segment_residence_time(base_config, profile, segment_index, inlet_reactor)
    numerics = NumericsConfig(
        time = base_config.numerics.time,
        solver = base_config.numerics.solver,
        space = base_config.numerics.space,
        residence_time = rt,
    )
    runtime = _build_segment_runtime(base_config.runtime, segment_case_path)

    return Config(;
        reactor = reactor,
        models = models,
        numerics = numerics,
        runtime = runtime,
    )
end

function _to_cgs_density(density::AbstractVector{<:Real}, unit_system::Symbol)
    rho = Float64.(density)
    if unit_system == :SI
        return convert_density_si_to_cgs(rho)
    end
    return rho
end

function _total_number_density_from_cgs_density(rho_cgs::Vector{Float64}, molecular_weights::Vector{Float64})
    return sum(rho_cgs .* AVOGADRO ./ molecular_weights)
end

function _extract_segment_endpoint_reactor(base_config::Config,
        segment_result::SimulationResult)
    if !segment_result.success
        throw(ArgumentError("Cannot extract endpoint reactor from an unsuccessful segment result."))
    end
    if isempty(segment_result.time)
        throw(ArgumentError("Cannot extract endpoint reactor: segment result contains no time samples."))
    end

    species = base_config.reactor.composition.species
    molecular_weights = get_molecular_weights(species)

    rho_endpoint = segment_result.species_densities[:, end]
    rho_endpoint_cgs = _to_cgs_density(rho_endpoint, base_config.runtime.unit_system)
    mole_fractions = mass_densities_to_mole_fractions(rho_endpoint_cgs, molecular_weights)

    n_total_cgs = _total_number_density_from_cgs_density(rho_endpoint_cgs, molecular_weights)
    n_total = base_config.runtime.unit_system == :SI ?
              convert_number_density_cgs_to_si(n_total_cgs) : n_total_cgs

    tt = segment_result.temperatures.tt[end]
    tv = segment_result.temperatures.tv[end]
    te = segment_result.temperatures.te[end]

    return ReactorConfig(
        composition = ReactorComposition(
            species,
            mole_fractions,
            n_total,
        ),
        thermal = ReactorThermalState(tt, tv, te, te),
    )
end

function _build_chain_diagnostics(base_config::Config,
        profile::AxialChainProfile,
        marching::AxialMarchingConfig,
        segment_end_reactors::Vector{ReactorConfig})
    diagnostics = Dict{String, Any}(
        "handoff_mode" => String(marching.handoff_mode),
        "termination_mode" => String(marching.termination_mode),
        "tee_policy" => String(marching.tee_policy),
        "override_tt_K" => marching.override_tt_K,
        "override_tv_K" => marching.override_tv_K,
    )

    if !isempty(profile.diagnostics)
        diagnostics["input_profile_diagnostics"] = Dict{String, Vector{Float64}}(
            key => copy(values) for (key, values) in pairs(profile.diagnostics)
        )
    end

    simulated_n_total = [reactor.composition.total_number_density for reactor in segment_end_reactors]
    diagnostics["simulated_total_number_density"] = simulated_n_total

    if haskey(profile.diagnostics, "n_total_m3")
        profile_n_total = if base_config.runtime.unit_system == :SI
            copy(profile.diagnostics["n_total_m3"])
        else
            profile.diagnostics["n_total_m3"] .* 1e-6
        end
        diagnostics["input_total_number_density"] = profile_n_total
        diagnostics["total_number_density_abs_error"] = abs.(simulated_n_total .- profile_n_total)
    end

    return diagnostics
end

"""
$(SIGNATURES)

Solve a steady axial chain of CSTR segments using the TERRA 0D solver.

This MVP implementation supports:
- `handoff_mode = :reinitialize`
- `termination_mode = :final_time`

# Arguments
- `config::Config`: Base TERRA config for solver controls and inlet state
- `profile::AxialChainProfile`: Axial chain profile

# Keyword Arguments
- `marching::AxialMarchingConfig`: Axial marching controls

# Returns
- `ChainSimulationResult`
"""
function solve_terra_chain_steady(config::Config,
        profile::AxialChainProfile;
        marching::AxialMarchingConfig = AxialMarchingConfig())
    validate_axial_chain_profile(profile)
    validate_axial_marching_config(marching)

    # Chain marching repeatedly reinitializes the Fortran API; keep MPI
    # finalization disabled to avoid shutting down MPI mid-process.
    try
        if !is_terra_loaded()
            load_terra_library!()
        end
        set_api_finalize_mpi_wrapper(false)
    catch
        # Older libraries may not expose this control; proceed best-effort.
    end

    n_segments = length(profile.z_m)
    segment_end_reactors = [config.reactor for _ in 1:n_segments]
    segment_success = fill(false, n_segments)
    segment_messages = fill("Not executed.", n_segments)
    segment_results = marching.store_segment_histories ?
                      [_failed_simulation_result("Not executed.") for _ in 1:n_segments] :
                      nothing

    inlet_reactor = config.reactor

    for k in 1:n_segments
        local_result = nothing
        try
            segment_config = _build_chain_segment_config(config, profile, k, inlet_reactor, marching)
            initialize_terra(segment_config, segment_config.runtime.case_path)
            local_result = solve_terra_0d(segment_config;
                residence_time = segment_config.numerics.residence_time,
                use_residence_time = true)
        catch e
            segment_messages[k] = "Segment setup/integration threw exception: $(e)"
            diagnostics = _build_chain_diagnostics(config, profile, marching, segment_end_reactors)
            return ChainSimulationResult(
                copy(profile.z_m),
                copy(profile.dx_m),
                copy(profile.te_K),
                copy(profile.u_neutral_m_s),
                copy(profile.u_ion_m_s),
                segment_end_reactors,
                segment_success,
                segment_messages,
                segment_results,
                diagnostics,
                false,
                k,
                "Chain failed at segment $k during setup/integration.",
            )
        end

        if segment_results !== nothing
            segment_results[k] = local_result
        end

        segment_success[k] = local_result.success
        segment_messages[k] = local_result.message

        if !local_result.success
            diagnostics = _build_chain_diagnostics(config, profile, marching, segment_end_reactors)
            return ChainSimulationResult(
                copy(profile.z_m),
                copy(profile.dx_m),
                copy(profile.te_K),
                copy(profile.u_neutral_m_s),
                copy(profile.u_ion_m_s),
                segment_end_reactors,
                segment_success,
                segment_messages,
                segment_results,
                diagnostics,
                false,
                k,
                "Chain failed at segment $k: $(local_result.message)",
            )
        end

        endpoint_reactor = _extract_segment_endpoint_reactor(config, local_result)
        segment_end_reactors[k] = endpoint_reactor
        inlet_reactor = endpoint_reactor
    end

    diagnostics = _build_chain_diagnostics(config, profile, marching, segment_end_reactors)
    return ChainSimulationResult(
        copy(profile.z_m),
        copy(profile.dx_m),
        copy(profile.te_K),
        copy(profile.u_neutral_m_s),
        copy(profile.u_ion_m_s),
        segment_end_reactors,
        segment_success,
        segment_messages,
        segment_results,
        diagnostics,
        true,
        nothing,
        "Chain integration completed successfully.",
    )
end

"""
$(SIGNATURES)

Load a profile from disk and solve the steady axial chain of CSTR segments.
"""
function solve_terra_chain_steady(config::Config,
        profile_path::AbstractString;
        marching::AxialMarchingConfig = AxialMarchingConfig())
    profile = load_chain_profile(profile_path)
    return solve_terra_chain_steady(config, profile; marching = marching)
end
