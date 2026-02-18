function _build_residence_time_continuity_indices(layout::ApiLayout)
    neutral_indices = Int[]
    ion_indices = Int[]

    idx = 1

    # Vibrational STS entries (if any)
    if layout.n_eq_vib > 0
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0 || layout.ih[isp] != 2
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 0
                    continue
                end
                mv_isp_iex = layout.mv[isp, iex]
                for _ivx in 0:mv_isp_iex
                    if isp != layout.esp
                        if layout.ie[isp] == 0
                            push!(neutral_indices, idx)
                        else
                            push!(ion_indices, idx)
                        end
                    end
                    idx += 1
                end
            end
        end
    end

    # Electronic STS entries
    if layout.is_elec_sts
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 1
                    continue
                end
                if isp != layout.esp
                    if layout.ie[isp] == 0
                        push!(neutral_indices, idx)
                    else
                        push!(ion_indices, idx)
                    end
                end
                idx += 1
            end
        end
    end

    # Species-density entries for non-electronic-specific species (except charge-balanced electrons)
    @inbounds for isp in 1:(layout.nsp)
        if layout.ies[isp] == 1
            continue
        end
        if layout.get_electron_density_by_charge_balance && isp == layout.esp
            continue
        end
        # Ignore explicit electrons for now (when charge balance is off).
        if isp != layout.esp
            if layout.ie[isp] == 0
                push!(neutral_indices, idx)
            else
                push!(ion_indices, idx)
            end
        end
        idx += 1
    end

    @assert idx==layout.mom_range.start "Internal error: continuity index mismatch while building residence-time indices"

    return (neutral_indices = neutral_indices, ion_indices = ion_indices)
end

function _build_residence_time_energy_indices(layout::ApiLayout, is_isothermal_teex::Bool)
    indices = collect(layout.energy_range)
    if is_isothermal_teex && layout.idx_eeex != 0
        filter!(i -> i != layout.idx_eeex, indices)
    end
    return indices
end

@inline function _apply_residence_time_term!(du::Vector{Float64}, u::Vector{Float64}, rt)
    if rt === nothing
        return nothing
    end
    u_in = rt.u_in
    @inbounds for i in rt.neutral_indices
        du[i] += (u_in[i] - u[i]) * rt.inv_tau_neutral
    end
    @inbounds for i in rt.ion_indices
        du[i] += (u_in[i] - u[i]) * rt.inv_tau_ion
    end
    @inbounds for i in rt.energy_indices
        du[i] += (u_in[i] - u[i]) * rt.inv_tau_energy
    end
    return nothing
end

function _prepare_residence_time_data(layout::ApiLayout, config::Config,
        u0::Vector{Float64}, rt_cfg::ResidenceTimeConfig)
    inv_tau_neutral = rt_cfg.U_neutral / rt_cfg.L
    inv_tau_ion = rt_cfg.U_ion / rt_cfg.L

    continuity = _build_residence_time_continuity_indices(layout)
    neutral_indices = continuity.neutral_indices
    ion_indices = continuity.ion_indices
    energy_indices = _build_residence_time_energy_indices(
        layout, config.models.physics.is_isothermal_teex)

    u_in = if rt_cfg.inlet_reactor === nothing
        copy(u0)
    else
        inlet_reactor = rt_cfg.inlet_reactor
        if inlet_reactor.composition.species != config.reactor.composition.species
            error("ResidenceTimeConfig inlet_reactor species must match the initialized TERRA species ordering.")
        end
        inlet_config = Config(;
            reactor = inlet_reactor,
            models = config.models,
            numerics = config.numerics,
            runtime = config.runtime)

        inlet_state = config_to_initial_state(inlet_config)
        energy_scalar_in = config.models.physics.is_isothermal_teex ? inlet_state.rho_rem :
                           inlet_state.rho_energy
        rho_ex_in = layout.is_elec_sts ? inlet_state.rho_ex : nothing
        rho_eeex_in = layout.eex_noneq ? inlet_state.rho_eeex : nothing
        rho_erot_in = layout.rot_noneq ? 0.0 : nothing
        rho_evib_in = layout.vib_noneq ? inlet_state.rho_evib : nothing

        pack_state_vector(layout, inlet_state.rho_sp, energy_scalar_in;
            rho_ex = rho_ex_in,
            rho_u = layout.nd >= 1 ? 0.0 : nothing,
            rho_v = layout.nd >= 2 ? 0.0 : nothing,
            rho_w = layout.nd >= 3 ? 0.0 : nothing,
            rho_eeex = rho_eeex_in,
            rho_erot = rho_erot_in,
            rho_evib = rho_evib_in)
    end

    U_energy = if rt_cfg.U_energy === nothing
        rho_n = 0.0
        @inbounds for i in neutral_indices
            rho_n += u_in[i]
        end
        rho_i = 0.0
        @inbounds for i in ion_indices
            rho_i += u_in[i]
        end
        rho_total = rho_n + rho_i
        if !(rho_total > 0.0) || !isfinite(rho_total)
            error("ResidenceTimeConfig: cannot infer U_energy from inlet state (rho_n+rho_i=$rho_total).")
        end
        (rho_n * rt_cfg.U_neutral + rho_i * rt_cfg.U_ion) / rho_total
    else
        rt_cfg.U_energy
    end

    inv_tau_energy = U_energy / rt_cfg.L

    return (
        u_in = u_in,
        neutral_indices = neutral_indices,
        ion_indices = ion_indices,
        energy_indices = energy_indices,
        inv_tau_neutral = inv_tau_neutral,
        inv_tau_ion = inv_tau_ion,
        inv_tau_energy = inv_tau_energy
    )
end

function _prepare_residence_time_data(layout::ApiLayout, config::TERRAConfig,
        u0::Vector{Float64}, rt_cfg::ResidenceTimeConfig)
    return _prepare_residence_time_data(layout, to_config(config), u0, rt_cfg)
end

@inline function _resolve_residence_time(
        residence_time::Union{Nothing, ResidenceTimeConfig},
        use_residence_time::Union{Nothing, Bool})
    default_enabled = residence_time !== nothing && residence_time.enabled
    enabled = use_residence_time === nothing ? default_enabled : use_residence_time
    if !enabled
        return nothing
    end
    return residence_time === nothing ? ResidenceTimeConfig() : residence_time
end

mutable struct NativeRampLimiter{T <: Real}
    base_dt::T
    understep_dt::T
    history_steps::Int
end

"""
$(SIGNATURES)

Apply the native TERRA ramp limiter update to the integrator.

# Arguments
- `integrator`: DifferentialEquations.jl integrator being stepped

# Returns
- `Nothing`
"""
function (lim::NativeRampLimiter)(integrator)
    if lim.history_steps <= 0
        return u_modified!(integrator, false)
    end

    steps_done = integrator.stats.naccept

    if steps_done < lim.history_steps
        set_proposed_dt!(integrator, lim.understep_dt)
    elseif steps_done == lim.history_steps && integrator.dt < lim.base_dt
        set_proposed_dt!(integrator, lim.base_dt)
    end

    u_modified!(integrator, false)
end

"""
$(SIGNATURES)

Condition function for the native ramp callback (always triggers).

# Arguments
- `u`: Current state vector (unused)
- `t`: Current simulation time (unused)
- `integrator`: DifferentialEquations.jl integrator (unused)

# Returns
- `true`
"""
native_ramp_condition(u, t, integrator) = true

"""
$(SIGNATURES)

Initializer for the native ramp callback.

# Arguments
- `cb`: Callback instance containing the ramp limiter
- `u`: Current state vector (unused)
- `t`: Current simulation time (unused)
- `integrator`: DifferentialEquations.jl integrator instance

# Returns
- `Nothing`
"""
function native_ramp_initialize(cb, u, t, integrator)
    if cb.affect!.history_steps > 0
        set_proposed_dt!(integrator, cb.affect!.understep_dt)
    end
    u_modified!(integrator, false)
end

"""
$(SIGNATURES)

Create the native TERRA step-size ramp callback.

# Arguments
- `initial_dt`: Baseline time step used after the ramp is complete
- `understep_ratio`: Optional ratio controlling the initial under-stepping
- `history_steps`: Optional number of accepted steps to maintain the ramp

# Returns
- `DiscreteCallback`: Callback implementing the ramp behaviour
"""
function native_ramp_callback(initial_dt; understep_ratio = inv(128), history_steps = 5)
    understep_ratio <= 0 && error("understep_ratio must be positive")
    understep_dt = min(initial_dt, initial_dt * understep_ratio)
    affect! = NativeRampLimiter(initial_dt, understep_dt, history_steps)
    DiscreteCallback(native_ramp_condition, affect!;
        initialize = native_ramp_initialize,
        save_positions = (false, false))
end

@inline function _is_electron_species(name::AbstractString, molecular_weight::Real)
    species = uppercase(strip(name))
    if species in ("E", "E-", "E+", "ELEC", "ELECTRON")
        return true
    end
    return molecular_weight < 1.0e-2
end

"""
$(SIGNATURES)

Compute the mixture pressure using translational and electron temperatures.

# Arguments
- `rho_sp::AbstractVector{<:Real}`: Species mass densities in CGS units
- `gas_constants::AbstractVector{<:Real}`: Species-specific gas constants
- `species_names::AbstractVector{<:AbstractString}`: Species identifiers used to detect electrons
- `molecular_weights::AbstractVector{<:Real}`: Species molecular weights
- `tt::Real`: Translational (heavy-particle) temperature
- `te::Real`: Electron temperature applied to electron species

# Returns
- `Float64`: Mixture pressure consistent with TERRA conventions
"""
