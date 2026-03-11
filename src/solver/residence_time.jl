function _build_residence_time_continuity_indices(layout::ApiLayout,
        species_names::AbstractVector{<:AbstractString})
    length(species_names) == layout.nsp || throw(ArgumentError(
        "Residence-time species ordering length ($(length(species_names))) does not match layout.nsp ($(layout.nsp))."
    ))

    species_order = String[String(species_names[isp]) for isp in 1:layout.nsp if isp != layout.esp]
    species_indices = Dict(name => Int[] for name in species_order)

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
                        push!(species_indices[String(species_names[isp])], idx)
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
                    push!(species_indices[String(species_names[isp])], idx)
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
        if isp != layout.esp
            push!(species_indices[String(species_names[isp])], idx)
        end
        idx += 1
    end

    @assert idx==layout.mom_range.start "Internal error: continuity index mismatch while building residence-time indices"

    return (species_order = species_order, species_indices = species_indices)
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
    @inbounds for name in rt.species_order
        inv_tau = rt.inv_tau_species[name]
        for i in rt.species_indices[name]
            du[i] += (u_in[i] - u[i]) * inv_tau
        end
    end
    @inbounds for i in rt.energy_indices
        du[i] += (u_in[i] - u[i]) * rt.inv_tau_energy
    end
    return nothing
end

function _validate_residence_time_species(rt_cfg::ResidenceTimeConfig,
        species_order::Vector{String})
    missing_species = String[name for name in species_order if !haskey(rt_cfg.U_species, name)]
    extra_species = sort!(String[
        name for name in keys(rt_cfg.U_species) if !(name in species_order)
    ])

    if !isempty(missing_species) || !isempty(extra_species)
        throw(ArgumentError(
            "ResidenceTimeConfig U_species keys must match the non-electron TERRA species ordering. " *
            "Missing: $(isempty(missing_species) ? "none" : join(missing_species, ", ")); " *
            "Extra: $(isempty(extra_species) ? "none" : join(extra_species, ", "))."
        ))
    end
end

function _prepare_residence_time_data(layout::ApiLayout, config::Config,
        u0::Vector{Float64}, rt_cfg::ResidenceTimeConfig)
    continuity = _build_residence_time_continuity_indices(
        layout, config.reactor.composition.species)
    species_order = continuity.species_order
    species_indices = continuity.species_indices
    _validate_residence_time_species(rt_cfg, species_order)

    inv_tau_species = Dict(
        name => rt_cfg.U_species[name] / rt_cfg.L for name in species_order
    )
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

    inlet_state = unpack_state_vector(u_in, layout)

    U_energy = if rt_cfg.U_energy === nothing
        rho_weighted_u = 0.0
        rho_total = 0.0
        @inbounds for isp in 1:layout.nsp
            if isp == layout.esp
                continue
            end
            species_name = String(config.reactor.composition.species[isp])
            rho_species = inlet_state.rho_sp[isp]
            rho_weighted_u += rho_species * rt_cfg.U_species[species_name]
            rho_total += rho_species
        end
        if !(rho_total > 0.0) || !isfinite(rho_total)
            error("ResidenceTimeConfig: cannot infer U_energy from inlet state (rho_total=$rho_total).")
        end
        rho_weighted_u / rho_total
    else
        rt_cfg.U_energy
    end

    inv_tau_energy = U_energy / rt_cfg.L

    return (
        u_in = u_in,
        species_order = species_order,
        species_indices = species_indices,
        energy_indices = energy_indices,
        inv_tau_species = inv_tau_species,
        inv_tau_energy = inv_tau_energy
    )
end

@inline function _resolve_residence_time(
        residence_time::Union{Nothing, ResidenceTimeConfig},
        use_residence_time::Union{Nothing, Bool})
    default_enabled = residence_time !== nothing && residence_time.enabled
    enabled = use_residence_time === nothing ? default_enabled : use_residence_time
    if !enabled
        return nothing
    end
    residence_time === nothing && throw(ArgumentError(
        "Residence-time forcing was requested, but no ResidenceTimeConfig was provided."
    ))
    return residence_time
end
