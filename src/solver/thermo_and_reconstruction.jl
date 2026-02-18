function compute_mixture_pressure(rho_sp::AbstractVector{<:Real},
        gas_constants::AbstractVector{<:Real},
        species_names::AbstractVector{<:AbstractString},
        molecular_weights::AbstractVector{<:Real},
        tt::Real, te::Real)
    @assert length(rho_sp) == length(gas_constants) == length(species_names) ==
            length(molecular_weights)

    pressure = 0.0
    @inbounds for i in eachindex(rho_sp, gas_constants, species_names, molecular_weights)
        rho_i = rho_sp[i]
        rho_i <= 0 && continue
        T_i = _is_electron_species(species_names[i], molecular_weights[i]) ? te : tt
        pressure += rho_i * gas_constants[i] * T_i
    end
    return pressure
end

"""
$(SIGNATURES)

Convert total energy density into enthalpy density using the ideal-gas closure.

# Arguments
- `rho_etot::Float64`: Total energy density (CGS units)
- `rho_sp::AbstractVector{<:Real}`: Species mass densities
- `gas_constants::AbstractVector{<:Real}`: Species-specific gas constants
- `species_names::AbstractVector{<:AbstractString}`: Species names used to detect electrons
- `molecular_weights::AbstractVector{<:Real}`: Species molecular weights
- `tt::Real`: Translational temperature applied to heavy species
- `te::Real`: Electron temperature applied to electron species

# Returns
- `Tuple{Float64, Float64}`: Enthalpy density and corresponding pressure
"""
function enthalpy_from_energy(rho_etot::Float64, rho_sp::AbstractVector{<:Real},
        gas_constants::AbstractVector{<:Real}, species_names::AbstractVector{<:AbstractString},
        molecular_weights::AbstractVector{<:Real},
        tt::Real, te::Real)
    pressure = compute_mixture_pressure(
        rho_sp, gas_constants, species_names, molecular_weights, tt, te)
    return rho_etot + pressure, pressure
end

@inline function _compact_nonnegative_ok(u::AbstractVector{<:Real}, layout::ApiLayout)
    if any(@view(u[layout.vib_range]) .< 0.0)
        return false
    end
    if any(@view(u[layout.elec_range]) .< 0.0)
        return false
    end
    if any(@view(u[layout.sp_range]) .< 0.0)
        return false
    end
    if u[layout.idx_etot] < 0.0
        return false
    end
    if layout.idx_eeex != 0 && u[layout.idx_eeex] < 0.0
        return false
    end
    if layout.idx_erot != 0 && u[layout.idx_erot] < 0.0
        return false
    end
    if layout.idx_evib != 0 && u[layout.idx_evib] < 0.0
        return false
    end
    return true
end

@inline function _compact_density_nonnegative_ok(u::AbstractVector{<:Real}, layout::ApiLayout)
    if any(@view(u[layout.vib_range]) .< 0.0)
        return false
    end
    if any(@view(u[layout.elec_range]) .< 0.0)
        return false
    end
    if any(@view(u[layout.sp_range]) .< 0.0)
        return false
    end
    return true
end

@inline function _compact_energy_nonnegative_ok(u::AbstractVector{<:Real}, layout::ApiLayout)
    if u[layout.idx_etot] < 0.0
        return false
    end
    if layout.idx_eeex != 0 && u[layout.idx_eeex] < 0.0
        return false
    end
    if layout.idx_erot != 0 && u[layout.idx_erot] < 0.0
        return false
    end
    if layout.idx_evib != 0 && u[layout.idx_evib] < 0.0
        return false
    end
    return true
end

@inline function _compact_project_density_nonnegative!(u_work::Vector{Float64},
        u::AbstractVector{<:Real},
        layout::ApiLayout)
    copyto!(u_work, u)

    density_stop = layout.mom_range.start - 1
    if density_stop <= 0
        return true
    end

    # Preserve total density for the continuity block. This avoids creating an
    # unphysical (rho, etot) pair by simply clamping negatives to zero, which can
    # drive temperatures negative/undefined inside the Fortran RHS.
    rho_target = 0.0
    rho_clamped = 0.0
    @inbounds for i in 1:density_stop
        ui_raw = Float64(u_work[i])
        rho_target += Float64(u[i])
        if ui_raw < 0.0
            ui_raw = 0.0
            u_work[i] = 0.0
        end
        rho_clamped += ui_raw
    end

    if !(rho_target > 0.0) || !isfinite(rho_target)
        return false
    end
    if !(rho_clamped > 0.0) || !isfinite(rho_clamped)
        return false
    end

    scale = rho_target / rho_clamped
    if !(isfinite(scale)) || scale <= 0.0
        return false
    end

    if scale != 1.0
        @inbounds for i in 1:density_stop
            u_work[i] *= scale
        end
    end

    return true
end

@inline function _compact_apply_shampine_positivity!(du::AbstractVector{<:Real},
        u::AbstractVector{<:Real},
        layout::ApiLayout)
    density_stop = layout.mom_range.start - 1
    if density_stop <= 0
        return nothing
    end
    @inbounds for i in 1:density_stop
        if u[i] < 0.0
            du[i] = max(Float64(du[i]), 0.0)
        end
    end
    return nothing
end

function _reconstruct_rho_sp_rho_ex_from_compact!(rho_sp::Vector{Float64},
        rho_ex::Union{Nothing, Matrix{Float64}},
        u::Vector{Float64},
        layout::ApiLayout)
    fill!(rho_sp, 0.0)
    if rho_ex !== nothing
        fill!(rho_ex, 0.0)
    end

    idx = 1

    # Vibrational STS densities (contribute to species totals)
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
                    val = u[idx]
                    rho_sp[isp] += val
                    idx += 1
                end
            end
        end
    end

    # Electronic STS densities (contribute to species totals and rho_ex buffer)
    if layout.is_elec_sts
        @assert rho_ex !== nothing
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 1
                    continue
                end
                val = u[idx]
                (rho_ex::Matrix{Float64})[iex, isp] = val
                rho_sp[isp] += val
                idx += 1
            end
        end
    end

    # Species densities for non-electronic-specific species (except charge-balanced electrons)
    @inbounds for isp in 1:(layout.nsp)
        if layout.ies[isp] == 1
            continue
        end
        if layout.get_electron_density_by_charge_balance && isp == layout.esp
            continue
        end
        val = u[idx]
        rho_sp[isp] = val
        idx += 1
    end

    @assert idx==layout.mom_range.start "Internal error: continuity reconstruction index mismatch"

    # Reconstruct electron mass density from charge balance if requested.
    if layout.get_electron_density_by_charge_balance
        esp = layout.esp
        if esp >= 1 && esp <= layout.nsp
            s = 0.0
            @inbounds for isp in 1:(layout.nsp)
                if isp == esp
                    continue
                end
                if layout.ie[isp] != 0
                    s += layout.ie[isp] * rho_sp[isp] / layout.spwt[isp]
                end
            end
            rho_sp[esp] = layout.spwt[esp] * s
        end
    end

    return nothing
end

function _reconstruct_drho_sp_from_compact!(drho_sp::Vector{Float64},
        dy::AbstractVector{<:Real},
        layout::ApiLayout)
    fill!(drho_sp, 0.0)

    idx = 1

    # Vibrational STS derivatives (contribute to species totals)
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
                    drho_sp[isp] += Float64(dy[idx])
                    idx += 1
                end
            end
        end
    end

    # Electronic STS derivatives (contribute to species totals)
    if layout.is_elec_sts
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 1
                    continue
                end
                drho_sp[isp] += Float64(dy[idx])
                idx += 1
            end
        end
    end

    # Species-density derivatives for non-electronic-specific species (except charge-balanced electrons)
    @inbounds for isp in 1:(layout.nsp)
        if layout.ies[isp] == 1
            continue
        end
        if layout.get_electron_density_by_charge_balance && isp == layout.esp
            continue
        end
        drho_sp[isp] = Float64(dy[idx])
        idx += 1
    end

    @assert idx==layout.mom_range.start "Internal error: continuity derivative reconstruction index mismatch"

    return nothing
end

@inline function _compact_drho_e_from_drho_sp(drho_sp::Vector{Float64}, layout::ApiLayout)
    esp = layout.esp
    if !(esp >= 1 && esp <= layout.nsp)
        return 0.0
    end
    if !layout.get_electron_density_by_charge_balance
        return drho_sp[esp]
    end

    s = 0.0
    @inbounds for isp in 1:(layout.nsp)
        if isp == esp
            continue
        end
        zi = layout.ie[isp]
        if zi != 0
            s += zi * drho_sp[isp] / layout.spwt[isp]
        end
    end
    return layout.spwt[esp] * s
end

@inline _config_electron_temperature(config::TERRAConfig) = config.temperatures.Te
@inline _config_electron_temperature(config::Config) = config.reactor.thermal.Te
@inline _config_is_isothermal_teex(config::TERRAConfig) = config.physics.is_isothermal_teex
@inline _config_is_isothermal_teex(config::Config) = config.models.physics.is_isothermal_teex

function _compact_isothermal_fill_fortran_y_work!(y_work::Vector{Float64},
        rho_sp::Vector{Float64},
        rho_ex::Union{Nothing, Matrix{Float64}},
        u::Vector{Float64},
        p,
        layout::ApiLayout)
    config = p.config
    teex_vec = hasproperty(p, :teex_const_vec) ? p.teex_const_vec : nothing
    teex_const = hasproperty(p, :teex_const) ? p.teex_const :
                 _config_electron_temperature(config)

    _reconstruct_rho_sp_rho_ex_from_compact!(rho_sp, rho_ex, u, layout)

    rho_evib = u[layout.idx_evib]
    tvib = calculate_vibrational_temperature_wrapper(rho_evib, rho_sp;
        rho_ex = rho_ex,
        tex = teex_vec)
    rho_eeex_eff = calculate_electron_electronic_energy_wrapper(teex_const, tvib, rho_sp)

    electron_idx = layout.esp
    gas_const_e = (electron_idx >= 1 && electron_idx <= layout.nsp) ?
                  Float64(p.gas_constants[electron_idx]) : 0.0
    electron_enthalpy = (electron_idx >= 1 && electron_idx <= layout.nsp) ?
                        rho_sp[electron_idx] * gas_const_e * teex_const : 0.0

    rho_rem = u[layout.idx_etot]
    rho_enthalpy_total = rho_rem + rho_eeex_eff + electron_enthalpy
    rho_erot = layout.idx_erot == 0 ? 0.0 : Float64(u[layout.idx_erot])
    conversion = energy_from_enthalpy_isothermal_teex_wrapper(rho_enthalpy_total, rho_sp, teex_const;
        rho_ex = rho_ex,
        rho_eeex = rho_eeex_eff,
        rho_erot = rho_erot,
        rho_evib = rho_evib)
    rho_etot = conversion.rho_etot

    copyto!(y_work, u)
    y_work[layout.idx_eeex] = rho_eeex_eff
    y_work[layout.idx_etot] = rho_etot

    temps_raw = calculate_temperatures_wrapper(rho_sp, rho_etot;
        rho_ex = rho_ex,
        rho_erot = layout.idx_erot == 0 ? nothing : rho_erot,
        rho_eeex = rho_eeex_eff,
        rho_evib = rho_evib)
    temps = (; temps_raw..., teex = teex_const)

    return (
        teex_const = teex_const,
        rho_eeex = rho_eeex_eff,
        rho_evib = rho_evib,
        temps = temps,
        rho_etot = rho_etot,
        pressure = conversion.pressure
    )
end
