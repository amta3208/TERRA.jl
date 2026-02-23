"""
$(SIGNATURES)

Pack state components into the `y` vector ordering used by Fortran `rhs_api`.
"""
function pack_state_vector(layout::ApiLayout,
        rho_sp::AbstractVector{<:Real},
        rho_etot::Real;
        rho_ex::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        rho_vx::Union{AbstractArray{<:Real, 3}, Nothing} = nothing,
        rho_u::Union{Real, Nothing} = nothing,
        rho_v::Union{Real, Nothing} = nothing,
        rho_w::Union{Real, Nothing} = nothing,
        rho_erot::Union{Real, Nothing} = nothing,
        rho_eeex::Union{Real, Nothing} = nothing,
        rho_evib::Union{Real, Nothing} = nothing)
    y = Vector{Float64}(undef, layout.neq)
    pack_state_vector!(y, layout, rho_sp, rho_etot;
        rho_ex = rho_ex,
        rho_vx = rho_vx,
        rho_u = rho_u,
        rho_v = rho_v,
        rho_w = rho_w,
        rho_erot = rho_erot,
        rho_eeex = rho_eeex,
        rho_evib = rho_evib)
    return y
end

"""
$(SIGNATURES)

In-place variant of `pack_state_vector`.
"""
function pack_state_vector!(y::Vector{Float64},
        layout::ApiLayout,
        rho_sp::AbstractVector{<:Real},
        rho_etot::Real;
        rho_ex::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        rho_vx::Union{AbstractArray{<:Real, 3}, Nothing} = nothing,
        rho_u::Union{Real, Nothing} = nothing,
        rho_v::Union{Real, Nothing} = nothing,
        rho_w::Union{Real, Nothing} = nothing,
        rho_erot::Union{Real, Nothing} = nothing,
        rho_eeex::Union{Real, Nothing} = nothing,
        rho_evib::Union{Real, Nothing} = nothing)
    if length(y) != layout.neq
        throw(DimensionMismatch("y length ($(length(y))) does not match layout.neq ($(layout.neq))"))
    end
    if length(rho_sp) != layout.nsp
        throw(DimensionMismatch("rho_sp length ($(length(rho_sp))) must match layout.nsp ($(layout.nsp))"))
    end
    if layout.is_elec_sts && rho_ex === nothing
        throw(ArgumentError("rho_ex must be provided when electronic STS is active (layout.is_elec_sts=true)."))
    end
    if layout.is_vib_sts && rho_vx === nothing
        throw(ArgumentError("rho_vx must be provided when vibrational STS is active (layout.is_vib_sts=true)."))
    end
    if layout.rot_noneq && rho_erot === nothing
        throw(ArgumentError("rho_erot must be provided when rot_noneq=true."))
    end
    if layout.eex_noneq && rho_eeex === nothing
        throw(ArgumentError("rho_eeex must be provided when eex_noneq=true."))
    end
    if layout.vib_noneq && rho_evib === nothing
        throw(ArgumentError("rho_evib must be provided when vib_noneq=true."))
    end
    if layout.nd >= 1 && rho_u === nothing
        throw(ArgumentError("rho_u must be provided when nd>=1."))
    end
    if layout.nd >= 2 && rho_v === nothing
        throw(ArgumentError("rho_v must be provided when nd>=2."))
    end
    if layout.nd >= 3 && rho_w === nothing
        throw(ArgumentError("rho_w must be provided when nd>=3."))
    end

    fill!(y, 0.0)
    idx = 1

    # Vibrational states
    if layout.n_eq_vib > 0
        @assert rho_vx !== nothing
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0 || layout.ih[isp] != 2
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 0
                    continue
                end
                mv_isp_iex = layout.mv[isp, iex]
                for ivx in 0:mv_isp_iex
                    y[idx] = Float64((rho_vx::AbstractArray{<:Real, 3})[ivx + 1, iex, isp])
                    idx += 1
                end
            end
        end
    end

    # Electronic states
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
                y[idx] = Float64((rho_ex::AbstractMatrix{<:Real})[iex, isp])
                idx += 1
            end
        end
    end

    # Species densities for non-electronic-specific species
    @inbounds for isp in 1:(layout.nsp)
        if layout.ies[isp] == 1
            continue
        end
        if layout.get_electron_density_by_charge_balance && isp == layout.esp
            continue
        end
        y[idx] = Float64(rho_sp[isp])
        idx += 1
    end

    # Momentum (if any)
    if layout.nd >= 1
        y[idx] = Float64(rho_u::Real)
        idx += 1
    end
    if layout.nd >= 2
        y[idx] = Float64(rho_v::Real)
        idx += 1
    end
    if layout.nd >= 3
        y[idx] = Float64(rho_w::Real)
        idx += 1
    end

    # Energies
    y[idx] = Float64(rho_etot)
    idx += 1
    if layout.eex_noneq
        y[idx] = Float64(rho_eeex::Real)
        idx += 1
    end
    if layout.rot_noneq
        y[idx] = Float64(rho_erot::Real)
        idx += 1
    end
    if layout.vib_noneq
        y[idx] = Float64(rho_evib::Real)
        idx += 1
    end

    @assert idx==layout.neq + 1 "Internal error: state packing length mismatch"
    return nothing
end

"""
$(SIGNATURES)

Unpack a compact `y` vector into component arrays.

This returns `rho_sp` for all active species (including charge-balanced electrons
when enabled), plus `rho_ex`/`rho_vx` when the corresponding STS modes are active.
"""
function unpack_state_vector(y::AbstractVector{<:Real}, layout::ApiLayout)
    if length(y) != layout.neq
        throw(DimensionMismatch("y length ($(length(y))) does not match layout.neq ($(layout.neq))"))
    end

    rho_sp = zeros(Float64, layout.nsp)
    rho_ex = layout.is_elec_sts ? zeros(Float64, layout.mnex, layout.nsp) : nothing
    rho_vx = layout.n_eq_vib > 0 ?
             zeros(Float64, layout.mnv + 1, layout.mmnex, layout.nsp) : nothing

    idx = 1
    rho_total = 0.0

    # Vibrational states first
    if layout.n_eq_vib > 0
        @assert rho_vx !== nothing
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0 || layout.ih[isp] != 2
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 0
                    continue
                end
                mv_isp_iex = layout.mv[isp, iex]
                for ivx in 0:mv_isp_iex
                    val = Float64(y[idx])
                    (rho_vx::Array{Float64, 3})[ivx + 1, iex, isp] = val
                    rho_sp[isp] += val
                    rho_total += val
                    idx += 1
                end
            end
        end
    end

    # Electronic states next
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
                val = Float64(y[idx])
                (rho_ex::Matrix{Float64})[iex, isp] = val
                rho_sp[isp] += val
                rho_total += val
                idx += 1
            end
        end
    end

    # Species densities
    @inbounds for isp in 1:(layout.nsp)
        if layout.ies[isp] == 1
            continue
        end
        if layout.get_electron_density_by_charge_balance && isp == layout.esp
            continue
        end
        val = Float64(y[idx])
        rho_sp[isp] = val
        rho_total += val
        idx += 1
    end

    # Momentum
    rho_u = 0.0
    rho_v = 0.0
    rho_w = 0.0
    if layout.nd >= 1
        rho_u = Float64(y[idx])
        idx += 1
    end
    if layout.nd >= 2
        rho_v = Float64(y[idx])
        idx += 1
    end
    if layout.nd >= 3
        rho_w = Float64(y[idx])
        idx += 1
    end

    # Energies
    rho_etot = Float64(y[idx])
    idx += 1
    rho_eeex = 0.0
    rho_erot = 0.0
    rho_evib = 0.0
    if layout.eex_noneq
        rho_eeex = Float64(y[idx])
        idx += 1
    end
    if layout.rot_noneq
        rho_erot = Float64(y[idx])
        idx += 1
    end
    if layout.vib_noneq
        rho_evib = Float64(y[idx])
        idx += 1
    end

    @assert idx==layout.neq + 1 "Internal error: state unpacking length mismatch"

    # Reconstruct electron density from charge balance if requested (matches Fortran convention).
    if layout.get_electron_density_by_charge_balance && layout.esp >= 1 &&
       layout.esp <= layout.nsp && rho_total > 0
        spgam = zeros(Float64, layout.nsp)
        @inbounds for isp in 1:(layout.nsp)
            if rho_sp[isp] == 0.0
                continue
            end
            spgam[isp] = rho_sp[isp] / rho_total / (AVOGADRO * layout.spwt[isp])
        end
        spgam_e = 0.0
        @inbounds for isp in 1:(layout.nsp)
            if isp == layout.esp
                continue
            end
            spgam_e += layout.ie[isp] * spgam[isp]
        end
        rho_sp[layout.esp] = spgam_e * rho_total * AVOGADRO * layout.spwt[layout.esp]
    end

    return (
        rho_sp = rho_sp,
        rho_etot = rho_etot,
        rho_ex = rho_ex,
        rho_vx = rho_vx,
        rho_u = rho_u,
        rho_v = rho_v,
        rho_w = rho_w,
        rho_erot = rho_erot,
        rho_eeex = rho_eeex,
        rho_evib = rho_evib
    )
end
