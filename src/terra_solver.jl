"""
# TERRA Solver Module

This module provides the high-level interface for running TERRA simulations
from Julia, hiding the complexity of the Fortran interface and providing
a clean, Julia-native API.
"""

# Debug: RHS call counter
const _TERRA_ODE_DBG_CALLS = Ref(0)

"""
Compact API layout information for the Fortran `rhs_api` / `calculate_rhs_api` state vectors.

This describes the *compact* ordering used internally by TERRA:
`[rho_vib_states, rho_elec_states, rho_species, rho_u..., rho_etot, rho_eeex?, rho_erot?, rho_evib?]`.
"""
struct ApiLayout
    layout_version::Int
    mnsp::Int
    mnex::Int
    mmnex::Int
    mnv::Int
    mneq::Int

    nsp::Int
    nd::Int
    neq::Int
    esp::Int

    get_electron_density_by_charge_balance::Bool
    eex_noneq::Bool
    rot_noneq::Bool
    vib_noneq::Bool
    is_isothermal::Bool
    is_isothermal_teex::Bool
    is_elec_sts::Bool
    is_vib_sts::Bool

    n_eq_vib::Int
    n_eq_elec::Int
    n_eq_sp::Int
    n_eq_mom::Int
    n_eq_energy::Int

    ih::Vector{Int}
    ie::Vector{Int}
    ies::Vector{Int}
    mex::Vector{Int}
    ivs::Matrix{Int}
    mv::Matrix{Int}
    spwt::Vector{Float64}

    vib_range::UnitRange{Int}
    elec_range::UnitRange{Int}
    sp_range::UnitRange{Int}
    mom_range::UnitRange{Int}
    energy_range::UnitRange{Int}

    idx_etot::Int
    idx_eeex::Int
    idx_erot::Int
    idx_evib::Int
end

function ApiLayout(layout::NamedTuple)
    ih = Int.(layout.ih)
    ie = Int.(layout.ie)
    ies = Int.(layout.ies)
    mex = Int.(layout.mex)
    ivs = Int.(layout.ivs)
    mv = Int.(layout.mv)
    spwt = Vector{Float64}(layout.spwt)

    n_eq_vib = Int(layout.n_eq_vib)
    n_eq_elec = Int(layout.n_eq_elec)
    n_eq_sp = Int(layout.n_eq_sp)
    n_eq_mom = Int(layout.n_eq_mom)
    n_eq_energy = Int(layout.n_eq_energy)
    neq = Int(layout.neq)

    vib_start = 1
    vib_stop = n_eq_vib
    elec_start = vib_stop + 1
    elec_stop = vib_stop + n_eq_elec
    sp_start = elec_stop + 1
    sp_stop = elec_stop + n_eq_sp
    mom_start = sp_stop + 1
    mom_stop = sp_stop + n_eq_mom
    energy_start = mom_stop + 1
    energy_stop = neq

    idx_etot = energy_start
    idx_eeex = (layout.eex_noneq == 1) ? (idx_etot + 1) : 0
    idx_erot = (layout.rot_noneq == 1) ? (idx_etot + 1 + (layout.eex_noneq == 1 ? 1 : 0)) :
               0
    idx_evib = (layout.vib_noneq == 1) ?
               (idx_etot + 1 + (layout.eex_noneq == 1 ? 1 : 0) +
                (layout.rot_noneq == 1 ? 1 : 0)) : 0

    return ApiLayout(
        Int(layout.layout_version),
        Int(layout.mnsp),
        Int(layout.mnex),
        Int(layout.mmnex),
        Int(layout.mnv),
        Int(layout.mneq),
        Int(layout.nsp),
        Int(layout.nd),
        neq,
        Int(layout.esp),
        layout.get_electron_density_by_charge_balance == 1,
        layout.eex_noneq == 1,
        layout.rot_noneq == 1,
        layout.vib_noneq == 1,
        layout.is_isothermal == 1,
        layout.is_isothermal_teex == 1,
        layout.is_elec_sts == 1,
        layout.is_vib_sts == 1,
        n_eq_vib,
        n_eq_elec,
        n_eq_sp,
        n_eq_mom,
        n_eq_energy,
        ih,
        ie,
        ies,
        mex,
        ivs,
        mv,
        spwt,
        vib_start:vib_stop,
        elec_start:elec_stop,
        sp_start:sp_stop,
        mom_start:mom_stop,
        energy_start:energy_stop,
        idx_etot,
        idx_eeex,
        idx_erot,
        idx_evib
    )
end

"""
$(SIGNATURES)

Query the Fortran API for the current `y`/`dy` layout.
"""
function get_api_layout()
    return ApiLayout(get_api_layout_wrapper())
end

"""
$(SIGNATURES)

Initialize the TERRA system.

This function must be called before any TERRA calculations can be performed.
It sets up the Fortran library, initializes internal data structures, and
prepares the system for simulation. The TERRA shared library is obtained from
the `TERRA_LIB_PATH` environment variable when it is not already loaded.

# Arguments
- `config::TERRAConfig`: Configuration for initialization
- `case_path::String`: Case directory path (optional, defaults to config.case_path)

# Returns
- `true` if initialization successful

# Throws
- `ErrorException` if initialization fails
"""
function initialize_terra(config::TERRAConfig, case_path::String = config.case_path)
    try
        # Ensure the shared library is loaded
        if !is_terra_loaded()
            load_terra_library!()
            @debug "TERRA library loaded successfully via TERRA_LIB_PATH"
        end

        # If the Fortran API is already initialized, finalize it so we can
        # reinitialize using the updated configuration and input files.
        try
            if is_api_initialized_wrapper()
                @debug "Finalizing existing TERRA API before re-initialization"
                finalize_api_wrapper()
            end
        catch e
            @warn "Unable to query/finalize existing TERRA API state" exception=e
        end

        # Generate input files for TERRA based on the provided configuration
        generate_input_files(config, case_path)
        @debug "TERRA input files generated" case_path=case_path

        # Initialize the API - get dimensions from Fortran
        result = initialize_api_wrapper(case_path = case_path)
        num_species = result.num_species
        num_dimensions = result.num_dimensions

        # Hard consistency check: configured species count must match TERRA setup
        if num_species != length(config.species)
            error("Configured species count ($(length(config.species))) does not match TERRA setup ($num_species). " *
                  "Ensure configuration and generated input match the TERRA database.")
        end

        @debug "TERRA initialized successfully" num_species=num_species num_dimensions=num_dimensions
        # Fetch and log runtime flags from TERRA for verification
        try
            flags = get_runtime_flags()
            @debug "TERRA runtime flags" ev_relax_set=flags.ev_relax_set vib_noneq=flags.vib_noneq eex_noneq=flags.eex_noneq rot_noneq=flags.rot_noneq bfe=flags.consider_elec_bfe bbh=flags.consider_elec_bbh bfh=flags.consider_elec_bfh bbe=flags.consider_elec_bbe
        catch e
            @warn "Unable to read TERRA runtime flags (rebuild library to enable)" exception=e
        end
        return true

    catch e
        @error "Failed to initialize TERRA" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Initialize TERRA using nested `Config`.
"""
function initialize_terra(config::Config, case_path::String = config.runtime.case_path)
    return initialize_terra(to_legacy_config(config), case_path)
end

"""
$(SIGNATURES)

Finalize the TERRA system and clean up resources.

This function should be called when TERRA is no longer needed to properly
clean up memory and resources.

# Arguments
- None

# Returns
- `Nothing`
"""
function finalize_terra()
    try
        # Call the wrapper finalization
        finalize_api_wrapper()

        # Close the library
        close_terra_library()

        @info "TERRA finalized successfully"
    catch e
        @error "Error during TERRA finalization" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Check if TERRA is properly initialized.

# Arguments
- None

# Returns
- `true` if TERRA is initialized and ready for use
"""
function is_terra_initialized()
    try
        return is_terra_loaded() && is_api_initialized_wrapper()
    catch
        return false
    end
end

"""
$(SIGNATURES)

Convert TERRAConfig to initial state vectors for TERRA.

# Arguments
- `config::TERRAConfig`: Configuration object

# Returns
- Named tuple with initial state vectors in CGS units (as required by TERRA)
"""
function config_to_initial_state(config::TERRAConfig)
    # Get molecular weights
    molecular_weights = get_molecular_weights(config.species)

    # Ensure we're working in CGS - convert if needed
    config_cgs = config.unit_system == :CGS ? config : convert_config_units(config, :CGS)

    # Convert mole fractions to mass densities (CGS units)
    mass_densities_cgs = mole_fractions_to_mass_densities(
        config_cgs.mole_fractions, molecular_weights, config_cgs.total_number_density
    )

    gas_constants_full = get_species_gas_constants_wrapper()
    gas_constants = gas_constants_full[1:length(molecular_weights)]

    # Calculate electronic state populations using TERRA's Boltzmann distribution
    # Use appropriate temperatures for each mode
    initial_electronic_states = set_electronic_boltzmann_wrapper(
        mass_densities_cgs,
        # Initialize electronic-state-resolved populations using the
        # electron-electronic temperature (TEE). Species that are not
        # electronically resolved do not contribute here (their mex = 0),
        # so using TEE only affects resolved species as intended.
        config_cgs.temperatures.Tee,  # Electronic-state populations use TEE
        config_cgs.temperatures.Tt,  # Rotational temperature (use translational as proxy)
        config_cgs.temperatures.Tv   # Vibrational temperature
    )

    # Calculate electron-electronic energy using TERRA's method.
    # IMPORTANT: Species without electronically resolved states contribute to
    # the electron-electronic mode; their energy must be initialized using the
    # electron temperature (TE), not TEE. TERRA internally excludes species with
    # resolved electronic states from this mode, so passing TE here is correct.
    initial_electron_electronic_energy = calculate_electron_electronic_energy_wrapper(
        config_cgs.temperatures.Te, config_cgs.temperatures.Tv, mass_densities_cgs
    )

    # Initialize vibrational state populations using TERRA's Boltzmann
    # distribution (for potential diagnostics), but do NOT use this to infer
    # STS. The active STS setting comes from the database/setup, which the
    # wrapper cannot reliably query here. To avoid a mismatch with TERRA,
    # always treat vibrational energy in mode form and do not include rho_vx
    # in Etot.
    initial_vibrational_states = set_vibrational_boltzmann_wrapper(
        initial_electronic_states,
        config_cgs.temperatures.Te,
        config_cgs.temperatures.Tt,
        config_cgs.temperatures.Tv
    )

    # Mode-level vibrational energy at Tv. Provide per-species electronic
    # temperature proxy as TE uniformly; species with resolved electronic
    # states are already accounted through rho_ex above.
    initial_vibrational_energy = calculate_vibrational_energy_wrapper(
        config_cgs.temperatures.Tv, mass_densities_cgs;
        rho_ex = initial_electronic_states,
        tex = fill(config_cgs.temperatures.Te, length(mass_densities_cgs))
    )

    # Total energy: do not include rho_vx to avoid double counting when TERRA
    # is not in vibrational STS mode.
    initial_total_energy = calculate_total_energy_wrapper(
        config_cgs.temperatures.Tt, mass_densities_cgs;
        rho_ex = initial_electronic_states,
        u = 0.0, v = 0.0, w = 0.0,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy
    )

    has_elec_sts = has_electronic_sts_wrapper()
    has_vib_sts = has_vibrational_sts_wrapper()
    rho_ex_arg = has_elec_sts ? initial_electronic_states : nothing
    rho_vx_arg = has_vib_sts ? initial_vibrational_states : nothing

    initial_temperatures = calculate_temperatures_wrapper(
        mass_densities_cgs, initial_total_energy;
        rho_ex = rho_ex_arg,
        rho_vx = rho_vx_arg,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy)

    initial_enthalpy, initial_pressure = enthalpy_from_energy(
        initial_total_energy,
        mass_densities_cgs,
        gas_constants,
        config.species,
        molecular_weights,
        initial_temperatures.tt,
        initial_temperatures.teex)

    electron_index = findfirst(
        i -> _is_electron_species(config.species[i], molecular_weights[i]),
        eachindex(config.species))
    electron_enthalpy = electron_index === nothing ? 0.0 :
                        mass_densities_cgs[electron_index] * gas_constants[electron_index] *
                        config_cgs.temperatures.Te

    initial_remainder_energy = initial_enthalpy - initial_electron_electronic_energy -
                               electron_enthalpy

    return (
        rho_sp = mass_densities_cgs,
        rho_etot = initial_enthalpy,
        rho_energy = initial_total_energy,
        rho_rem = initial_remainder_energy,
        rho_ex = initial_electronic_states,
        rho_vx = nothing,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy,
        number_density = config_cgs.total_number_density,
        molecular_weights = molecular_weights,
        teex_const = config_cgs.temperatures.Te,
        gas_constants = gas_constants,
        pressure = initial_pressure,
        initial_temperatures = initial_temperatures
    )
end

"""
$(SIGNATURES)

Convert nested `Config` to initial state vectors for TERRA.
"""
function config_to_initial_state(config::Config)
    return config_to_initial_state(to_legacy_config(config))
end

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

function _prepare_residence_time_data(layout::ApiLayout, config::TERRAConfig,
        u0::Vector{Float64}, rt_cfg::ResidenceTimeConfig)
    inv_tau_neutral = rt_cfg.U_neutral / rt_cfg.L
    inv_tau_ion = rt_cfg.U_ion / rt_cfg.L

    continuity = _build_residence_time_continuity_indices(layout)
    neutral_indices = continuity.neutral_indices
    ion_indices = continuity.ion_indices
    energy_indices = _build_residence_time_energy_indices(
        layout, config.physics.is_isothermal_teex)

    u_in = if rt_cfg.inlet_config === nothing
        copy(u0)
    else
        inlet_config = rt_cfg.inlet_config
        if inlet_config.species != config.species
            error("ResidenceTimeConfig inlet_config species must match the initialized TERRA species ordering.")
        end
        if inlet_config.physics.is_isothermal_teex != config.physics.is_isothermal_teex
            error("ResidenceTimeConfig inlet_config isothermal_teex flag must match the simulation config.")
        end

        inlet_state = config_to_initial_state(inlet_config)
        energy_scalar_in = config.physics.is_isothermal_teex ? inlet_state.rho_rem :
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

function _compact_isothermal_fill_fortran_y_work!(y_work::Vector{Float64},
        rho_sp::Vector{Float64},
        rho_ex::Union{Nothing, Matrix{Float64}},
        u::Vector{Float64},
        p,
        layout::ApiLayout)
    config = p.config
    teex_vec = hasproperty(p, :teex_const_vec) ? p.teex_const_vec : nothing
    teex_const = hasproperty(p, :teex_const) ? p.teex_const : config.temperatures.Te

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

"""
$(SIGNATURES)

ODE system that calls the Fortran `rhs_api` directly (via `calculate_rhs_api`).

Assumptions about the state `u`:
- Non-isothermal (`config.physics.is_isothermal_teex == false`): `u` matches the Fortran `rhs_api` `y` layout and
  stores total energy density at `layout.idx_etot`.
- Isothermal Teex (`== true`): `u[layout.idx_etot]` stores the remainder energy
  `rho_rem = rho_h - rho_eeex - rho_e * R_e * Teex`. The `rho_eeex` slot is treated
  as a dummy (its derivative is forced to zero). The Fortran call uses a working vector with `rho_eeex` recomputed from
  the prescribed Teex, and with `rho_etot` reconstructed from `rho_rem` by converting total enthalpy back to energy.
"""
function terra_ode_system!(du::Vector{Float64}, u::Vector{Float64}, p, t::Float64)
    layout = p.layout

    if any(!isfinite, u)
        fill!(du, 0.0)
        return nothing
    end
    if !_compact_energy_nonnegative_ok(u, layout)
        fill!(du, 0.0)
        return nothing
    end

    # Enforce a nonnegative evaluation state for the Fortran RHS. This keeps the RHS
    # informative for stiffness detection (autoswitch) while protecting the Fortran
    # routines from invalid intermediate stage states.
    needs_clamp = !_compact_density_nonnegative_ok(u, layout)
    u_eval = u
    if needs_clamp
        u_work = if hasproperty(p, :work_u)
            p.work_u
        elseif hasproperty(p, :work_y)
            p.work_y
        else
            nothing
        end
        if u_work === nothing
            u_work = similar(u)
        end
        ok = _compact_project_density_nonnegative!(u_work, u, layout)
        if !ok
            fill!(du, 0.0)
            return nothing
        end
        u_eval = u_work
    end

    config = hasproperty(p, :config) ? p.config : nothing
    is_isothermal = config !== nothing && config.physics.is_isothermal_teex

    if !is_isothermal
        calculate_rhs_api_wrapper!(du, u_eval)
        if hasproperty(p, :residence_time)
            _apply_residence_time_term!(du, u_eval, p.residence_time)
        end
        if needs_clamp
            _compact_apply_shampine_positivity!(du, u, layout)
        end
        return nothing
    end

    if !layout.eex_noneq || layout.idx_eeex == 0
        error("Isothermal Teex requires eex_noneq=1 (electron-electronic energy equation present) in the API layout.")
    end
    if layout.idx_evib == 0
        error("Isothermal Teex requires vib_noneq=1 (vibrational energy present) so Tvib can be reconstructed.")
    end
    if config === nothing
        error("Isothermal Teex requires a config in the parameter tuple.")
    end
    teex_const = hasproperty(p, :teex_const) ? p.teex_const : config.temperatures.Te
    teex_vec = hasproperty(p, :teex_const_vec) ? p.teex_const_vec : nothing

    calculate_rhs_api_isothermal_teex_wrapper!(du, u_eval, teex_const;
        tex = teex_vec)

    if hasproperty(p, :residence_time)
        _apply_residence_time_term!(du, u_eval, p.residence_time)
    end
    if needs_clamp
        _compact_apply_shampine_positivity!(du, u, layout)
    end

    return nothing
end

function _electronic_state_print_metadata(rho_ex_initial::AbstractMatrix{<:Real},
        n_species::Int;
        ex_tol::Float64 = 1e-80)
    electronic_state_counts = zeros(Int, n_species)
    has_electronic_states = falses(n_species)

    @inbounds for isp in 1:n_species
        col = @view rho_ex_initial[:, isp]
        last_idx = findlast(x -> abs(x) > ex_tol, col)
        if last_idx === nothing
            electronic_state_counts[isp] = 0
            has_electronic_states[isp] = false
        else
            electronic_state_counts[isp] = last_idx
            has_electronic_states[isp] = true
        end
    end

    return (electronic_state_counts = electronic_state_counts,
        has_electronic_states = has_electronic_states)
end

function _print_terra_integration_output(t::Float64,
        dt::Float64,
        rho_sp::Vector{Float64},
        molecular_weights::Vector{Float64},
        temps,
        rho_etot::Float64;
        rho_ex::Union{Nothing, Matrix{Float64}} = nothing,
        rho_eeex::Union{Nothing, Float64} = nothing,
        rho_evib::Union{Nothing, Float64} = nothing,
        has_electronic_states::Union{Nothing, AbstractVector{Bool}} = nothing,
        electronic_state_counts::Union{Nothing, AbstractVector{Int}} = nothing)
    # Mass fraction sum error
    rho_total = sum(rho_sp)
    ysum = 0.0
    @inbounds for val in rho_sp
        ysum += val / rho_total
    end
    yerr = abs(ysum - 1.0)
    println(@sprintf(" Ytot,err   = % .3E", yerr))

    # Relative enthalpy change (%) at current Tt vs reconstructed energy
    Ecomp = calculate_total_energy_wrapper(temps.tt, rho_sp;
        rho_ex = rho_ex,
        rho_eeex = rho_eeex,
        rho_evib = rho_evib)
    dEnth = 100.0 * (Ecomp - rho_etot) / rho_etot
    println(@sprintf(" dEnth (%%)  = % .5E", dEnth))

    t_us = t * 1e6
    dt_us = dt * 1e6
    println(@sprintf(" time       = % .2E mu-s ", t_us))
    println(@sprintf(" dt         = % .2E mu-s", dt_us))
    println(@sprintf(" T(t,e,r,v) = % .3E % .3E % .3E % .3E K",
        Float64(temps.tt), Float64(temps.teex), Float64(temps.trot), Float64(temps.tvib)))

    # Mole fractions
    denom = 0.0
    @inbounds for i in eachindex(rho_sp, molecular_weights)
        denom += rho_sp[i] / molecular_weights[i]
    end
    xbuf = IOBuffer()
    print(xbuf, " X          =")
    @inbounds for i in eachindex(rho_sp, molecular_weights)
        xi = (rho_sp[i] / molecular_weights[i]) / denom
        print(xbuf, @sprintf(" % .3E", xi))
    end
    println(String(take!(xbuf)))

    tex = hasproperty(temps, :tex) ? temps.tex : nothing
    if rho_ex !== nothing && has_electronic_states !== nothing &&
       electronic_state_counts !== nothing && tex !== nothing
        n_species = length(molecular_weights)
        @inbounds for isp in 1:n_species
            if !has_electronic_states[isp]
                continue
            end
            mex = min(electronic_state_counts[isp], size(rho_ex, 1))
            if mex <= 0
                continue
            end
            if isp < 10
                println(@sprintf(" Tex(%d)     = % .3E", isp, Float64(tex[isp])))
            else
                println(@sprintf(" Tex(%d)    = % .3E", isp, Float64(tex[isp])))
            end

            states_to_show = min(mex, 7)
            nbuf = IOBuffer()
            if isp < 10
                print(nbuf, @sprintf(" n(%d,1:%d)   =", isp, states_to_show))
            else
                print(nbuf, @sprintf(" n(%d,1:%d)  =", isp, states_to_show))
            end

            for iex in 1:states_to_show
                nval = Float64(rho_ex[iex, isp]) / molecular_weights[isp] * AVOGADRO
                print(nbuf, @sprintf(" % .3E", nval))
            end
            println(String(take!(nbuf)))
        end
    end

    println()
    return nothing
end

"""
$(SIGNATURES)

Integrate the 0D system over time using the Fortran `rhs_api` layout.

This uses `terra_ode_system!` and constructs a `y` vector that matches the
ordering returned by `get_api_layout()`.

Notes:
- Vibrational STS is not yet supported in this wrapper (vibrational *mode* energy is supported).
- For isothermal Teex cases, the stored energy slot in `u` is `rho_rem` (legacy semantics), and
  the working `y` passed to Fortran is reconstructed each call.
- Optional residence-time (CSTR) terms can be enabled via `residence_time`.
"""
function integrate_0d_system(config::TERRAConfig, initial_state;
        residence_time::Union{Nothing, ResidenceTimeConfig} = nothing)
    dt = config.time_params.dt
    tlim = config.time_params.tlim

    @info "Setting up ODE integration" tlim=tlim

    native_outputs_requested = config.write_native_outputs
    print_integration_output = config.print_integration_output
    outputs_opened = false

    layout = get_api_layout()
    if layout.neq <= 0
        error("Invalid API layout: layout.neq=$(layout.neq). Ensure the Fortran API is initialized.")
    end
    if layout.nd != 0
        error("integrate_0d_system currently supports nd=0 only (got nd=$(layout.nd)).")
    end
    if layout.n_eq_vib > 0 || layout.is_vib_sts
        error("Vibrational STS is not yet supported (layout.n_eq_vib=$(layout.n_eq_vib)).")
    end

    n_species = layout.nsp
    is_isothermal = config.physics.is_isothermal_teex

    electronic_state_counts = zeros(Int, n_species)
    has_electronic_states = falses(n_species)
    if print_integration_output && layout.is_elec_sts
        meta = _electronic_state_print_metadata(initial_state.rho_ex, n_species)
        electronic_state_counts = meta.electronic_state_counts
        has_electronic_states = meta.has_electronic_states
    end

    # Initial state vector (rhs_api layout, CGS units)
    energy_scalar0 = is_isothermal ? initial_state.rho_rem : initial_state.rho_energy
    rho_ex0 = layout.is_elec_sts ? initial_state.rho_ex : nothing
    rho_eeex0 = layout.eex_noneq ? initial_state.rho_eeex : nothing
    rho_erot0 = layout.rot_noneq ? 0.0 : nothing
    rho_evib0 = layout.vib_noneq ? initial_state.rho_evib : nothing

    u0 = pack_state_vector(layout, initial_state.rho_sp, energy_scalar0;
        rho_ex = rho_ex0,
        rho_u = layout.nd >= 1 ? 0.0 : nothing,
        rho_v = layout.nd >= 2 ? 0.0 : nothing,
        rho_w = layout.nd >= 3 ? 0.0 : nothing,
        rho_eeex = rho_eeex0,
        rho_erot = rho_erot0,
        rho_evib = rho_evib0)

    @info "Initial state vector created" length_u0=length(u0) neq=layout.neq n_species=n_species

    tspan = (0.0, tlim)

    molecular_weights = get_molecular_weights(config.species)
    species_names = config.species
    gas_constants = initial_state.gas_constants
    teex_const_vec = fill(config.temperatures.Te, n_species)

    # Preallocate RHS work buffers
    work_y = similar(u0)
    work_u = similar(u0)
    work_rho_sp = zeros(Float64, layout.nsp)
    work_rho_ex = layout.is_elec_sts ? zeros(Float64, layout.mnex, layout.nsp) : nothing

    residence_time_data = (residence_time === nothing || !residence_time.enabled) ? nothing :
                          _prepare_residence_time_data(layout, config, u0, residence_time)

    p = (
        layout = layout,
        config = config,
        molecular_weights = molecular_weights,
        species = species_names,
        gas_constants = gas_constants,
        teex_const = initial_state.teex_const,
        teex_const_vec = teex_const_vec,
        work_u = work_u,
        residence_time = residence_time_data
    )

    prob = ODEProblem(terra_ode_system!, u0, tspan, p)
    ramp_callback = native_ramp_callback(dt; understep_ratio = inv(128), history_steps = 5)

    @info "ODE problem created, starting integration..."

    # Component-wise tolerances (CGS). The compact state mixes densities (g/cm^3)
    # and energies (erg/cm^3). A scalar abstol tends to make the energy equation
    # far more restrictive than the density equations. For isothermal Teex runs,
    # loosen the energy absolute tolerance based on a target temperature
    # resolution:
    #   abstol_E ≈ (3/2) kB n_total ΔT  [erg/cm^3]
    local_reltol = is_isothermal ? 1e-7 : 1e-8
    local_abstol_density = 1e-10
    local_abstol_vec = is_isothermal ? fill(local_abstol_density, layout.neq) : nothing
    if is_isothermal
        kB_erg_per_K = 1.380650524e-16
        n_total = initial_state.number_density
        rho_total = sum(initial_state.rho_sp)
        local_abstol_density = max(1e-20, 1e-12 * rho_total)

        deltaT = let s = strip(get(ENV, "TERRA_ISO_ENERGY_DELTAT_K", ""))
            isempty(s) ? 0.5 : parse(Float64, s)
        end
        energy_abstol = max(1.5 * kB_erg_per_K * n_total * deltaT, 1e-8)

        @inbounds begin
            # Density-like blocks (g/cm^3)
            local_abstol_vec[layout.vib_range] .= local_abstol_density
            local_abstol_vec[layout.elec_range] .= local_abstol_density
            local_abstol_vec[layout.sp_range] .= local_abstol_density
            local_abstol_vec[layout.mom_range] .= local_abstol_density

            # Energy-like block (erg/cm^3)
            local_abstol_vec[layout.energy_range] .= energy_abstol

            # In isothermal Teex mode, the rho_eeex slot is effectively fixed
            # (du=0); avoid letting it influence step control.
            if layout.idx_eeex != 0
                local_abstol_vec[layout.idx_eeex] = max(energy_abstol, 1.0)
            end
        end
    end

    try
        if native_outputs_requested && !outputs_opened
            open_api_output_files_wrapper()
            outputs_opened = true
        end

        sol = solve(prob;
            alg_hints = [:stiff],
            dt = dt,
            callback = ramp_callback,
            isoutofdomain = (u, _p, _t) -> (!_compact_nonnegative_ok(u, layout) || any(!isfinite, u)),
            reltol = local_reltol,
            abstol = (local_abstol_vec === nothing ? local_abstol_density :
                      local_abstol_vec),
            saveat = range(0.0, tlim; length = 100),
            save_everystep = false
        )

        @info "ODE integration completed" retcode=sol.retcode

        # Native output snapshots (optional)
        if outputs_opened
            first_dt = length(sol.t) >= 2 ? (sol.t[2] - sol.t[1]) : dt
            for (i, t) in enumerate(sol.t)
                local_dt = i == 1 ? first_dt : (t - sol.t[i - 1])
                if is_isothermal
                    _compact_isothermal_fill_fortran_y_work!(
                        work_y, work_rho_sp, work_rho_ex, sol.u[i], p, layout)
                    write_api_outputs_wrapper(
                        i - 1, t, local_dt, work_y; dist = 0.0, dx = 0.0)
                else
                    write_api_outputs_wrapper(
                        i - 1, t, local_dt, sol.u[i]; dist = 0.0, dx = 0.0)
                end
            end
        end

        # Extract results at output times
        n_times = length(sol.t)
        time_points = collect(sol.t)

        species_densities = zeros(n_species, n_times)
        temperatures_tt = Vector{Float64}(undef, n_times)
        temperatures_te = Vector{Float64}(undef, n_times)
        temperatures_tv = Vector{Float64}(undef, n_times)
        total_energies = Vector{Float64}(undef, n_times)

        rho_ex_arg = layout.is_elec_sts ? work_rho_ex : nothing
        first_dt = length(sol.t) >= 2 ? (sol.t[2] - sol.t[1]) : dt

        for i in 1:n_times
            t = sol.t[i]
            local_dt = i == 1 ? first_dt : (t - sol.t[i - 1])
            ui = sol.u[i]

            if is_isothermal
                iso = _compact_isothermal_fill_fortran_y_work!(
                    work_y, work_rho_sp, rho_ex_arg, ui, p, layout)
                species_densities[:, i] = work_rho_sp
                total_energies[i] = iso.rho_etot

                temps = iso.temps
                temperatures_tt[i] = temps.tt
                temperatures_te[i] = temps.teex
                temperatures_tv[i] = temps.tvib

                if print_integration_output
                    _print_terra_integration_output(
                        t, local_dt, work_rho_sp, molecular_weights, temps,
                        iso.rho_etot;
                        rho_ex = rho_ex_arg,
                        rho_eeex = iso.rho_eeex,
                        rho_evib = iso.rho_evib,
                        has_electronic_states = has_electronic_states,
                        electronic_state_counts = electronic_state_counts)
                end
            else
                _reconstruct_rho_sp_rho_ex_from_compact!(
                    work_rho_sp, rho_ex_arg, ui, layout)
                species_densities[:, i] = work_rho_sp

                rho_etot = Float64(ui[layout.idx_etot])
                total_energies[i] = rho_etot

                temps = calculate_temperatures_wrapper(work_rho_sp, rho_etot;
                    rho_ex = rho_ex_arg,
                    rho_erot = layout.idx_erot == 0 ? nothing :
                               Float64(ui[layout.idx_erot]),
                    rho_eeex = layout.idx_eeex == 0 ? nothing :
                               Float64(ui[layout.idx_eeex]),
                    rho_evib = layout.idx_evib == 0 ? nothing :
                               Float64(ui[layout.idx_evib]))
                temperatures_tt[i] = temps.tt
                temperatures_te[i] = temps.teex
                temperatures_tv[i] = temps.tvib

                if print_integration_output
                    _print_terra_integration_output(
                        t, local_dt, work_rho_sp, molecular_weights, temps,
                        rho_etot;
                        rho_ex = rho_ex_arg,
                        rho_eeex = layout.idx_eeex == 0 ? nothing :
                                   Float64(ui[layout.idx_eeex]),
                        rho_evib = layout.idx_evib == 0 ? nothing :
                                   Float64(ui[layout.idx_evib]),
                        has_electronic_states = has_electronic_states,
                        electronic_state_counts = electronic_state_counts)
                end
            end
        end

        # Convert results back to SI units if needed
        if config.unit_system == :SI
            species_densities_si = zeros(size(species_densities))
            for i in axes(species_densities, 2)
                species_densities_si[:, i] = convert_density_cgs_to_si(species_densities[
                    :, i])
            end
            total_energies_si = [convert_energy_density_cgs_to_si(e)
                                 for e in total_energies]
        else
            species_densities_si = species_densities
            total_energies_si = total_energies
        end

        temperatures = (tt = temperatures_tt, te = temperatures_te, tv = temperatures_tv)

        rc = sol.retcode
        success = rc isa Symbol ? (rc in (:Success, :Terminated)) :
                  (occursin("Success", string(rc)) || occursin("Terminated", string(rc)))
        message = success ? "ODE integration completed successfully" :
                  "ODE integration terminated: $(rc)"

        return TERRAResults(
            time_points,
            species_densities_si,
            temperatures,
            total_energies_si,
            nothing,
            success,
            message
        )

    catch e
        @error "ODE integration failed" exception=e
        return TERRAResults(
            [0.0],
            reshape(initial_state.rho_sp, :, 1),
            (tt = [config.temperatures.Tt], te = [config.temperatures.Te],
                tv = [config.temperatures.Tv],
                tee = [config.temperatures.Te]),
            [initial_state.rho_energy],
            nothing,
            false,
            "ODE integration failed: $(string(e))"
        )
    finally
        if outputs_opened
            try
                close_api_output_files_wrapper()
            catch close_err
                @warn "Failed to close native TERRA outputs after integration" exception=close_err
            end
            outputs_opened = false
        end
    end
end

"""
$(SIGNATURES)

Integrate the 0D system using nested `Config`.
"""
function integrate_0d_system(config::Config, initial_state;
        residence_time::Union{Nothing, ResidenceTimeConfig} = config.numerics.residence_time,
        use_residence_time::Union{Nothing, Bool} = nothing)
    effective_residence_time = _resolve_residence_time(residence_time, use_residence_time)
    return integrate_0d_system(
        to_legacy_config(config),
        initial_state;
        residence_time = effective_residence_time)
end

"""
$(SIGNATURES)

Solve a 0D TERRA simulation.

This is the main high-level interface for running TERRA simulations.
It handles all the complexity of data conversion, Fortran interfacing,
and result processing.

# Arguments
- `config::TERRAConfig`: Configuration for the simulation

# Returns
- `TERRAResults`: Results of the simulation

# Throws
- `ErrorException` if TERRA not initialized or simulation fails
"""
function solve_terra_0d(config::TERRAConfig;
        residence_time::Union{Nothing, ResidenceTimeConfig} = nothing)
    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        @info "Starting TERRA 0D simulation" species=config.species

        # Convert configuration to initial conditions (SI to CGS)
        initial_state = config_to_initial_state(config)

        # Run the time integration
        results = integrate_0d_system(config, initial_state; residence_time = residence_time)

        @info "TERRA simulation completed successfully"
        return results

    catch e
        @error "TERRA simulation failed" exception=e
        return TERRAResults(
            Float64[], zeros(0, 0), (;), Float64[], nothing, false,
            "Simulation failed: $(string(e))"
        )
    end
end

"""
$(SIGNATURES)

Solve a 0D TERRA simulation using nested `Config`.
"""
function solve_terra_0d(config::Config;
        residence_time::Union{Nothing, ResidenceTimeConfig} = config.numerics.residence_time,
        use_residence_time::Union{Nothing, Bool} = nothing)
    effective_residence_time = _resolve_residence_time(residence_time, use_residence_time)
    return solve_terra_0d(
        to_legacy_config(config);
        residence_time = effective_residence_time)
end

"""
$(SIGNATURES)

Run the 0D Nitrogen Te=10eV example case.

This function provides a convenient way to run the reference test case
that matches the TERRA example in `/terra/examples/0D_Nitrogen_Te_10eV`.
Requires the `TERRA_LIB_PATH` environment variable to point to the TERRA shared library.

# Arguments
- `case_path::String`: Case directory path (optional, creates temp directory if not provided)

# Returns
- `TERRAResults`: Results of the simulation

# Example
```julia
results = nitrogen_10ev_example()
```
"""
function nitrogen_10ev_example(case_path::String = mktempdir();
        isothermal::Bool = false)
    # Create configuration for the example case
    config = nitrogen_10ev_config(; isothermal = isothermal)

    # Update config with case path
    config_with_path = TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = config.time_params,
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = config.validate_species_against_terra,
        print_source_terms = config.print_source_terms,
        write_native_outputs = config.write_native_outputs,
        print_integration_output = config.print_integration_output
    )

    @info "Running 0D Nitrogen Te=10eV example case"
    @info "Configuration" species=config_with_path.species mole_fractions=config_with_path.mole_fractions
    @info "Temperatures" Tt=config_with_path.temperatures.Tt Te=config_with_path.temperatures.Te
    @info "Time parameters" dt=config_with_path.time_params.dt tlim=config_with_path.time_params.tlim
    @info "Case path" case_path=case_path

    try
        # Initialize TERRA with config
        initialize_terra(config_with_path, case_path)

        # Run simulation
        results = solve_terra_0d(config_with_path)

        if results.success
            @info "Example simulation completed successfully"
            @info "Final conditions" time=results.time[end]

            # Print final species densities
            for (i, species) in enumerate(config_with_path.species)
                final_density = results.species_densities[i, end]
                unit_str = config_with_path.unit_system == :SI ? "kg/m³" : "g/cm³"
                @info "Final density" species=species density=final_density unit=unit_str
            end

            @info "Final temperatures" Tt=results.temperatures.tt[end] Tv=results.temperatures.tv[end] Te=results.temperatures.te[end]
        else
            @error "Example simulation failed" message=results.message
        end

        return results

    finally
        # Clean up TERRA resources
        try
            finalize_terra()
        catch e
            @warn "Error during cleanup" exception=e
        end
    end
end

"""
$(SIGNATURES)

Validate simulation results for physical consistency.

# Arguments
- `results::TERRAResults`: Simulation results to validate

# Returns
- `true` if results pass validation, `false` otherwise
"""
function validate_results(results::TERRAResults)
    if !results.success
        @warn "Simulation was not successful"
        return false
    end

    # Check for negative densities
    if any(results.species_densities .< 0)
        @warn "Negative species densities found"
        return false
    end

    # Check for NaN or Inf values
    if any(isnan.(results.species_densities)) || any(isinf.(results.species_densities))
        @warn "NaN or Inf values found in species densities"
        return false
    end

    if any(isnan.(results.temperatures.tt)) || any(isinf.(results.temperatures.tt))
        @warn "NaN or Inf values found in temperatures"
        return false
    end

    # Check temperature ranges (should be positive and reasonable)
    if any(results.temperatures.tt .<= 0) || any(results.temperatures.tt .> 1e6)
        @warn "Unreasonable translational temperatures found"
        return false
    end

    if any(results.temperatures.te .<= 0) || any(results.temperatures.te .> 1e6)
        @warn "Unreasonable electron temperatures found"
        return false
    end

    @info "Results validation passed"
    return true
end

"""
$(SIGNATURES)

Save TERRA results to file.

# Arguments
- `results::TERRAResults`: Results to save
- `filename::String`: Output filename (CSV format)

# Returns
- `true` if save successful
"""
function save_results(results::TERRAResults, filename::String)
    try
        # Prepare data for CSV output
        n_times = length(results.time)
        n_species = size(results.species_densities, 1)

        # Create header
        header = ["time", "total_energy", "T_trans", "T_electron", "T_vib"]
        for i in 1:n_species
            push!(header, "species_$(i)_density")
        end

        # Create data matrix
        data = zeros(n_times, length(header))
        data[:, 1] = results.time
        data[:, 2] = results.total_energy
        data[:, 3] = results.temperatures.tt
        data[:, 4] = results.temperatures.te
        data[:, 5] = results.temperatures.tv

        for i in 1:n_species
            data[:, 5 + i] = results.species_densities[i, :]
        end

        # Write to file
        open(filename, "w") do io
            # Write header
            println(io, join(header, ","))

            # Write data
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end

        @info "Results saved successfully" filename=filename
        return true

    catch e
        @error "Failed to save results" filename=filename exception=e
        return false
    end
end
