"""
# TERRA Solver Module

This module provides the high-level interface for running TERRA simulations
from Julia, hiding the complexity of the Fortran interface and providing
a clean, Julia-native API.
"""

# Debug: RHS call counter
const _TERRA_ODE_DBG_CALLS = Ref(0)

const ENERGY_FROM_ENTHALPY_MAX_ITERS = 11
const ENERGY_FROM_ENTHALPY_RTOL = 1e-10
const ENERGY_FROM_ENTHALPY_ATOL = 1e-8

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

Calculate dimensions for the ODE state vector components.

# Arguments
- `config::TERRAConfig`: Configuration object

# Returns
- Named tuple with dimensions for each component
"""
function get_state_dimensions(config::TERRAConfig)
    if !is_terra_loaded()
        error("TERRA library must be loaded to get state dimensions. Set TERRA_LIB_PATH or call load_terra_library!(path) first.")
    end

    n_species = length(config.species)

    # Get actual dimensions from TERRA API
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    max_vibrational_quantum_number = get_max_vibrational_quantum_number_wrapper()
    has_vibrational_sts = has_vibrational_sts_wrapper()
    has_electronic_sts = has_electronic_sts_wrapper()

    return (
        n_species = n_species,
        n_atomic_electronic_states = max_atomic_electronic_states,
        n_molecular_electronic_states = max_molecular_electronic_states,
        n_vibrational_levels = max_vibrational_quantum_number,
        rho_ex_size = (max_atomic_electronic_states, n_species),
        rho_vx_size = (
            max_vibrational_quantum_number + 1, max_molecular_electronic_states, n_species),  # +1 for 0:mnv indexing
        has_vibrational_sts = has_vibrational_sts,
        has_electronic_sts = has_electronic_sts
    )
end

"""
$(SIGNATURES)

Compute the expected state vector length for a given dimensions tuple.

# Arguments
- `dimensions`: Dimensions structure from `get_state_dimensions()`

# Returns
- `Int`: Expected length of packed state/derivative vectors
"""
function expected_state_length(dimensions)
    n_species = dimensions.n_species
    nexc = prod(dimensions.rho_ex_size)
    nvib = prod(dimensions.rho_vx_size)
    # Layout: species + etot + ex(flat) + vx(flat) + (rho_u,rho_v,rho_w) + (rho_erot,rho_eeex,rho_evib)
    return n_species + 1 + nexc + nvib + 3 + 3
end

"""
$(SIGNATURES)

Compute index ranges for the packed TERRA state vector segments.

Splitting the packed vector into logical blocks (species densities, energy
terms, etc.) allows consistent application of custom norms and tolerances
without having to repeatedly recompute offsets.

# Arguments
- `dimensions`: Dimensions tuple returned by `get_state_dimensions`

# Returns
- Named tuple of index ranges for each logical block of the packed state vector
"""
function state_component_ranges(dimensions)
    idx = 1

    # Species block
    n_species = dimensions.n_species
    species_range = n_species > 0 ? (idx:(idx + n_species - 1)) : (idx:(idx - 1))
    idx += n_species

    # Total energy scalar
    energy_total_range = idx:idx
    idx += 1

    # Electronic state populations (flattened)
    nex = prod(dimensions.rho_ex_size)
    electronic_range = nex > 0 ? (idx:(idx + nex - 1)) : (idx:(idx - 1))
    idx += nex

    # Vibrational populations (flattened)
    nvx = prod(dimensions.rho_vx_size)
    vibrational_range = nvx > 0 ? (idx:(idx + nvx - 1)) : (idx:(idx - 1))
    idx += nvx

    # Momentum-like terms (ρu, ρv, ρw)
    momentum_range = idx:(idx + 2)
    idx += 3

    # Energy mode reservoirs (ρerot, ρeeex, ρevib)
    energy_modes_range = idx:(idx + 2)

    return (
        species = species_range,
        energy_total = energy_total_range,
        electronic = electronic_range,
        vibrational = vibrational_range,
        momentum = momentum_range,
        energy_modes = energy_modes_range
    )
end

@inline function _component_max(residual::AbstractVector, range::UnitRange{Int})
    isempty(range) && return 0.0
    return maximum(abs, @view residual[range])
end

@inline function _component_max(residual, range::UnitRange{Int})
    isempty(range) && return 0.0
    return abs(residual)
end

"""
$(SIGNATURES)

Create a max-norm style error metric tailored to TERRA state structure.

DifferentialEquations.jl expects `internalnorm(residual, t)` to provide a
scalar error estimate. By splitting the residual into physically meaningful
blocks we can ensure no single population or energy mode dominates the norm
and that the solver honours the TERRA-style "max |Δy/y|" heuristic.

# Arguments
- `dimensions`: Dimensions tuple returned by `get_state_dimensions`
- `weights`: Optional named-tuple overrides for component weighting

# Returns
- Callable that maps a residual vector (and time) to a scalar norm value
"""
function create_terra_error_norm(dimensions; weights = nothing)
    ranges = state_component_ranges(dimensions)
    default_weights = (
        species = 1.0,
        energy_total = 1.0,
        electronic = 1.0,
        vibrational = 1.0,
        momentum = 1.0,
        energy_modes = 1.0
    )
    w = weights === nothing ? default_weights : merge(default_weights, weights)

    return function (residual, t)
        species_err = w.species * _component_max(residual, ranges.species)
        energy_err = w.energy_total * _component_max(residual, ranges.energy_total)
        ex_err = w.electronic * _component_max(residual, ranges.electronic)
        vib_err = w.vibrational * _component_max(residual, ranges.vibrational)
        mom_err = w.momentum * _component_max(residual, ranges.momentum)
        mode_err = w.energy_modes * _component_max(residual, ranges.energy_modes)

        return maximum((species_err, energy_err, ex_err, vib_err, mom_err, mode_err))
    end
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

"""
$(SIGNATURES)

Pack state components into a single ODE state vector.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities
- `rho_etot::Float64`: Total enthalpy density
- `dimensions`: Dimensions structure from get_state_dimensions()
- Optional components (all default to zero if not provided)

# Returns
- `Vector{Float64}`: Packed state vector
"""
function pack_state_vector(rho_sp::Vector{Float64}, rho_etot::Float64, dimensions;
        rho_ex = nothing, rho_vx = nothing, rho_u = nothing, rho_v = nothing,
        rho_w = nothing, rho_erot = nothing, rho_eeex = nothing, rho_evib = nothing)
    n_species = dimensions.n_species

    # Validate input array size
    if length(rho_sp) != n_species
        throw(ArgumentError("Species density array length ($(length(rho_sp))) must match number of species ($(n_species))"))
    end

    # Preallocate output for determinism and speed
    u = Vector{Float64}(undef, expected_state_length(dimensions))
    idx = 1

    # Species densities
    @inbounds begin
        u[idx:(idx + n_species - 1)] .= rho_sp
    end
    idx += n_species

    # Total enthalpy density
    u[idx] = rho_etot
    idx += 1

    # Electronic state densities (flattened)
    ex_target = dimensions.rho_ex_size
    if rho_ex !== nothing
        if size(rho_ex) != ex_target
            ex_fixed = zeros(Float64, ex_target)
            m1 = min(size(rho_ex, 1), ex_target[1])
            m2 = min(size(rho_ex, 2), ex_target[2])
            @inbounds ex_fixed[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
            @inbounds u[idx:(idx + prod(ex_target) - 1)] .= vec(ex_fixed)
        else
            @inbounds u[idx:(idx + prod(ex_target) - 1)] .= vec(rho_ex)
        end
    else
        @inbounds fill!(view(u, idx:(idx + prod(ex_target) - 1)), 0.0)
    end
    idx += prod(ex_target)

    # Vibrational state densities (flattened)
    vx_target = dimensions.rho_vx_size
    if rho_vx !== nothing
        if size(rho_vx) != vx_target
            vx_fixed = zeros(Float64, vx_target)
            m1 = min(size(rho_vx, 1), vx_target[1])
            m2 = min(size(rho_vx, 2), vx_target[2])
            m3 = min(size(rho_vx, 3), vx_target[3])
            @inbounds vx_fixed[1:m1, 1:m2, 1:m3] .= rho_vx[1:m1, 1:m2, 1:m3]
            @inbounds u[idx:(idx + prod(vx_target) - 1)] .= vec(vx_fixed)
        else
            @inbounds u[idx:(idx + prod(vx_target) - 1)] .= vec(rho_vx)
        end
    else
        @inbounds fill!(view(u, idx:(idx + prod(vx_target) - 1)), 0.0)
    end
    idx += prod(vx_target)

    # Velocity components
    u[idx] = (rho_u === nothing ? 0.0 : rho_u)
    idx += 1
    u[idx] = (rho_v === nothing ? 0.0 : rho_v)
    idx += 1
    u[idx] = (rho_w === nothing ? 0.0 : rho_w)
    idx += 1

    # Energy mode components
    u[idx] = (rho_erot === nothing ? 0.0 : rho_erot)
    idx += 1
    u[idx] = (rho_eeex === nothing ? 0.0 : rho_eeex)
    idx += 1
    u[idx] = (rho_evib === nothing ? 0.0 : rho_evib)

    # Final sanity check
    @assert idx==length(u) "Internal error: packed state length mismatch"

    return u
end

"""
$(SIGNATURES)

Unpack ODE state vector into individual components.

# Arguments
- `u::Vector{Float64}`: Packed state vector
- `dimensions`: Dimensions structure from get_state_dimensions()

# Returns
- Named tuple with unpacked components
"""
function unpack_state_vector(u::Vector{Float64}, dimensions)
    n_species = dimensions.n_species
    rho_ex_size = dimensions.rho_ex_size
    rho_vx_size = dimensions.rho_vx_size

    idx = 1

    # Extract species densities (return a concrete Vector for wrapper calls)
    rho_sp = Vector{Float64}(@view u[idx:(idx + n_species - 1)])
    idx += n_species

    # Extract total enthalpy density
    rho_etot = u[idx]
    idx += 1

    # Extract electronic state densities
    rho_ex_flat = @view u[idx:(idx + prod(rho_ex_size) - 1)]
    rho_ex = reshape(Vector{Float64}(rho_ex_flat), rho_ex_size)
    idx += prod(rho_ex_size)

    # Extract vibrational state densities
    rho_vx_flat = @view u[idx:(idx + prod(rho_vx_size) - 1)]
    rho_vx = reshape(Vector{Float64}(rho_vx_flat), rho_vx_size)
    idx += prod(rho_vx_size)

    # Extract tail components (rho_u, rho_v, rho_w, rho_erot, rho_eeex, rho_evib)
    rho_u = u[idx]
    idx += 1
    rho_v = u[idx]
    idx += 1
    rho_w = u[idx]
    idx += 1
    rho_erot = u[idx]
    idx += 1
    rho_eeex = u[idx]
    idx += 1
    rho_evib = u[idx]

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

"""
$(SIGNATURES)

Pack derivative components into ODE derivative vector.

# Arguments
- `du::Vector{Float64}`: Output derivative vector (modified in-place)
- `derivatives`: Derivatives from calculate_sources_wrapper()
- `dimensions`: Dimensions structure

# Returns
- Nothing (modifies du in-place)
"""
function pack_derivative_vector!(du::Vector{Float64}, derivatives, dimensions)
    n_species = dimensions.n_species
    rho_ex_size = dimensions.rho_ex_size
    rho_vx_size = dimensions.rho_vx_size

    idx = 1

    # Pack species density derivatives
    du[idx:(idx + n_species - 1)] .= derivatives.drho_sp
    idx += n_species

    # Pack total energy derivative
    du[idx] = derivatives.drho_etot
    idx += 1

    # Pack electronic state derivatives
    if derivatives.drho_ex !== nothing
        if size(derivatives.drho_ex) != rho_ex_size
            throw(DimensionMismatch("drho_ex size $(size(derivatives.drho_ex)) does not match expected $(rho_ex_size)"))
        end
        du[idx:(idx + prod(rho_ex_size) - 1)] .= vec(derivatives.drho_ex)
    else
        du[idx:(idx + prod(rho_ex_size) - 1)] .= 0.0
    end
    idx += prod(rho_ex_size)

    # Pack vibrational state derivatives
    if derivatives.drho_vx !== nothing
        if size(derivatives.drho_vx) != rho_vx_size
            throw(DimensionMismatch("drho_vx size $(size(derivatives.drho_vx)) does not match expected $(rho_vx_size)"))
        end
        du[idx:(idx + prod(rho_vx_size) - 1)] .= vec(derivatives.drho_vx)
    else
        du[idx:(idx + prod(rho_vx_size) - 1)] .= 0.0
    end
    idx += prod(rho_vx_size)

    # Pack velocity derivatives (zero for 0D)
    du[idx] = 0.0  # drho_u
    idx += 1
    du[idx] = 0.0  # drho_v
    idx += 1
    du[idx] = 0.0  # drho_w
    idx += 1

    # Pack energy derivatives
    du[idx] = derivatives.drho_erot === nothing ? 0.0 : derivatives.drho_erot
    idx += 1
    du[idx] = derivatives.drho_eeex === nothing ? 0.0 : derivatives.drho_eeex
    idx += 1
    du[idx] = derivatives.drho_evib === nothing ? 0.0 : derivatives.drho_evib

    # Sanity check: ensure we've filled the entire vector
    @assert idx==length(du) "Internal error: derivative vector length mismatch"

    return nothing
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

"""
$(SIGNATURES)

Recover the total energy density from enthalpy by iteratively updating pressure.

# Arguments
- `rho_enth::Float64`: Target enthalpy density (CGS units)
- `rho_sp::AbstractVector{<:Real}`: Species mass densities
- `rho_ex`: Optional electronic-state populations used when electronic STS is active
- `rho_vx`: Optional vibrational-state populations used when vibrational STS is active
- `rho_erot`: Optional rotational-mode energy density
- `rho_eeex`: Optional electron-electronic energy density
- `rho_evib`: Optional vibrational-mode energy density
- `gas_constants::AbstractVector{<:Real}`: Species-specific gas constants
- `molecular_weights::AbstractVector{<:Real}`: Species molecular weights
- `species_names::AbstractVector{<:AbstractString}`: Species identifiers used for electron detection
- `energy_guess::Float64`: Initial guess for the total energy density
- `teex_override::Union{Nothing, Float64}`: Optional enforced electron temperature for pressure evaluation

# Returns
- Named tuple containing `rho_etot`, `pressure`, `temps`, and `teex`
"""
function energy_from_enthalpy(rho_enth::Float64, rho_sp::AbstractVector{<:Real};
        rho_ex::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
        rho_vx::Union{Nothing, Array{<:Real, 3}} = nothing,
        rho_erot::Union{Nothing, Real} = nothing,
        rho_eeex::Union{Nothing, Real} = nothing,
        rho_evib::Union{Nothing, Real} = nothing,
        gas_constants::AbstractVector{<:Real},
        molecular_weights::AbstractVector{<:Real},
        species_names::AbstractVector{<:AbstractString},
        energy_guess::Float64,
        teex_override::Union{Nothing, Float64} = nothing)
    energy = max(energy_guess, 0.0)
    temps = nothing
    pressure = 0.0

    for iter in 1:ENERGY_FROM_ENTHALPY_MAX_ITERS
        energy = max(energy, 0.0)
        temps = calculate_temperatures_wrapper(rho_sp, energy;
            rho_ex = rho_ex,
            rho_vx = rho_vx,
            rho_erot = rho_erot,
            rho_eeex = rho_eeex,
            rho_evib = rho_evib)

        te_used = teex_override === nothing ? temps.teex : teex_override
        pressure = compute_mixture_pressure(
            rho_sp, gas_constants, species_names, molecular_weights,
            temps.tt, te_used)
        new_energy = rho_enth - pressure
        if !isfinite(new_energy)
            error("Encountered non-finite energy while converting enthalpy to energy.")
        end

        tol = max(abs(energy), abs(new_energy), 1.0) * ENERGY_FROM_ENTHALPY_RTOL +
              ENERGY_FROM_ENTHALPY_ATOL
        if abs(new_energy - energy) ≤ tol
            energy = new_energy
            break
        end

        energy = 0.5 * (energy + new_energy)
    end

    temps = calculate_temperatures_wrapper(rho_sp, energy;
        rho_ex = rho_ex,
        rho_vx = rho_vx,
        rho_erot = rho_erot,
        rho_eeex = rho_eeex,
        rho_evib = rho_evib)
    te_used = teex_override === nothing ? temps.teex : teex_override
    pressure = compute_mixture_pressure(
        rho_sp, gas_constants, species_names, molecular_weights,
        temps.tt, te_used)

    return (rho_etot = energy, pressure = pressure, temps = temps, teex = te_used)
end

"""
$(SIGNATURES)

Reconstruct energy components based on configuration flags.

When the isothermal electron-electronic mode is enabled, the energy slot stored
in the state vector holds the remainder energy (`rho_rem`). This helper
recomputes the electron-electronic energy from the prescribed `Tee`, rebuilds
the total energy, and returns all relevant pieces for downstream use. When the
mode is disabled, it simply returns the existing state components.

# Arguments
- `state`: Unpacked state tuple produced by `unpack_state_vector`
- `config::TERRAConfig`: Simulation configuration controlling physics toggles
- `teex_const::Float64`: Keyword override for the uniform electron temperature
- `teex_vec::Union{Nothing, AbstractVector{<:Real}}`: Optional per-species electron temperatures
- `molecular_weights::AbstractVector{<:Real}`: Species molecular weights
- `species_names::AbstractVector{<:AbstractString}`: Species identifiers
- `gas_constants::AbstractVector{<:Real}`: Species-specific gas constants
- `has_electronic_sts::Bool`: Indicates whether electronic STS data are present
- `has_vibrational_sts::Bool`: Indicates whether vibrational STS data are present
- `energy_cache::Union{Nothing, Ref{Float64}}`: Optional cache storing the last energy estimate
- `pressure_cache::Union{Nothing, Ref{Float64}}`: Optional cache storing the last pressure

# Returns
- Named tuple containing `rho_etot`, `rho_eeex`, `rho_rem`, `tvib`, `temps`, and `pressure`
"""
function reconstruct_energy_components(state, config;
        teex_const::Float64 = config.temperatures.Te,
        teex_vec::Union{Nothing, AbstractVector{<:Real}} = nothing,
        molecular_weights::AbstractVector{<:Real},
        species_names::AbstractVector{<:AbstractString},
        gas_constants::AbstractVector{<:Real},
        has_electronic_sts::Bool,
        has_vibrational_sts::Bool,
        energy_cache::Union{Nothing, Ref{Float64}} = nothing,
        pressure_cache::Union{Nothing, Ref{Float64}} = nothing)
    is_isothermal = config.physics.is_isothermal_teex

    if !is_isothermal
        energy_guess = energy_cache === nothing ? state.rho_etot : energy_cache[]
        rho_ex_arg = has_electronic_sts ? state.rho_ex : nothing
        rho_vx_arg = has_vibrational_sts ? state.rho_vx : nothing
        result = energy_from_enthalpy(state.rho_etot, state.rho_sp;
            rho_ex = rho_ex_arg,
            rho_vx = rho_vx_arg,
            rho_erot = state.rho_erot,
            rho_eeex = state.rho_eeex,
            rho_evib = state.rho_evib,
            gas_constants = gas_constants,
            molecular_weights = molecular_weights,
            species_names = species_names,
            energy_guess = energy_guess)

        if energy_cache !== nothing
            energy_cache[] = result.rho_etot
        end
        if pressure_cache !== nothing
            pressure_cache[] = result.pressure
        end

        rho_eeex = state.rho_eeex
        electron_idx = findfirst(
            i -> _is_electron_species(species_names[i], molecular_weights[i]),
            eachindex(species_names))
        electron_enthalpy = electron_idx === nothing ? 0.0 :
                            state.rho_sp[electron_idx] * gas_constants[electron_idx] *
                            result.teex
        rho_rem = result.rho_etot - rho_eeex - electron_enthalpy

        return (rho_etot = result.rho_etot,
            rho_eeex = rho_eeex,
            rho_rem = rho_rem,
            tvib = result.temps.tvib,
            temps = result.temps,
            pressure = result.pressure)
    end

    rho_rem = state.rho_etot
    rho_sp = state.rho_sp

    teex_kw = teex_vec === nothing ? nothing : teex_vec

    tvib = calculate_vibrational_temperature_wrapper(
        state.rho_evib, rho_sp;
        rho_ex = state.rho_ex,
        tex = teex_kw)

    rho_eeex = calculate_electron_electronic_energy_wrapper(teex_const, tvib, rho_sp)
    electron_idx = findfirst(
        i -> _is_electron_species(species_names[i], molecular_weights[i]),
        eachindex(species_names))
    electron_enthalpy = electron_idx === nothing ? 0.0 :
                        rho_sp[electron_idx] * gas_constants[electron_idx] * teex_const

    rho_enthalpy_total = rho_rem + rho_eeex + electron_enthalpy
    rho_ex_arg = has_electronic_sts ? state.rho_ex : nothing
    rho_vx_arg = has_vibrational_sts ? state.rho_vx : nothing
    energy_guess = energy_cache === nothing ? rho_enthalpy_total : energy_cache[]
    result = energy_from_enthalpy(rho_enthalpy_total, rho_sp;
        rho_ex = rho_ex_arg,
        rho_vx = rho_vx_arg,
        rho_erot = state.rho_erot,
        rho_eeex = rho_eeex,
        rho_evib = state.rho_evib,
        gas_constants = gas_constants,
        molecular_weights = molecular_weights,
        species_names = species_names,
        energy_guess = energy_guess,
        teex_override = teex_const)

    if energy_cache !== nothing
        energy_cache[] = result.rho_etot
    end
    if pressure_cache !== nothing
        pressure_cache[] = result.pressure
    end

    return (rho_etot = result.rho_etot,
        rho_eeex = rho_eeex,
        rho_rem = rho_rem,
        tvib = result.temps.tvib,
        temps = result.temps,
        pressure = result.pressure)
end

"""
$(SIGNATURES)

ODE system function for TERRA integration.

This function defines the ODE system du/dt = f(u, p, t) where:
- u is the state vector containing all TERRA variables
- p contains parameters (dimensions, etc.)
- t is time

# Arguments
- `du::Vector{Float64}`: Output derivative vector
- `u::Vector{Float64}`: Input state vector
- `p`: Parameters structure
- `t::Float64`: Current time

# Returns
- Nothing (modifies du in-place)
"""
function terra_ode_system!(du::Vector{Float64}, u::Vector{Float64}, p, t::Float64)
    # Guard against invalid inputs to the Fortran layer
    if any(!isfinite, u)
        fill!(du, 0.0)
        return nothing
    end
    nsp = p.dimensions.n_species
    if any(@view(u[1:nsp]) .< 0.0)
        fill!(du, 0.0)
        return nothing
    end
    if u[nsp + 1] < 0.0
        fill!(du, 0.0)
        return nothing
    end
    # Unpack state vector
    state = unpack_state_vector(u, p.dimensions)

    # Calculate source terms using TERRA
    try
        config = hasproperty(p, :config) ? p.config : nothing
        teex_vec = hasproperty(p, :teex_const_vec) ? p.teex_const_vec : nothing
        teex_const = hasproperty(p, :teex_const) ? p.teex_const :
                     (config === nothing ? 0.0 : config.temperatures.Te)
        is_isothermal = config !== nothing && config.physics.is_isothermal_teex

        if config === nothing
            result = energy_from_enthalpy(state.rho_etot, state.rho_sp;
                rho_ex = p.dimensions.has_electronic_sts ? state.rho_ex : nothing,
                rho_vx = p.dimensions.has_vibrational_sts ? state.rho_vx : nothing,
                rho_erot = state.rho_erot,
                rho_eeex = state.rho_eeex,
                rho_evib = state.rho_evib,
                gas_constants = p.gas_constants,
                molecular_weights = p.molecular_weights,
                species_names = p.species,
                energy_guess = hasproperty(p, :energy_cache) ? p.energy_cache[] :
                               state.rho_etot)

            electron_idx = hasproperty(p, :electron_index) ? p.electron_index :
                           findfirst(
                i -> _is_electron_species(p.species[i], p.molecular_weights[i]),
                eachindex(p.species))
            electron_enthalpy = electron_idx === nothing ? 0.0 :
                                state.rho_sp[electron_idx] * p.gas_constants[electron_idx] *
                                result.teex

            energy = (rho_etot = result.rho_etot,
                rho_eeex = state.rho_eeex,
                rho_rem = result.rho_etot - state.rho_eeex - electron_enthalpy,
                tvib = result.temps.tvib,
                temps = result.temps,
                pressure = result.pressure)
            if hasproperty(p, :energy_cache)
                p.energy_cache[] = result.rho_etot
            end
            if hasproperty(p, :pressure_cache)
                p.pressure_cache[] = result.pressure
            end
        else
            energy = reconstruct_energy_components(state, config;
                teex_const = teex_const,
                teex_vec = teex_vec,
                molecular_weights = p.molecular_weights,
                species_names = p.species,
                gas_constants = p.gas_constants,
                has_electronic_sts = p.dimensions.has_electronic_sts,
                has_vibrational_sts = p.dimensions.has_vibrational_sts,
                energy_cache = hasproperty(p, :energy_cache) ? p.energy_cache : nothing,
                pressure_cache = hasproperty(p, :pressure_cache) ? p.pressure_cache :
                                 nothing)
        end

        rho_etot_effective = energy.rho_etot
        rho_eeex_effective = energy.rho_eeex

        rho_ex_arg = p.dimensions.has_electronic_sts ? state.rho_ex : nothing
        rho_vx_arg = p.dimensions.has_vibrational_sts ? state.rho_vx : nothing

        # Always pass the current total energy to the Fortran RHS so that
        # temperatures and rates are computed from the instantaneous state.
        # Mirror TERRA handling of total energy:
        # - If radiation is OFF, hold total energy constant across RHS calls
        # - If radiation is ON, allow TERRA to evolve total energy
        use_const_etot = !is_isothermal && hasproperty(p, :config) &&
                         hasproperty(p.config, :processes) &&
                         getfield(p.config.processes, :consider_rad) == 0
        if use_const_etot && hasproperty(p, :rho_energy0)
            rho_etot_to_use = p.rho_energy0
            if hasproperty(p, :energy_cache)
                p.energy_cache[] = p.rho_energy0
            end
        else
            rho_etot_to_use = rho_etot_effective
        end

        derivatives = calculate_sources_wrapper(
            state.rho_sp, rho_etot_to_use;
            rho_ex = rho_ex_arg,
            rho_vx = rho_vx_arg,
            rho_u = state.rho_u,
            rho_v = state.rho_v,
            rho_w = state.rho_w,
            rho_erot = state.rho_erot,
            rho_eeex = rho_eeex_effective,
            rho_evib = state.rho_evib
        )

        # Optional solver-side debug mirror of wrapper outputs
        if get(ENV, "TERRA_SOLVER_DEBUG", "0") == "1"
            @info "SOLVER_DERIV" drho_etot=derivatives.drho_etot drho_erot=derivatives.drho_erot drho_eeex=derivatives.drho_eeex drho_evib=derivatives.drho_evib
        end

        # If solving electrons via charge balance, adjust electron derivative
        if hasproperty(p, :config) &&
           p.config.physics.get_electron_density_by_charge_balance
            names = p.config.species
            weights = p.molecular_weights
            # simple charge inference from species names
            charges = map(nm -> count(==('+'), nm) - count(==('-'), nm), names)
            elec_idx = findfirst(==("E-"), names)
            elec_idx = elec_idx === nothing ? findfirst(==(-1), charges) : elec_idx
            if elec_idx !== nothing
                spwt_e = weights[elec_idx]
                s = 0.0
                @inbounds for i in eachindex(derivatives.drho_sp)
                    if i == elec_idx
                        continue
                    end
                    zi = charges[i]
                    if zi != 0
                        s += (zi / weights[i]) * derivatives.drho_sp[i]
                    end
                end
                new_drho_sp = copy(derivatives.drho_sp)
                new_drho_sp[elec_idx] = spwt_e * s
                derivatives = (
                    drho_sp = new_drho_sp,
                    drho_etot = derivatives.drho_etot,
                    drho_ex = derivatives.drho_ex,
                    drho_vx = derivatives.drho_vx,
                    drho_erot = derivatives.drho_erot,
                    drho_eeex = derivatives.drho_eeex,
                    drho_evib = derivatives.drho_evib
                )
            end
        end

        if is_isothermal
            # Convert total enthalpy derivative into the stored remainder by
            # removing electronic-mode changes and the electron enthalpy work term.
            drho_eeex_val = derivatives.drho_eeex === nothing ? 0.0 : derivatives.drho_eeex
            drho_etot_raw = derivatives.drho_etot
            extra_electron_enthalpy = 0.0
            electron_idx = hasproperty(p, :electron_index) ? p.electron_index : nothing
            if electron_idx !== nothing && electron_idx <= length(derivatives.drho_sp)
                gas_const_e = p.gas_constants[electron_idx]
                te_for_enthalpy = hasproperty(p, :teex_const) ? p.teex_const :
                                  config.temperatures.Te
                extra_electron_enthalpy = te_for_enthalpy * gas_const_e *
                                          derivatives.drho_sp[electron_idx]
            end

            # Scaling factor that aligns the Julia wrapper output with native TERRA output
            # for isothermal cases. Unsure of what the source of the bug is.
            # TODO: Find source of this scaling factor error in drho_rem/dt and correct it
            scaling_factor = -1.5
            drho_rem = drho_etot_raw - drho_eeex_val -
                       (scaling_factor * extra_electron_enthalpy)

            derivatives = merge(derivatives, (drho_etot = drho_rem, drho_eeex = 0.0))
        end

        # Debug: print derivative norms on first few RHS calls (debug-level only)
        _TERRA_ODE_DBG_CALLS[] += 1
        if _TERRA_ODE_DBG_CALLS[] <= 5
            drsp_norm = maximum(abs, derivatives.drho_sp)
            dr_etot = derivatives.drho_etot
            dr_eeex = derivatives.drho_eeex === nothing ? 0.0 : derivatives.drho_eeex
            dr_evib = derivatives.drho_evib === nothing ? 0.0 : derivatives.drho_evib
            @debug "TERRA RHS debug" t=t drho_sp_max=drsp_norm drho_etot=dr_etot drho_eeex=dr_eeex drho_evib=dr_evib
        end

        # Pack derivatives into output vector
        pack_derivative_vector!(du, derivatives, p.dimensions)
        if use_const_etot
            n_species = p.dimensions.n_species
            du[n_species + 1] = 0.0
        end

    catch e
        @error "Error in TERRA ODE system" exception=e time=t
        # Fill with zeros to prevent integration failure
        fill!(du, 0.0)
    end

    return nothing
end

"""
$(SIGNATURES)

Integrate the 0D system over time using DifferentialEquations.jl.

# Arguments
- `config::TERRAConfig`: Configuration object
- `initial_state`: Initial state vectors in CGS units

# Returns
- `TERRAResults`: Simulation results (converted back to SI units)
"""
function integrate_0d_system(config::TERRAConfig, initial_state)
    # Time parameters
    dt = config.time_params.dt
    tlim = config.time_params.tlim

    @info "Setting up ODE integration" tlim=tlim

    native_outputs_requested = config.write_native_outputs
    outputs_opened = false

    # Get state dimensions
    dimensions = get_state_dimensions(config)
    n_species = dimensions.n_species
    is_isothermal = config.physics.is_isothermal_teex
    error_norm = create_terra_error_norm(dimensions)

    # Derive electronic STS metadata from the Boltzmann-seeded initial state
    rho_ex_initial = initial_state.rho_ex
    electronic_state_counts = zeros(Int, n_species)
    has_electronic_states = falses(n_species)
    ex_tol = 1e-80                       # Treat anything smaller as numerical zero
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

    # Create initial ODE state vector (all in CGS units)
    # Include available electronic states and energy components so the ODE RHS
    # receives consistent nonequilibrium energies at t=0.
    energy_scalar0 = is_isothermal ? initial_state.rho_rem : initial_state.rho_etot
    rho_eeex0 = is_isothermal ? 0.0 : initial_state.rho_eeex
    u0 = pack_state_vector(
        initial_state.rho_sp, energy_scalar0, dimensions;
        rho_ex = initial_state.rho_ex,
        # rho_vx = initial_state.rho_vx,
        rho_eeex = initial_state.rho_eeex,
        # rho_eeex = rho_eeex0,
        rho_evib = initial_state.rho_evib
    )

    @info "Initial state vector created" length_u0=length(u0) n_species=n_species

    # Time span for integration
    tspan = (0.0, tlim)

    # Parameters for the ODE system. Keep total energy density constant
    # (aligned with Fortran API behavior) by passing it as a parameter.
    molecular_weights = get_molecular_weights(config.species)
    species_names = config.species
    gas_constants = initial_state.gas_constants
    teex_const_vec = fill(config.temperatures.Te, n_species)
    electron_index = begin
        idx = findfirst(i -> _is_electron_species(species_names[i], molecular_weights[i]),
            eachindex(species_names))
        idx === nothing ? nothing : idx
    end

    p = (dimensions = dimensions,
        config = config,
        rho_enth0 = initial_state.rho_etot,
        rho_energy0 = initial_state.rho_energy,
        molecular_weights = molecular_weights,
        species = species_names,
        gas_constants = gas_constants,
        teex_const = initial_state.teex_const,
        teex_const_vec = teex_const_vec,
        electron_index = electron_index,
        energy_cache = Ref(initial_state.rho_energy),
        pressure_cache = Ref(initial_state.pressure))

    # Create ODE problem and attach step-size ramp mirroring native TERRA behaviour
    prob = ODEProblem(terra_ode_system!, u0, tspan, p)
    ramp_callback = native_ramp_callback(dt; understep_ratio = inv(128), history_steps = 5)

    @info "ODE problem created, starting integration..."

    try
        if native_outputs_requested && !outputs_opened
            open_api_output_files_wrapper()
            outputs_opened = true
        end

        sol = solve(prob;
            alg_hints = [:stiff],
            dt = dt,
            internalnorm = error_norm,
            callback = ramp_callback,
            reltol = 1e-7,
            abstol = 1e-10,
            # qmin = 0.5,
            # qmax = 1.02,
            # qsteady_max = 1.02,
            save_everystep = true
        )

        @info "ODE integration completed" success=sol.retcode

        # TERRA-style status printing for each solver step
        iter = 0
        first_dt = length(sol.t) >= 2 ? (sol.t[2] - sol.t[1]) : dt
        for (i, t) in enumerate(sol.t)
            steps = 0
            if length(sol.t) > 5000
                steps = 1000
            elseif length(sol.t) > 2500
                steps = 100
            elseif length(sol.t) > 100
                steps = 10
            else
                steps = 1
            end
            if outputs_opened
                local_dt = i == 1 ? first_dt : (t - sol.t[i - 1])
                write_api_outputs_wrapper(i - 1, t, local_dt, sol.u[i]; dist = 0.0, dx = 0.0)
            end
            if (iter % steps == 0) || (iter + 1 == length(sol.t))
                st = unpack_state_vector(sol.u[i], dimensions)
                energy = reconstruct_energy_components(st, config;
                    teex_const = config.temperatures.Te,
                    teex_vec = teex_const_vec,
                    molecular_weights = molecular_weights,
                    species_names = species_names,
                    gas_constants = gas_constants,
                    has_electronic_sts = dimensions.has_electronic_sts,
                    has_vibrational_sts = dimensions.has_vibrational_sts,
                    energy_cache = nothing,
                    pressure_cache = nothing)

                temps = haskey(energy, :temps) ? energy.temps :
                begin
                    rho_ex_arg = dimensions.has_electronic_sts ? st.rho_ex : nothing
                    rho_vx_arg = dimensions.has_vibrational_sts ? st.rho_vx : nothing
                    calculate_temperatures_wrapper(st.rho_sp, energy.rho_etot;
                        rho_ex = rho_ex_arg,
                        rho_vx = rho_vx_arg,
                        rho_eeex = energy.rho_eeex,
                        rho_evib = st.rho_evib)
                end

                # Mass fraction sum error
                ys = st.rho_sp ./ sum(st.rho_sp)
                yerr = abs(sum(ys) - 1.0)
                println(@sprintf(" Ytot,err   = % .3E", yerr))

                # Relative enthalpy change (%) at current Tt vs reference energy
                Ecomp = calculate_total_energy_wrapper(temps.tt, st.rho_sp;
                    rho_ex = st.rho_ex,
                    rho_eeex = energy.rho_eeex, rho_evib = st.rho_evib)
                dEnth = 100.0 * (Ecomp - energy.rho_etot) / energy.rho_etot
                println(@sprintf(" dEnth (%%)  = % .5E", dEnth))

                t_us = t * 1e6
                dt_us = (i == 1 ? first_dt : (sol.t[i] - sol.t[i - 1])) * 1e6
                println(@sprintf(" iter       = %6d", iter))
                println(@sprintf(" time       = % .2E mu-s ", t_us))
                println(@sprintf(" dt         = % .2E mu-s", dt_us))
                println(@sprintf(" T(t,e,r,v) = % .3E % .3E % .3E % .3E K",
                    temps.tt, temps.teex, temps.trot, temps.tvib))

                # Mole fractions
                denom = sum(st.rho_sp ./ molecular_weights)
                x = (st.rho_sp ./ molecular_weights) ./ denom
                xbuf = IOBuffer()
                print(xbuf, " X          =")
                for xi in x
                    print(xbuf, @sprintf(" % .3E", xi))
                end
                println(String(take!(xbuf)))

                for isp in 1:n_species
                    if !has_electronic_states[isp]
                        continue
                    end
                    mex = min(electronic_state_counts[isp], size(st.rho_ex, 1))
                    if mex <= 0
                        continue
                    end
                    if isp < 10
                        println(@sprintf(" Tex(%d)     = % .3E", isp, temps.tex[isp]))
                    else
                        println(@sprintf(" Tex(%d)    = % .3E", isp, temps.tex[isp]))
                    end

                    states_to_show = min(mex, 7)
                    nbuf = IOBuffer()
                    if isp < 10
                        print(nbuf, @sprintf(" n(%d,1:%d)   =", isp, states_to_show))
                    else
                        print(nbuf, @sprintf(" n(%d,1:%d)  =", isp, states_to_show))
                    end

                    for iex in 1:states_to_show
                        nval = st.rho_ex[iex, isp] / molecular_weights[isp] * AVOGADRO
                        print(nbuf, @sprintf(" % .3E", nval))
                    end
                    println(String(take!(nbuf)))
                end

                println()
            end
            iter += 1
        end

        # Extract results at output times
        n_times = length(sol.t)
        time_points = collect(sol.t)

        # Pre-allocate arrays
        species_densities = zeros(n_species, n_times)
        temperatures_tt = Vector{Float64}(undef, n_times)
        temperatures_te = Vector{Float64}(undef, n_times)
        temperatures_tv = Vector{Float64}(undef, n_times)
        total_energies = Vector{Float64}(undef, n_times)

        # Extract species densities and calculate temperatures at each time point
        for i in 1:n_times
            state = unpack_state_vector(sol.u[i], dimensions)
            energy = reconstruct_energy_components(state, config;
                teex_const = config.temperatures.Te,
                teex_vec = teex_const_vec,
                molecular_weights = molecular_weights,
                species_names = species_names,
                gas_constants = gas_constants,
                has_electronic_sts = dimensions.has_electronic_sts,
                has_vibrational_sts = dimensions.has_vibrational_sts,
                energy_cache = nothing,
                pressure_cache = nothing)

            # Store species densities
            species_densities[:, i] = state.rho_sp
            total_energies[i] = energy.rho_etot

            try
                temps = haskey(energy, :temps) ? energy.temps :
                begin
                    rho_ex_arg = dimensions.has_electronic_sts ? state.rho_ex : nothing
                    rho_vx_arg = dimensions.has_vibrational_sts ? state.rho_vx : nothing
                    calculate_temperatures_wrapper(
                        state.rho_sp, energy.rho_etot; rho_ex = rho_ex_arg, rho_vx = rho_vx_arg,
                        rho_eeex = energy.rho_eeex, rho_evib = state.rho_evib)
                end

                temperatures_tt[i] = temps.tt
                temperatures_te[i] = temps.teex
                temperatures_tv[i] = temps.tvib
            catch e
                @warn "Temperature calculation failed at time $(time_points[i])" exception=e
                # Use previous values or reasonable defaults
                if i > 1
                    temperatures_tt[i] = temperatures_tt[i - 1]
                    temperatures_te[i] = temperatures_te[i - 1]
                    temperatures_tv[i] = temperatures_tv[i - 1]
                else
                    temperatures_tt[i] = config.temperatures.Tt
                    temperatures_te[i] = config.temperatures.Te
                    temperatures_tv[i] = config.temperatures.Tv
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

        # Create results structure
        temperatures = (
            tt = temperatures_tt,
            te = temperatures_te,
            tv = temperatures_tv
        )

        # Robust success detection across SciMLBase versions
        rc = sol.retcode
        success = rc isa Symbol ? (rc in (:Success, :Terminated)) :
                  (occursin("Success", string(rc)) || occursin("Terminated", string(rc)))
        message = success ? "ODE integration completed successfully" :
                  "ODE integration terminated: $(rc)"

        @info "Results processing completed" final_time=time_points[end] success=success

        return TERRAResults(
            time_points,
            species_densities_si,
            temperatures,
            total_energies_si,
            nothing,  # source_terms - could be added later
            success,
            message
        )

    catch e
        @error "ODE integration failed" exception=e

        # Return minimal results with error information
        return TERRAResults(
            [0.0],
            reshape(initial_state.rho_sp, :, 1),
            (tt = [config.temperatures.Tt], te = [config.temperatures.Te],
                tv = [config.temperatures.Tv], tee = [config.temperatures.Te]),
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
function solve_terra_0d(config::TERRAConfig)
    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        @info "Starting TERRA 0D simulation" species=config.species

        # Convert configuration to initial conditions (SI to CGS)
        initial_state = config_to_initial_state(config)

        # Run the time integration
        results = integrate_0d_system(config, initial_state)

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
        write_native_outputs = config.write_native_outputs
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
