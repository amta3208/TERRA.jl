"""
# Fortran Wrapper Module

This module provides low-level Julia bindings to the TERRA Fortran library using `ccall`.
It handles direct interfacing with the Fortran API functions defined in `interface.f90`.
"""

# TERRA library state
const TERRA_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const LOADED_TERRA_LIB_PATH = Ref{String}("")
const TERRA_INITIALIZED = Ref{Bool}(false)
const DEBUG_WRAPPER = Ref{Bool}(false)
const TERRA_CASE_PATH = Ref{String}("")
const TERRA_OUTPUTS_OPEN = Ref{Bool}(false)

# Environment variable used to locate the shared library when not provided explicitly
const TERRA_ENV_VAR_NAME = "TERRA_LIB_PATH"

"""
$(SIGNATURES)

Resolve the TERRA shared library path from the `$(TERRA_ENV_VAR_NAME)` environment variable.

# Returns
- `String`: Validated absolute path to the TERRA shared library

# Throws
- `ErrorException`: If the environment variable is unset or points to a missing file
"""
function resolve_terra_library_path()
    path_raw = String(strip(get(ENV, TERRA_ENV_VAR_NAME, "")))
    path = expanduser(path_raw)
    if isempty(path)
        error("TERRA library path not provided. Set $(TERRA_ENV_VAR_NAME) in the environment before running TERRA.")
    elseif !isfile(path)
        # If a relative path is provided, interpret it relative to the package root.
        # This makes `TERRA_LIB_PATH=terra/source/libterra.so` work under `Pkg.test`,
        # where the working directory is a temporary environment directory.
        candidate = if isabspath(path) || isempty(PACKAGE_ROOT)
            ""
        else
            abspath(joinpath(PACKAGE_ROOT, path))
        end
        if !isempty(candidate) && isfile(candidate)
            path = candidate
        else
            error("TERRA library file not found: $path\n" *
                  (isempty(candidate) ? "" : "Also tried: $candidate\n") *
                  "Set $(TERRA_ENV_VAR_NAME) to the full path of the TERRA shared library.")
        end
    end
    return path
end

"""
$(SIGNATURES)

Set the path to the TERRA shared library and load it.

# Arguments
- `path::String`: Path to the TERRA shared library file

# Throws
- `ErrorException`: If the library cannot be loaded
"""
function load_terra_library!(path::String)
    # Close existing handle if open
    if TERRA_HANDLE[] != C_NULL
        Libdl.dlclose(TERRA_HANDLE[])
        TERRA_HANDLE[] = C_NULL
    end

    # Validate path exists
    if !isfile(path)
        error("TERRA library file not found: $path\n" *
              "Set environment variable $(TERRA_ENV_VAR_NAME) to the full path of the TERRA shared library.")
    end

    # Open new library
    try
        TERRA_HANDLE[] = Libdl.dlopen(path)
        LOADED_TERRA_LIB_PATH[] = path  # Store the path for ccall usage
        # Reset initialization state for the freshly loaded library
        try
            TERRA_INITIALIZED[] = is_api_initialized_wrapper()
        catch
            TERRA_INITIALIZED[] = false
        end
    catch e
        error("Failed to load TERRA library from $path: $(e.msg)")
    end
end

"""
$(SIGNATURES)

Load the TERRA shared library using the `$(TERRA_ENV_VAR_NAME)` environment variable.

Throws an error if the environment variable is not set or does not point to an existing file.
"""
function load_terra_library!()
    path = resolve_terra_library_path()
    return load_terra_library!(path)
end

"""
$(SIGNATURES)

Enable or disable verbose wrapper debug logging.
"""
function set_wrapper_debug!(flag::Bool)
    DEBUG_WRAPPER[] = flag
end

"""
$(SIGNATURES)

Get runtime setup flags from TERRA (for verification).
"""
function get_runtime_flags()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    ev_relax_set = ccall((:get_ev_relax_set, get_terra_lib_path()), Int32, ())
    vib_noneq = ccall((:get_vib_noneq, get_terra_lib_path()), Int32, ())
    eex_noneq = ccall((:get_eex_noneq, get_terra_lib_path()), Int32, ())
    rot_noneq = ccall((:get_rot_noneq, get_terra_lib_path()), Int32, ())
    has_elec_bfe = ccall((:get_consider_elec_bfe, get_terra_lib_path()), Int32, ())
    has_elec_bbh = ccall((:get_consider_elec_bbh, get_terra_lib_path()), Int32, ())
    has_elec_bfh = ccall((:get_consider_elec_bfh, get_terra_lib_path()), Int32, ())
    has_elec_bbe = ccall((:get_consider_elec_bbe, get_terra_lib_path()), Int32, ())

    return (
        ev_relax_set = Int(ev_relax_set),
        vib_noneq = Int(vib_noneq),
        eex_noneq = Int(eex_noneq),
        rot_noneq = Int(rot_noneq),
        consider_elec_bfe = Int(has_elec_bfe),
        consider_elec_bbh = Int(has_elec_bbh),
        consider_elec_bfh = Int(has_elec_bfh),
        consider_elec_bbe = Int(has_elec_bbe)
    )
end

"""
$(SIGNATURES)

Get the API layout for the `y`/`dy` vectors used by the Fortran `rhs_api`.

This calls the Fortran `get_api_layout()` routine and returns sizes, flags, and
species/state metadata needed to construct an ODE state in the `rhs_api` ordering on the Julia side.
"""
function get_api_layout_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    # Allocate scalar outputs
    layout_version = Ref{Int32}(0)
    mnsp_out = Ref{Int32}(0)
    mnex_out = Ref{Int32}(0)
    mmnex_out = Ref{Int32}(0)
    mnv_out = Ref{Int32}(0)
    mneq_out = Ref{Int32}(0)
    nsp_out = Ref{Int32}(0)
    nd_out = Ref{Int32}(0)
    neq_out = Ref{Int32}(0)
    esp_out = Ref{Int32}(0)
    get_electron_density_by_charge_balance_out = Ref{Int32}(0)
    eex_noneq_out = Ref{Int32}(0)
    rot_noneq_out = Ref{Int32}(0)
    vib_noneq_out = Ref{Int32}(0)
    is_isothermal_out = Ref{Int32}(0)
    is_isothermal_teex_out = Ref{Int32}(0)
    is_elec_sts_out = Ref{Int32}(0)
    is_vib_sts_out = Ref{Int32}(0)
    n_eq_vib_out = Ref{Int32}(0)
    n_eq_elec_out = Ref{Int32}(0)
    n_eq_sp_out = Ref{Int32}(0)
    n_eq_mom_out = Ref{Int32}(0)
    n_eq_energy_out = Ref{Int32}(0)

    # Allocate array outputs using library maxima
    max_species = Int(get_max_number_of_species_wrapper())
    max_atomic_electronic_states = Int(get_max_number_of_atomic_electronic_states_wrapper())
    max_molecular_electronic_states = Int(get_max_number_of_molecular_electronic_states_wrapper())

    ih_full = zeros(Int32, max_species)
    ie_full = zeros(Int32, max_species)
    ies_full = zeros(Int32, max_species)
    mex_full = zeros(Int32, max_species)
    ivs_full = zeros(Int32, max_atomic_electronic_states, max_species)
    mv_full = zeros(Int32, max_species, max_molecular_electronic_states)
    spwt_full = zeros(Float64, max_species)

    ccall((:get_api_layout, get_terra_lib_path()), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
            Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
            Ref{Int32},
            Ref{Int32}, Ref{Int32}, Ref{Int32},
            Ref{Int32}, Ref{Int32},
            Ref{Int32}, Ref{Int32},
            Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        layout_version, mnsp_out, mnex_out, mmnex_out, mnv_out, mneq_out,
        nsp_out, nd_out, neq_out, esp_out,
        get_electron_density_by_charge_balance_out,
        eex_noneq_out, rot_noneq_out, vib_noneq_out,
        is_isothermal_out, is_isothermal_teex_out,
        is_elec_sts_out, is_vib_sts_out,
        n_eq_vib_out, n_eq_elec_out, n_eq_sp_out, n_eq_mom_out, n_eq_energy_out,
        ih_full, ie_full, ies_full, mex_full, ivs_full, mv_full, spwt_full)

    if layout_version[] == 0
        error("get_api_layout() returned layout_version=0 (Fortran API not initialized). Call initialize_api_wrapper() first.")
    end

    nsp = Int(nsp_out[])
    ih = copy(@view ih_full[1:nsp])
    ie = copy(@view ie_full[1:nsp])
    ies = copy(@view ies_full[1:nsp])
    mex = copy(@view mex_full[1:nsp])
    ivs = copy(@view ivs_full[:, 1:nsp])
    mv = copy(@view mv_full[1:nsp, :])
    spwt = copy(@view spwt_full[1:nsp])

    return (
        layout_version = Int(layout_version[]),
        mnsp = Int(mnsp_out[]),
        mnex = Int(mnex_out[]),
        mmnex = Int(mmnex_out[]),
        mnv = Int(mnv_out[]),
        mneq = Int(mneq_out[]),
        nsp = nsp,
        nd = Int(nd_out[]),
        neq = Int(neq_out[]),
        esp = Int(esp_out[]),
        get_electron_density_by_charge_balance = Int(get_electron_density_by_charge_balance_out[]),
        eex_noneq = Int(eex_noneq_out[]),
        rot_noneq = Int(rot_noneq_out[]),
        vib_noneq = Int(vib_noneq_out[]),
        is_isothermal = Int(is_isothermal_out[]),
        is_isothermal_teex = Int(is_isothermal_teex_out[]),
        is_elec_sts = Int(is_elec_sts_out[]),
        is_vib_sts = Int(is_vib_sts_out[]),
        n_eq_vib = Int(n_eq_vib_out[]),
        n_eq_elec = Int(n_eq_elec_out[]),
        n_eq_sp = Int(n_eq_sp_out[]),
        n_eq_mom = Int(n_eq_mom_out[]),
        n_eq_energy = Int(n_eq_energy_out[]),
        ih = ih,
        ie = ie,
        ies = ies,
        mex = mex,
        ivs = ivs,
        mv = mv,
        spwt = spwt
    )
end

"""
$(SIGNATURES)

Get the path to the loaded TERRA library.

# Returns
- `String`: Path to the loaded library

# Throws
- `ErrorException`: If no library is loaded
"""
function get_terra_lib_path()
    if TERRA_HANDLE[] == C_NULL
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    return LOADED_TERRA_LIB_PATH[]
end

"""
$(SIGNATURES)

Check if the TERRA library is currently loaded.

# Returns
- `Bool`: True if library is loaded, false otherwise
"""
function is_terra_loaded()
    return TERRA_HANDLE[] != C_NULL
end

"""
$(SIGNATURES)

Close the TERRA library and free resources.
"""
function close_terra_library()
    if TERRA_HANDLE[] != C_NULL
        Libdl.dlclose(TERRA_HANDLE[])
        TERRA_HANDLE[] = C_NULL
        LOADED_TERRA_LIB_PATH[] = ""  # Clear the path
        TERRA_INITIALIZED[] = false
    end
end

"""
$(SIGNATURES)

Initialize the TERRA API system.

# Arguments
- `case_path::String`: Path to directory containing input/ subdirectory (default: current directory)

# Returns
- `NamedTuple`: Contains `num_species` and `num_dimensions` as determined by TERRA from input files

# Throws
- `ErrorException`: If case_path doesn't exist, input file is missing, or Fortran call fails
"""
function initialize_api_wrapper(; case_path::String = pwd())
    # Reconcile Julia/Fortran state first to avoid mismatches in tests
    try
        TERRA_INITIALIZED[] = is_api_initialized_wrapper()
    catch
        # ignore if library not loaded yet
        TERRA_INITIALIZED[] = false
    end

    # Always validate inputs and prepare filesystem, even if Fortran is already initialized.
    # This preserves input validation semantics and directory creation guarantees.
    # Validate inputs
    if !isdir(case_path)
        error("Case path does not exist: $case_path")
    end

    input_file = joinpath(case_path, "input", "prob_setup.inp")
    if !isfile(input_file)
        error("Required input file not found: $input_file")
    end

    TERRA_CASE_PATH[] = case_path
    TERRA_OUTPUTS_OPEN[] = false

    # Ensure output directory structure exists
    output_dir = joinpath(case_path, "output")
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Ensure required output subdirectories exist
    sources_dir = joinpath(output_dir, "sources")
    states_dir = joinpath(output_dir, "states")

    if !isdir(sources_dir)
        mkpath(sources_dir)
    end

    if !isdir(states_dir)
        mkpath(states_dir)
    end

    # If Fortran already initialized, skip reinitialization but still return dims
    if TERRA_INITIALIZED[]
        num_species = get_number_of_active_species_wrapper()
        num_dimensions = get_number_of_dimensions_wrapper()
        return (num_species = num_species, num_dimensions = num_dimensions)
    end

    # Store current directory and change to case path
    original_dir = pwd()

    try
        cd(case_path)

        # Create references for the Fortran call (as output parameters)
        # Fortran interface expects (num_species, num_dimensions) order
        num_species_ref = Ref{Int32}(0)
        num_dimensions_ref = Ref{Int32}(0)

        # Call Fortran function with correct parameter order
        ccall((:initialize_api, get_terra_lib_path()), Cvoid,
            (Ref{Int32}, Ref{Int32}),
            num_species_ref, num_dimensions_ref)

        # If we get here, the call succeeded
        TERRA_INITIALIZED[] = true

        # Return the values determined by Fortran
        return (num_species = num_species_ref[], num_dimensions = num_dimensions_ref[])

    catch e
        # Re-throw with more context about the failure
        if isa(e, Base.SystemError) || isa(e, ErrorException)
            error("TERRA initialization failed in directory $case_path: $(e.msg)")
        else
            rethrow(e)
        end
    finally
        # Always restore original directory
        cd(original_dir)
    end
end

"""
$(SIGNATURES)

Finalize the TERRA API system and clean up resources.
"""
function finalize_api_wrapper()
    # If library isn't loaded, there's nothing to do
    if !is_terra_loaded()
        return nothing
    end
    if TERRA_OUTPUTS_OPEN[]
        try
            close_api_output_files_wrapper()
        catch e
            @warn "Failed to close TERRA output files during finalize" exception=e
        end
    end
    # Query Fortran-side state and finalize quietly if needed
    TERRA_INITIALIZED[] = is_api_initialized_wrapper()
    if TERRA_INITIALIZED[]
        ccall((:finalize_api, get_terra_lib_path()), Cvoid, ())
        TERRA_INITIALIZED[] = false
    end
    TERRA_CASE_PATH[] = ""
    TERRA_OUTPUTS_OPEN[] = false
    return nothing
end

"""
$(SIGNATURES)

Control whether finalize_api() will call MPI_Finalize on the Fortran side.
Default is true; tests may disable it to allow reinitialization in one process.
"""
function set_api_finalize_mpi_wrapper(enable::Bool)
    if !is_terra_loaded()
        return nothing
    end
    flag = enable ? Int32(1) : Int32(0)
    ccall((:set_api_finalize_mpi, get_terra_lib_path()), Cvoid, (Int32,), flag)
    return nothing
end

"""
$(SIGNATURES)

Get the maximum number of species supported by TERRA.

# Returns
- `Int32`: Maximum number of species

# Throws
- `ErrorException`: If TERRA library is not loaded
"""
function get_max_number_of_species_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    return ccall((:get_max_number_of_species, get_terra_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get the active number of species (nsp) from TERRA.

# Returns
- `Int32`: Active species count
"""
function get_number_of_active_species_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    return ccall((:get_number_of_species, get_terra_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get the maximum number of electronic states per atomic species supported by TERRA.

# Returns
- `Int32`: Maximum number of atomic electronic states

# Throws
- `ErrorException`: If TERRA library is not loaded
"""
function get_max_number_of_atomic_electronic_states_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    return ccall(
        (:get_max_number_of_atomic_electronic_states, get_terra_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get the maximum number of electronic states per molecular species supported by TERRA.

# Returns
- `Int32`: Maximum number of molecular electronic states

# Throws
- `ErrorException`: If TERRA library is not loaded
"""
function get_max_number_of_molecular_electronic_states_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    return ccall(
        (:get_max_number_of_molecular_electronic_states, get_terra_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Get the maximum vibrational quantum number supported by TERRA.

# Returns
- `Int32`: Maximum vibrational quantum number (mnv from Fortran parameters)

# Throws
- `ErrorException`: If TERRA library is not loaded
"""
function get_max_vibrational_quantum_number_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    return ccall((:get_max_vibrational_quantum_number, get_terra_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Return whether the TERRA Fortran API reports itself initialized.

# Returns
- `Bool`: True if Fortran side is initialized, false otherwise
"""
function is_api_initialized_wrapper()
    # If library isn't loaded, treat as not initialized instead of erroring
    if !is_terra_loaded()
        return false
    end
    v = ccall((:is_api_initialized, get_terra_lib_path()), Int32, ())
    return v != 0
end

"""
$(SIGNATURES)

Get the number of spatial dimensions (`nd`) from TERRA.

# Returns
- `Int32`: Number of spatial dimensions
"""
function get_number_of_dimensions_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    return ccall((:get_number_of_dimensions, get_terra_lib_path()), Int32, ())
end

"""
$(SIGNATURES)

Report whether the current TERRA setup includes vibrational state-to-state data.
"""
function has_vibrational_sts_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    return ccall((:get_has_vibrational_sts, get_terra_lib_path()), Int32, ()) != 0
end

"""
$(SIGNATURES)

Report whether the current TERRA setup includes electronic state-to-state data.
"""
function has_electronic_sts_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    return ccall((:get_has_electronic_sts, get_terra_lib_path()), Int32, ()) != 0
end

"""
$(SIGNATURES)

Get species names from TERRA.

# Returns
- `Vector{String}`: Array of species names

# Throws
- `ErrorException`: If TERRA library is not loaded
"""
function get_species_names_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    # Query dimensions from Fortran
    max_species = get_max_number_of_species_wrapper()
    active_species = get_number_of_active_species_wrapper()
    name_length = ccall((:get_species_name_length, get_terra_lib_path()), Int32, ())

    # Allocate buffer for species names (exact size expected by Fortran)
    names_buffer = zeros(UInt8, name_length * max_species)

    # Call Fortran subroutine
    ccall((:get_species_names, get_terra_lib_path()), Cvoid, (Ptr{UInt8},), names_buffer)

    # Convert fixed-size, null-terminated blocks to Julia strings
    species_names = Vector{String}(undef, active_species)
    offset = 0
    @inbounds for i in 1:active_species
        block = @view names_buffer[(offset + 1):(offset + name_length)]
        # Find first null terminator within the block
        z = findfirst(==(0), block)
        last = z === nothing ? name_length : (z - 1)
        species_names[i] = String(block[1:last]) |> strip
        offset += name_length
    end
    return species_names
end

"""
$(SIGNATURES)

Retrieve per-species gas constants from TERRA.

Values are returned in CGS units (erg/(g·K)) to maintain consistency with the
rest of the wrapper.
"""
function get_species_gas_constants_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    max_species = get_max_number_of_species_wrapper()
    buffer = zeros(Float64, max_species)

    ccall((:get_species_gas_constants, get_terra_lib_path()), Cvoid,
        (Ptr{Float64},), buffer)

    n_active = get_number_of_active_species_wrapper()
    if n_active <= 0
        return Float64[]
    end

    # Fortran returns J/(kg·K); convert to erg/(g·K)
    return buffer[1:n_active] .* 1.0e4
end

"""
$(SIGNATURES)

Calculate nonequilibrium source terms.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_etot::Float64`: Total energy density
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `rho_vx::Array{Float64,3}`: Vibrational state densities (optional)
- `rho_erot::Float64`: Rotational energy density (optional)
- `rho_eeex::Float64`: Electron-electronic energy density (optional)
- `rho_evib::Float64`: Vibrational energy density (optional)

# Returns
- Tuple of derivative arrays corresponding to input arrays

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
"""
function calculate_sources_wrapper(rho_sp::Vector{Float64},
        rho_etot::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        rho_u::Union{Float64, Nothing} = nothing,
        rho_v::Union{Float64, Nothing} = nothing,
        rho_w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    # Validate required optional energies/velocities based on runtime flags
    flags = nothing
    nd = 0
    try
        flags = get_runtime_flags()
        nd = Int(get_number_of_dimensions_wrapper())
    catch
        flags = nothing
        nd = 0
    end
    if flags !== nothing
        if flags.vib_noneq == 1 && rho_evib === nothing
            throw(ArgumentError("rho_evib must be provided when vib_noneq=1"))
        end
        if flags.eex_noneq == 1 && rho_eeex === nothing
            throw(ArgumentError("rho_eeex must be provided when eex_noneq=1"))
        end
        if flags.rot_noneq == 1 && rho_erot === nothing
            throw(ArgumentError("rho_erot must be provided when rot_noneq=1"))
        end
    end
    if nd >= 1 && rho_u === nothing
        throw(ArgumentError("rho_u must be provided when nd >= 1"))
    end
    if nd >= 2 && rho_v === nothing
        throw(ArgumentError("rho_v must be provided when nd >= 2"))
    end
    if nd >= 3 && rho_w === nothing
        throw(ArgumentError("rho_w must be provided when nd >= 3"))
    end

    has_vib_sts = false
    has_elec_sts = false
    try
        has_vib_sts = has_vibrational_sts_wrapper()
        has_elec_sts = has_electronic_sts_wrapper()
    catch
        has_vib_sts = false
        has_elec_sts = false
    end
    if has_elec_sts && rho_ex === nothing
        throw(ArgumentError("rho_ex must be provided when electronic STS is active"))
    end
    if has_vib_sts && rho_vx === nothing
        throw(ArgumentError("rho_vx must be provided when vibrational STS is active"))
    end

    # Query maxima and build full-size buffers expected by Fortran
    max_species = get_max_number_of_species_wrapper()
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    max_vibrational_quantum_number = get_max_vibrational_quantum_number_wrapper()

    nsp = length(rho_sp)
    if nsp > max_species
        throw(ArgumentError("rho_sp length ($nsp) exceeds library maximum species ($max_species)"))
    end

    # Input species densities: pad to mnsp
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    # Optional inputs: ensure full-size shape for Fortran
    rho_ex_full = nothing
    if rho_ex !== nothing
        if size(rho_ex, 1) > max_atomic_electronic_states || size(rho_ex, 2) > max_species
            throw(ArgumentError("rho_ex size $(size(rho_ex)) exceeds library maxima ($(max_atomic_electronic_states), $(max_species))"))
        end
        rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
        m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
        m2 = min(size(rho_ex, 2), max_species)
        @inbounds (rho_ex_full::Matrix{Float64})[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
    end

    rho_vx_full = nothing
    if rho_vx !== nothing
        if size(rho_vx, 1) > (max_vibrational_quantum_number + 1) ||
           size(rho_vx, 2) > max_molecular_electronic_states ||
           size(rho_vx, 3) > max_species
            throw(ArgumentError("rho_vx size $(size(rho_vx)) exceeds library maxima ($(max_vibrational_quantum_number + 1), $(max_molecular_electronic_states), $(max_species))"))
        end
        rho_vx_full = zeros(Float64, max_vibrational_quantum_number + 1,
            max_molecular_electronic_states, max_species)
        m1 = min(size(rho_vx, 1), max_vibrational_quantum_number + 1)
        m2 = min(size(rho_vx, 2), max_molecular_electronic_states)
        m3 = min(size(rho_vx, 3), max_species)
        @inbounds (rho_vx_full::Array{Float64, 3})[1:m1, 1:m2, 1:m3] .= rho_vx[
            1:m1, 1:m2, 1:m3]
    end

    # Outputs: full buffers for Fortran
    drho_sp_full = zeros(Float64, max_species)
    drho_etot = Ref{Float64}(0.0)
    drho_ex_full = rho_ex !== nothing ?
                   zeros(Float64, max_atomic_electronic_states, max_species) : nothing
    drho_vx_full = rho_vx !== nothing ?
                   zeros(Float64, max_vibrational_quantum_number + 1,
        max_molecular_electronic_states, max_species) : nothing
    # Always request scalar energy-mode derivatives to avoid positional ambiguity
    drho_erot_ref = Ref{Float64}(0.0)
    drho_eeex_ref = Ref{Float64}(0.0)
    drho_evib_ref = Ref{Float64}(0.0)

    # Call Fortran subroutine with proper optional argument handling
    try
        if DEBUG_WRAPPER[]
            flags_dbg = try
                get_runtime_flags()
            catch
                nothing
            end
            @info "WRAP_IN calculate_sources" nsp=length(rho_sp) rho_etot=rho_etot has_ex=(rho_ex !==
                                                                                           nothing) has_vx=(rho_vx !==
                                                                                                            nothing) has_u=(rho_u !==
                                                                                                                            nothing) has_v=(rho_v !==
                                                                                                                                            nothing) has_w=(rho_w !==
                                                                                                                                                            nothing) has_erot=(rho_erot !==
                                                                                                                                                                               nothing) has_eeex=(rho_eeex !==
                                                                                                                                                                                                  nothing) has_evib=(rho_evib !==
                                                                                                                                                                                                                     nothing) flags=flags_dbg
        end
        ccall((:calculate_nonequilibrium_sources, get_terra_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_sp
                Ptr{Cvoid},                                      # rho_ex (optional)
                Ptr{Cvoid},                                      # rho_vx (optional)
                Ptr{Cvoid},                                      # rho_u (optional)
                Ptr{Cvoid},                                      # rho_v (optional)
                Ptr{Cvoid},                                      # rho_w (optional)
                Ref{Float64},                                    # rho_etot
                Ptr{Cvoid},                                      # rho_erot (optional)
                Ptr{Cvoid},                                      # rho_eeex (optional)
                Ptr{Cvoid},                                      # rho_evib (optional)
                Ptr{Float64},                                    # drho_sp
                Ptr{Cvoid},                                      # drho_ex (optional)
                Ptr{Cvoid},                                      # drho_vx (optional)
                Ref{Float64},                                    # drho_etot
                Ptr{Cvoid},                                      # drho_erot (optional)
                Ptr{Cvoid},                                      # drho_eeex (optional)
                Ptr{Cvoid}),                                     # drho_evib (optional)
            rho_sp_full,
            rho_ex_full !== nothing ? (rho_ex_full::Matrix{Float64}) : C_NULL,
            rho_vx_full !== nothing ? (rho_vx_full::Array{Float64, 3}) : C_NULL,
            rho_u !== nothing ? Ref{Float64}(rho_u) : C_NULL,
            rho_v !== nothing ? Ref{Float64}(rho_v) : C_NULL,
            rho_w !== nothing ? Ref{Float64}(rho_w) : C_NULL,
            rho_etot,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,
            drho_sp_full,
            drho_ex_full !== nothing ? (drho_ex_full::Matrix{Float64}) : C_NULL,
            drho_vx_full !== nothing ? (drho_vx_full::Array{Float64, 3}) : C_NULL,
            drho_etot,
            drho_erot_ref,
            drho_eeex_ref,
            drho_evib_ref)
    catch e
        error("Failed to calculate source terms: $(e)")
    end

    # Trim outputs to active species to match solver dimensions
    drho_sp = copy(@view(drho_sp_full[1:nsp]))
    drho_ex_out = nothing
    if rho_ex !== nothing
        # shape: (max_atomic_electronic_states, nsp)
        drho_ex_out = copy(@view((drho_ex_full::Matrix{Float64})[:, 1:nsp]))
    end
    drho_vx_out = nothing
    if rho_vx !== nothing
        # shape: (max_vibrational_quantum_number+1, max_molecular_electronic_states, nsp)
        drho_vx_out = copy(@view((drho_vx_full::Array{Float64, 3})[:, :, 1:nsp]))
    end

    if DEBUG_WRAPPER[]
        @info "WRAP_OUT calculate_sources" drho_etot=drho_etot[] drho_erot=drho_erot_ref[] drho_eeex=drho_eeex_ref[] drho_evib=drho_evib_ref[] nsp=nsp has_drho_ex=(drho_ex_full !==
                                                                                                                                                                    nothing) has_drho_vx=(drho_vx_full !==
                                                                                                                                                                                          nothing)
    end

    return (drho_sp = drho_sp,
        drho_etot = drho_etot[],
        drho_ex = drho_ex_out,
        drho_vx = drho_vx_out,
        drho_erot = drho_erot_ref[],
        drho_eeex = drho_eeex_ref[],
        drho_evib = drho_evib_ref[])
end

"""
$(SIGNATURES)

Calculate temperatures from thermodynamic state.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities
- `rho_etot::Float64`: Total energy density
- Additional optional energy components

# Returns
- Named tuple with temperatures (tt, trot, teex, tvib, tex, tvx)

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
"""
function calculate_temperatures_wrapper(rho_sp::Vector{Float64},
        rho_etot::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        rho_u::Union{Float64, Nothing} = nothing,
        rho_v::Union{Float64, Nothing} = nothing,
        rho_w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    # Input validation
    if any(rho_sp .< 0)
        error("Negative species densities found in temperature calculation")
    end
    if isnan(rho_etot) || isinf(rho_etot)
        error("Invalid total energy value: $rho_etot")
    end

    # Prepare output variables
    tt = Ref{Float64}(0.0)
    trot = Ref{Float64}(0.0)
    teex = Ref{Float64}(0.0)
    tvib = Ref{Float64}(0.0)

    # Prepare arrays for species-specific temperatures and full-size buffers
    max_species = get_max_number_of_species_wrapper()
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    max_vibrational_quantum_number = get_max_vibrational_quantum_number_wrapper()

    nsp = length(rho_sp)
    if nsp > max_species
        throw(ArgumentError("rho_sp length ($nsp) exceeds library maximum species ($max_species)"))
    end
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    has_vib_sts = false
    has_elec_sts = false
    try
        has_vib_sts = has_vibrational_sts_wrapper()
        has_elec_sts = has_electronic_sts_wrapper()
    catch
        has_vib_sts = false
        has_elec_sts = false
    end
    if has_elec_sts && rho_ex === nothing
        throw(ArgumentError("rho_ex must be provided when electronic STS is active"))
    end
    if has_vib_sts && rho_vx === nothing
        throw(ArgumentError("rho_vx must be provided when vibrational STS is active"))
    end

    # Optional inputs resized to full maxima
    rho_ex_full = nothing
    if rho_ex !== nothing
        if size(rho_ex, 1) > max_atomic_electronic_states || size(rho_ex, 2) > max_species
            throw(ArgumentError("rho_ex size $(size(rho_ex)) exceeds library maxima ($(max_atomic_electronic_states), $(max_species))"))
        end
        rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
        m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
        m2 = min(size(rho_ex, 2), max_species)
        @inbounds (rho_ex_full::Matrix{Float64})[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
    end
    rho_vx_full = nothing
    if rho_vx !== nothing
        if size(rho_vx, 1) > (max_vibrational_quantum_number + 1) ||
           size(rho_vx, 2) > max_molecular_electronic_states ||
           size(rho_vx, 3) > max_species
            throw(ArgumentError("rho_vx size $(size(rho_vx)) exceeds library maxima ($(max_vibrational_quantum_number + 1), $(max_molecular_electronic_states), $(max_species))"))
        end
        rho_vx_full = zeros(Float64, max_vibrational_quantum_number + 1,
            max_molecular_electronic_states, max_species)
        m1 = min(size(rho_vx, 1), max_vibrational_quantum_number + 1)
        m2 = min(size(rho_vx, 2), max_molecular_electronic_states)
        m3 = min(size(rho_vx, 3), max_species)
        @inbounds (rho_vx_full::Array{Float64, 3})[1:m1, 1:m2, 1:m3] .= rho_vx[
            1:m1, 1:m2, 1:m3]
    end

    tex = zeros(Float64, max_species)
    tvx = zeros(Float64, max_molecular_electronic_states, max_species)

    # Guard against missing required optional energies based on runtime flags
    # to avoid Fortran-side fatal errors. Only catch errors from the flag query
    # itself; do not swallow our validation exceptions.
    flags = nothing
    try
        flags = get_runtime_flags()
    catch
        flags = nothing
    end
    if flags !== nothing
        if flags.eex_noneq == 1 && rho_eeex === nothing
            throw(ArgumentError("rho_eeex must be provided when eex_noneq=1"))
        end
        if flags.rot_noneq == 1 && rho_erot === nothing
            throw(ArgumentError("rho_erot must be provided when rot_noneq=1"))
        end
        if flags.vib_noneq == 1 && rho_evib === nothing
            throw(ArgumentError("rho_evib must be provided when vib_noneq=1"))
        end
    end

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_temperatures, get_terra_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_sp
                Ptr{Cvoid},                                      # rho_ex (optional)
                Ptr{Cvoid},                                      # rho_vx (optional)
                Ptr{Cvoid},                                      # rho_u (optional)
                Ptr{Cvoid},                                      # rho_v (optional)
                Ptr{Cvoid},                                      # rho_w (optional)
                Ref{Float64},                                    # rho_etot
                Ptr{Cvoid},                                      # rho_erot (optional)
                Ptr{Cvoid},                                      # rho_eeex (optional)
                Ptr{Cvoid},                                      # rho_evib (optional)
                Ref{Float64},                                    # tt
                Ref{Float64},                                    # trot
                Ref{Float64},                                    # teex
                Ref{Float64},                                    # tvib
                Ptr{Float64},                                    # tex
                Ptr{Float64}),                                   # tvx
            rho_sp_full,
            rho_ex_full !== nothing ? (rho_ex_full::Matrix{Float64}) : C_NULL,
            rho_vx_full !== nothing ? (rho_vx_full::Array{Float64, 3}) : C_NULL,
            rho_u !== nothing ? Ref{Float64}(rho_u) : C_NULL,
            rho_v !== nothing ? Ref{Float64}(rho_v) : C_NULL,
            rho_w !== nothing ? Ref{Float64}(rho_w) : C_NULL,
            rho_etot,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL,
            tt,
            trot,
            teex,
            tvib,
            tex,
            tvx)
    catch e
        error("Failed to calculate temperatures: $(e)")
    end

    return (tt = tt[], trot = trot[], teex = teex[], tvib = tvib[],
        tex = tex, tvx = tvx)
end

"""
$(SIGNATURES)

Calculate total energy from state variables.

# Arguments
- `tt::Float64`: Translational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `rho_vx::Array{Float64,3}`: Vibrational state densities (optional)
- `u::Float64`: x-velocity component (optional)
- `v::Float64`: y-velocity component (optional)
- `w::Float64`: z-velocity component (optional)
- `rho_erot::Float64`: Rotational energy density (optional)
- `rho_eeex::Float64`: Electron-electronic energy density (optional)
- `rho_evib::Float64`: Vibrational energy density (optional)

# Returns
- `Float64`: Total energy density

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
"""
function calculate_total_energy_wrapper(tt::Float64,
        rho_sp::Vector{Float64};
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_vx::Union{Array{Float64, 3}, Nothing} = nothing,
        u::Union{Float64, Nothing} = nothing,
        v::Union{Float64, Nothing} = nothing,
        w::Union{Float64, Nothing} = nothing,
        rho_erot::Union{Float64, Nothing} = nothing,
        rho_eeex::Union{Float64, Nothing} = nothing,
        rho_evib::Union{Float64, Nothing} = nothing)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    rho_etot = Ref{Float64}(0.0)

    # Validate required optional energies/velocities based on runtime flags
    # and number of spatial dimensions to avoid Fortran-side aborts.
    flags = nothing
    nd = 0
    try
        flags = get_runtime_flags()
        nd = Int(get_number_of_dimensions_wrapper())
    catch
        flags = nothing
        nd = 0
    end
    if flags !== nothing
        if flags.vib_noneq == 1 && rho_evib === nothing
            throw(ArgumentError("rho_evib must be provided when vib_noneq=1"))
        end
        if flags.eex_noneq == 1 && rho_eeex === nothing
            throw(ArgumentError("rho_eeex must be provided when eex_noneq=1"))
        end
        if flags.rot_noneq == 1 && rho_erot === nothing
            throw(ArgumentError("rho_erot must be provided when rot_noneq=1"))
        end
    end
    if nd >= 1 && u === nothing
        throw(ArgumentError("u must be provided when nd >= 1"))
    end
    if nd >= 2 && v === nothing
        throw(ArgumentError("v must be provided when nd >= 2"))
    end
    if nd >= 3 && w === nothing
        throw(ArgumentError("w must be provided when nd >= 3"))
    end
    # Fortran expects full-size arrays using maxima
    max_species = get_max_number_of_species_wrapper()
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    max_vibrational_quantum_number = get_max_vibrational_quantum_number_wrapper()

    nsp = length(rho_sp)
    if nsp > max_species
        throw(ArgumentError("rho_sp length ($nsp) exceeds library maximum species ($max_species)"))
    end
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    rho_ex_full = nothing
    if rho_ex !== nothing
        if size(rho_ex, 1) > max_atomic_electronic_states || size(rho_ex, 2) > max_species
            throw(ArgumentError("rho_ex size $(size(rho_ex)) exceeds library maxima ($(max_atomic_electronic_states), $(max_species))"))
        end
        rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
        m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
        m2 = min(size(rho_ex, 2), max_species)
        @inbounds (rho_ex_full::Matrix{Float64})[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
    end
    rho_vx_full = nothing
    if rho_vx !== nothing
        if size(rho_vx, 1) > (max_vibrational_quantum_number + 1) ||
           size(rho_vx, 2) > max_molecular_electronic_states ||
           size(rho_vx, 3) > max_species
            throw(ArgumentError("rho_vx size $(size(rho_vx)) exceeds library maxima ($(max_vibrational_quantum_number + 1), $(max_molecular_electronic_states), $(max_species))"))
        end
        rho_vx_full = zeros(Float64, max_vibrational_quantum_number + 1,
            max_molecular_electronic_states, max_species)
        m1 = min(size(rho_vx, 1), max_vibrational_quantum_number + 1)
        m2 = min(size(rho_vx, 2), max_molecular_electronic_states)
        m3 = min(size(rho_vx, 3), max_species)
        @inbounds (rho_vx_full::Array{Float64, 3})[1:m1, 1:m2, 1:m3] .= rho_vx[
            1:m1, 1:m2, 1:m3]
    end

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_total_energy, get_terra_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_etot (output)
                Ref{Float64},                                    # tt
                Ptr{Float64},                                    # rho_sp
                Ptr{Cvoid},                                      # rho_ex (optional)
                Ptr{Cvoid},                                      # rho_vx (optional)
                Ptr{Cvoid},                                      # u (optional)
                Ptr{Cvoid},                                      # v (optional)
                Ptr{Cvoid},                                      # w (optional)
                Ptr{Cvoid},                                      # rho_erot (optional)
                Ptr{Cvoid},                                      # rho_eeex (optional)
                Ptr{Cvoid}),                                     # rho_evib (optional)
            rho_etot,
            tt,
            rho_sp_full,
            rho_ex_full !== nothing ? (rho_ex_full::Matrix{Float64}) : C_NULL,
            rho_vx_full !== nothing ? (rho_vx_full::Array{Float64, 3}) : C_NULL,
            u !== nothing ? Ref{Float64}(u) : C_NULL,
            v !== nothing ? Ref{Float64}(v) : C_NULL,
            w !== nothing ? Ref{Float64}(w) : C_NULL,
            rho_erot !== nothing ? Ref{Float64}(rho_erot) : C_NULL,
            rho_eeex !== nothing ? Ref{Float64}(rho_eeex) : C_NULL,
            rho_evib !== nothing ? Ref{Float64}(rho_evib) : C_NULL)
    catch e
        error("Failed to calculate total energy: $(e)")
    end

    return rho_etot[]
end

"""
$(SIGNATURES)

Calculate vibrational energy from state variables.

# Arguments
- `tvib::Float64`: Vibrational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `rho_ex::Matrix{Float64}`: Electronic state densities (optional)
- `tex::Vector{Float64}`: Electronic temperatures per species (optional)
- `teex::Float64`: Electron-electronic temperature (optional)

# Returns
- `Float64`: Vibrational energy density

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
"""
function calculate_vibrational_energy_wrapper(tvib::Float64,
        rho_sp::Vector{Float64};
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        tex::Union{Vector{Float64}, Nothing} = nothing,
        teex::Union{Float64, Nothing} = nothing)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    if tvib <= 0.0
        error("Invalid vibrational temperature: $tvib")
    end

    rho_evib = Ref{Float64}(0.0)
    # Fortran expects full-size arrays using maxima
    max_species = get_max_number_of_species_wrapper()
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    nsp = length(rho_sp)
    if nsp > max_species
        throw(ArgumentError("rho_sp length ($nsp) exceeds library maximum species ($max_species)"))
    end
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    rho_ex_full = nothing
    if rho_ex !== nothing
        # Here MNEX refers to molecular electronic states in vibrational energy context
        rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
        m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
        m2 = min(size(rho_ex, 2), max_species)
        @inbounds (rho_ex_full::Matrix{Float64})[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
    end

    tex_full = nothing
    if tex !== nothing
        tex_full = zeros(Float64, max_species)
        m = min(length(tex), max_species)
        @inbounds (tex_full::Vector{Float64})[1:m] .= tex[1:m]
    end

    # Call Fortran subroutine with proper optional argument handling
    try
        ccall((:calculate_vibrational_energy, get_terra_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_evib (output)
                Ref{Float64},                                    # tvib
                Ptr{Float64},                                    # rho_sp
                Ptr{Cvoid},                                      # rho_ex (optional)
                Ptr{Cvoid},                                      # tex (optional)
                Ptr{Cvoid}),                                     # teex (optional)
            rho_evib,
            tvib,
            rho_sp_full,
            rho_ex_full !== nothing ? (rho_ex_full::Matrix{Float64}) : C_NULL,
            tex_full !== nothing ? (tex_full::Vector{Float64}) : C_NULL,
            teex !== nothing ? Ref{Float64}(teex) : C_NULL)
    catch e
        error("Failed to calculate vibrational energy: $(e)")
    end

    return rho_evib[]
end

"""
$(SIGNATURES)

Calculate vibrational temperature from vibrational energy density and species densities.

# Arguments
- `rho_evib::Float64`: Vibrational energy density (erg/cm^3)
- `rho_sp::Vector{Float64}`: Species mass densities (g/cm^3)
- `rho_ex::Union{Matrix{Float64}, Nothing}`: Optional electronic state densities (layout: mnex × nsp)
- `tex::Union{Vector{Float64}, Nothing}`: Optional species electronic temperatures (K)

# Returns
- `Float64`: Vibrational temperature (K)

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
- `ArgumentError`: If inputs are invalid or dimensions exceed library maxima
"""
function calculate_vibrational_temperature_wrapper(rho_evib::Float64,
        rho_sp::Vector{Float64};
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        tex::Union{Vector{Float64}, Nothing} = nothing)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    if length(rho_sp) == 0
        throw(ArgumentError("Species density array cannot be empty"))
    end
    if any(!isfinite, rho_sp)
        throw(ArgumentError("Species densities must be finite"))
    end
    if !isfinite(rho_evib)
        throw(ArgumentError("Vibrational energy density must be finite"))
    end

    # Prepare full-size inputs for Fortran
    max_species = get_max_number_of_species_wrapper()
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()

    nsp = length(rho_sp)
    if nsp > max_species
        throw(ArgumentError("rho_sp length ($nsp) exceeds library maximum species ($max_species)"))
    end
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    rho_ex_full = nothing
    if rho_ex !== nothing
        if size(rho_ex, 1) > max_atomic_electronic_states || size(rho_ex, 2) > max_species
            throw(ArgumentError("rho_ex size $(size(rho_ex)) exceeds library maxima ($(max_atomic_electronic_states), $(max_species))"))
        end
        rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
        m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
        m2 = min(size(rho_ex, 2), max_species)
        @inbounds (rho_ex_full::Matrix{Float64})[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
    end

    tex_full = nothing
    if tex !== nothing
        if length(tex) > max_species
            throw(ArgumentError("tex length ($(length(tex))) exceeds library maximum species ($max_species)"))
        end
        tex_full = zeros(Float64, max_species)
        @inbounds tex_full[1:length(tex)] .= tex
    end

    # Output
    tvib_ref = Ref{Float64}(0.0)

    # Call Fortran subroutine
    try
        ccall((:calculate_vibrational_temperature, get_terra_lib_path()), Cvoid,
            (Ref{Float64},                                    # tvib (output)
                Ref{Float64},                                    # rho_evib
                Ptr{Float64},                                    # rho_sp
                Ptr{Cvoid},                                      # rho_ex (optional)
                Ptr{Cvoid}),                                     # tex (optional)
            tvib_ref,
            rho_evib,
            rho_sp_full,
            rho_ex_full !== nothing ? (rho_ex_full::Matrix{Float64}) : C_NULL,
            tex_full !== nothing ? (tex_full::Vector{Float64}) : C_NULL)
    catch e
        error("Failed to calculate vibrational temperature: $(e)")
    end

    return tvib_ref[]
end

"""
$(SIGNATURES)

Calculate electron-electronic energy from state variables.

# Arguments
- `teex::Float64`: Electron-electronic temperature (K)
- `tvib::Float64`: Vibrational temperature (K)
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)

# Returns
- `Float64`: Electron-electronic energy density

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
- `ArgumentError`: If teex ≤ 0, tvib ≤ 0, or arrays are invalid
"""
function calculate_electron_electronic_energy_wrapper(teex::Float64,
        tvib::Float64, rho_sp::Vector{Float64})
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    if teex <= 0.0
        throw(ArgumentError("Electron-electronic temperature must be positive, got: $teex"))
    end

    if tvib <= 0.0
        throw(ArgumentError("Vibrational temperature must be positive, got: $tvib"))
    end

    if length(rho_sp) == 0
        throw(ArgumentError("Species density array cannot be empty"))
    end

    rho_eeex = Ref{Float64}(0.0)

    # Fortran expects rho_sp of length mnsp
    max_species = get_max_number_of_species_wrapper()
    nsp = length(rho_sp)
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    # Call Fortran subroutine
    try
        ccall((:calculate_electron_electronic_energy, get_terra_lib_path()), Cvoid,
            (Ref{Float64},                                    # rho_eeex (output)
                Ref{Float64},                                    # teex
                Ref{Float64},                                    # tvib
                Ptr{Float64}),                                   # rho_sp
            rho_eeex,
            teex,
            tvib,
            rho_sp_full)
    catch e
        error("Failed to calculate electron-electronic energy: $(e)")
    end

    return rho_eeex[]
end

"""
$(SIGNATURES)

Set electronic state densities to Boltzmann distribution.

# Arguments
- `rho_sp::Vector{Float64}`: Species densities (mass/volume)
- `tex::Float64`: Electronic temperature (K)
- `trot::Float64`: Rotational temperature (K)
- `tvib::Float64`: Vibrational temperature (K)

# Returns
- `Matrix{Float64}`: Electronic state densities in Boltzmann distribution

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
- `ArgumentError`: If temperatures ≤ 0 or arrays are invalid
"""
function set_electronic_boltzmann_wrapper(rho_sp::Vector{Float64},
        tex::Float64, trot::Float64, tvib::Float64)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end

    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end

    if tex <= 0.0 || trot <= 0.0 || tvib <= 0.0
        throw(ArgumentError("All temperatures must be positive. Got: tex=$tex, trot=$trot, tvib=$tvib"))
    end

    if length(rho_sp) == 0
        throw(ArgumentError("Species density array cannot be empty"))
    end

    # Get dimensions for electronic state array and pad rho_sp to mnsp
    max_species = get_max_number_of_species_wrapper()
    if length(rho_sp) > max_species
        throw(ArgumentError("rho_sp length ($(length(rho_sp))) exceeds library maximum species ($max_species)"))
    end
    # Get max atomic electronic states for rho_ex array
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()

    # Prepare output array
    rho_ex = zeros(Float64, max_atomic_electronic_states, max_species)

    # Fortran expects rho_sp of length mnsp
    nsp = length(rho_sp)
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    # Call Fortran subroutine
    try
        ccall((:set_electronic_boltzmann, get_terra_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_ex (output)
                Ptr{Float64},                                    # rho_sp
                Ref{Float64},                                    # tex
                Ref{Float64},                                    # trot
                Ref{Float64}),                                  # tvib
            rho_ex,
            rho_sp_full,
            tex,
            trot,
            tvib)
    catch e
        error("Failed to set electronic Boltzmann distribution: $(e)")
    end

    return rho_ex
end

"""
$(SIGNATURES)

Set vibrational state densities to a Boltzmann distribution given rho_ex.

# Arguments
- `rho_ex::Matrix{Float64}`: Electronic state densities (mass/volume)
- `tex::Float64`: Electronic temperature (K)
- `trot::Float64`: Rotational temperature (K)
- `tvib::Float64`: Vibrational temperature (K)

# Returns
- `Array{Float64,3}`: Vibrational state densities (0:mnv, mmnex, mnsp)

# Throws
- `ErrorException`: If TERRA library is not loaded or not initialized
- `ArgumentError`: If temperatures ≤ 0 or arrays are invalid
"""
function set_vibrational_boltzmann_wrapper(rho_ex::Matrix{Float64},
        tex::Float64, trot::Float64, tvib::Float64)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end
    if tex <= 0.0 || trot <= 0.0 || tvib <= 0.0
        throw(ArgumentError("All temperatures must be positive. Got: tex=$tex, trot=$trot, tvib=$tvib"))
    end

    # Dimensions
    max_species = get_max_number_of_species_wrapper()
    max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_electronic_states = get_max_number_of_molecular_electronic_states_wrapper()
    max_vibrational_quantum_number = get_max_vibrational_quantum_number_wrapper()

    # Prepare input/output: ensure rho_ex has full shape [mnex, mnsp]
    if size(rho_ex, 1) > max_atomic_electronic_states || size(rho_ex, 2) > max_species
        throw(ArgumentError("rho_ex size $(size(rho_ex)) exceeds library maxima ($(max_atomic_electronic_states), $(max_species))"))
    end
    rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
    m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
    m2 = min(size(rho_ex, 2), max_species)
    @inbounds rho_ex_full[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]

    rho_vx = zeros(Float64, max_vibrational_quantum_number + 1,
        max_molecular_electronic_states, max_species)

    # Call Fortran subroutine
    try
        ccall((:set_vibrational_boltzmann, get_terra_lib_path()), Cvoid,
            (Ptr{Float64},                                    # rho_vx (output)
                Ptr{Float64},                                    # rho_ex
                Ref{Float64},                                    # tex
                Ref{Float64},                                    # trot
                Ref{Float64}),                                   # tvib
            rho_vx,
            rho_ex_full,
            tex,
            trot,
            tvib)
    catch e
        error("Failed to set vibrational Boltzmann distribution: $(e)")
    end

    return rho_vx
end

"""
$(SIGNATURES)

Open native TERRA Tecplot-style output files for the current API session.

Requires a successful call to `initialize_api_wrapper` so the case path is tracked.
"""
function open_api_output_files_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Call load_terra_library! before opening outputs.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() before opening outputs.")
    end
    if isempty(TERRA_CASE_PATH[])
        error("Case path not set. Initialize the API with a valid case_path before opening outputs.")
    end
    if TERRA_OUTPUTS_OPEN[]
        return nothing
    end

    case_path = TERRA_CASE_PATH[]
    try
        cd(case_path) do
            ccall((:open_api_output_files, get_terra_lib_path()), Cvoid, ())
        end
        TERRA_OUTPUTS_OPEN[] = true
    catch e
        error("Failed to open TERRA output files: $(e)")
    end

    return nothing
end

"""
$(SIGNATURES)

Write a snapshot of the current API state to the native TERRA output files, using
the `y` layout expected by `rhs_api`.
"""
function write_api_outputs_wrapper(istep::Integer, time::Real, dt::Real,
        state::AbstractVector{<:Real}; dist::Real = 0.0, dx::Real = 0.0)
    if !TERRA_OUTPUTS_OPEN[]
        error("TERRA output files are not open. Call open_api_output_files_wrapper() first.")
    end
    if isempty(TERRA_CASE_PATH[])
        error("Case path not set. Cannot write TERRA outputs without an initialized case path.")
    end

    state_vec = state isa Vector{Float64} ? state : Vector{Float64}(state)
    istep32 = Int32(istep)
    neq32 = Int32(length(state_vec))
    time64 = Float64(time)
    dt64 = Float64(dt)
    dist64 = Float64(dist)
    dx64 = Float64(dx)
    case_path = TERRA_CASE_PATH[]

    GC.@preserve state_vec begin
        ptr = Base.unsafe_convert(Ptr{Float64}, state_vec)
        cd(case_path) do
            ccall((:write_api_outputs, get_terra_lib_path()), Cvoid,
                (Int32, Float64, Float64, Float64, Float64, Int32, Ptr{Float64}),
                istep32, time64, dt64, dist64, dx64, neq32, ptr)
        end
    end

    return nothing
end

"""
$(SIGNATURES)

Close the native TERRA output files opened through the API wrappers.
"""
function close_api_output_files_wrapper()
    if !TERRA_OUTPUTS_OPEN[]
        return nothing
    end
    if !is_terra_loaded()
        TERRA_OUTPUTS_OPEN[] = false
        return nothing
    end
    if isempty(TERRA_CASE_PATH[])
        TERRA_OUTPUTS_OPEN[] = false
        return nothing
    end

    case_path = TERRA_CASE_PATH[]
    try
        cd(case_path) do
            ccall((:close_api_output_files, get_terra_lib_path()), Cvoid, ())
        end
    catch e
        error("Failed to close TERRA output files: $(e)")
    finally
        TERRA_OUTPUTS_OPEN[] = false
    end

    return nothing
end

"""
$(SIGNATURES)

Compute `dy = rhs_api(y)` for a `y` vector using the `rhs_api` layout.

`y` and `dy` must use the ordering returned by `get_api_layout_wrapper()`.
"""
function calculate_rhs_api_wrapper!(dy::Vector{Float64}, y::Vector{Float64})
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end
    if length(dy) != length(y)
        throw(ArgumentError("dy length ($(length(dy))) must match y length ($(length(y)))."))
    end

    neq32 = Int32(length(y))
    GC.@preserve y dy begin
        y_ptr = Base.unsafe_convert(Ptr{Float64}, y)
        dy_ptr = Base.unsafe_convert(Ptr{Float64}, dy)
        ccall((:calculate_rhs_api, get_terra_lib_path()), Cvoid,
            (Int32, Ptr{Float64}, Ptr{Float64}),
            neq32, y_ptr, dy_ptr)
    end

    return nothing
end

"""
$(SIGNATURES)

Convert total enthalpy density to total energy density for 0D isothermal Teex cases.

This calls the Fortran API routine `energy_from_enthalpy_isothermal_teex_api`, which performs
the closed-form enthalpy→energy inversion used by the isothermal Teex API RHS.

# Returns
- Named tuple with `rho_etot`, `pressure`, and `tt` (all CGS units; `tt` in K).
"""
function energy_from_enthalpy_isothermal_teex_wrapper(rho_enth::Float64,
        rho_sp::Vector{Float64},
        teex_const::Float64;
        rho_ex::Union{Matrix{Float64}, Nothing} = nothing,
        rho_erot::Float64 = 0.0,
        rho_eeex::Float64 = 0.0,
        rho_evib::Float64 = 0.0)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end
    if !isfinite(rho_enth)
        throw(ArgumentError("rho_enth must be finite (got $rho_enth)."))
    end
    if any(rho_sp .< 0)
        throw(ArgumentError("rho_sp contains negative densities."))
    end
    if !isfinite(teex_const) || teex_const <= 0.0
        throw(ArgumentError("teex_const must be finite and positive (got $teex_const)."))
    end

    # Fortran expects rho_sp of length mnsp
    max_species = get_max_number_of_species_wrapper()
    nsp = length(rho_sp)
    if nsp > max_species
        throw(ArgumentError("rho_sp length ($nsp) exceeds library maximum species ($max_species)."))
    end
    rho_sp_full = zeros(Float64, max_species)
    @inbounds rho_sp_full[1:nsp] .= rho_sp

    # Optional rho_ex padded to [mnex, mnsp]
    rho_ex_full = nothing
    if rho_ex !== nothing
        max_atomic_electronic_states = get_max_number_of_atomic_electronic_states_wrapper()
        if size(rho_ex, 1) > max_atomic_electronic_states || size(rho_ex, 2) > max_species
            throw(ArgumentError("rho_ex size $(size(rho_ex)) exceeds library maxima ($(max_atomic_electronic_states), $(max_species))."))
        end
        rho_ex_full = zeros(Float64, max_atomic_electronic_states, max_species)
        m1 = min(size(rho_ex, 1), max_atomic_electronic_states)
        m2 = min(size(rho_ex, 2), max_species)
        @inbounds (rho_ex_full::Matrix{Float64})[1:m1, 1:m2] .= rho_ex[1:m1, 1:m2]
    end

    rho_etot = Ref{Float64}(0.0)
    pres = Ref{Float64}(0.0)
    tt = Ref{Float64}(0.0)

    GC.@preserve rho_sp_full rho_ex_full begin
        rho_ex_ptr = rho_ex_full === nothing ? C_NULL :
                     Ptr{Cvoid}(Base.unsafe_convert(Ptr{Float64}, rho_ex_full))

        ccall((:energy_from_enthalpy_isothermal_teex_api, get_terra_lib_path()), Cvoid,
            (Float64, Ptr{Float64}, Ptr{Cvoid}, Float64, Float64, Float64, Float64,
                Ref{Float64}, Ref{Float64}, Ref{Float64}),
            rho_enth, rho_sp_full, rho_ex_ptr,
            rho_erot, rho_eeex, rho_evib, teex_const,
            rho_etot, pres, tt)
    end

    return (rho_etot = rho_etot[], pressure = pres[], tt = tt[])
end

"""
$(SIGNATURES)

Compute `du = rhs_api(u)` for a `u` vector in isothermal Teex mode, using the `rhs_api` layout.

This calls the Fortran API routine `calculate_rhs_api_isothermal_teex`, which expects:
- Ordering from `get_api_layout_wrapper()`
- `u[idx_etot]` stores the legacy enthalpy remainder `rho_rem`
- `u[idx_eeex]` is treated as a dummy; `du[idx_eeex]` is forced to zero

Optional inputs:
- `tex`: per-species electronic temperatures passed to the Tvib inversion (length `nsp`)
"""
function calculate_rhs_api_isothermal_teex_wrapper!(du::Vector{Float64}, u::Vector{Float64},
        teex_const::Float64;
        tex = nothing)
    if !is_terra_loaded()
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() first.")
    end
    if length(du) != length(u)
        throw(ArgumentError("du length ($(length(du))) must match u length ($(length(u)))."))
    end
    if !isfinite(teex_const) || teex_const <= 0.0
        throw(ArgumentError("teex_const must be finite and positive (got $teex_const)."))
    end

    neq32 = Int32(length(u))

    GC.@preserve u du tex begin
        u_ptr = Base.unsafe_convert(Ptr{Float64}, u)
        du_ptr = Base.unsafe_convert(Ptr{Float64}, du)
        tex_ptr = tex === nothing ? C_NULL :
                  Ptr{Cvoid}(Base.unsafe_convert(Ptr{Float64}, tex))

        ccall((:calculate_rhs_api_isothermal_teex, get_terra_lib_path()), Cvoid,
            (Int32, Ptr{Float64}, Ptr{Float64}, Float64, Ptr{Cvoid}, Ptr{Cvoid}),
            neq32, u_ptr, du_ptr, teex_const, C_NULL, tex_ptr)
    end

    return nothing
end
