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
