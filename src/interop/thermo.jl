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
