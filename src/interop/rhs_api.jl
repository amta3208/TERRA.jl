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
