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
    is_isothermal = config !== nothing && _config_is_isothermal_teex(config)

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
    teex_const = hasproperty(p, :teex_const) ? p.teex_const :
                 _config_electron_temperature(config)
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
