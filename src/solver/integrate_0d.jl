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
function integrate_0d_system(config::Config, initial_state;
        residence_time::Union{Nothing, ResidenceTimeConfig} = config.numerics.residence_time,
        use_residence_time::Union{Nothing, Bool} = nothing)
    dt = config.numerics.time.dt
    tlim = config.numerics.time.duration
    solver_cfg = config.numerics.solver

    @info "Setting up ODE integration" tlim=tlim

    native_outputs_requested = config.runtime.write_native_outputs
    print_integration_output = config.runtime.print_integration_output
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
    is_isothermal = config.models.physics.is_isothermal_teex

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

    molecular_weights = get_molecular_weights(config.reactor.composition.species)
    species_names = config.reactor.composition.species
    gas_constants = initial_state.gas_constants
    teex_const_vec = fill(config.reactor.thermal.Te, n_species)

    # Preallocate RHS work buffers
    work_y = similar(u0)
    work_u = similar(u0)
    work_rho_sp = zeros(Float64, layout.nsp)
    work_rho_ex = layout.is_elec_sts ? zeros(Float64, layout.mnex, layout.nsp) : nothing

    effective_residence_time = _resolve_residence_time(residence_time, use_residence_time)
    residence_time_data = (effective_residence_time === nothing ||
                           !effective_residence_time.enabled) ? nothing :
                          _prepare_residence_time_data(
        layout, config, u0, effective_residence_time)

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
    ramp_callback = native_ramp_callback(dt;
        understep_ratio = solver_cfg.ramp_understep_ratio,
        history_steps = solver_cfg.ramp_history_steps)

    @info "ODE problem created, starting integration..."

    # Component-wise tolerances (CGS). The compact state mixes densities (g/cm^3)
    # and energies (erg/cm^3). A scalar abstol tends to make the energy equation
    # far more restrictive than the density equations. For isothermal Teex runs,
    # loosen the energy absolute tolerance based on a target temperature
    # resolution:
    #   abstol_E ≈ (3/2) kB n_total ΔT  [erg/cm^3]
    local_reltol = is_isothermal ? max(1e-7, solver_cfg.reltol) : solver_cfg.reltol
    local_abstol_density = solver_cfg.abstol_density
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
            saveat = range(0.0, tlim; length = solver_cfg.saveat_count),
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
        if config.runtime.unit_system == :SI
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
            (tt = [config.reactor.thermal.Tt], te = [config.reactor.thermal.Te],
                tv = [config.reactor.thermal.Tv],
                tee = [config.reactor.thermal.Te]),
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

function integrate_0d_system(config::TERRAConfig, initial_state;
        residence_time::Union{Nothing, ResidenceTimeConfig} = nothing,
        use_residence_time::Union{Nothing, Bool} = nothing)
    return integrate_0d_system(
        to_config(config),
        initial_state;
        residence_time = residence_time,
        use_residence_time = use_residence_time)
end
