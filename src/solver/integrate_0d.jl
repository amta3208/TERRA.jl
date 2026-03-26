"""
$(SIGNATURES)

Integrate the 0D system over time using the Fortran `rhs_api` layout.

This uses `terra_ode_system!` and constructs a `y` vector that matches the
ordering returned by `get_api_layout()`.

Notes:
- Vibrational STS is not yet supported in this wrapper (vibrational *mode* energy is supported).
- For isothermal Teex cases, the stored energy slot in `u` is `rho_rem` (legacy semantics), and
  the working `y` passed to Fortran is reconstructed each call.
- Wrapper-managed additive source terms are prepared from `config.sources` or an explicit
  `sources` override.
"""
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

const STANDALONE_0D_BANNER = "\n" * "="^12 * " TERRA 0D Simulation " * "="^12
const STANDALONE_0D_FOOTER = "="^(length(STANDALONE_0D_BANNER) - 1)

function _integration_presentation_emits_banner(presentation::Symbol)
    if presentation == :standalone_0d
        return true
    elseif presentation == :chain_segment
        return false
    end
    throw(ArgumentError("Unsupported integration presentation: :$presentation"))
end

function _integration_start_message(presentation::Symbol)
    if _integration_presentation_emits_banner(presentation)
        return string(STANDALONE_0D_BANNER, "\n", "starting ODE integration...")
    end
    return "starting ODE integration..."
end

function _integration_completion_message(presentation::Symbol, message::AbstractString)
    if _integration_presentation_emits_banner(presentation)
        return string(message, "\n", STANDALONE_0D_FOOTER)
    end
    return String(message)
end

function _integration_completion_console_visibility(presentation::Symbol, success::Bool)
    if presentation == :standalone_0d
        return :minimal
    elseif presentation == :chain_segment
        return success ? :verbose : :minimal
    end
    throw(ArgumentError("Unsupported integration presentation: :$presentation"))
end

function integration_progress_condition(u, t, integrator)
    reporter = integrator.p.progress_reporter
    reporter === nothing && return false
    reporter.next_fraction <= 1.0 || return false
    reporter.tlim > 0.0 || return false

    fraction = Float64(t) / reporter.tlim
    return fraction + 1e-12 >= reporter.next_fraction
end

function integration_progress_affect!(integrator)
    reporter = integrator.p.progress_reporter
    reporter === nothing && return u_modified!(integrator, false)

    _report_progress!(integrator.p.config.runtime, reporter, integrator.t)

    fraction = reporter.tlim > 0.0 ?
               clamp(Float64(integrator.t) / reporter.tlim, 0.0,
                     1.0) : 1.0
    while reporter.next_fraction <= 1.0 && fraction + 1e-12 >= reporter.next_fraction
        reporter.next_fraction += reporter.fraction_step
    end
    if reporter.next_fraction > 1.0
        reporter.next_fraction = Inf
    end

    u_modified!(integrator, false)
end

function integration_progress_callback()
    return DiscreteCallback(integration_progress_condition, integration_progress_affect!;
                            save_positions = (false, false))
end

function _extract_reactor_state_cache(species::AbstractVector{<:AbstractString},
                                      p,
                                      layout::ApiLayout,
                                      u_final::Vector{Float64},
                                      is_isothermal::Bool)
    rho_sp = zeros(Float64, layout.nsp)
    rho_ex = layout.is_elec_sts ? zeros(Float64, layout.mnex, layout.nsp) : nothing

    if is_isothermal
        y_work = similar(u_final)
        _compact_isothermal_fill_fortran_y_work!(y_work, rho_sp, rho_ex, u_final, p, layout)
    else
        _reconstruct_rho_sp_rho_ex_from_compact!(rho_sp, rho_ex, u_final, layout)
    end

    return ReactorStateCache(;
                             species = species,
                             rho_sp_cgs = rho_sp,
                             rho_ex_cgs = rho_ex,)
end

function integrate_0d_system(config::Config, initial_state;
                             sources::Union{Nothing, SourceTermsConfig} = config.sources)
    _validate_direct_wall_loss_usage(sources)
    results, _ = _integrate_0d_system(config, initial_state;
                                      sources = sources,
                                      presentation = :standalone_0d)
    return results
end

function _integrate_0d_system(config::Config, initial_state;
                              sources::Union{Nothing, SourceTermsConfig} = config.sources,
                              wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                              inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing,
                              presentation::Symbol = :standalone_0d)
    runtime = config.runtime
    dt = config.numerics.time.dt
    tlim = config.numerics.time.duration
    solver_cfg = config.numerics.solver

    _log_run_event(runtime, :info, "Preparing ODE integration";
                   console = :never,
                   :tlim => tlim,
                   :dt => dt,
                   :saveat_count => solver_cfg.saveat_count)

    native_outputs_requested = runtime.write_native_state_files
    integration_detail_requested = runtime.logging.integration_detail_mode != :off
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
    if integration_detail_requested && layout.is_elec_sts
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

    _log_run_event(runtime, :info, "Initial state vector prepared";
                   console = :never,
                   :length_u0 => length(u0),
                   :neq => layout.neq,
                   :n_species => n_species)

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

    prepared_sources = _prepare_source_terms_data(layout, config, u0, sources;
                                                  wall_inputs = wall_inputs,
                                                  inlet_state_cache = inlet_state_cache)
    result_metadata = _source_terms_result_metadata(prepared_sources)

    progress_reporter = _progress_reporter(runtime, tlim)

    p = (layout = layout,
         config = config,
         molecular_weights = molecular_weights,
         species = species_names,
         gas_constants = gas_constants,
         teex_const = initial_state.teex_const,
         teex_const_vec = teex_const_vec,
         work_u = work_u,
         sources = prepared_sources,
         progress_reporter = progress_reporter)

    prob = ODEProblem(terra_ode_system!, u0, tspan, p)
    ramp_callback = native_ramp_callback(dt;
                                         understep_ratio = solver_cfg.ramp_understep_ratio,
                                         history_steps = solver_cfg.ramp_history_steps)
    callback = progress_reporter === nothing ? ramp_callback :
               CallbackSet(ramp_callback, integration_progress_callback())

    _log_run_event(runtime, :info, _integration_start_message(presentation);
                   console = :minimal,
                   :reltol => solver_cfg.reltol,
                   :abstol_density => solver_cfg.abstol_density)

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
                    callback = callback,
                    isoutofdomain = (u, _p, _t) -> (!_compact_nonnegative_ok(u, layout) ||
                                                    any(!isfinite, u)),
                    reltol = local_reltol,
                    abstol = (local_abstol_vec === nothing ? local_abstol_density :
                              local_abstol_vec),
                    saveat = range(0.0, tlim; length = solver_cfg.saveat_count),
                    save_everystep = false)

        # Native output snapshots (optional)
        if outputs_opened
            first_dt = length(sol.t) >= 2 ? (sol.t[2] - sol.t[1]) : dt
            for (i, t) in enumerate(sol.t)
                local_dt = i == 1 ? first_dt : (t - sol.t[i - 1])
                if is_isothermal
                    _compact_isothermal_fill_fortran_y_work!(work_y, work_rho_sp,
                                                             work_rho_ex, sol.u[i], p,
                                                             layout)
                    write_api_outputs_wrapper(i - 1, t, local_dt, work_y; dist = 0.0,
                                              dx = 0.0)
                else
                    write_api_outputs_wrapper(i - 1, t, local_dt, sol.u[i]; dist = 0.0,
                                              dx = 0.0)
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
                iso = _compact_isothermal_fill_fortran_y_work!(work_y, work_rho_sp,
                                                               rho_ex_arg, ui, p, layout)
                species_densities[:, i] = work_rho_sp
                total_energies[i] = iso.rho_etot

                temps = iso.temps
                temperatures_tt[i] = temps.tt
                temperatures_te[i] = temps.teex
                temperatures_tv[i] = temps.tvib

                if integration_detail_requested
                    detail_text = _format_terra_integration_output(t, local_dt, work_rho_sp,
                                                                   molecular_weights, temps,
                                                                   iso.rho_etot;
                                                                   rho_ex = rho_ex_arg,
                                                                   rho_eeex = iso.rho_eeex,
                                                                   rho_evib = iso.rho_evib,
                                                                   has_electronic_states = has_electronic_states,
                                                                   electronic_state_counts = electronic_state_counts)
                    _emit_integration_detail(runtime, detail_text)
                end
            else
                _reconstruct_rho_sp_rho_ex_from_compact!(work_rho_sp, rho_ex_arg, ui,
                                                         layout)
                species_densities[:, i] = work_rho_sp

                rho_etot = Float64(ui[layout.idx_etot])
                total_energies[i] = rho_etot

                temps = calculate_temperatures_wrapper(work_rho_sp, rho_etot;
                                                       rho_ex = rho_ex_arg,
                                                       rho_erot = layout.idx_erot == 0 ?
                                                                  nothing :
                                                                  Float64(ui[layout.idx_erot]),
                                                       rho_eeex = layout.idx_eeex == 0 ?
                                                                  nothing :
                                                                  Float64(ui[layout.idx_eeex]),
                                                       rho_evib = layout.idx_evib == 0 ?
                                                                  nothing :
                                                                  Float64(ui[layout.idx_evib]))
                temperatures_tt[i] = temps.tt
                temperatures_te[i] = temps.teex
                temperatures_tv[i] = temps.tvib

                if integration_detail_requested
                    detail_text = _format_terra_integration_output(t, local_dt, work_rho_sp,
                                                                   molecular_weights, temps,
                                                                   rho_etot;
                                                                   rho_ex = rho_ex_arg,
                                                                   rho_eeex = layout.idx_eeex ==
                                                                              0 ?
                                                                              nothing :
                                                                              Float64(ui[layout.idx_eeex]),
                                                                   rho_evib = layout.idx_evib ==
                                                                              0 ?
                                                                              nothing :
                                                                              Float64(ui[layout.idx_evib]),
                                                                   has_electronic_states = has_electronic_states,
                                                                   electronic_state_counts = electronic_state_counts)
                    _emit_integration_detail(runtime, detail_text)
                end
            end
        end

        # Convert results back to SI units if needed
        if config.runtime.unit_system == :SI
            species_densities_si = zeros(size(species_densities))
            for i in axes(species_densities, 2)
                species_densities_si[:, i] = convert_density_cgs_to_si(species_densities[:,
                                                                                         i])
            end
            total_energies_si = [convert_energy_density_cgs_to_si(e)
                                 for e in total_energies]
        else
            species_densities_si = species_densities
            total_energies_si = total_energies
        end

        final_state_cache = _extract_reactor_state_cache(species_names, p, layout,
                                                         sol.u[end], is_isothermal)

        rc = sol.retcode
        success = rc isa Symbol ? (rc in (:Success, :Terminated)) :
                  (occursin("Success", string(rc)) || occursin("Terminated", string(rc)))
        message = success ? "success!" :
                  "ODE integration terminated: $(rc)"
        completion_message = _integration_completion_message(presentation, message)
        _log_run_event(runtime, success ? :info : :warn, completion_message;
                       console = _integration_completion_console_visibility(presentation,
                                                                            success),
                       :retcode => sol.retcode,
                       :saved_points => n_times)
        frames = Vector{ReactorFrame}(undef, n_times)
        for i in 1:n_times
            frame_temps = (tt = temperatures_tt[i], te = temperatures_te[i],
                           tv = temperatures_tv[i])
            frame_source_terms = _source_terms_frame_snapshot(sol.u[i], prepared_sources,
                                                              config.runtime.unit_system)
            frames[i] = ReactorFrame(;
                                     t = time_points[i],
                                     species_densities = species_densities_si[:, i],
                                     temperatures = frame_temps,
                                     total_energy = total_energies_si[i],
                                     source_terms = frame_source_terms)
        end

        return ReactorResult(;
                             t = time_points,
                             frames = frames,
                             source_terms = nothing,
                             success = success,
                             message = message,
                             metadata = result_metadata,), final_state_cache

    catch e
        _log_run_exception(runtime, :error,
                           _integration_completion_message(presentation,
                                                           "ODE integration failed"),
                           e;
                           console = :minimal)
        fallback_species = config.runtime.unit_system == :SI ?
                           convert_density_cgs_to_si(initial_state.rho_sp) :
                           copy(initial_state.rho_sp)
        fallback_energy = config.runtime.unit_system == :SI ?
                          convert_energy_density_cgs_to_si(initial_state.rho_energy) :
                          initial_state.rho_energy
        fallback_frame = ReactorFrame(;
                                      t = 0.0,
                                      species_densities = fallback_species,
                                      temperatures = (tt = config.reactor.thermal.Tt,
                                                      te = config.reactor.thermal.Te,
                                                      tv = config.reactor.thermal.Tv,
                                                      tee = config.reactor.thermal.Te),
                                      total_energy = fallback_energy,
                                      source_terms = nothing)
        return ReactorResult(;
                             t = [0.0],
                             frames = [fallback_frame],
                             source_terms = nothing,
                             success = false,
                             message = "ODE integration failed: $(string(e))",
                             metadata = result_metadata,), nothing
    finally
        if outputs_opened
            try
                close_api_output_files_wrapper()
            catch close_err
                _log_run_exception(runtime, :warn,
                                   "Failed to close native TERRA outputs after integration",
                                   close_err;
                                   console = :never)
            end
            outputs_opened = false
        end
    end
end
