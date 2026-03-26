const CHAIN_SEGMENT_DIRNAME = "chain_segments"

function _segment_case_path(base_case_path::AbstractString, segment_index::Integer)
    return normpath(joinpath(base_case_path, CHAIN_SEGMENT_DIRNAME,
                             @sprintf("segment_%04d", segment_index)))
end

function _failed_simulation_result(message::AbstractString)
    return ReactorResult(;
                         t = Float64[],
                         frames = ReactorFrame[],
                         source_terms = nothing,
                         success = false,
                         message = String(message))
end

function _build_segment_wall_inputs(profile::AxialChainProfile,
                                    segment_index::Integer,
                                    wall_cfg::Union{Nothing, WallLossConfig} = nothing)
    wall_values = _segment_wall_profile_values(profile, segment_index, wall_cfg)
    wall_values === nothing && return nothing
    return SegmentWallInputs(; wall_values...)
end

function _non_electron_species_names(config::Config)
    species = config.reactor.composition.species
    molecular_weights = get_molecular_weights(species)
    return String[species[i]
                  for i in eachindex(species)
                  if !_is_electron_species(species[i], molecular_weights[i])]
end

function _validate_profile_species(profile::AxialChainProfile, config::Config)
    expected_species = _non_electron_species_names(config)
    profile_species = sort!(collect(keys(profile.species_u_m_s)))
    expected_sorted = sort!(copy(expected_species))

    if profile_species != expected_sorted
        missing_species = [name for name in expected_sorted if !(name in profile_species)]
        extra_species = [name for name in profile_species if !(name in expected_sorted)]
        throw(ArgumentError("AxialChainProfile species_u_m_s keys must match the non-electron reactor species. " *
                            "Missing: $(isempty(missing_species) ? "none" : join(missing_species, ", ")); " *
                            "Extra: $(isempty(extra_species) ? "none" : join(extra_species, ", "))."))
    end
end

function _validate_profile_inlet(profile::AxialChainProfile, config::Config)
    inlet_species = profile.inlet.composition.species
    expected_species = config.reactor.composition.species
    inlet_species == expected_species ||
        throw(ArgumentError("AxialChainProfile inlet composition species must match the initialized TERRA species ordering. " *
                            "Expected: $(join(expected_species, ", ")); got: $(join(inlet_species, ", "))."))
    profile.inlet.source_compact_index == 1 ||
        throw(ArgumentError("AxialChainProfile inlet.source_compact_index must be 1 for segment-1 initialization."))
end

function _build_chain_cells(profile::AxialChainProfile,
                            compact_to_source_index::AbstractVector{<:Integer},
                            segment_end_reactors::AbstractVector{<:ReactorConfig},
                            segment_success::AbstractVector{Bool},
                            segment_messages::AbstractVector{<:AbstractString},
                            segment_results::AbstractVector{<:ReactorResult})
    n_segments = length(profile.z_m)

    length(compact_to_source_index) == n_segments ||
        throw(ArgumentError("Chain cell/source mapping length must match profile segment count."))
    length(segment_end_reactors) == n_segments ||
        throw(ArgumentError("Segment endpoint reactor count must match profile segment count."))
    length(segment_success) == n_segments ||
        throw(ArgumentError("Segment success count must match profile segment count."))
    length(segment_messages) == n_segments ||
        throw(ArgumentError("Segment message count must match profile segment count."))
    length(segment_results) == n_segments ||
        throw(ArgumentError("Segment result count must match profile segment count."))

    cells = Vector{ChainCellResult}(undef, n_segments)
    for k in 1:n_segments
        cells[k] = ChainCellResult(;
                                   compact_cell_index = k,
                                   source_cell_index = compact_to_source_index[k],
                                   z_m = profile.z_m[k],
                                   dx_m = profile.dx_m[k],
                                   te_K = profile.te_K[k],
                                   species_u_m_s = _profile_species_velocity_at_segment(profile,
                                                                                        k),
                                   reactor = segment_results[k],
                                   endpoint_reactor = segment_end_reactors[k],
                                   success = segment_success[k],
                                   message = segment_messages[k],)
    end
    return cells
end

function _build_chain_models(base_models::ModelConfig,
                             marching::AxialMarchingConfig)
    physics = base_models.physics
    chain_physics = PhysicsConfig(bbh_model = physics.bbh_model,
                                  esc_model = physics.esc_model,
                                  ar_et_model = physics.ar_et_model,
                                  eex_noneq = physics.eex_noneq,
                                  ev_relax_set = physics.ev_relax_set,
                                  et_relax_set = physics.et_relax_set,
                                  radiation_length = physics.radiation_length,
                                  get_electron_density_by_charge_balance = physics.get_electron_density_by_charge_balance,
                                  min_sts_frac = physics.min_sts_frac,
                                  is_isothermal_teex = marching.is_isothermal_teex,
                                  energy_loss_per_eii = physics.energy_loss_per_eii)
    return ModelConfig(; physics = chain_physics, processes = base_models.processes)
end

@inline function _resolve_segment_tee(inlet_reactor::ReactorConfig,
                                      te_profile::Float64,
                                      marching::AxialMarchingConfig,
                                      segment_index::Integer)
    if segment_index == 1
        return inlet_reactor.thermal.Tee
    end
    if marching.handoff_mode == :full_state
        return te_profile
    end
    return inlet_reactor.thermal.Tee
end

function _build_segment_reactor(inlet_reactor::ReactorConfig,
                                te_profile::Float64,
                                marching::AxialMarchingConfig,
                                segment_index::Integer)
    tt = marching.override_tt_K === nothing ? inlet_reactor.thermal.Tt :
         marching.override_tt_K
    tv = marching.override_tv_K === nothing ? inlet_reactor.thermal.Tv :
         marching.override_tv_K
    tee = _resolve_segment_tee(inlet_reactor, te_profile, marching, segment_index)

    thermal = ReactorThermalState(; Tt = tt, Tv = tv, Tee = tee, Te = te_profile)
    return ReactorConfig(; composition = inlet_reactor.composition, thermal = thermal)
end

function _build_segment_sources(base_config::Config,
                                profile::AxialChainProfile,
                                segment_index::Integer,
                                inlet_reactor::ReactorConfig)
    base_rt = base_config.sources.residence_time
    u_energy = base_rt === nothing ? nothing : base_rt.U_energy
    species_u = _profile_species_velocity_at_segment(profile, segment_index)

    rt = ResidenceTimeConfig(enabled = true,
                             L = profile.dx_m[segment_index],
                             U_species = species_u,
                             U_energy = u_energy,
                             inlet_reactor = inlet_reactor)
    return SourceTermsConfig(;
                             residence_time = rt,
                             wall_losses = base_config.sources.wall_losses,)
end

function _build_segment_runtime(base_runtime::RuntimeConfig,
                                segment_case_path::AbstractString)
    return with_runtime(base_runtime;
                        case_path = segment_case_path,
                        logging = _segment_logging(base_runtime, segment_case_path))
end

function _build_chain_segment_config(base_config::Config,
                                     profile::AxialChainProfile,
                                     segment_index::Integer,
                                     inlet_reactor::ReactorConfig,
                                     marching::AxialMarchingConfig)
    segment_case_path = _segment_case_path(base_config.runtime.case_path, segment_index)
    mkpath(segment_case_path)

    reactor = _build_segment_reactor(inlet_reactor, profile.te_K[segment_index], marching,
                                     segment_index)
    models = _build_chain_models(base_config.models, marching)
    sources = _build_segment_sources(base_config, profile, segment_index, inlet_reactor)
    numerics = NumericsConfig(time = base_config.numerics.time,
                              solver = base_config.numerics.solver,
                              space = base_config.numerics.space)
    runtime = _build_segment_runtime(base_config.runtime, segment_case_path)

    return Config(;
                  reactor = reactor,
                  models = models,
                  sources = sources,
                  numerics = numerics,
                  runtime = runtime,)
end

function _to_cgs_density(density::AbstractVector{<:Real}, unit_system::Symbol)
    rho = Float64.(density)
    if unit_system == :SI
        return convert_density_si_to_cgs(rho)
    end
    return rho
end

function _total_number_density_from_cgs_density(rho_cgs::Vector{Float64},
                                                molecular_weights::Vector{Float64})
    return sum(rho_cgs .* AVOGADRO ./ molecular_weights)
end

function _extract_segment_endpoint_reactor(base_config::Config,
                                           segment_result::ReactorResult)
    if !segment_result.success
        throw(ArgumentError("Cannot extract endpoint reactor from an unsuccessful segment result."))
    end
    if isempty(segment_result.t)
        throw(ArgumentError("Cannot extract endpoint reactor: segment result contains no time samples."))
    end
    endpoint_frame = segment_result.frames[end]

    species = base_config.reactor.composition.species
    molecular_weights = get_molecular_weights(species)

    rho_endpoint = endpoint_frame.species_densities
    rho_endpoint_cgs = _to_cgs_density(rho_endpoint, base_config.runtime.unit_system)
    mole_fractions = mass_densities_to_mole_fractions(rho_endpoint_cgs, molecular_weights)

    n_total_cgs = _total_number_density_from_cgs_density(rho_endpoint_cgs,
                                                         molecular_weights)
    n_total = base_config.runtime.unit_system == :SI ?
              convert_number_density_cgs_to_si(n_total_cgs) : n_total_cgs

    tt = endpoint_frame.temperatures.tt
    tv = endpoint_frame.temperatures.tv
    te = endpoint_frame.temperatures.te

    return ReactorConfig(composition = ReactorComposition(species,
                                                          mole_fractions,
                                                          n_total),
                         thermal = ReactorThermalState(tt, tv, te, te))
end

function _build_chain_diagnostics(base_config::Config,
                                  profile::AxialChainProfile,
                                  marching::AxialMarchingConfig,
                                  segment_end_reactors::Vector{ReactorConfig},
                                  segment_state_cache_used::AbstractVector{Bool},
                                  full_state_handoff_supported::Union{Nothing, Bool})
    diagnostics = Dict{String, Any}("handoff_mode" => String(marching.handoff_mode),
                                    "termination_mode" => String(marching.termination_mode),
                                    "override_tt_K" => marching.override_tt_K,
                                    "override_tv_K" => marching.override_tv_K,
                                    "segment_rho_ex_handoff_applied" => Bool[flag
                                                                             for flag in segment_state_cache_used])

    if marching.handoff_mode == :full_state
        diagnostics["full_state_rho_ex_handoff_supported"] = full_state_handoff_supported ===
                                                             true
        if full_state_handoff_supported === false
            diagnostics["full_state_rho_ex_handoff_reason"] = "electronic_sts_inactive"
        end
    end

    if !isempty(profile.diagnostics)
        diagnostics["input_profile_diagnostics"] = Dict{String, Vector{Float64}}(key => copy(values)
                                                                                 for (key, values) in pairs(profile.diagnostics))
    end

    if profile.wall_profile !== nothing
        diagnostics["wall_profile"] = _chain_wall_profile_to_dict(profile.wall_profile)
    end

    simulated_n_total = [reactor.composition.total_number_density
                         for reactor in segment_end_reactors]
    diagnostics["simulated_total_number_density"] = simulated_n_total

    if haskey(profile.diagnostics, "n_total_m3")
        profile_n_total = if base_config.runtime.unit_system == :SI
            copy(profile.diagnostics["n_total_m3"])
        else
            profile.diagnostics["n_total_m3"] .* 1e-6
        end
        diagnostics["input_total_number_density"] = profile_n_total
        diagnostics["total_number_density_abs_error"] = abs.(simulated_n_total .-
                                                             profile_n_total)
    end

    return diagnostics
end

"""
$(SIGNATURES)

Solve a steady axial chain of CSTR segments using the TERRA 0D solver.

This implementation supports:
- `handoff_mode = :reinitialize`
- `handoff_mode = :full_state`
- `termination_mode = :final_time`

# Arguments
- `config::Config`: Base TERRA config for models, sources, solver controls, and inlet state
- `profile::AxialChainProfile`: Axial chain profile

# Keyword Arguments
- `marching::AxialMarchingConfig`: Axial marching controls

# Returns
- `ChainSimulationResult`
"""
function solve_terra_chain_steady(config::Config,
                                  profile::AxialChainProfile;
                                  marching::AxialMarchingConfig = AxialMarchingConfig())
    validate_axial_chain_profile(profile)
    validate_axial_marching_config(marching)
    _validate_profile_species(profile, config)
    _validate_profile_inlet(profile, config)
    chain_runtime = config.runtime

    return with_active_runtime_for_logging(chain_runtime) do
        try
            if !is_terra_loaded()
                load_terra_library!()
            end
            set_api_finalize_mpi_wrapper(false)
        catch
            # Older libraries may not expose this control; proceed best-effort.
        end

        n_segments = length(profile.z_m)
        compact_to_source_index = _resolve_compact_to_source_index(profile)
        header_entry = ChainHeaderEntry(config, profile, marching, compact_to_source_index)
        _prepare_chain_logging(chain_runtime)
        _emit_chain_summary(chain_runtime, header_entry)
        _emit_chain_detail(chain_runtime, header_entry)

        segment_end_reactors = [config.reactor for _ in 1:n_segments]
        segment_success = fill(false, n_segments)
        segment_messages = fill("Not executed.", n_segments)
        segment_results = [_failed_simulation_result("Not executed.") for _ in 1:n_segments]
        segment_state_cache_used = fill(false, n_segments)

        inlet_reactor = _build_profile_inlet_reactor(profile, config.runtime.unit_system)
        upstream_state_cache = nothing
        full_state_handoff_supported = nothing

        for k in 1:n_segments
            local_result = nothing
            local_state_cache = nothing
            segment_case_path = _segment_case_path(config.runtime.case_path, k)
            segment_entry = ChainSegmentEntry(config, profile, compact_to_source_index, k,
                                              segment_case_path,
                                              segment_results[k])
            _emit_chain_summary(chain_runtime, segment_entry)

            try
                segment_config = _build_chain_segment_config(config, profile, k, inlet_reactor,
                                                             marching)
                wall_inputs = _build_segment_wall_inputs(profile, k,
                                                         segment_config.sources.wall_losses)
                initialize_terra(segment_config, segment_config.runtime.case_path;
                                 lifecycle_console = :never,
                                 preserve_active_runtime = false)
                if marching.handoff_mode == :full_state &&
                   full_state_handoff_supported === nothing
                    full_state_handoff_supported = has_electronic_sts_wrapper()
                end
                requested_state_cache = (marching.handoff_mode == :full_state && k > 1) ?
                                        upstream_state_cache : nothing
                segment_state_cache_used[k] = requested_state_cache !== nothing &&
                                              full_state_handoff_supported === true &&
                                              requested_state_cache.rho_ex_cgs !== nothing
                local_result, local_state_cache = _solve_terra_0d_internal(segment_config;
                                                                           wall_inputs = wall_inputs,
                                                                           state_cache = requested_state_cache,
                                                                           presentation = :chain_segment)
            catch e
                segment_messages[k] = "Segment setup/integration threw exception: $(e)"
                segment_results[k] = _failed_simulation_result(segment_messages[k])
                _emit_chain_detail(chain_runtime,
                                   ChainSegmentEntry(config, profile,
                                                     compact_to_source_index, k,
                                                     segment_case_path,
                                                     segment_results[k];
                                                     state_cache_used = segment_state_cache_used[k]))
                diagnostics = _build_chain_diagnostics(config, profile, marching,
                                                       segment_end_reactors,
                                                       segment_state_cache_used,
                                                       full_state_handoff_supported)
                cells = _build_chain_cells(profile, compact_to_source_index,
                                           segment_end_reactors, segment_success,
                                           segment_messages, segment_results)
                metadata = _chain_profile_metadata(profile, diagnostics,
                                                   compact_to_source_index)
                chain_result = ChainSimulationResult(;
                                                     cells = cells,
                                                     metadata = metadata,
                                                     success = false,
                                                     failed_cell = k,
                                                     message = "Chain failed at segment $k during setup/integration.",)
                result_entry = ChainResultEntry(chain_result)
                _emit_chain_detail(chain_runtime, result_entry)
                _emit_chain_summary(chain_runtime, result_entry)
                return chain_result
            end

            segment_results[k] = local_result
            segment_success[k] = local_result.success
            segment_messages[k] = local_result.message

            if !local_result.success
                _emit_chain_detail(chain_runtime,
                                   ChainSegmentEntry(config, profile,
                                                     compact_to_source_index, k,
                                                     segment_case_path, local_result;
                                                     state_cache_used = segment_state_cache_used[k]))
                diagnostics = _build_chain_diagnostics(config, profile, marching,
                                                       segment_end_reactors,
                                                       segment_state_cache_used,
                                                       full_state_handoff_supported)
                cells = _build_chain_cells(profile, compact_to_source_index,
                                           segment_end_reactors, segment_success,
                                           segment_messages, segment_results)
                metadata = _chain_profile_metadata(profile, diagnostics,
                                                   compact_to_source_index)
                chain_result = ChainSimulationResult(;
                                                     cells = cells,
                                                     metadata = metadata,
                                                     success = false,
                                                     failed_cell = k,
                                                     message = "Chain failed at segment $k: $(local_result.message)",)
                result_entry = ChainResultEntry(chain_result)
                _emit_chain_detail(chain_runtime, result_entry)
                _emit_chain_summary(chain_runtime, result_entry)
                return chain_result
            end

            endpoint_reactor = _extract_segment_endpoint_reactor(config, local_result)
            segment_end_reactors[k] = endpoint_reactor
            _emit_chain_detail(chain_runtime,
                               ChainSegmentEntry(config, profile,
                                                 compact_to_source_index, k,
                                                 segment_case_path, local_result;
                                                 endpoint_reactor = endpoint_reactor,
                                                 state_cache_used = segment_state_cache_used[k]))
            inlet_reactor = endpoint_reactor
            upstream_state_cache = local_state_cache
        end

        diagnostics = _build_chain_diagnostics(config, profile, marching,
                                               segment_end_reactors,
                                               segment_state_cache_used,
                                               full_state_handoff_supported)
        cells = _build_chain_cells(profile, compact_to_source_index,
                                   segment_end_reactors, segment_success,
                                   segment_messages, segment_results)
        metadata = _chain_profile_metadata(profile, diagnostics, compact_to_source_index)
        chain_result = ChainSimulationResult(;
                                             cells = cells,
                                             metadata = metadata,
                                             success = true,
                                             failed_cell = nothing,
                                             message = "chain success!",)
        result_entry = ChainResultEntry(chain_result)
        _emit_chain_detail(chain_runtime, result_entry)
        _emit_chain_summary(chain_runtime, result_entry)
        return chain_result
    end
end

"""
$(SIGNATURES)

Load a profile from disk and solve the steady axial chain of CSTR segments.
"""
function solve_terra_chain_steady(config::Config,
                                  profile_path::AbstractString;
                                  marching::AxialMarchingConfig = AxialMarchingConfig())
    profile = load_chain_profile(profile_path)
    return solve_terra_chain_steady(config, profile; marching = marching)
end
