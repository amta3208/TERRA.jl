const CHAIN_SEGMENT_DIRNAME = "chain_segments"

_chain_segment_case_path(base_case_path::AbstractString, segment_index::Integer) =
    normpath(joinpath(base_case_path, CHAIN_SEGMENT_DIRNAME,
                      @sprintf("segment_%04d", segment_index)))

function _validate_chain_inputs(config::Config,
                                profile::AxialChainProfile,
                                marching::AxialMarchingConfig)
    validate_axial_chain_profile(profile)
    validate_axial_marching_config(marching)

    species = config.reactor.composition.species
    molecular_weights = get_molecular_weights(species)
    expected_species = String[species[i]
                              for i in eachindex(species)
                              if !_is_electron_species(species[i], molecular_weights[i])]
    profile_species = sort!(collect(keys(profile.species_u_m_s)))
    expected_sorted = sort!(copy(expected_species))

    if profile_species != expected_sorted
        missing_species = [name for name in expected_sorted if !(name in profile_species)]
        extra_species = [name for name in profile_species if !(name in expected_sorted)]
        throw(ArgumentError("AxialChainProfile species_u_m_s keys must match the non-electron reactor species. " *
                            "Missing: $(isempty(missing_species) ? "none" : join(missing_species, ", ")); " *
                            "Extra: $(isempty(extra_species) ? "none" : join(extra_species, ", "))."))
    end

    inlet_species = profile.inlet.composition.species
    inlet_species == species ||
        throw(ArgumentError("AxialChainProfile inlet composition species must match the initialized TERRA species ordering. " *
                            "Expected: $(join(species, ", ")); got: $(join(inlet_species, ", "))."))
    profile.inlet.source_compact_index == 1 ||
        throw(ArgumentError("AxialChainProfile inlet.source_compact_index must be 1 for segment-1 initialization."))
    return nothing
end

function _segment_wall_inputs(profile::AxialChainProfile,
                              segment_index::Integer,
                              wall_cfg::Union{Nothing, WallLossConfig} = nothing)
    wall_values = _segment_wall_profile_values(profile, segment_index, wall_cfg)
    wall_values === nothing && return nothing
    return SegmentWallInputs(; wall_values...)
end

function _chain_segment_config(base_config::Config,
                               profile::AxialChainProfile,
                               segment_index::Integer,
                               inlet_reactor::ReactorConfig,
                               marching::AxialMarchingConfig)
    segment_case_path = _chain_segment_case_path(base_config.runtime.case_path, segment_index)
    mkpath(segment_case_path)

    te_profile = profile.te_K[segment_index]
    tt = marching.override_tt_K === nothing ? inlet_reactor.thermal.Tt :
         marching.override_tt_K
    tv = marching.override_tv_K === nothing ? inlet_reactor.thermal.Tv :
         marching.override_tv_K
    tee = segment_tee(marching.handoff_policy, segment_index, te_profile, inlet_reactor)

    chain_physics = PhysicsConfig(bbh_model = base_config.models.physics.bbh_model,
                                  esc_model = base_config.models.physics.esc_model,
                                  ar_et_model = base_config.models.physics.ar_et_model,
                                  eex_noneq = base_config.models.physics.eex_noneq,
                                  ev_relax_set = base_config.models.physics.ev_relax_set,
                                  et_relax_set = base_config.models.physics.et_relax_set,
                                  radiation_length = base_config.models.physics.radiation_length,
                                  get_electron_density_by_charge_balance = base_config.models.physics.get_electron_density_by_charge_balance,
                                  min_sts_frac = base_config.models.physics.min_sts_frac,
                                  is_isothermal_teex = marching.is_isothermal_teex,
                                  energy_loss_per_eii = base_config.models.physics.energy_loss_per_eii)

    base_rt = base_config.sources.residence_time
    residence_time = ResidenceTimeConfig(L = profile.dx_m[segment_index],
                                         U_species = _profile_species_velocity_at_segment(profile,
                                                                                          segment_index),
                                         U_energy = base_rt === nothing ? nothing :
                                                    base_rt.U_energy,
                                         inlet_reactor = inlet_reactor)

    segment_log_dir = if base_config.runtime.logging.log_dir === nothing
        nothing
    else
        segment_suffix = relpath(segment_case_path, base_config.runtime.case_path)
        normpath(joinpath(log_dir(base_config.runtime), segment_suffix))
    end

    return Config(;
                  reactor = ReactorConfig(;
                                          composition = inlet_reactor.composition,
                                          thermal = ReactorThermalState(; Tt = tt,
                                                                        Tv = tv,
                                                                        Tee = tee,
                                                                        Te = te_profile)),
                  models = ModelConfig(; physics = chain_physics,
                                       processes = base_config.models.processes),
                  sources = SourceTermsConfig(;
                                              residence_time = residence_time,
                                              wall_losses = base_config.sources.wall_losses,),
                  numerics = base_config.numerics,
                  runtime = with_runtime(base_config.runtime;
                                         case_path = segment_case_path,
                                         logging = with_logging(base_config.runtime.logging;
                                                                native_stream_mode = _logging_mode_file_only(base_config.runtime.logging.native_stream_mode),
                                                                integration_detail_mode = _logging_mode_file_only(base_config.runtime.logging.integration_detail_mode),
                                                                chain_detail_mode = :off,
                                                                log_dir = segment_log_dir)),)
end

function _segment_endpoint_reactor(base_config::Config,
                                   segment_result::ReactorResult)
    segment_result.success ||
        throw(ArgumentError("Cannot extract endpoint reactor from an unsuccessful segment result."))
    isempty(segment_result.t) &&
        throw(ArgumentError("Cannot extract endpoint reactor: segment result contains no time samples."))

    endpoint_frame = segment_result.frames[end]
    species = base_config.reactor.composition.species
    molecular_weights = get_molecular_weights(species)
    rho_endpoint = Float64.(endpoint_frame.species_densities)
    rho_endpoint_cgs = base_config.runtime.unit_system == :SI ?
                       convert_density_si_to_cgs(rho_endpoint) : rho_endpoint
    mole_fractions = mass_densities_to_mole_fractions(rho_endpoint_cgs, molecular_weights)
    n_total_cgs = mass_densities_to_total_number_density(rho_endpoint_cgs,
                                                         molecular_weights)
    n_total = base_config.runtime.unit_system == :SI ?
              convert_number_density_cgs_to_si(n_total_cgs) : n_total_cgs

    return ReactorConfig(;
                         composition = ReactorComposition(species,
                                                          mole_fractions,
                                                          n_total),
                         thermal = ReactorThermalState(endpoint_frame.temperatures.tt,
                                                       endpoint_frame.temperatures.tv,
                                                       endpoint_frame.temperatures.te,
                                                       endpoint_frame.temperatures.te))
end

"""
$(SIGNATURES)

Solve a steady axial chain of CSTR segments using the TERRA 0D solver.

This implementation supports:
- `ReinitializeHandoff()`
- `FullStateHandoff()`
- `FinalTimeTermination()`

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
    _validate_chain_inputs(config, profile, marching)
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
        prepare!(CHAIN_LOG, chain_runtime)
        emit!(RUN_LOG, chain_runtime, header_entry)
        emit!(CHAIN_LOG, chain_runtime, header_entry)

        segment_end_reactors = [config.reactor for _ in 1:n_segments]
        segment_success = fill(false, n_segments)
        segment_messages = fill("Not executed.", n_segments)
        segment_results = [_failed_reactor_result("Not executed.") for _ in 1:n_segments]
        segment_state_cache_used = fill(false, n_segments)

        inlet_reactor = _profile_inlet_reactor(profile, config.runtime.unit_system)
        upstream_state_cache = nothing
        full_state_handoff_supported = nothing

        for k in 1:n_segments
            local_result = nothing
            local_state_cache = nothing
            segment_case_path = _chain_segment_case_path(config.runtime.case_path, k)
            emit!(RUN_LOG, chain_runtime,
                  ChainSegmentEntry(config, profile, compact_to_source_index, k,
                                    segment_case_path, segment_results[k]))

            try
                segment_config = _chain_segment_config(config, profile, k, inlet_reactor,
                                                       marching)
                wall_inputs = _segment_wall_inputs(profile, k,
                                                   segment_config.sources.wall_losses)
                initialize_terra(segment_config, segment_config.runtime.case_path;
                                 lifecycle_console = :never,
                                 preserve_active_runtime = false)
                if requires_state_cache_handoff(marching.handoff_policy) &&
                   full_state_handoff_supported === nothing
                    full_state_handoff_supported = has_electronic_sts_wrapper()
                end
                requested_state_cache = requested_handoff_state_cache(marching.handoff_policy,
                                                                      k,
                                                                      upstream_state_cache)
                segment_state_cache_used[k] = state_cache_handoff_used(marching.handoff_policy,
                                                                       requested_state_cache,
                                                                       full_state_handoff_supported)
                local_result, local_state_cache = _solve_terra_0d_internal(segment_config;
                                                                           wall_inputs = wall_inputs,
                                                                           state_cache = requested_state_cache,
                                                                           presentation = CHAIN_SEGMENT_PRESENTATION)
            catch e
                segment_messages[k] = "Segment setup/integration threw exception: $(e)"
                segment_results[k] = _failed_reactor_result(segment_messages[k])
                emit!(CHAIN_LOG, chain_runtime,
                      ChainSegmentEntry(config, profile,
                                        compact_to_source_index, k,
                                        segment_case_path,
                                        segment_results[k];
                                        state_cache_used = segment_state_cache_used[k]))
                chain_result = _chain_result(config, profile, marching,
                                             compact_to_source_index,
                                             segment_end_reactors,
                                             segment_success,
                                             segment_messages,
                                             segment_results,
                                             segment_state_cache_used,
                                             full_state_handoff_supported;
                                             success = false,
                                             failed_cell = k,
                                             message = "Chain failed at segment $k during setup/integration.",)
                result_entry = ChainResultEntry(chain_result)
                emit!(CHAIN_LOG, chain_runtime, result_entry)
                emit!(RUN_LOG, chain_runtime, result_entry)
                return chain_result
            end

            segment_results[k] = local_result
            segment_success[k] = local_result.success
            segment_messages[k] = local_result.message

            if !local_result.success
                emit!(CHAIN_LOG, chain_runtime,
                      ChainSegmentEntry(config, profile,
                                        compact_to_source_index, k,
                                        segment_case_path, local_result;
                                        state_cache_used = segment_state_cache_used[k]))
                chain_result = _chain_result(config, profile, marching,
                                             compact_to_source_index,
                                             segment_end_reactors,
                                             segment_success,
                                             segment_messages,
                                             segment_results,
                                             segment_state_cache_used,
                                             full_state_handoff_supported;
                                             success = false,
                                             failed_cell = k,
                                             message = "Chain failed at segment $k: $(local_result.message)",)
                result_entry = ChainResultEntry(chain_result)
                emit!(CHAIN_LOG, chain_runtime, result_entry)
                emit!(RUN_LOG, chain_runtime, result_entry)
                return chain_result
            end

            endpoint_reactor = _segment_endpoint_reactor(config, local_result)
            segment_end_reactors[k] = endpoint_reactor
            emit!(CHAIN_LOG, chain_runtime,
                  ChainSegmentEntry(config, profile,
                                    compact_to_source_index, k,
                                    segment_case_path, local_result;
                                    endpoint_reactor = endpoint_reactor,
                                    state_cache_used = segment_state_cache_used[k]))
            inlet_reactor = endpoint_reactor
            upstream_state_cache = local_state_cache
        end

        chain_result = _chain_result(config, profile, marching,
                                     compact_to_source_index,
                                     segment_end_reactors,
                                     segment_success,
                                     segment_messages,
                                     segment_results,
                                     segment_state_cache_used,
                                     full_state_handoff_supported;
                                     success = true,
                                     failed_cell = nothing,
                                     message = "chain success!",)
        result_entry = ChainResultEntry(chain_result)
        emit!(CHAIN_LOG, chain_runtime, result_entry)
        emit!(RUN_LOG, chain_runtime, result_entry)
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
