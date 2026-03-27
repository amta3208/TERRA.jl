function _chain_cells(profile::AxialChainProfile,
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

function _chain_diagnostics(base_config::Config,
                            profile::AxialChainProfile,
                            marching::AxialMarchingConfig,
                            segment_end_reactors::AbstractVector{<:ReactorConfig},
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

function _chain_result(base_config::Config,
                       profile::AxialChainProfile,
                       marching::AxialMarchingConfig,
                       compact_to_source_index::AbstractVector{<:Integer},
                       segment_end_reactors::AbstractVector{<:ReactorConfig},
                       segment_success::AbstractVector{Bool},
                       segment_messages::AbstractVector{<:AbstractString},
                       segment_results::AbstractVector{<:ReactorResult},
                       segment_state_cache_used::AbstractVector{Bool},
                       full_state_handoff_supported::Union{Nothing, Bool};
                       success::Bool,
                       failed_cell::Union{Nothing, Integer},
                       message::AbstractString)
    diagnostics = _chain_diagnostics(base_config, profile, marching,
                                     segment_end_reactors,
                                     segment_state_cache_used,
                                     full_state_handoff_supported)
    cells = _chain_cells(profile, compact_to_source_index,
                         segment_end_reactors, segment_success,
                         segment_messages, segment_results)
    metadata = _chain_profile_metadata(profile, diagnostics, compact_to_source_index)
    return ChainSimulationResult(;
                                 cells = cells,
                                 metadata = metadata,
                                 success = success,
                                 failed_cell = failed_cell,
                                 message = message,)
end
