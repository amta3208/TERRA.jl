@inline function _prepare_source_terms_data(layout::ApiLayout, config::Config,
                                            u0::Vector{Float64}, ::Nothing;
                                            wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                                            inlet_state_cache::Union{Nothing,
                                                                     ReactorStateCache} = nothing)
    return (residence_time = nothing, wall_losses = nothing)
end

function _prepare_source_terms_data(layout::ApiLayout, config::Config,
                                    u0::Vector{Float64}, sources::SourceTermsConfig;
                                    wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                                    inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    rt_cfg = sources.residence_time
    residence_time = (rt_cfg === nothing || !rt_cfg.enabled) ? nothing :
                     _prepare_residence_time_data(layout, config, u0, rt_cfg;
                                                  inlet_state_cache = inlet_state_cache)
    wall_cfg = sources.wall_losses
    wall_losses = (wall_cfg === nothing || !wall_cfg.enabled) ? nothing :
                  _prepare_wall_loss_data(layout, config, wall_cfg;
                                          wall_inputs = wall_inputs)
    return (residence_time = residence_time, wall_losses = wall_losses)
end

@inline function _apply_source_terms!(du::Vector{Float64}, u::Vector{Float64}, sources)
    sources === nothing && return nothing
    if hasproperty(sources, :residence_time)
        _apply_residence_time_term!(du, u, sources.residence_time)
    end
    if hasproperty(sources, :wall_losses)
        _apply_wall_loss_term!(du, u, sources.wall_losses)
    end
    return nothing
end

function _source_terms_frame_snapshot(u::Vector{Float64}, sources, unit_system::Symbol)
    sources === nothing && return nothing

    wall_snapshot = if hasproperty(sources, :wall_losses)
        _wall_loss_source_snapshot(u, sources.wall_losses, unit_system)
    else
        nothing
    end

    wall_snapshot === nothing && return nothing
    return (wall_losses = wall_snapshot,)
end

function _source_terms_result_metadata(sources)
    metadata = Dict{String, Any}()
    sources === nothing && return metadata

    if hasproperty(sources, :wall_losses)
        wall_metadata = _wall_loss_metadata_dict(sources.wall_losses)
        wall_metadata === nothing || (metadata["wall_losses"] = wall_metadata)
    end

    return metadata
end
