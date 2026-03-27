struct PreparedSources
    residence_time::Union{Nothing, PreparedResidenceTimeSource}
    wall_losses::Union{Nothing, PreparedWallLossData}
end

PreparedSources() = PreparedSources(nothing, nothing)

function _prepare_sources(layout::ApiLayout, config::Config, u0::Vector{Float64}, ::Nothing;
                          wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                          inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    return PreparedSources()
end

function _prepare_sources(layout::ApiLayout, config::Config,
                          u0::Vector{Float64}, sources::SourceTermsConfig;
                          wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                          inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    residence_time = _prepare_residence_time_source(layout, config, u0,
                                                    sources.residence_time;
                                                    inlet_state_cache = inlet_state_cache)
    wall_losses = _prepare_wall_losses(layout, config, sources.wall_losses;
                                       wall_inputs = wall_inputs)
    return PreparedSources(residence_time, wall_losses)
end

@inline _apply_sources!(du::Vector{Float64}, u::Vector{Float64}, ::Nothing) = nothing

function _apply_sources!(du::Vector{Float64}, u::Vector{Float64}, sources::PreparedSources)
    _apply_residence_time!(du, u, sources.residence_time)
    _apply_wall_losses!(du, u, sources.wall_losses)
    return nothing
end

function _source_frame_snapshot(u::Vector{Float64},
                                sources::PreparedSources,
                                unit_system::Symbol)
    wall_snapshot = _wall_loss_frame_snapshot(u, sources.wall_losses, unit_system)
    wall_snapshot === nothing && return nothing
    return (wall_losses = wall_snapshot,)
end

function _source_result_metadata(sources::PreparedSources)
    metadata = Dict{String, Any}()
    wall_metadata = _wall_loss_metadata(sources.wall_losses)
    wall_metadata === nothing || (metadata["wall_losses"] = wall_metadata)
    return metadata
end
