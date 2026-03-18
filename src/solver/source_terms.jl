@inline function _prepare_source_terms_data(layout::ApiLayout, config::Config,
        u0::Vector{Float64}, ::Nothing;
        wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
        inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    return (residence_time = nothing, wall_losses = nothing)
end

function _prepare_source_terms_data(layout::ApiLayout, config::Config,
        u0::Vector{Float64}, sources::SourceTermsConfig;
        wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
        inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    rt_cfg = sources.residence_time
    residence_time = (rt_cfg === nothing || !rt_cfg.enabled) ? nothing :
                     _prepare_residence_time_data(
        layout, config, u0, rt_cfg;
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
