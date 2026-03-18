@inline function _prepare_source_terms_data(layout::ApiLayout, config::Config,
        u0::Vector{Float64}, ::Nothing;
        inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    return (residence_time = nothing,)
end

function _prepare_source_terms_data(layout::ApiLayout, config::Config,
        u0::Vector{Float64}, sources::SourceTermsConfig;
        inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    rt_cfg = sources.residence_time
    residence_time = (rt_cfg === nothing || !rt_cfg.enabled) ? nothing :
                     _prepare_residence_time_data(
        layout, config, u0, rt_cfg;
        inlet_state_cache = inlet_state_cache)
    return (residence_time = residence_time,)
end

@inline function _apply_source_terms!(du::Vector{Float64}, u::Vector{Float64}, sources)
    sources === nothing && return nothing
    if hasproperty(sources, :residence_time)
        _apply_residence_time_term!(du, u, sources.residence_time)
    end
    return nothing
end
