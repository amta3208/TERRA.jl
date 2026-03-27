abstract type AbstractPreparedSource end

struct PreparedSources
    operators::Vector{AbstractPreparedSource}
end

PreparedSources() = PreparedSources(AbstractPreparedSource[])
PreparedSources(operators::AbstractPreparedSource...) =
    PreparedSources(AbstractPreparedSource[operator for operator in operators])

function source_operator(sources::PreparedSources, ::Type{T}) where {T <: AbstractPreparedSource}
    for operator in sources.operators
        operator isa T && return operator
    end
    return nothing
end

prepare_source(layout::ApiLayout, config::Config, u0::Vector{Float64}, ::Nothing;
               wall_inputs = nothing,
               inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing) = nothing

function _prepare_sources(layout::ApiLayout, config::Config, u0::Vector{Float64}, ::Nothing;
                          wall_inputs = nothing,
                          inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    return PreparedSources()
end

function _prepare_sources(layout::ApiLayout, config::Config,
                          u0::Vector{Float64}, sources::SourceTermsConfig;
                          wall_inputs = nothing,
                          inlet_state_cache::Union{Nothing, ReactorStateCache} = nothing)
    operators = AbstractPreparedSource[]

    residence_time = prepare_source(layout, config, u0, sources.residence_time;
                                    inlet_state_cache = inlet_state_cache)
    residence_time === nothing || push!(operators, residence_time)

    wall_losses = prepare_source(layout, config, u0, sources.wall_losses;
                                 wall_inputs = wall_inputs)
    wall_losses === nothing || push!(operators, wall_losses)

    return PreparedSources(operators)
end

@inline _apply_sources!(du::Vector{Float64}, u::Vector{Float64}, ::Nothing) = nothing

function _apply_sources!(du::Vector{Float64}, u::Vector{Float64}, sources::PreparedSources)
    for operator in sources.operators
        apply_source!(du, u, operator)
    end
    return nothing
end

function _source_snapshot_namedtuple(entries::Vector{Pair{Symbol, Any}})
    keys = Tuple(first(entry) for entry in entries)
    values = Tuple(last(entry) for entry in entries)
    return NamedTuple{keys}(values)
end

function _source_frame_snapshot(u::Vector{Float64},
                                sources::PreparedSources,
                                unit_system::Symbol)
    entries = Pair{Symbol, Any}[]
    for operator in sources.operators
        snapshot = source_frame_snapshot(u, operator, unit_system)
        snapshot === nothing || push!(entries, source_name(operator) => snapshot)
    end
    isempty(entries) && return nothing
    return _source_snapshot_namedtuple(entries)
end

function _source_result_metadata(sources::PreparedSources)
    metadata = Dict{String, Any}()
    for operator in sources.operators
        operator_metadata = source_result_metadata(operator)
        operator_metadata === nothing ||
            (metadata[String(source_name(operator))] = operator_metadata)
    end
    return metadata
end
