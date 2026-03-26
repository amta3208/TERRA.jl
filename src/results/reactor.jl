"""
$(SIGNATURES)

Time-frame result for a single reactor solve.

# Fields
- `t::Float64`: Frame timestamp
- `species_densities::Vector{Float64}`: Species mass densities at this frame
- `temperatures::NamedTuple`: Temperature state at this frame
- `total_energy::Float64`: Total energy density at this frame
- `source_terms::Union{NamedTuple, Nothing}`: Optional source-term snapshot
- `diagnostics::Dict{String, Any}`: Optional frame-level diagnostics
"""
struct ReactorFrame
    t::Float64
    species_densities::Vector{Float64}
    temperatures::NamedTuple
    total_energy::Float64
    source_terms::Union{NamedTuple, Nothing}
    diagnostics::Dict{String, Any}

    function ReactorFrame(; t::Real,
                          species_densities,
                          temperatures::NamedTuple,
                          total_energy::Real,
                          source_terms::Union{NamedTuple, Nothing} = nothing,
                          diagnostics::AbstractDict = Dict{String, Any}())
        species_densities_vec = Float64.(species_densities)
        diagnostics_dict = Dict{String, Any}(String(k) => v
                                             for (k, v) in pairs(diagnostics))

        return new(Float64(t),
                   species_densities_vec,
                   temperatures,
                   Float64(total_energy),
                   source_terms,
                   diagnostics_dict)
    end
end

"""
$(SIGNATURES)

Time-history container for one reactor (one chain cell).

# Fields
- `t::Vector{Float64}`: Saved times
- `frames::Vector{ReactorFrame}`: Saved per-time reactor frames
- `success::Bool`: Reactor success flag
- `message::String`: Reactor status message
- `source_terms::Union{NamedTuple, Nothing}`: Optional source history payload
- `metadata::Dict{String, Any}`: Optional reactor metadata

# Notes
- Supports HallThruster-like slicing: `reactor[i]` returns a one-frame `ReactorResult`.
"""
struct ReactorResult
    t::Vector{Float64}
    frames::Vector{ReactorFrame}
    success::Bool
    message::String
    source_terms::Union{NamedTuple, Nothing}
    metadata::Dict{String, Any}

    function ReactorResult(; t,
                           frames,
                           success::Bool = true,
                           message::AbstractString = "",
                           source_terms::Union{NamedTuple, Nothing} = nothing,
                           metadata::AbstractDict = Dict{String, Any}())
        t_vec = Float64.(t)
        frames_vec = ReactorFrame[frame for frame in frames]
        length(frames_vec) == length(t_vec) ||
            throw(ArgumentError("ReactorResult: `t` and `frames` must have identical lengths."))

        metadata_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(metadata))

        return new(t_vec,
                   frames_vec,
                   success,
                   String(message),
                   source_terms,
                   metadata_dict)
    end
end

Base.firstindex(result::ReactorResult) = 1
Base.lastindex(result::ReactorResult) = length(result.frames)

function Base.getindex(result::ReactorResult, frame::Integer)
    return ReactorResult(;
                         t = [result.t[frame]],
                         frames = [result.frames[frame]],
                         success = result.success,
                         message = result.message,
                         source_terms = result.source_terms,
                         metadata = copy(result.metadata))
end

function Base.getindex(result::ReactorResult, frames::AbstractVector{<:Integer})
    return ReactorResult(; t = result.t[frames],
                         frames = result.frames[frames],
                         success = result.success,
                         message = result.message,
                         source_terms = result.source_terms,
                         metadata = copy(result.metadata))
end
