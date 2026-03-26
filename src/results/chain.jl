"""
$(SIGNATURES)

Metadata for chain-level result organization and segmentation provenance.
"""
struct ChainMetadata
    schema_version::String
    generator::Dict{String, Any}
    selection::Dict{String, Any}
    source_snapshot::Union{Nothing, Dict{String, Any}}
    diagnostics::Dict{String, Any}
    compact_to_source_index::Vector{Int}
    original_point_count::Union{Nothing, Int}
    retained_point_count::Int

    function ChainMetadata(; schema_version::AbstractString = "terra_chain_profile_v4",
                           generator::AbstractDict = Dict{String, Any}(),
                           selection::AbstractDict = Dict{String, Any}(),
                           source_snapshot::Union{Nothing, AbstractDict} = nothing,
                           diagnostics::AbstractDict = Dict{String, Any}(),
                           compact_to_source_index::AbstractVector{<:Integer} = Int[],
                           original_point_count::Union{Nothing, Integer} = nothing,
                           retained_point_count::Union{Nothing, Integer} = nothing)
        compact_to_source = Int.(compact_to_source_index)
        retained_points = retained_point_count === nothing ? length(compact_to_source) :
                          Int(retained_point_count)
        retained_points >= 0 ||
            throw(ArgumentError("ChainMetadata: retained_point_count must be non-negative."))
        if !isempty(compact_to_source) && length(compact_to_source) != retained_points
            throw(ArgumentError("ChainMetadata: compact_to_source_index length must match retained_point_count."))
        end

        original_points = original_point_count === nothing ? nothing :
                          Int(original_point_count)
        if original_points !== nothing && original_points < retained_points
            throw(ArgumentError("ChainMetadata: original_point_count must be >= retained_point_count."))
        end

        generator_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(generator))
        selection_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(selection))
        source_snapshot_dict = source_snapshot === nothing ? nothing :
                               Dict{String, Any}(String(k) => v
                                                 for (k, v) in pairs(source_snapshot))
        diagnostics_dict = Dict{String, Any}(String(k) => v
                                             for (k, v) in pairs(diagnostics))

        return new(String(schema_version),
                   generator_dict,
                   selection_dict,
                   source_snapshot_dict,
                   diagnostics_dict,
                   compact_to_source,
                   original_points,
                   retained_points)
    end
end

"""
$(SIGNATURES)

Per-cell chain result, including cell-aligned profile inputs and reactor time history.
"""
struct ChainCellResult
    compact_cell_index::Int
    source_cell_index::Int
    z_m::Float64
    dx_m::Float64
    te_K::Float64
    species_u_m_s::Dict{String, Float64}
    reactor::ReactorResult
    endpoint_reactor::Union{Nothing, ReactorConfig}
    success::Bool
    message::String

    function ChainCellResult(; compact_cell_index::Integer,
                             source_cell_index::Integer = compact_cell_index,
                             z_m::Real,
                             dx_m::Real,
                             te_K::Real,
                             species_u_m_s::AbstractDict,
                             reactor::ReactorResult,
                             endpoint_reactor::Union{Nothing, ReactorConfig} = nothing,
                             success::Bool = reactor.success,
                             message::AbstractString = reactor.message)
        compact_idx = Int(compact_cell_index)
        source_idx = Int(source_cell_index)
        compact_idx >= 1 ||
            throw(ArgumentError("ChainCellResult: compact_cell_index must be >= 1."))
        source_idx >= 1 ||
            throw(ArgumentError("ChainCellResult: source_cell_index must be >= 1."))

        z_val = Float64(z_m)
        dx_val = Float64(dx_m)
        te_val = Float64(te_K)
        species_u_dict = Dict{String, Float64}()
        for (name, value) in pairs(species_u_m_s)
            name_str = String(name)
            isempty(strip(name_str)) &&
                throw(ArgumentError("ChainCellResult: species_u_m_s keys must be non-empty."))
            value_f64 = Float64(value)
            isfinite(value_f64) && value_f64 > 0 ||
                throw(ArgumentError("ChainCellResult: species_u_m_s[$name_str] must be finite and positive."))
            species_u_dict[name_str] = value_f64
        end

        isfinite(z_val) || throw(ArgumentError("ChainCellResult: z_m must be finite."))
        isfinite(dx_val) && dx_val > 0 ||
            throw(ArgumentError("ChainCellResult: dx_m must be finite and positive."))
        isfinite(te_val) && te_val > 0 ||
            throw(ArgumentError("ChainCellResult: te_K must be finite and positive."))
        isempty(species_u_dict) &&
            throw(ArgumentError("ChainCellResult: species_u_m_s must contain at least one species."))

        return new(compact_idx,
                   source_idx,
                   z_val,
                   dx_val,
                   te_val,
                   species_u_dict,
                   reactor,
                   endpoint_reactor,
                   success,
                   String(message))
    end
end

"""
$(SIGNATURES)

Results container for axial-marching chain-of-CSTR simulations.

# Fields
- `cells::Vector{ChainCellResult}`: Solved chain cells (compact retained indexing)
- `metadata::ChainMetadata`: Chain-level metadata, diagnostics, and index mapping
- `success::Bool`: Overall chain success flag
- `failed_cell::Union{Nothing, Int}`: First failing retained-cell index, if any
- `message::String`: Overall status message
"""
struct ChainSimulationResult
    cells::Vector{ChainCellResult}
    metadata::ChainMetadata
    success::Bool
    failed_cell::Union{Nothing, Int}
    message::String

    function ChainSimulationResult(cells::Vector{ChainCellResult},
                                   metadata::ChainMetadata,
                                   success::Bool,
                                   failed_cell::Union{Nothing, Integer},
                                   message::AbstractString)
        failed_cell_val = failed_cell === nothing ? nothing : Int(failed_cell)
        if failed_cell_val !== nothing &&
           (failed_cell_val < 1 || failed_cell_val > length(cells))
            throw(ArgumentError("ChainSimulationResult: failed_cell must be `nothing` or in 1:length(cells)."))
        end
        return new(cells, metadata, success, failed_cell_val, String(message))
    end
end

function ChainSimulationResult(;
                               cells,
                               metadata::ChainMetadata = ChainMetadata(compact_to_source_index = [cell.source_cell_index
                                                                                                  for cell in cells],
                                                                       retained_point_count = length(cells)),
                               success::Bool = all(cell.success for cell in cells),
                               failed_cell::Union{Nothing, Integer} = nothing,
                               message::AbstractString = success ?
                                                         "Chain integration completed successfully." :
                                                         "Chain integration failed.")
    cells_vec = ChainCellResult[cell for cell in cells]
    return ChainSimulationResult(cells_vec, metadata, success, failed_cell, message)
end

Base.firstindex(chain::ChainSimulationResult) = 1
Base.lastindex(chain::ChainSimulationResult) = length(chain.cells)

function _slice_chain_metadata(metadata::ChainMetadata, cells::Vector{ChainCellResult})
    return ChainMetadata(; schema_version = metadata.schema_version,
                         generator = copy(metadata.generator),
                         selection = copy(metadata.selection),
                         source_snapshot = metadata.source_snapshot === nothing ? nothing :
                                           copy(metadata.source_snapshot),
                         diagnostics = copy(metadata.diagnostics),
                         compact_to_source_index = [cell.source_cell_index
                                                    for cell in cells],
                         original_point_count = metadata.original_point_count,
                         retained_point_count = length(cells))
end

function _slice_chain_result(chain::ChainSimulationResult,
                             indices::AbstractVector{<:Integer})
    cells = chain.cells[indices]
    metadata = _slice_chain_metadata(chain.metadata, cells)
    failed_cell = findfirst(!, [cell.success for cell in cells])
    success = failed_cell === nothing
    message = success ? chain.message : cells[failed_cell].message
    return ChainSimulationResult(cells, metadata, success, failed_cell, message)
end

function Base.getindex(chain::ChainSimulationResult, cell::Integer)
    return _slice_chain_result(chain, [cell])
end

function Base.getindex(chain::ChainSimulationResult, cells::AbstractVector{<:Integer})
    return _slice_chain_result(chain, cells)
end
