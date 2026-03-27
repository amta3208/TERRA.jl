# -----------------------------------------------------------------------------
# Reactor CSV persistence
# -----------------------------------------------------------------------------

"""
$(SIGNATURES)

Save TERRA results to file.

# Arguments
- `results::ReactorResult`: Results to save
- `filename::String`: Output filename (CSV format)

# Returns
- `true` if save successful
"""
function save_results(results::ReactorResult, filename::String)
    try
        species_densities = species_density_matrix(results)
        temperatures = temperature_history(results)
        total_energy = total_energy_history(results)

        n_times = length(results.t)
        n_species = size(species_densities, 1)

        header = ["time", "total_energy", "T_trans", "T_electron", "T_vib"]
        for i in 1:n_species
            push!(header, "species_$(i)_density")
        end

        data = zeros(n_times, length(header))
        data[:, 1] = results.t
        data[:, 2] = total_energy
        data[:, 3] = temperatures.tt
        data[:, 4] = temperatures.te
        data[:, 5] = temperatures.tv

        for i in 1:n_species
            data[:, 5 + i] = species_densities[i, :]
        end

        open(filename, "w") do io
            println(io, join(header, ","))
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end

        @info "Results saved successfully" filename=filename
        return true
    catch e
        @error "Failed to save results" filename=filename exception=e
        return false
    end
end

# -----------------------------------------------------------------------------
# Shared JSON coercion and validation
# -----------------------------------------------------------------------------

const CHAIN_RESULTS_SCHEMA_VERSION = "terra_chain_results_v3"

function json_value(value)
    if value === nothing
        return nothing
    elseif value isa Bool || value isa Real
        return value
    elseif value isa AbstractString
        return String(value)
    elseif value isa Symbol
        return String(value)
    elseif value isa NamedTuple
        dict = Dict{String, Any}()
        for (key, entry) in pairs(value)
            dict[String(key)] = json_value(entry)
        end
        return dict
    elseif value isa AbstractDict
        dict = Dict{String, Any}()
        for (key, entry) in pairs(value)
            dict[String(key)] = json_value(entry)
        end
        return dict
    elseif value isa AbstractVector
        return [json_value(entry) for entry in value]
    end
    return string(value)
end

function read_json_value(value)
    if value isa AbstractDict
        dict = Dict{String, Any}()
        for (key, entry) in pairs(value)
            dict[String(key)] = read_json_value(entry)
        end
        return dict
    elseif value isa AbstractVector
        return Any[read_json_value(entry) for entry in value]
    end
    return value
end

function _wall_species_model_metadata(raw::AbstractDict)
    model = Dict{String, Any}(String(key) => value for (key, value) in pairs(raw))
    haskey(model, "model_type") ||
        throw(ArgumentError("Wall-loss species model metadata is missing required `model_type`."))
    return model
end

function _reactor_metadata(metadata::AbstractDict)
    result = Dict{String, Any}(String(key) => value for (key, value) in pairs(metadata))
    wall_losses = get(result, "wall_losses", nothing)
    wall_losses === nothing && return result
    wall_losses isa AbstractDict ||
        throw(ArgumentError("Expected `wall_losses` in reactor metadata to be an object."))

    wall_metadata = Dict{String, Any}(String(key) => value for (key, value) in pairs(wall_losses))
    species_models = get(wall_metadata, "species_models", nothing)
    if species_models !== nothing
        species_models isa AbstractDict ||
            throw(ArgumentError("Expected `species_models` in reactor metadata wall losses to be an object."))
        wall_metadata["species_models"] = Dict{String, Any}(String(reactant) =>
                                                            begin
                                                                raw_model isa AbstractDict ||
                                                                    throw(ArgumentError("Expected wall-loss metadata for reactant `$reactant` to be an object."))
                                                                _wall_species_model_metadata(raw_model)
                                                            end
                                                            for (reactant, raw_model) in pairs(species_models))
    end

    result["wall_losses"] = wall_metadata
    return result
end

function required_array(container::AbstractDict, key::AbstractString,
                        context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    value isa AbstractVector ||
        throw(ArgumentError("Expected `$key` in $context to be an array."))
    return value
end

function required_dict(container::AbstractDict, key::AbstractString,
                       context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    value isa AbstractDict ||
        throw(ArgumentError("Expected `$key` in $context to be an object."))
    return value
end

function required_bool(container::AbstractDict, key::AbstractString,
                       context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    value isa Bool || throw(ArgumentError("Expected `$key` in $context to be a Bool."))
    return value
end

function required_string(container::AbstractDict, key::AbstractString,
                         context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    value isa AbstractString ||
        throw(ArgumentError("Expected `$key` in $context to be a String."))
    return String(value)
end

function int_value(value, key::AbstractString, context::AbstractString)
    if value isa Integer
        return Int(value)
    end
    if value isa Real && isfinite(value) && isinteger(value)
        return Int(round(value))
    end
    throw(ArgumentError("Expected `$key` in $context to be an integer."))
end

function required_int(container::AbstractDict, key::AbstractString,
                      context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    return int_value(container[key], key, context)
end

function optional_int(value, key::AbstractString, context::AbstractString)
    value === nothing && return nothing
    return int_value(value, key, context)
end

function float_array(container::AbstractDict, key::AbstractString,
                     context::AbstractString)
    values = required_array(container, key, context)
    result = Float64[]
    sizehint!(result, length(values))
    for (i, value) in pairs(values)
        value isa Real ||
            throw(ArgumentError("Expected numeric values in `$context.$key`; got $(typeof(value)) at index $i."))
        push!(result, Float64(value))
    end
    return result
end

function float_dict(container::AbstractDict, key::AbstractString,
                    context::AbstractString)
    raw = required_dict(container, key, context)
    result = Dict{String, Float64}()
    for (name, value) in pairs(raw)
        value isa Real ||
            throw(ArgumentError("Expected numeric values in `$context.$key`; got $(typeof(value)) for key $(name)."))
        result[String(name)] = Float64(value)
    end
    return result
end

function int_array(container::AbstractDict, key::AbstractString,
                   context::AbstractString)
    values = required_array(container, key, context)
    result = Int[]
    sizehint!(result, length(values))
    for (i, value) in pairs(values)
        if value isa Integer
            push!(result, Int(value))
        elseif value isa Real && isfinite(value) && isinteger(value)
            push!(result, Int(round(value)))
        else
            throw(ArgumentError("Expected integer values in `$context.$key`; got $(typeof(value)) at index $i."))
        end
    end
    return result
end

function namedtuple_value(raw::AbstractDict)
    keys_sorted = sort!(collect(String(key) for key in keys(raw)))
    syms = Symbol.(keys_sorted)
    vals = Tuple(read_json_value(raw[key]) for key in keys_sorted)
    return NamedTuple{Tuple(syms)}(vals)
end

function numeric_namedtuple_value(raw::AbstractDict, context::AbstractString)
    keys_sorted = sort!(collect(String(key) for key in keys(raw)))
    values = Vector{Float64}(undef, length(keys_sorted))
    for (i, key) in pairs(keys_sorted)
        value = raw[key]
        value isa Real ||
            throw(ArgumentError("Expected numeric values in `$context`; got $(typeof(value)) for key `$key`."))
        values[i] = Float64(value)
    end
    syms = Symbol.(keys_sorted)
    return NamedTuple{Tuple(syms)}(Tuple(values))
end

# -----------------------------------------------------------------------------
# Chain JSON codecs
# -----------------------------------------------------------------------------

reactor_frame_dict(frame::ReactorFrame) = Dict{String, Any}(
    "t" => frame.t,
    "species_densities" => copy(frame.species_densities),
    "temperatures" => json_value(frame.temperatures),
    "total_energy" => frame.total_energy,
    "source_terms" => frame.source_terms === nothing ? nothing : json_value(frame.source_terms),
    "diagnostics" => json_value(frame.diagnostics),
)

function reactor_frame(raw::AbstractDict)
    source_terms_raw = get(raw, "source_terms", nothing)
    if source_terms_raw !== nothing && !isa(source_terms_raw, AbstractDict)
        throw(ArgumentError("Expected `source_terms` in reactor frame to be an object or null."))
    end
    diagnostics_raw = get(raw, "diagnostics", Dict{String, Any}())
    diagnostics_raw isa AbstractDict ||
        throw(ArgumentError("Expected `diagnostics` in reactor frame to be an object."))

    return ReactorFrame(; t = begin
                            t = get(raw, "t", nothing)
                            t isa Real ||
                                throw(ArgumentError("Expected `t` in reactor frame to be numeric."))
                            Float64(t)
                        end,
                        species_densities = float_array(raw, "species_densities",
                                                        "reactor frame"),
                        temperatures = numeric_namedtuple_value(required_dict(raw,
                                                                              "temperatures",
                                                                              "reactor frame"),
                                                                "reactor frame.temperatures"),
                        total_energy = begin
                            e = get(raw, "total_energy", nothing)
                            e isa Real ||
                                throw(ArgumentError("Expected `total_energy` in reactor frame to be numeric."))
                            Float64(e)
                        end,
                        source_terms = source_terms_raw === nothing ? nothing :
                                       namedtuple_value(source_terms_raw),
                        diagnostics = read_json_value(diagnostics_raw))
end

reactor_result_dict(result::ReactorResult) = Dict{String, Any}(
    "t" => copy(result.t),
    "frames" => [reactor_frame_dict(frame) for frame in result.frames],
    "success" => result.success,
    "message" => result.message,
    "source_terms" => result.source_terms === nothing ? nothing : json_value(result.source_terms),
    "metadata" => json_value(result.metadata),
)

function reactor_result(raw::AbstractDict)
    frames_raw = required_array(raw, "frames", "reactor result")
    frames = Vector{ReactorFrame}(undef, length(frames_raw))
    for (i, frame_raw) in pairs(frames_raw)
        frame_raw isa AbstractDict ||
            throw(ArgumentError("Expected `frames[$i]` in reactor result to be an object."))
        frames[i] = reactor_frame(frame_raw)
    end

    source_terms_raw = get(raw, "source_terms", nothing)
    if source_terms_raw !== nothing && !isa(source_terms_raw, AbstractDict)
        throw(ArgumentError("Expected `source_terms` in reactor result to be an object or null."))
    end
    metadata_raw = get(raw, "metadata", Dict{String, Any}())
    metadata_raw isa AbstractDict ||
        throw(ArgumentError("Expected `metadata` in reactor result to be an object."))

    return ReactorResult(; t = float_array(raw, "t", "reactor result"),
                         frames = frames,
                         success = required_bool(raw, "success", "reactor result"),
                         message = required_string(raw, "message", "reactor result"),
                         source_terms = source_terms_raw === nothing ? nothing :
                                        namedtuple_value(source_terms_raw),
                         metadata = _reactor_metadata(read_json_value(metadata_raw)))
end

reactor_config_dict(config::ReactorConfig) = Dict{String, Any}(
    "composition" => Dict{String, Any}(
        "species" => copy(config.composition.species),
        "mole_fractions" => copy(config.composition.mole_fractions),
        "total_number_density" => config.composition.total_number_density,
    ),
    "thermal" => Dict{String, Any}(
        "Tt" => config.thermal.Tt,
        "Tv" => config.thermal.Tv,
        "Tee" => config.thermal.Tee,
        "Te" => config.thermal.Te,
    ),
)

function reactor_config(raw::AbstractDict)
    composition_raw = required_dict(raw, "composition", "reactor config")
    thermal_raw = required_dict(raw, "thermal", "reactor config")

    species_raw = required_array(composition_raw, "species", "reactor config.composition")
    species = String[]
    sizehint!(species, length(species_raw))
    for (i, value) in pairs(species_raw)
        value isa AbstractString ||
            throw(ArgumentError("Expected string values in `reactor config.composition.species`; got $(typeof(value)) at index $i."))
        push!(species, String(value))
    end

    return ReactorConfig(; composition = ReactorComposition(;
                                                            species = species,
                                                            mole_fractions = float_array(composition_raw,
                                                                                         "mole_fractions",
                                                                                         "reactor config.composition"),
                                                            total_number_density = begin
                                                                n = get(composition_raw,
                                                                        "total_number_density",
                                                                        nothing)
                                                                n isa Real ||
                                                                    throw(ArgumentError("Expected `total_number_density` in reactor config.composition to be numeric."))
                                                                Float64(n)
                                                            end),
                         thermal = ReactorThermalState(; Tt = begin
                                                            val = get(thermal_raw, "Tt",
                                                                      nothing)
                                                            val isa Real ||
                                                                throw(ArgumentError("Expected `Tt` in reactor config.thermal to be numeric."))
                                                            Float64(val)
                                                        end,
                                                        Tv = begin
                                                            val = get(thermal_raw, "Tv",
                                                                      nothing)
                                                            val isa Real ||
                                                                throw(ArgumentError("Expected `Tv` in reactor config.thermal to be numeric."))
                                                            Float64(val)
                                                        end,
                                                        Tee = begin
                                                            val = get(thermal_raw, "Tee",
                                                                      nothing)
                                                            val isa Real ||
                                                                throw(ArgumentError("Expected `Tee` in reactor config.thermal to be numeric."))
                                                            Float64(val)
                                                        end,
                                                        Te = begin
                                                            val = get(thermal_raw, "Te",
                                                                      nothing)
                                                            val isa Real ||
                                                                throw(ArgumentError("Expected `Te` in reactor config.thermal to be numeric."))
                                                            Float64(val)
                                                        end))
end

chain_metadata_dict(metadata::ChainMetadata) = Dict{String, Any}(
    "schema_version" => metadata.schema_version,
    "generator" => json_value(metadata.generator),
    "selection" => json_value(metadata.selection),
    "source_snapshot" => metadata.source_snapshot === nothing ? nothing :
                         json_value(metadata.source_snapshot),
    "diagnostics" => json_value(metadata.diagnostics),
    "compact_to_source_index" => copy(metadata.compact_to_source_index),
    "original_point_count" => metadata.original_point_count,
    "retained_point_count" => metadata.retained_point_count,
)

function chain_metadata(raw::AbstractDict)
    source_snapshot_raw = get(raw, "source_snapshot", nothing)
    if source_snapshot_raw !== nothing && !isa(source_snapshot_raw, AbstractDict)
        throw(ArgumentError("Expected `source_snapshot` in chain metadata to be an object or null."))
    end

    return ChainMetadata(; schema_version = required_string(raw,
                                                            "schema_version",
                                                            "chain metadata"),
                         generator = read_json_value(required_dict(raw,
                                                                  "generator",
                                                                  "chain metadata")),
                         selection = read_json_value(required_dict(raw,
                                                                  "selection",
                                                                  "chain metadata")),
                         source_snapshot = source_snapshot_raw === nothing ? nothing :
                                           read_json_value(source_snapshot_raw),
                         diagnostics = read_json_value(required_dict(raw,
                                                                    "diagnostics",
                                                                    "chain metadata")),
                         compact_to_source_index = int_array(raw,
                                                             "compact_to_source_index",
                                                             "chain metadata"),
                         original_point_count = optional_int(get(raw,
                                                                "original_point_count",
                                                                nothing),
                                                            "original_point_count",
                                                            "chain metadata"),
                         retained_point_count = required_int(raw,
                                                             "retained_point_count",
                                                             "chain metadata"))
end

chain_cell_dict(cell::ChainCellResult) = Dict{String, Any}(
    "compact_cell_index" => cell.compact_cell_index,
    "source_cell_index" => cell.source_cell_index,
    "z_m" => cell.z_m,
    "dx_m" => cell.dx_m,
    "te_K" => cell.te_K,
    "species_u_m_s" => Dict{String, Float64}(name => velocity
                                             for (name, velocity) in pairs(cell.species_u_m_s)),
    "reactor" => reactor_result_dict(cell.reactor),
    "endpoint_reactor" => cell.endpoint_reactor === nothing ? nothing :
                          reactor_config_dict(cell.endpoint_reactor),
    "success" => cell.success,
    "message" => cell.message,
)

function chain_cell(raw::AbstractDict)
    endpoint_reactor_raw = get(raw, "endpoint_reactor", nothing)
    if endpoint_reactor_raw !== nothing && !isa(endpoint_reactor_raw, AbstractDict)
        throw(ArgumentError("Expected `endpoint_reactor` in chain cell to be an object or null."))
    end

    return ChainCellResult(; compact_cell_index = required_int(raw,
                                                               "compact_cell_index",
                                                               "chain cell"),
                           source_cell_index = required_int(raw,
                                                            "source_cell_index",
                                                            "chain cell"),
                           z_m = begin
                               z = get(raw, "z_m", nothing)
                               z isa Real ||
                                   throw(ArgumentError("Expected `z_m` in chain cell to be numeric."))
                               Float64(z)
                           end,
                           dx_m = begin
                               dx = get(raw, "dx_m", nothing)
                               dx isa Real ||
                                   throw(ArgumentError("Expected `dx_m` in chain cell to be numeric."))
                               Float64(dx)
                           end,
                           te_K = begin
                               te = get(raw, "te_K", nothing)
                               te isa Real ||
                                   throw(ArgumentError("Expected `te_K` in chain cell to be numeric."))
                               Float64(te)
                           end,
                           species_u_m_s = float_dict(raw, "species_u_m_s", "chain cell"),
                           reactor = reactor_result(required_dict(raw,
                                                                  "reactor",
                                                                  "chain cell")),
                           endpoint_reactor = endpoint_reactor_raw === nothing ? nothing :
                                              reactor_config(endpoint_reactor_raw),
                           success = required_bool(raw, "success", "chain cell"),
                           message = required_string(raw, "message", "chain cell"))
end

"""
$(SIGNATURES)

Save chain simulation results to a JSON file.

# Arguments
- `results::ChainSimulationResult`: Chain results to save
- `filename::String`: Output JSON filename

# Returns
- `true` if save successful
"""
function save_results(results::ChainSimulationResult, filename::String)
    try
        output_dir = dirname(filename)
        output_dir != "." && !isempty(output_dir) && mkpath(output_dir)

        payload = Dict{String, Any}(
            "schema_version" => CHAIN_RESULTS_SCHEMA_VERSION,
            "success" => results.success,
            "failed_cell" => results.failed_cell,
            "message" => results.message,
            "metadata" => chain_metadata_dict(results.metadata),
            "cells" => [chain_cell_dict(cell) for cell in results.cells],
        )

        open(filename, "w") do io
            JSON.print(io, payload, 2)
            println(io)
        end

        @info "Chain results saved successfully" filename=filename cells=length(results.cells)
        return true
    catch e
        @error "Failed to save chain results" filename=filename exception=e
        return false
    end
end

"""
$(SIGNATURES)

Load chain simulation results from a JSON file generated by `save_results`.

# Arguments
- `filename::String`: Input JSON filename

# Returns
- `ChainSimulationResult`
"""
function load_results_chain(filename::String)
    if !isfile(filename)
        throw(ArgumentError("Chain results file does not exist: $filename"))
    end

    raw = JSON.parsefile(filename)
    raw isa AbstractDict ||
        throw(ArgumentError("Chain results root must be a JSON object."))

    schema_version = get(raw, "schema_version", nothing)
    schema_version == CHAIN_RESULTS_SCHEMA_VERSION ||
        throw(ArgumentError("Unsupported chain results schema version: $(schema_version). Expected `$CHAIN_RESULTS_SCHEMA_VERSION`."))

    metadata = chain_metadata(required_dict(raw, "metadata", "chain results"))
    cells_raw = required_array(raw, "cells", "chain results")
    cells = Vector{ChainCellResult}(undef, length(cells_raw))
    for (i, cell_raw) in pairs(cells_raw)
        cell_raw isa AbstractDict ||
            throw(ArgumentError("Expected `cells[$i]` in chain results to be an object."))
        cells[i] = chain_cell(cell_raw)
    end

    if metadata.retained_point_count != length(cells)
        throw(ArgumentError("Chain metadata retained_point_count ($(metadata.retained_point_count)) does not match saved cell count ($(length(cells)))."))
    end

    failed_cell = optional_int(get(raw, "failed_cell", nothing),
                               "failed_cell",
                               "chain results")

    return ChainSimulationResult(cells,
                                 metadata,
                                 required_bool(raw, "success", "chain results"),
                                 failed_cell,
                                 required_string(raw, "message", "chain results"))
end
