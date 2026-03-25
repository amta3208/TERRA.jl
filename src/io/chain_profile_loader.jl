const CHAIN_PROFILE_SCHEMA_VERSION = "terra_chain_profile_v4"

function _read_chain_profile_json(path::AbstractString)
    if !isfile(path)
        throw(ArgumentError("Chain profile file does not exist: $path"))
    end
    return JSON.parsefile(path)
end

function _require_dict_field(container::AbstractDict, key::AbstractString,
                             context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    isa(value, AbstractDict) ||
        throw(ArgumentError("Expected `$key` in $context to be a JSON object."))
    return value
end

function _require_array_field(container::AbstractDict, key::AbstractString,
                              context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    isa(value, Vector) ||
        throw(ArgumentError("Expected `$key` in $context to be a JSON array."))
    return value
end

function _require_number_field(container::AbstractDict, key::AbstractString,
                               context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    value isa Real || throw(ArgumentError("Expected `$key` in $context to be numeric."))
    return Float64(value)
end

function _require_float_array_field(container::AbstractDict, key::AbstractString,
                                    context::AbstractString)
    values = _require_array_field(container, key, context)
    result = Float64[]
    sizehint!(result, length(values))
    for (i, value) in pairs(values)
        if value isa Real
            push!(result, Float64(value))
        else
            throw(ArgumentError("Expected numeric values in `$context.$key`; got $(typeof(value)) at index $i."))
        end
    end
    return result
end

function _optional_float_array_field(container::AbstractDict, key::AbstractString,
                                     context::AbstractString)
    haskey(container, key) || return nothing
    return _require_float_array_field(container, key, context)
end

function _require_string_array_field(container::AbstractDict, key::AbstractString,
                                     context::AbstractString)
    values = _require_array_field(container, key, context)
    result = String[]
    sizehint!(result, length(values))
    for (i, value) in pairs(values)
        value isa AbstractString ||
            throw(ArgumentError("Expected string values in `$context.$key`; got $(typeof(value)) at index $i."))
        push!(result, String(value))
    end
    return result
end

function _coerce_optional_diagnostics(raw::Union{Nothing, AbstractDict})
    diagnostics = Dict{String, Vector{Float64}}()
    raw === nothing && return diagnostics

    for (key, value) in pairs(raw)
        key_str = String(key)
        isa(value, Vector) ||
            throw(ArgumentError("Diagnostic `$key_str` must be a JSON array."))
        values_vec = Float64[]
        sizehint!(values_vec, length(value))
        for (i, entry) in pairs(value)
            if entry isa Real
                push!(values_vec, Float64(entry))
            else
                throw(ArgumentError("Diagnostic `$key_str` must contain numeric values (failed at index $i)."))
            end
        end
        diagnostics[key_str] = values_vec
    end

    return diagnostics
end

function _require_float_array_dict_field(container::AbstractDict,
                                         key::AbstractString,
                                         context::AbstractString)
    raw = _require_dict_field(container, key, context)
    result = Dict{String, Vector{Float64}}()
    for (name, values) in pairs(raw)
        name_str = String(name)
        isa(values, Vector) ||
            throw(ArgumentError("Expected `$context.$key.$name_str` to be a JSON array."))
        values_vec = Float64[]
        sizehint!(values_vec, length(values))
        for (i, value) in pairs(values)
            if value isa Real
                push!(values_vec, Float64(value))
            else
                throw(ArgumentError("Expected numeric values in `$context.$key.$name_str`; got $(typeof(value)) at index $i."))
            end
        end
        result[name_str] = values_vec
    end
    return result
end

function _build_chain_profile_inlet(raw::AbstractDict)
    inlet = _require_dict_field(raw, "inlet", "chain profile root")
    composition = _require_dict_field(inlet, "composition", "chain profile inlet")
    thermal = _require_dict_field(inlet, "thermal", "chain profile inlet")

    source_compact_index = get(inlet, "source_compact_index", nothing)
    source_compact_index isa Integer ||
        (source_compact_index isa Real && isfinite(source_compact_index) &&
         isinteger(source_compact_index)) ||
        throw(ArgumentError("Expected `source_compact_index` in chain profile inlet to be an integer."))

    inlet_composition = ChainProfileInletComposition(;
                                                     species = _require_string_array_field(composition,
                                                                                           "species",
                                                                                           "inlet.composition"),
                                                     mole_fractions = _require_float_array_field(composition,
                                                                                                 "mole_fractions",
                                                                                                 "inlet.composition"),
                                                     total_number_density_m3 = _require_number_field(composition,
                                                                                                     "total_number_density_m3",
                                                                                                     "inlet.composition"),)
    inlet_thermal = ReactorThermalState(;
                                        Tt = _require_number_field(thermal, "Tt_K",
                                                                   "inlet.thermal"),
                                        Tv = _require_number_field(thermal, "Tv_K",
                                                                   "inlet.thermal"),
                                        Tee = _require_number_field(thermal, "Tee_K",
                                                                    "inlet.thermal"),
                                        Te = _require_number_field(thermal, "Te_K",
                                                                   "inlet.thermal"),)

    return ChainProfileInlet(;
                             composition = inlet_composition,
                             thermal = inlet_thermal,
                             source_compact_index = Int(round(source_compact_index)),)
end

function _build_chain_wall_profile(raw::AbstractDict)
    haskey(raw, "wall_profile") || return nothing
    wall_profile = _require_dict_field(raw, "wall_profile", "chain profile root")
    return ChainWallProfile(;
                            a_wall_over_v_m_inv = _require_float_array_field(wall_profile,
                                                                             "a_wall_over_v_m_inv",
                                                                             "wall_profile"),
                            channel_gap_m = _optional_float_array_field(wall_profile,
                                                                        "channel_gap_m",
                                                                        "wall_profile"),
                            wall_temperature_K = _optional_float_array_field(wall_profile,
                                                                             "wall_temperature_K",
                                                                             "wall_profile"),
                            ion_edge_to_center_ratio = _optional_float_array_field(wall_profile,
                                                                                   "ion_edge_to_center_ratio",
                                                                                   "wall_profile"),)
end

function _validate_chain_profile_schema(raw::AbstractDict)
    schema_version = get(raw, "schema_version", nothing)
    schema_version === CHAIN_PROFILE_SCHEMA_VERSION ||
        throw(ArgumentError("Unsupported chain profile schema version: $(schema_version). Expected `$CHAIN_PROFILE_SCHEMA_VERSION`."))

    generator = _require_dict_field(raw, "generator", "chain profile root")
    selection = _require_dict_field(raw, "selection", "chain profile root")
    profile = _require_dict_field(raw, "profile", "chain profile root")
    inlet = _require_dict_field(raw, "inlet", "chain profile root")
    inlet_composition = _require_dict_field(inlet, "composition", "chain profile inlet")
    inlet_thermal = _require_dict_field(inlet, "thermal", "chain profile inlet")

    for key in ("tool", "tool_version", "created_utc")
        haskey(generator, key) ||
            throw(ArgumentError("Missing required generator field `$key`."))
    end

    for key in ("average_start_time_s",
                "exported_species",
                "ion_velocity_policy",
                "u_ion_floor",
                "min_consecutive_positive",
                "trim_start_index",
                "trim_start_z_m",
                "trimmed_point_count",
                "original_point_count")
        haskey(selection, key) ||
            throw(ArgumentError("Missing required selection field `$key`."))
    end

    for key in ("z_m", "dx_m", "te_K", "species_u_m_s")
        haskey(profile, key) ||
            throw(ArgumentError("Missing required profile field `$key`."))
    end

    if haskey(raw, "wall_profile")
        _require_dict_field(raw, "wall_profile", "chain profile root")
    end

    for key in ("species", "mole_fractions", "total_number_density_m3")
        haskey(inlet_composition, key) ||
            throw(ArgumentError("Missing required inlet composition field `$key`."))
    end

    for key in ("Tt_K", "Tv_K", "Tee_K", "Te_K")
        haskey(inlet_thermal, key) ||
            throw(ArgumentError("Missing required inlet thermal field `$key`."))
    end
    haskey(inlet, "source_compact_index") ||
        throw(ArgumentError("Missing required inlet field `source_compact_index`."))

    return nothing
end

"""
$(SIGNATURES)

Load a chain profile JSON artifact into normalized `AxialChainProfile`.

# Arguments
- `path::AbstractString`: Path to `chain_profile_v4.json`

# Returns
- `AxialChainProfile`

# Throws
- `ArgumentError` on file, schema, or data validation failures
"""
function load_chain_profile(path::AbstractString)
    raw = _read_chain_profile_json(path)
    isa(raw, AbstractDict) ||
        throw(ArgumentError("Chain profile root must be a JSON object."))

    _validate_chain_profile_schema(raw)

    generator = _require_dict_field(raw, "generator", "chain profile root")
    selection = _require_dict_field(raw, "selection", "chain profile root")
    profile = _require_dict_field(raw, "profile", "chain profile root")

    diagnostics_raw = haskey(raw, "diagnostics") ? raw["diagnostics"] : nothing
    if diagnostics_raw !== nothing && !isa(diagnostics_raw, AbstractDict)
        throw(ArgumentError("`diagnostics` must be a JSON object when provided."))
    end

    source_snapshot_raw = haskey(raw, "source_snapshot") ? raw["source_snapshot"] : nothing
    if source_snapshot_raw !== nothing && !isa(source_snapshot_raw, AbstractDict)
        throw(ArgumentError("`source_snapshot` must be a JSON object when provided."))
    end

    chain_profile = AxialChainProfile(z_m = _require_float_array_field(profile, "z_m",
                                                                       "profile"),
                                      dx_m = _require_float_array_field(profile, "dx_m",
                                                                        "profile"),
                                      te_K = _require_float_array_field(profile, "te_K",
                                                                        "profile"),
                                      species_u_m_s = _require_float_array_dict_field(profile,
                                                                                      "species_u_m_s",
                                                                                      "profile"),
                                      wall_profile = _build_chain_wall_profile(raw),
                                      inlet = _build_chain_profile_inlet(raw),
                                      diagnostics = _coerce_optional_diagnostics(diagnostics_raw),
                                      generator = generator,
                                      selection = selection,
                                      schema_version = String(raw["schema_version"]),
                                      source_snapshot = source_snapshot_raw)

    validate_axial_chain_profile(chain_profile)
    return chain_profile
end
