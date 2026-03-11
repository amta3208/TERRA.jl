using Dates
using JSON

const SCHEMA_VERSION = "terra_chain_profile_v3"
const SCRIPT_VERSION = "0.4.0"
const EV_TO_K = 11604.518121550082
const DEFAULT_OUTPUT_BASENAME = "chain_profile_v3.json"

function usage()
    return """
    Usage:
      julia --project=tools/hallthruster_jl tools/hallthruster_jl/export_chain_profile.jl <hallthruster_input.json> <terra_case_path> --average_start_time=<seconds> [--write_source_snapshot=<true|false>] [--ion_velocity_policy=<trim_to_first_positive|strict_positive>] [--u_ion_floor=<float>] [--min_consecutive_positive=<int>]

    Notes:
      - The HallThruster JSON input must already contain an `output.average` block.
      - For JSON inputs, `--average_start_time` is recorded as metadata only; this script does not recompute the time average.
      - The output is written to `<terra_case_path>/input/chain_profile_v3.json`.
      - Default ion-velocity handling is `trim_to_first_positive` with `u_ion_floor=0.0` and `min_consecutive_positive=3`.
      - All non-electron neutral and ion species present in HallThruster are exported.
    """
end

function parse_bool(value::AbstractString)
    lowered = lowercase(strip(value))
    if lowered in ("true", "1", "yes", "y")
        return true
    elseif lowered in ("false", "0", "no", "n")
        return false
    end
    throw(ArgumentError("Invalid boolean value: $(value)"))
end

function parse_options(args::Vector{String})
    if any(arg -> arg in ("-h", "--help"), args)
        println(usage())
        return nothing
    end

    positionals = String[]
    options = Dict{String, String}()

    for arg in args
        if startswith(arg, "--")
            parts = split(arg[3:end], "="; limit = 2)
            length(parts) == 2 || throw(ArgumentError("Options must use --key=value syntax: $(arg)"))
            key, value = parts
            options[key] = value
        else
            push!(positionals, arg)
        end
    end

    length(positionals) == 2 || throw(ArgumentError("Expected 2 positional arguments.\n\n" * usage()))
    haskey(options, "average_start_time") ||
        throw(ArgumentError("Missing required option --average_start_time=...\n\n" * usage()))

    write_source_snapshot = haskey(options, "write_source_snapshot") ?
        parse_bool(options["write_source_snapshot"]) : true
    ion_velocity_policy = haskey(options, "ion_velocity_policy") ?
        options["ion_velocity_policy"] : "trim_to_first_positive"
    u_ion_floor = haskey(options, "u_ion_floor") ? parse(Float64, options["u_ion_floor"]) : 0.0
    min_consecutive_positive = haskey(options, "min_consecutive_positive") ?
        parse(Int, options["min_consecutive_positive"]) : 3

    return (
        hallthruster_input = abspath(positionals[1]),
        terra_case_path = abspath(positionals[2]),
        average_start_time_s = parse(Float64, options["average_start_time"]),
        write_source_snapshot = write_source_snapshot,
        ion_velocity_policy = ion_velocity_policy,
        u_ion_floor = u_ion_floor,
        min_consecutive_positive = min_consecutive_positive,
    )
end

function normalize_ion_velocity_policy(policy::AbstractString)
    allowed = ("trim_to_first_positive", "strict_positive")
    policy in allowed || throw(ArgumentError(
        "Unsupported ion velocity policy `$(policy)`. Allowed values: $(join(allowed, ", "))."
    ))
    return policy
end

function ensure_file_exists(path::AbstractString)
    isfile(path) || throw(ArgumentError("Input file does not exist: $(path)"))
    return path
end

function load_source_json(path::AbstractString)
    ensure_file_exists(path)
    return JSON.parsefile(path)
end

function require_dict(container::AbstractDict, key::AbstractString, context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$(key)` in $(context)."))
    value = container[key]
    isa(value, AbstractDict) || throw(ArgumentError("Expected `$(key)` in $(context) to be a JSON object."))
    return value
end

function require_array(container::AbstractDict, key::AbstractString, context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$(key)` in $(context)."))
    value = container[key]
    isa(value, Vector) || throw(ArgumentError("Expected `$(key)` in $(context) to be a JSON array."))
    return value
end

function require_number_array(container, key::AbstractString, context::AbstractString)
    values = require_array(container, key, context)
    return to_float_vector(values; field_name = "$(context).$(key)")
end

function require_number(container, key::AbstractString, context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$(key)` in $(context)."))
    value = container[key]
    value isa Real || throw(ArgumentError("Expected `$(key)` in $(context) to be numeric."))
    return Float64(value)
end

function require_string(container, key::AbstractString, context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$(key)` in $(context)."))
    value = container[key]
    value isa AbstractString || throw(ArgumentError("Expected `$(key)` in $(context) to be a string."))
    return String(value)
end

function to_float_vector(values::Vector; field_name::AbstractString)
    result = Float64[]
    sizehint!(result, length(values))
    for (i, value) in pairs(values)
        if value isa Real
            push!(result, Float64(value))
        else
            throw(ArgumentError("Expected numeric values in $(field_name); found $(typeof(value)) at index $(i)."))
        end
    end
    return result
end

function load_average_block(source::AbstractDict)
    output = require_dict(source, "output", "top-level JSON")
    retcode = get(output, "retcode", nothing)
    retcode == "success" || throw(ArgumentError("HallThruster output retcode must be `success`; got `$(retcode)`."))
    haskey(output, "average") || throw(ArgumentError(
        "HallThruster JSON is missing `output.average`. Re-export with `write_to_json(...; average_start_time=..., save_time_resolved=...)` so the average block is present."
    ))
    average = output["average"]
    isa(average, AbstractDict) || throw(ArgumentError("Expected `output.average` to be a JSON object."))
    return output, average
end

function sorted_joined_keys(dict::AbstractDict)
    return join(sort!(collect(keys(dict))), ", ")
end

function coerce_charge_state(value, context::AbstractString)
    if value isa Integer
        charge_state = Int(value)
    elseif value isa Real && isfinite(value) && isinteger(value)
        charge_state = Int(round(value))
    else
        throw(ArgumentError("Expected integer charge state in $(context); got $(typeof(value))."))
    end
    charge_state >= 1 || throw(ArgumentError("Charge state in $(context) must be >= 1."))
    return charge_state
end

function normalize_species_name(name::AbstractString)
    species = strip(String(name))
    lowercase(species) in ("e", "e-", "electron") && return "E-"
    return species
end

function normalize_propellant_species_name(name::AbstractString)
    token = replace(lowercase(strip(String(name))), r"[^a-z0-9]+" => "")
    token in ("n", "nitrogen", "atomicnitrogen") && return "N"
    token in ("n2", "molecularnitrogen", "dinitrogen") && return "N2"
    return strip(String(name))
end

function map_neutral_species_name(ht_species::AbstractString)
    terra_name = normalize_species_name(ht_species)
    terra_name == "E-" && throw(ArgumentError(
        "Electron species must not be exported as a neutral convective species."
    ))
    return terra_name
end

function map_ion_species_name(ht_species::AbstractString, charge_state::Integer)
    base = map_neutral_species_name(ht_species)
    charge_state == 1 && return "$(base)+"
    return "$(base)$(charge_state)+"
end

base_species_name(species::AbstractString) = replace(String(species), "+" => "")

function compute_point_aligned_dx(z_m::Vector{Float64})
    n = length(z_m)
    n >= 2 || throw(ArgumentError("At least two axial points are required to compute `dx_m`."))

    dx_m = Vector{Float64}(undef, n)
    dx_m[1] = z_m[2] - z_m[1]
    for i in 2:(n - 1)
        dx_m[i] = 0.5 * (z_m[i + 1] - z_m[i - 1])
    end
    dx_m[end] = z_m[end] - z_m[end - 1]
    return dx_m
end

function find_first_consecutive_positive_index(
    u_ion_m_s::Vector{Float64};
    floor::Float64 = 0.0,
    min_consecutive_positive::Int = 3,
)
    n = length(u_ion_m_s)
    min_consecutive_positive >= 1 || throw(ArgumentError("`min_consecutive_positive` must be >= 1."))
    n >= min_consecutive_positive || return nothing

    max_start = n - min_consecutive_positive + 1
    for i in 1:max_start
        all_positive = true
        @inbounds for j in i:(i + min_consecutive_positive - 1)
            value = u_ion_m_s[j]
            if !isfinite(value) || value <= floor
                all_positive = false
                break
            end
        end
        all_positive && return i
    end

    return nothing
end

function trim_arrays_from_index(arrays::Dict{String, Vector{Float64}}, start_idx::Int)
    trimmed = Dict{String, Vector{Float64}}()
    for (key, values) in arrays
        trimmed[key] = values[start_idx:end]
    end
    return trimmed
end

function validate_profile!(
    z_m::Vector{Float64},
    dx_m::Vector{Float64},
    te_K::Vector{Float64},
    species_u_m_s::Dict{String, Vector{Float64}},
)
    n = length(z_m)
    (n == length(dx_m) && n == length(te_K)) ||
        throw(ArgumentError("Profile arrays must all have identical lengths."))

    for (i, value) in pairs(z_m)
        isfinite(value) || throw(ArgumentError("`z_m[$(i)]` must be finite."))
    end
    for i in 2:n
        z_m[i] > z_m[i - 1] || throw(ArgumentError("`z_m` must be strictly increasing; failed at indices $(i - 1) and $(i)."))
    end

    for (name, values) in (("dx_m", dx_m), ("te_K", te_K))
        for (i, value) in pairs(values)
            isfinite(value) || throw(ArgumentError("`$(name)[$(i)]` must be finite."))
            value > 0.0 || throw(ArgumentError("`$(name)[$(i)]` must be strictly positive."))
        end
    end

    isempty(species_u_m_s) && throw(ArgumentError(
        "`species_u_m_s` must contain at least one exported species."
    ))
    for (species_name, velocities) in pairs(species_u_m_s)
        length(velocities) == n || throw(ArgumentError(
            "`species_u_m_s.$(species_name)` must have length $(n); got $(length(velocities))."
        ))
        for (i, value) in pairs(velocities)
            isfinite(value) || throw(ArgumentError(
                "`species_u_m_s.$(species_name)[$(i)]` must be finite."
            ))
            value > 0.0 || throw(ArgumentError(
                "`species_u_m_s.$(species_name)[$(i)]` must be strictly positive."
            ))
        end
    end
end

function validate_optional_diagnostics!(diagnostics::Dict{String, Any}, n::Int)
    for (name, values) in diagnostics
        isa(values, Vector{Float64}) || throw(ArgumentError("Diagnostic `$(name)` must be numeric."))
        length(values) == n || throw(ArgumentError("Diagnostic `$(name)` must have length $(n); got $(length(values))."))
        for (i, value) in pairs(values)
            isfinite(value) || throw(ArgumentError("Diagnostic `$(name)[$(i)]` must be finite."))
        end
    end
end

function load_propellant_thermal_map(source::AbstractDict)
    input = require_dict(source, "input", "top-level JSON")
    config = require_dict(input, "config", "input")
    propellants = require_array(config, "propellants", "input.config")
    isempty(propellants) && throw(ArgumentError(
        "Expected `input.config.propellants` to contain at least one propellant entry."
    ))

    thermal_map = Dict{String, NamedTuple{(:neutral_K, :ion_K), Tuple{Float64, Float64}}}()
    for (i, entry_raw) in pairs(propellants)
        entry_raw isa AbstractDict || throw(ArgumentError(
            "Expected `input.config.propellants[$i]` to be a JSON object."
        ))
        entry = entry_raw::AbstractDict
        species_name = normalize_propellant_species_name(
            require_string(entry, "gas", "input.config.propellants[$i]"))
        haskey(thermal_map, species_name) && throw(ArgumentError(
            "Duplicate propellant thermal entry for species `$(species_name)`."
        ))
        thermal_map[species_name] = (
            neutral_K = require_number(entry, "temperature_K", "input.config.propellants[$i]"),
            ion_K = require_number(entry, "ion_temperature_K", "input.config.propellants[$i]"),
        )
    end

    return thermal_map
end

function build_source_snapshot(
    source::AbstractDict,
    input_path::AbstractString,
    retcode::AbstractString,
)
    output = require_dict(source, "output", "top-level JSON")
    average = require_dict(output, "average", "output")
    available_neutrals = sort!(collect(keys(require_dict(average, "neutrals", "output.average"))))
    available_ions = sort!(collect(keys(require_dict(average, "ions", "output.average"))))

    summary = "Converted HallThruster JSON average block from $(basename(input_path))."

    snapshot = Dict{String, Any}(
        "enabled" => true,
        "source_type" => "hallthruster_json",
        "retcode" => retcode,
        "summary" => summary,
        "input_file" => input_path,
        "has_average_block" => haskey(output, "average"),
        "available_neutral_species" => available_neutrals,
        "available_ion_species" => available_ions,
    )

    if haskey(output, "error")
        snapshot["error"] = output["error"]
    end

    return snapshot
end

function _insert_species_series!(
    destination::Dict{String, Vector{Float64}},
    species_name::AbstractString,
    values::Vector{Float64},
    context::AbstractString,
)
    haskey(destination, species_name) && throw(ArgumentError(
        "Duplicate exported species `$(species_name)` while processing $(context)."
    ))
    destination[String(species_name)] = values
    return nothing
end

function collect_neutral_species_data(average::AbstractDict)
    neutrals = require_dict(average, "neutrals", "output.average")
    velocities = Dict{String, Vector{Float64}}()
    densities = Dict{String, Vector{Float64}}()

    for (ht_species, state_raw) in pairs(neutrals)
        state_raw isa AbstractDict || throw(ArgumentError(
            "Expected `output.average.neutrals.$(ht_species)` to be a JSON object."
        ))
        species_name = map_neutral_species_name(String(ht_species))
        state = state_raw::AbstractDict
        _insert_species_series!(
            velocities,
            species_name,
            require_number_array(state, "u", "output.average.neutrals.$(ht_species)"),
            "neutral species `$(ht_species)`",
        )
        _insert_species_series!(
            densities,
            species_name,
            require_number_array(state, "n", "output.average.neutrals.$(ht_species)"),
            "neutral species `$(ht_species)` density",
        )
    end

    return velocities, densities
end

function collect_ion_species_data(average::AbstractDict)
    ions = require_dict(average, "ions", "output.average")
    velocities = Dict{String, Vector{Float64}}()
    densities = Dict{String, Vector{Float64}}()
    ion_species = String[]
    charge_states = Dict{String, Int}()

    for (ht_species, species_states_raw) in pairs(ions)
        species_states_raw isa Vector || throw(ArgumentError(
            "Expected `output.average.ions.$(ht_species)` to be a JSON array."
        ))
        for (state_index, state_raw) in pairs(species_states_raw)
            state_raw isa AbstractDict || throw(ArgumentError(
                "Ion state entries for `$(ht_species)` must be JSON objects."
            ))
            state = state_raw::AbstractDict
            haskey(state, "Z") || throw(ArgumentError(
                "Missing `Z` in `output.average.ions.$(ht_species)[$(state_index)]`."
            ))
            charge_state = coerce_charge_state(
                state["Z"],
                "output.average.ions.$(ht_species)[$(state_index)].Z",
            )
            species_name = map_ion_species_name(String(ht_species), charge_state)
            push!(ion_species, species_name)
            charge_states[species_name] = charge_state
            _insert_species_series!(
                velocities,
                species_name,
                require_number_array(state, "u",
                    "output.average.ions.$(ht_species)[$(state_index)]"),
                "ion species `$(ht_species)` charge state $(charge_state)",
            )
            _insert_species_series!(
                densities,
                species_name,
                require_number_array(state, "n",
                    "output.average.ions.$(ht_species)[$(state_index)]"),
                "ion species `$(ht_species)` charge state $(charge_state) density",
            )
        end
    end

    return velocities, densities, sort!(ion_species), charge_states
end

function warn_on_charge_balance_mismatch(
    ne_m3::Float64,
    species_n_m3::Dict{String, Vector{Float64}},
    ion_species::Vector{String},
    ion_charge_states::Dict{String, Int};
    source_compact_index::Int,
)
    ion_charge_density = 0.0
    for species_name in ion_species
        ion_charge_density += ion_charge_states[species_name] *
                              species_n_m3[species_name][source_compact_index]
    end

    rel_error = abs(ion_charge_density - ne_m3) / max(ne_m3, eps(Float64))
    if rel_error > 1e-2
        @warn "HallThruster inlet charge balance mismatch at retained compact cell $(source_compact_index)." ne_m3 ion_charge_density rel_error
    end
end

function build_inlet_payload(
    source::AbstractDict,
    species_n_m3::Dict{String, Vector{Float64}},
    ne_m3::Vector{Float64},
    te_K::Vector{Float64},
    neutral_species::Vector{String},
    ion_species::Vector{String},
    ion_charge_states::Dict{String, Int},
)
    source_compact_index = 1
    propellant_thermal_map = load_propellant_thermal_map(source)
    warn_on_charge_balance_mismatch(
        ne_m3[source_compact_index],
        species_n_m3,
        ion_species,
        ion_charge_states;
        source_compact_index = source_compact_index,
    )

    heavy_species = vcat(neutral_species, ion_species)
    species_order = vcat(heavy_species, ["E-"])
    n_total_m3 = ne_m3[source_compact_index]
    heavy_n_total = 0.0
    heavy_temperature_numerator = 0.0

    for species_name in heavy_species
        density = species_n_m3[species_name][source_compact_index]
        n_total_m3 += density
        heavy_n_total += density

        base_species = base_species_name(species_name)
        haskey(propellant_thermal_map, base_species) || throw(ArgumentError(
            "Missing propellant thermal data for heavy species `$(base_species)`."
        ))
        temperatures = propellant_thermal_map[base_species]
        heavy_temperature_numerator += density * (
            occursin("+", species_name) ? temperatures.ion_K : temperatures.neutral_K
        )
    end

    heavy_n_total > 0.0 || throw(ArgumentError(
        "Cannot build inlet payload: total heavy-species number density is zero."
    ))
    n_total_m3 > 0.0 || throw(ArgumentError(
        "Cannot build inlet payload: total number density is zero."
    ))

    mole_fractions = Float64[
        species_n_m3[species_name][source_compact_index] / n_total_m3 for
        species_name in heavy_species
    ]
    push!(mole_fractions, ne_m3[source_compact_index] / n_total_m3)

    heavy_temperature_K = heavy_temperature_numerator / heavy_n_total
    electron_temperature_K = te_K[source_compact_index]

    return Dict{String, Any}(
        "composition" => Dict{String, Any}(
            "species" => species_order,
            "mole_fractions" => mole_fractions,
            "total_number_density_m3" => n_total_m3,
        ),
        "thermal" => Dict{String, Any}(
            "Tt_K" => heavy_temperature_K,
            "Tv_K" => heavy_temperature_K,
            "Tee_K" => heavy_temperature_K,
            "Te_K" => electron_temperature_K,
        ),
        "source_compact_index" => source_compact_index,
    )
end

function determine_trim_start_index(
    ion_velocity_arrays::Dict{String, Vector{Float64}},
    options,
)
    ion_velocity_policy = normalize_ion_velocity_policy(options.ion_velocity_policy)
    if ion_velocity_policy == "strict_positive" || isempty(ion_velocity_arrays)
        return 1
    end

    trim_indices = Int[]
    for species_name in sort!(collect(keys(ion_velocity_arrays)))
        idx = find_first_consecutive_positive_index(
            ion_velocity_arrays[species_name];
            floor = options.u_ion_floor,
            min_consecutive_positive = options.min_consecutive_positive,
        )
        idx === nothing && throw(ArgumentError(
            "Unable to find a valid trim point for ion species `$(species_name)` with policy `trim_to_first_positive` (u_ion_floor=$(options.u_ion_floor), min_consecutive_positive=$(options.min_consecutive_positive))."
        ))
        push!(trim_indices, idx)
    end

    return maximum(trim_indices)
end

function build_chain_profile(
    source::AbstractDict,
    options,
)
    output, average = load_average_block(source)

    z_m = require_number_array(average, "z", "output.average")
    neutral_velocities, neutral_densities = collect_neutral_species_data(average)
    ion_velocities, ion_densities, exported_ion_species, ion_charge_states = collect_ion_species_data(average)

    species_u_m_s = merge(copy(neutral_velocities), ion_velocities)
    species_n_m3 = merge(copy(neutral_densities), ion_densities)

    tev = require_number_array(average, "Tev", "output.average")
    ne_m3 = require_number_array(average, "ne", "output.average")
    channel_area_m2 = require_number_array(average, "channel_area", "output.average")
    potential_V = require_number_array(average, "potential", "output.average")
    electric_field_V_m = require_number_array(average, "E", "output.average")

    original_point_count = length(z_m)
    all_arrays = Dict{String, Vector{Float64}}(
        "z_m" => z_m,
        "tev_eV" => tev,
        "ne_m3" => ne_m3,
        "channel_area_m2" => channel_area_m2,
        "potential_V" => potential_V,
        "electric_field_V_m" => electric_field_V_m,
    )
    for (species_name, values) in pairs(species_u_m_s)
        all_arrays["species_u_m_s.$(species_name)"] = values
    end
    for (species_name, values) in pairs(species_n_m3)
        all_arrays["species_n_m3.$(species_name)"] = values
    end

    for (name, values) in all_arrays
        length(values) == original_point_count || throw(ArgumentError(
            "HallThruster average arrays must have identical lengths. `$(name)` has length $(length(values)), expected $(original_point_count)."
        ))
    end

    trim_start_index = determine_trim_start_index(ion_velocities, options)
    if trim_start_index > 1
        all_arrays = trim_arrays_from_index(all_arrays, trim_start_index)
    end

    z_m = all_arrays["z_m"]
    tev = all_arrays["tev_eV"]
    ne_m3 = all_arrays["ne_m3"]
    channel_area_m2 = all_arrays["channel_area_m2"]
    potential_V = all_arrays["potential_V"]
    electric_field_V_m = all_arrays["electric_field_V_m"]

    neutral_species = sort!(collect(keys(neutral_velocities)))
    species_u_m_s = Dict{String, Vector{Float64}}()
    for species_name in neutral_species
        species_u_m_s[species_name] = all_arrays["species_u_m_s.$(species_name)"]
    end
    for species_name in exported_ion_species
        species_u_m_s[species_name] = all_arrays["species_u_m_s.$(species_name)"]
    end

    species_n_m3 = Dict{String, Vector{Float64}}()
    for species_name in neutral_species
        species_n_m3[species_name] = all_arrays["species_n_m3.$(species_name)"]
    end
    for species_name in exported_ion_species
        species_n_m3[species_name] = all_arrays["species_n_m3.$(species_name)"]
    end

    te_K = tev .* EV_TO_K
    dx_m = compute_point_aligned_dx(z_m)

    validate_profile!(z_m, dx_m, te_K, species_u_m_s)

    n_total_m3 = copy(ne_m3)
    for values in values(species_n_m3)
        n_total_m3 .+= values
    end

    diagnostics = Dict{String, Any}(
        "ne_m3" => ne_m3,
        "n_total_m3" => n_total_m3,
        "channel_area_m2" => channel_area_m2,
        "potential_V" => potential_V,
        "electric_field_V_m" => electric_field_V_m,
    )
    validate_optional_diagnostics!(diagnostics, length(z_m))

    exported_species = sort!(collect(keys(species_u_m_s)))
    generator = Dict{String, Any}(
        "tool" => "hallthruster_jl_export_chain_profile",
        "tool_version" => SCRIPT_VERSION,
        "created_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
    )

    selection = Dict{String, Any}(
        "average_start_time_s" => options.average_start_time_s,
        "exported_species" => exported_species,
        "ion_velocity_policy" => normalize_ion_velocity_policy(options.ion_velocity_policy),
        "u_ion_floor" => options.u_ion_floor,
        "min_consecutive_positive" => options.min_consecutive_positive,
        "trim_start_index" => trim_start_index,
        "trim_start_z_m" => z_m[1],
        "trimmed_point_count" => trim_start_index - 1,
        "original_point_count" => original_point_count,
    )

    inlet = build_inlet_payload(
        source,
        species_n_m3,
        ne_m3,
        te_K,
        neutral_species,
        exported_ion_species,
        ion_charge_states,
    )

    profile = Dict{String, Any}(
        "z_m" => z_m,
        "dx_m" => dx_m,
        "te_K" => te_K,
        "species_u_m_s" => Dict{String, Any}(
            name => species_u_m_s[name] for name in exported_species
        ),
    )

    payload = Dict{String, Any}(
        "schema_version" => SCHEMA_VERSION,
        "generator" => generator,
        "selection" => selection,
        "inlet" => inlet,
        "profile" => profile,
        "diagnostics" => diagnostics,
    )

    if options.write_source_snapshot
        payload["source_snapshot"] = build_source_snapshot(source, options.hallthruster_input, output["retcode"])
    end

    return payload
end

function output_path_for_case(terra_case_path::AbstractString)
    return joinpath(terra_case_path, "input", DEFAULT_OUTPUT_BASENAME)
end

function write_chain_profile(payload::Dict{String, Any}, terra_case_path::AbstractString)
    output_path = output_path_for_case(terra_case_path)
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        JSON.print(io, payload, 2)
        write(io, '\n')
    end
    return output_path
end

function export_chain_profile(args::Vector{String} = ARGS)
    options = parse_options(args)
    options === nothing && return nothing
    options.average_start_time_s >= 0.0 || throw(ArgumentError("`--average_start_time` must be non-negative."))
    options.min_consecutive_positive >= 1 || throw(ArgumentError("`--min_consecutive_positive` must be >= 1."))

    source = load_source_json(options.hallthruster_input)
    payload = build_chain_profile(source, options)
    output_path = write_chain_profile(payload, options.terra_case_path)

    println("HallThruster chain profile export completed.")
    println("input_file: ", options.hallthruster_input)
    println("output_file: ", output_path)
    println("average_start_time_s: ", options.average_start_time_s)
    println("exported_species: ", join(payload["selection"]["exported_species"], ", "))
    println("ion_velocity_policy: ", options.ion_velocity_policy)
    println("u_ion_floor: ", options.u_ion_floor)
    println("min_consecutive_positive: ", options.min_consecutive_positive)
    println("trim_start_index: ", payload["selection"]["trim_start_index"])
    println("trimmed_point_count: ", payload["selection"]["trimmed_point_count"])
    println("write_source_snapshot: ", options.write_source_snapshot)

    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_chain_profile()
end
