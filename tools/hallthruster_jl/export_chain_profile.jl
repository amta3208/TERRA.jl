using Dates
using JSON

const SCHEMA_VERSION = "terra_chain_profile_v1"
const SCRIPT_VERSION = "0.2.0"
const EV_TO_K = 11604.518121550082
const DEFAULT_OUTPUT_BASENAME = "chain_profile_v1.json"

function usage()
    return """
    Usage:
      julia --project=tools/hallthruster_jl tools/hallthruster_jl/export_chain_profile.jl <hallthruster_input.json> <terra_case_path> --average_start_time=<seconds> --neutral_species=<symbol> --ion_species=<symbol> --ion_charge_state=<int> [--write_source_snapshot=<true|false>] [--ion_velocity_policy=<trim_to_first_positive|strict_positive>] [--u_ion_floor=<float>] [--min_consecutive_positive=<int>]

    Notes:
      - The HallThruster JSON input must already contain an `output.average` block.
      - For JSON inputs, `--average_start_time` is recorded as metadata only; this script does not recompute the time average.
      - The output is written to `<terra_case_path>/input/chain_profile_v1.json`.
      - Default ion-velocity handling is `trim_to_first_positive` with `u_ion_floor=0.0` and `min_consecutive_positive=3`.
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

    required_option_keys = (
        "average_start_time",
        "neutral_species",
        "ion_species",
        "ion_charge_state",
    )

    for key in required_option_keys
        haskey(options, key) || throw(ArgumentError("Missing required option --$(key)=...\n\n" * usage()))
    end

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
        neutral_species = options["neutral_species"],
        ion_species = options["ion_species"],
        ion_charge_state = parse(Int, options["ion_charge_state"]),
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

function select_neutral_state(average::AbstractDict, neutral_species::AbstractString)
    neutrals = require_dict(average, "neutrals", "output.average")
    haskey(neutrals, neutral_species) || throw(ArgumentError(
        "Neutral species `$(neutral_species)` not found. Available neutrals: $(sorted_joined_keys(neutrals))."
    ))
    state = neutrals[neutral_species]
    isa(state, AbstractDict) || throw(ArgumentError("Expected `output.average.neutrals.$(neutral_species)` to be a JSON object."))
    return state
end

function select_ion_state(average::AbstractDict, ion_species::AbstractString, ion_charge_state::Int)
    ions = require_dict(average, "ions", "output.average")
    haskey(ions, ion_species) || throw(ArgumentError(
        "Ion species `$(ion_species)` not found. Available ion species: $(sorted_joined_keys(ions))."
    ))
    species_states = ions[ion_species]
    isa(species_states, Vector) || throw(ArgumentError("Expected `output.average.ions.$(ion_species)` to be a JSON array."))

    for state in species_states
        isa(state, AbstractDict) || throw(ArgumentError("Ion state entries for `$(ion_species)` must be JSON objects."))
        state_charge = get(state, "Z", nothing)
        state_charge == ion_charge_state && return state
    end

    available_charge_states = sort!(Int[
        Int(state["Z"]) for state in species_states if isa(state, AbstractDict) && haskey(state, "Z") && state["Z"] isa Real
    ])
    throw(ArgumentError(
        "Charge state Z=$(ion_charge_state) not found for ion species `$(ion_species)`. Available charge states: $(join(available_charge_states, ", "))."
    ))
end

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

function count_nonpositive(values::Vector{Float64}; floor::Float64 = 0.0)
    count = 0
    for value in values
        value > floor || (count += 1)
    end
    return count
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
    u_neutral_m_s::Vector{Float64},
    u_ion_m_s::Vector{Float64},
)
    n = length(z_m)
    (
        n == length(dx_m) &&
        n == length(te_K) &&
        n == length(u_neutral_m_s) &&
        n == length(u_ion_m_s)
    ) ||
        throw(ArgumentError("Profile arrays must all have identical lengths."))

    for (i, value) in pairs(z_m)
        isfinite(value) || throw(ArgumentError("`z_m[$(i)]` must be finite."))
    end
    for i in 2:n
        z_m[i] > z_m[i - 1] || throw(ArgumentError("`z_m` must be strictly increasing; failed at indices $(i - 1) and $(i)."))
    end

    for (name, values) in (
        ("dx_m", dx_m),
        ("te_K", te_K),
        ("u_neutral_m_s", u_neutral_m_s),
    )
        for (i, value) in pairs(values)
            isfinite(value) || throw(ArgumentError("`$(name)[$(i)]` must be finite."))
            value > 0.0 || throw(ArgumentError("`$(name)[$(i)]` must be strictly positive."))
        end
    end

    for (i, value) in pairs(u_ion_m_s)
        isfinite(value) || throw(ArgumentError("`u_ion_m_s[$(i)]` must be finite."))
    end

    nonpositive_ion_count = count_nonpositive(u_ion_m_s)
    if nonpositive_ion_count > 0
        throw(ArgumentError(
            "`u_ion_m_s` must be strictly positive for the Phase 1 contract. Found $(nonpositive_ion_count) non-positive values; minimum = $(minimum(u_ion_m_s))."
        ))
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

function build_chain_profile(
    source::AbstractDict,
    options,
)
    output, average = load_average_block(source)

    z_m = require_number_array(average, "z", "output.average")
    neutral_state = select_neutral_state(average, options.neutral_species)
    ion_state = select_ion_state(average, options.ion_species, options.ion_charge_state)

    u_neutral_m_s = require_number_array(neutral_state, "u", "output.average.neutrals.$(options.neutral_species)")
    n_neutral_m3 = require_number_array(neutral_state, "n", "output.average.neutrals.$(options.neutral_species)")
    u_ion_m_s = require_number_array(ion_state, "u", "output.average.ions.$(options.ion_species)")
    n_ion_m3 = require_number_array(ion_state, "n", "output.average.ions.$(options.ion_species)")

    tev = require_number_array(average, "Tev", "output.average")
    ne_m3 = require_number_array(average, "ne", "output.average")
    channel_area_m2 = require_number_array(average, "channel_area", "output.average")
    potential_V = require_number_array(average, "potential", "output.average")
    electric_field_V_m = require_number_array(average, "E", "output.average")

    original_point_count = length(z_m)
    all_arrays = Dict{String, Vector{Float64}}(
        "z_m" => z_m,
        "u_neutral_m_s" => u_neutral_m_s,
        "n_neutral_m3" => n_neutral_m3,
        "u_ion_m_s" => u_ion_m_s,
        "n_ion_m3" => n_ion_m3,
        "tev_eV" => tev,
        "ne_m3" => ne_m3,
        "channel_area_m2" => channel_area_m2,
        "potential_V" => potential_V,
        "electric_field_V_m" => electric_field_V_m,
    )

    for (name, values) in all_arrays
        length(values) == original_point_count || throw(ArgumentError(
            "HallThruster average arrays must have identical lengths. `$(name)` has length $(length(values)), expected $(original_point_count)."
        ))
    end

    ion_velocity_policy = normalize_ion_velocity_policy(options.ion_velocity_policy)
    trim_start_index = 1

    if ion_velocity_policy == "trim_to_first_positive"
        idx = find_first_consecutive_positive_index(
            all_arrays["u_ion_m_s"];
            floor = options.u_ion_floor,
            min_consecutive_positive = options.min_consecutive_positive,
        )
        idx === nothing && throw(ArgumentError(
            "Unable to find a valid trim point for `u_ion_m_s` with policy `trim_to_first_positive` (u_ion_floor=$(options.u_ion_floor), min_consecutive_positive=$(options.min_consecutive_positive))."
        ))
        trim_start_index = idx
        all_arrays = trim_arrays_from_index(all_arrays, trim_start_index)
    elseif ion_velocity_policy == "strict_positive"
        trim_start_index = 1
    end

    z_m = all_arrays["z_m"]
    u_neutral_m_s = all_arrays["u_neutral_m_s"]
    n_neutral_m3 = all_arrays["n_neutral_m3"]
    u_ion_m_s = all_arrays["u_ion_m_s"]
    n_ion_m3 = all_arrays["n_ion_m3"]
    tev = all_arrays["tev_eV"]
    ne_m3 = all_arrays["ne_m3"]
    channel_area_m2 = all_arrays["channel_area_m2"]
    potential_V = all_arrays["potential_V"]
    electric_field_V_m = all_arrays["electric_field_V_m"]

    te_K = tev .* EV_TO_K
    dx_m = compute_point_aligned_dx(z_m)

    validate_profile!(z_m, dx_m, te_K, u_neutral_m_s, u_ion_m_s)

    diagnostics = Dict{String, Any}(
        "n_neutral_m3" => n_neutral_m3,
        "n_ion_m3" => n_ion_m3,
        "ne_m3" => ne_m3,
        "n_total_m3" => n_neutral_m3 .+ n_ion_m3 .+ ne_m3,
        "channel_area_m2" => channel_area_m2,
        "potential_V" => potential_V,
        "electric_field_V_m" => electric_field_V_m,
    )
    validate_optional_diagnostics!(diagnostics, length(z_m))

    generator = Dict{String, Any}(
        "tool" => "hallthruster_jl_export_chain_profile",
        "tool_version" => SCRIPT_VERSION,
        "created_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
    )

    selection = Dict{String, Any}(
        "neutral_species" => options.neutral_species,
        "ion_species" => options.ion_species,
        "ion_charge_state" => options.ion_charge_state,
        "average_start_time_s" => options.average_start_time_s,
        "ion_velocity_policy" => ion_velocity_policy,
        "u_ion_floor" => options.u_ion_floor,
        "min_consecutive_positive" => options.min_consecutive_positive,
        "trim_start_index" => trim_start_index,
        "trim_start_z_m" => z_m[1],
        "trimmed_point_count" => trim_start_index - 1,
        "original_point_count" => original_point_count,
    )

    profile = Dict{String, Any}(
        "z_m" => z_m,
        "dx_m" => dx_m,
        "te_K" => te_K,
        "u_neutral_m_s" => u_neutral_m_s,
        "u_ion_m_s" => u_ion_m_s,
    )

    payload = Dict{String, Any}(
        "schema_version" => SCHEMA_VERSION,
        "generator" => generator,
        "selection" => selection,
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
    println("neutral_species: ", options.neutral_species)
    println("ion_species: ", options.ion_species)
    println("ion_charge_state: ", options.ion_charge_state)
    println("average_start_time_s: ", options.average_start_time_s)
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
