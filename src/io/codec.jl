json_value(::Nothing) = nothing
json_value(value::Bool) = value
json_value(value::Real) = value
json_value(value::AbstractString) = String(value)
json_value(value::Symbol) = String(value)

function json_value(value::NamedTuple)
    dict = Dict{String, Any}()
    for (key, entry) in pairs(value)
        dict[String(key)] = json_value(entry)
    end
    return dict
end

function json_value(value::AbstractDict)
    dict = Dict{String, Any}()
    for (key, entry) in pairs(value)
        dict[String(key)] = json_value(entry)
    end
    return dict
end

json_value(value::AbstractVector) = [json_value(entry) for entry in value]
json_value(value) = string(value)

function read_json_value(value::AbstractDict)
    dict = Dict{String, Any}()
    for (key, entry) in pairs(value)
        dict[String(key)] = read_json_value(entry)
    end
    return dict
end

read_json_value(value::AbstractVector) = Any[read_json_value(entry) for entry in value]
read_json_value(value) = value

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

function required_number(container::AbstractDict, key::AbstractString,
                         context::AbstractString)
    haskey(container, key) || throw(ArgumentError("Missing `$key` in $context."))
    value = container[key]
    value isa Real || throw(ArgumentError("Expected `$key` in $context to be numeric."))
    return Float64(value)
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

optional_int(value, key::AbstractString, context::AbstractString) =
    value === nothing ? nothing : int_value(value, key, context)

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

optional_float_array(container::AbstractDict, key::AbstractString,
                     context::AbstractString) =
    haskey(container, key) ? float_array(container, key, context) : nothing

function string_array(container::AbstractDict, key::AbstractString,
                      context::AbstractString)
    values = required_array(container, key, context)
    result = String[]
    sizehint!(result, length(values))
    for (i, value) in pairs(values)
        value isa AbstractString ||
            throw(ArgumentError("Expected string values in `$context.$key`; got $(typeof(value)) at index $i."))
        push!(result, String(value))
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

function float_array_dict(container::AbstractDict,
                          key::AbstractString,
                          context::AbstractString)
    raw = required_dict(container, key, context)
    result = Dict{String, Vector{Float64}}()
    for (name, values) in pairs(raw)
        name_str = String(name)
        values isa AbstractVector ||
            throw(ArgumentError("Expected `$context.$key.$name_str` to be an array."))
        values_vec = Float64[]
        sizehint!(values_vec, length(values))
        for (i, value) in pairs(values)
            value isa Real ||
                throw(ArgumentError("Expected numeric values in `$context.$key.$name_str`; got $(typeof(value)) at index $i."))
            push!(values_vec, Float64(value))
        end
        result[name_str] = values_vec
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
