"""
$(SIGNATURES)

Residence-time (CSTR / flow-through) model controls.
"""
struct ResidenceTimeConfig
    L::Float64
    U_species::Dict{String, Float64}
    U_energy::Union{Nothing, Float64}
    inlet_reactor::Union{Nothing, ReactorConfig}

    function ResidenceTimeConfig(L, U_species, U_energy = nothing, inlet_reactor = nothing)
        if !isfinite(L) || L <= 0
            throw(ArgumentError("ResidenceTimeConfig: L must be finite and positive (got $L)."))
        end
        u_species_dict = Dict{String, Float64}()
        for (name, value) in pairs(U_species)
            name_str = String(name)
            isempty(strip(name_str)) &&
                throw(ArgumentError("ResidenceTimeConfig: species keys in U_species must be non-empty."))
            if !isfinite(value) || value <= 0
                throw(ArgumentError("ResidenceTimeConfig: U_species[$name_str] must be finite and positive (got $value)."))
            end
            u_species_dict[name_str] = Float64(value)
        end
        if U_energy !== nothing && (!isfinite(U_energy) || U_energy <= 0)
            throw(ArgumentError("ResidenceTimeConfig: U_energy must be finite and positive (got $U_energy)."))
        end

        return new(Float64(L),
                   u_species_dict,
                   U_energy === nothing ? nothing : Float64(U_energy),
                   _coerce_residence_time_inlet_reactor(inlet_reactor))
    end
end

function ResidenceTimeConfig(; L::Real = 1.0,
                             U_species::AbstractDict = Dict{String, Float64}(),
                             U_energy::Union{Nothing, Real} = nothing,
                             inlet_reactor = nothing)
    return ResidenceTimeConfig(L, U_species, U_energy, inlet_reactor)
end

function _coerce_residence_time_inlet_reactor(inlet)
    if inlet === nothing
        return nothing
    elseif inlet isa ReactorConfig
        return inlet
    end
    throw(ArgumentError("ResidenceTimeConfig: inlet_reactor must be `ReactorConfig` or `nothing`."))
end
