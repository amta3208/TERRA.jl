"""
$(SIGNATURES)

Residence-time (CSTR / flow-through) model controls.
"""
struct ResidenceTimeConfig
    enabled::Bool
    L::Float64
    U_species::Dict{String, Float64}
    U_energy::Union{Nothing, Float64}
    inlet_reactor::Union{Nothing, ReactorConfig}

    function ResidenceTimeConfig(enabled::Bool, L, U_species,
                                 U_energy = nothing, inlet_reactor = nothing)
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
        inlet_reactor_val = _coerce_residence_time_inlet_reactor(inlet_reactor)

        return new(enabled,
                   Float64(L),
                   u_species_dict,
                   U_energy === nothing ? nothing : Float64(U_energy),
                   inlet_reactor_val)
    end
end

function ResidenceTimeConfig(L, U_species, U_energy = nothing, inlet_config = nothing)
    ResidenceTimeConfig(true, L, U_species, U_energy, inlet_config)
end

function ResidenceTimeConfig(; enabled::Bool = true,
                             L::Real = 1.0,
                             U_species::AbstractDict = Dict{String, Float64}(),
                             U_energy::Union{Nothing, Real} = nothing,
                             inlet_reactor = nothing,
                             inlet_config = nothing)
    if inlet_reactor !== nothing && inlet_config !== nothing
        throw(ArgumentError("ResidenceTimeConfig: provide only one of `inlet_reactor` or legacy `inlet_config`."))
    end
    inlet_value = inlet_reactor === nothing ? inlet_config : inlet_reactor
    return ResidenceTimeConfig(enabled, L, U_species, U_energy, inlet_value)
end

function _coerce_residence_time_inlet_reactor(inlet)
    if inlet === nothing
        return nothing
    elseif inlet isa ReactorConfig
        return inlet
    else
        throw(ArgumentError("ResidenceTimeConfig: inlet_reactor/inlet_config must be `ReactorConfig`, `Config`, or `nothing`."))
    end
end

"""
$(SIGNATURES)

Per-species wall-interaction closure for wrapper-managed wall losses.
"""
struct SpeciesWallModel
    class::Symbol
    rate_model::Symbol
    parameters::Dict{String, Float64}
    products::Dict{String, Float64}

    function SpeciesWallModel(; class::Symbol,
                              rate_model::Symbol,
                              parameters::AbstractDict = Dict{String, Float64}(),
                              products::AbstractDict = Dict{String, Float64}())
        class in (:ion_neutralization, :neutral_recombination, :electronic_quench) ||
            throw(ArgumentError("SpeciesWallModel: class must be one of `:ion_neutralization`, `:neutral_recombination`, or `:electronic_quench`."))
        rate_model in (:bohm_gap, :ballistic_sticking, :constant) ||
            throw(ArgumentError("SpeciesWallModel: rate_model must be one of `:bohm_gap`, `:ballistic_sticking`, or `:constant`."))

        parameter_dict = Dict{String, Float64}()
        for (name, value) in pairs(parameters)
            name_str = String(name)
            isempty(strip(name_str)) &&
                throw(ArgumentError("SpeciesWallModel: parameter keys must be non-empty."))
            value_f64 = Float64(value)
            isfinite(value_f64) && value_f64 >= 0.0 ||
                throw(ArgumentError("SpeciesWallModel: parameters[$name_str] must be finite and nonnegative."))
            parameter_dict[name_str] = value_f64
        end

        product_dict = Dict{String, Float64}()
        for (name, value) in pairs(products)
            name_str = String(name)
            isempty(strip(name_str)) &&
                throw(ArgumentError("SpeciesWallModel: product keys must be non-empty."))
            value_f64 = Float64(value)
            isfinite(value_f64) && value_f64 >= 0.0 ||
                throw(ArgumentError("SpeciesWallModel: products[$name_str] must be finite and nonnegative."))
            product_dict[name_str] = value_f64
        end

        return new(class, rate_model, parameter_dict, product_dict)
    end
end

"""
$(SIGNATURES)

Wrapper-managed wall-loss configuration.
"""
struct WallLossConfig
    enabled::Bool
    use_ion_losses::Bool
    use_neutral_recombination::Bool
    use_electronic_quenching::Bool
    species_models::Dict{String, SpeciesWallModel}

    function WallLossConfig(; enabled::Bool = true,
                            use_ion_losses::Bool = true,
                            use_neutral_recombination::Bool = false,
                            use_electronic_quenching::Bool = false,
                            species_models::AbstractDict = Dict{String, SpeciesWallModel}())
        species_model_dict = Dict{String, SpeciesWallModel}()
        for (name, model) in pairs(species_models)
            name_str = String(name)
            isempty(strip(name_str)) &&
                throw(ArgumentError("WallLossConfig: species_models keys must be non-empty."))
            model isa SpeciesWallModel ||
                throw(ArgumentError("WallLossConfig: species_models[$name_str] must be a SpeciesWallModel."))
            species_model_dict[name_str] = model
        end

        return new(enabled,
                   use_ion_losses,
                   use_neutral_recombination,
                   use_electronic_quenching,
                   species_model_dict)
    end
end

"""
$(SIGNATURES)

Wrapper-managed additive source-term configuration.
"""
struct SourceTermsConfig
    residence_time::Union{Nothing, ResidenceTimeConfig}
    wall_losses::Union{Nothing, WallLossConfig}

    function SourceTermsConfig(;
                               residence_time::Union{Nothing, ResidenceTimeConfig} = nothing,
                               wall_losses::Union{Nothing, WallLossConfig} = nothing)
        return new(residence_time, wall_losses)
    end
end
