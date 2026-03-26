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
