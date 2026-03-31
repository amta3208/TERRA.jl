"""
$(SIGNATURES)

Abstract root for per-species wall-interaction closures.
"""
abstract type SpeciesWallModel end

@inline function _coerce_wall_parameter(value, field_name::AbstractString)
    value_f64 = Float64(value)
    isfinite(value_f64) && value_f64 >= 0.0 ||
        throw(ArgumentError("$(field_name) must be finite and nonnegative."))
    return value_f64
end

function _coerce_wall_products(products::AbstractDict)
    product_dict = Dict{String, Float64}()
    for (name, value) in pairs(products)
        name_str = String(name)
        isempty(strip(name_str)) &&
            throw(ArgumentError("Wall model product keys must be non-empty."))
        product_dict[name_str] = _coerce_wall_parameter(value,
                                                        "Wall model product `$name_str`")
    end
    return product_dict
end

"""
$(SIGNATURES)

Ion wall-loss closure using a Bohm-gap style rate.
"""
struct IonNeutralizationWallModel <: SpeciesWallModel
    bohm_scale::Float64
    products::Dict{String, Float64}

    function IonNeutralizationWallModel(;
                                        bohm_scale::Real = 1.0,
                                        products::AbstractDict = Dict{String, Float64}())
        return new(_coerce_wall_parameter(bohm_scale, "IonNeutralizationWallModel.bohm_scale"),
                   _coerce_wall_products(products))
    end
end

"""
$(SIGNATURES)

Neutral wall-loss closure using a ballistic sticking rate.
"""
struct BallisticNeutralRecombinationWallModel <: SpeciesWallModel
    gamma::Float64
    products::Dict{String, Float64}

    function BallisticNeutralRecombinationWallModel(;
                                                    gamma::Real,
                                                    products::AbstractDict = Dict{String, Float64}())
        return new(_coerce_wall_parameter(gamma,
                                          "BallisticNeutralRecombinationWallModel.gamma"),
                   _coerce_wall_products(products))
    end
end

"""
$(SIGNATURES)

Neutral wall-loss closure using a constant rate.
"""
struct ConstantNeutralRecombinationWallModel <: SpeciesWallModel
    k_wall_1_s::Float64
    products::Dict{String, Float64}

    function ConstantNeutralRecombinationWallModel(;
                                                   k_wall_1_s::Real,
                                                   products::AbstractDict = Dict{String, Float64}())
        return new(_coerce_wall_parameter(k_wall_1_s,
                                          "ConstantNeutralRecombinationWallModel.k_wall_1_s"),
                   _coerce_wall_products(products))
    end
end

@inline _wall_model_products(model::SpeciesWallModel) = model.products
@inline _wall_model_type_name(::IonNeutralizationWallModel) = "IonNeutralizationWallModel"
@inline _wall_model_type_name(::BallisticNeutralRecombinationWallModel) =
    "BallisticNeutralRecombinationWallModel"
@inline _wall_model_type_name(::ConstantNeutralRecombinationWallModel) =
    "ConstantNeutralRecombinationWallModel"

@inline _wall_parameter_dict(model::IonNeutralizationWallModel) = Dict("bohm_scale" => model.bohm_scale)
@inline _wall_parameter_dict(model::BallisticNeutralRecombinationWallModel) = Dict("gamma" => model.gamma)
@inline _wall_parameter_dict(model::ConstantNeutralRecombinationWallModel) =
    Dict("k_wall_1_s" => model.k_wall_1_s)

function _coerce_wall_species_models(species_models::AbstractDict)
    species_model_dict = Dict{String, SpeciesWallModel}()
    for (name, model) in pairs(species_models)
        name_str = String(name)
        isempty(strip(name_str)) &&
            throw(ArgumentError("WallLossConfig: species_models keys must be non-empty."))
        model isa SpeciesWallModel ||
            throw(ArgumentError("WallLossConfig: species_models[$name_str] must be a SpeciesWallModel."))
        species_model_dict[name_str] = model
    end
    return species_model_dict
end

"""
$(SIGNATURES)

Wrapper-managed wall-loss configuration.
"""
struct WallLossConfig
    species_models::Dict{String, SpeciesWallModel}

    function WallLossConfig(; species_models::AbstractDict = Dict{String, SpeciesWallModel}())
        species_model_dict = _coerce_wall_species_models(species_models)
        return new(species_model_dict)
    end
end
