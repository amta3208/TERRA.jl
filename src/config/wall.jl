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

function _coerce_compat_wall_parameters(parameters::AbstractDict)
    parameter_dict = Dict{String, Float64}()
    for (name, value) in pairs(parameters)
        name_str = String(name)
        isempty(strip(name_str)) &&
            throw(ArgumentError("SpeciesWallModel: parameter keys must be non-empty."))
        parameter_dict[name_str] = _coerce_wall_parameter(value,
                                                          "SpeciesWallModel parameter `$name_str`")
    end
    return parameter_dict
end

function _compat_parameter_value(parameters::Dict{String, Float64},
                                 field_name::AbstractString;
                                 required::Bool,
                                 default::Union{Nothing, Float64} = nothing)
    allowed = Set([String(field_name)])
    extra = sort!(String[name for name in keys(parameters) if !(name in allowed)])
    isempty(extra) ||
        throw(ArgumentError("SpeciesWallModel: unexpected parameters for `$field_name`: $(join(extra, ", "))."))

    if haskey(parameters, field_name)
        return parameters[field_name]
    end
    required &&
        throw(ArgumentError("SpeciesWallModel: missing required parameter `$field_name`."))
    return default
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

@inline _wall_model_class(::IonNeutralizationWallModel) = :ion_neutralization
@inline _wall_model_class(::BallisticNeutralRecombinationWallModel) = :neutral_recombination
@inline _wall_model_class(::ConstantNeutralRecombinationWallModel) = :neutral_recombination

@inline _wall_model_rate_model(::IonNeutralizationWallModel) = :bohm_gap
@inline _wall_model_rate_model(::BallisticNeutralRecombinationWallModel) = :ballistic_sticking
@inline _wall_model_rate_model(::ConstantNeutralRecombinationWallModel) = :constant

@inline _wall_model_products(model::SpeciesWallModel) = model.products

@inline _wall_parameter_dict(model::IonNeutralizationWallModel) = Dict("bohm_scale" => model.bohm_scale)
@inline _wall_parameter_dict(model::BallisticNeutralRecombinationWallModel) = Dict("gamma" => model.gamma)
@inline _wall_parameter_dict(model::ConstantNeutralRecombinationWallModel) =
    Dict("k_wall_1_s" => model.k_wall_1_s)

"""
$(SIGNATURES)

Compatibility constructor for the legacy symbol-driven wall model API.
"""
function SpeciesWallModel(; class::Symbol,
                          rate_model::Symbol,
                          parameters::AbstractDict = Dict{String, Float64}(),
                          products::AbstractDict = Dict{String, Float64}())
    Base.depwarn("`SpeciesWallModel(; class=..., rate_model=...)` is deprecated; construct a concrete wall model type instead.",
                 :SpeciesWallModel)

    parameter_dict = _coerce_compat_wall_parameters(parameters)

    if class == :ion_neutralization
        rate_model == :bohm_gap ||
            throw(ArgumentError("SpeciesWallModel: `:ion_neutralization` supports only `:bohm_gap`."))
        bohm_scale = _compat_parameter_value(parameter_dict, "bohm_scale";
                                             required = false,
                                             default = 1.0)
        return IonNeutralizationWallModel(; bohm_scale = bohm_scale, products = products)
    end

    if class == :neutral_recombination
        if rate_model == :ballistic_sticking
            gamma = _compat_parameter_value(parameter_dict, "gamma"; required = true)
            return BallisticNeutralRecombinationWallModel(; gamma = gamma,
                                                          products = products)
        elseif rate_model == :constant
            k_wall_1_s = _compat_parameter_value(parameter_dict, "k_wall_1_s";
                                                 required = true)
            return ConstantNeutralRecombinationWallModel(; k_wall_1_s = k_wall_1_s,
                                                         products = products)
        end
        throw(ArgumentError("SpeciesWallModel: `:neutral_recombination` supports only `:ballistic_sticking` or `:constant`."))
    end

    if class == :electronic_quench
        throw(ArgumentError("SpeciesWallModel: `:electronic_quench` is not supported in phase 4."))
    end

    throw(ArgumentError("SpeciesWallModel: class must be one of `:ion_neutralization` or `:neutral_recombination`."))
end

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

function _validate_wall_loss_compat_flags(species_models::Dict{String, SpeciesWallModel};
                                          use_ion_losses::Union{Nothing, Bool} = nothing,
                                          use_neutral_recombination::Union{Nothing, Bool} = nothing,
                                          use_electronic_quenching::Union{Nothing, Bool} = nothing)
    if use_ion_losses !== nothing ||
       use_neutral_recombination !== nothing ||
       use_electronic_quenching !== nothing
        Base.depwarn("`WallLossConfig(use_ion_losses=..., use_neutral_recombination=..., use_electronic_quenching=...)` is deprecated; select the desired concrete wall model types in `species_models` instead.",
                     :WallLossConfig)
    end

    ion_models_present = any(model -> model isa IonNeutralizationWallModel,
                             values(species_models))
    neutral_models_present = any(model -> model isa BallisticNeutralRecombinationWallModel ||
                                          model isa ConstantNeutralRecombinationWallModel,
                                 values(species_models))

    use_ion_losses === false && ion_models_present &&
        throw(ArgumentError("WallLossConfig: `use_ion_losses = false` is incompatible with ion wall models in `species_models`."))
    use_neutral_recombination === false && neutral_models_present &&
        throw(ArgumentError("WallLossConfig: `use_neutral_recombination = false` is incompatible with neutral recombination wall models in `species_models`."))
    use_electronic_quenching === true &&
        throw(ArgumentError("WallLossConfig: `:electronic_quench` is not supported in phase 4."))

    return nothing
end

"""
$(SIGNATURES)

Wrapper-managed wall-loss configuration.
"""
struct WallLossConfig
    enabled::Bool
    species_models::Dict{String, SpeciesWallModel}

    function WallLossConfig(; enabled::Bool = true,
                            species_models::AbstractDict = Dict{String, SpeciesWallModel}(),
                            use_ion_losses::Union{Nothing, Bool} = nothing,
                            use_neutral_recombination::Union{Nothing, Bool} = nothing,
                            use_electronic_quenching::Union{Nothing, Bool} = nothing)
        species_model_dict = _coerce_wall_species_models(species_models)
        _validate_wall_loss_compat_flags(species_model_dict;
                                         use_ion_losses = use_ion_losses,
                                         use_neutral_recombination = use_neutral_recombination,
                                         use_electronic_quenching = use_electronic_quenching)
        return new(enabled, species_model_dict)
    end
end
