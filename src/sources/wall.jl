struct SegmentWallInputs
    a_wall_over_v_m_inv::Float64
    channel_gap_m::Union{Nothing, Float64}
    wall_temperature_K::Union{Nothing, Float64}
    ion_edge_to_center_ratio::Union{Nothing, Float64}

    function SegmentWallInputs(;
                               a_wall_over_v_m_inv::Real,
                               channel_gap_m::Union{Nothing, Real} = nothing,
                               wall_temperature_K::Union{Nothing, Real} = nothing,
                               ion_edge_to_center_ratio::Union{Nothing, Real} = nothing)
        function _coerce_optional_scalar(value, field_name::AbstractString)
            value === nothing && return nothing
            value_f64 = Float64(value)
            isfinite(value_f64) && value_f64 > 0.0 ||
                throw(ArgumentError("SegmentWallInputs: $field_name must be finite and strictly positive when provided."))
            return value_f64
        end

        a_wall_over_v_val = Float64(a_wall_over_v_m_inv)
        isfinite(a_wall_over_v_val) && a_wall_over_v_val > 0.0 ||
            throw(ArgumentError("SegmentWallInputs: a_wall_over_v_m_inv must be finite and strictly positive."))

        return new(a_wall_over_v_val,
                   _coerce_optional_scalar(channel_gap_m, "channel_gap_m"),
                   _coerce_optional_scalar(wall_temperature_K, "wall_temperature_K"),
                   _coerce_optional_scalar(ion_edge_to_center_ratio,
                                           "ion_edge_to_center_ratio"))
    end
end

const BOLTZMANN_J_PER_K = 1.380649e-23
const DEFAULT_ION_EDGE_TO_CENTER_RATIO = 0.86 / sqrt(3.0)

struct WallLossSpeciesIndexData
    species_order::Vector{String}
    species_positions::Dict{String, Int}
    all_indices::Dict{String, Vector{Int}}
    ground_indices::Dict{String, Int}
    density_indices::Dict{String, Int}
    charge_states::Dict{String, Int}
    is_electronic_resolved::Dict{String, Bool}
end

struct PreparedWallReactionData
    reactant::String
    reactant_indices::Vector{Int}
    reactant_ground_index::Int
    product_indices::Dict{String, Int}
    product_branching::Dict{String, Float64}
    reactant_molecular_weight::Float64
    product_molecular_weights::Dict{String, Float64}
    charge_state::Int
end

abstract type PreparedWallSpeciesModel end

struct PreparedWallReaction{M <: SpeciesWallModel} <: PreparedWallSpeciesModel
    reaction::PreparedWallReactionData
    model::M
end

struct PreparedWallLossData
    wall_config::WallLossConfig
    wall_inputs::SegmentWallInputs
    species_index_data::WallLossSpeciesIndexData
    models::Vector{PreparedWallSpeciesModel}
    segment_tt_K::Float64
    segment_te_K::Float64
end

@inline _reaction(model::PreparedWallReaction) = model.reaction
@inline _config_wall_model(model::PreparedWallReaction) = model.model

@inline function _wall_losses_enabled(sources::Union{Nothing, SourceTermsConfig})
    return sources !== nothing &&
           sources.wall_losses !== nothing &&
           sources.wall_losses.enabled
end

function _validate_direct_wall_loss_usage(sources::Union{Nothing, SourceTermsConfig})
    _wall_losses_enabled(sources) || return nothing
    throw(ArgumentError("WallLossConfig is currently supported only for profile-driven chain runs. " *
                        "Use `solve_terra_chain_steady` with a `terra_chain_profile_v4` profile containing `wall_profile`; " *
                        "direct 0D solves do not accept segment wall inputs."))
end

function _build_wall_loss_index_data(layout::ApiLayout,
                                     species_names::AbstractVector{<:AbstractString})
    length(species_names) == layout.nsp ||
        throw(ArgumentError("Wall-loss species ordering length ($(length(species_names))) does not match layout.nsp ($(layout.nsp))."))

    species_order = String[String(species_names[isp])
                           for isp in 1:(layout.nsp) if isp != layout.esp]
    species_positions = Dict{String, Int}(name => i for (i, name) in pairs(species_order))
    all_indices = Dict(name => Int[] for name in species_order)
    ground_indices = Dict{String, Int}()
    density_indices = Dict{String, Int}()
    charge_states = Dict{String, Int}()
    is_electronic_resolved = Dict{String, Bool}()

    @inbounds for isp in 1:(layout.nsp)
        isp == layout.esp && continue
        name = String(species_names[isp])
        charge_states[name] = layout.ie[isp]
        is_electronic_resolved[name] = layout.ies[isp] == 1
    end

    idx = 1
    if layout.n_eq_vib > 0
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0 || layout.ih[isp] != 2
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 0
                    continue
                end
                mv_isp_iex = layout.mv[isp, iex]
                for _ivx in 0:mv_isp_iex
                    if isp != layout.esp
                        push!(all_indices[String(species_names[isp])], idx)
                    end
                    idx += 1
                end
            end
        end
    end

    if layout.is_elec_sts
        @inbounds for isp in 1:(layout.nsp)
            if layout.ies[isp] == 0
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 1
                    continue
                end
                if isp != layout.esp
                    name = String(species_names[isp])
                    push!(all_indices[name], idx)
                    if iex == 1 && !haskey(ground_indices, name)
                        ground_indices[name] = idx
                    end
                end
                idx += 1
            end
        end
    end

    @inbounds for isp in 1:(layout.nsp)
        if layout.ies[isp] == 1
            continue
        end
        if layout.get_electron_density_by_charge_balance && isp == layout.esp
            continue
        end
        if isp != layout.esp
            name = String(species_names[isp])
            push!(all_indices[name], idx)
            density_indices[name] = idx
            ground_indices[name] = idx
        end
        idx += 1
    end

    @assert idx == layout.mom_range.start "Internal error: wall-loss continuity index mismatch."

    return WallLossSpeciesIndexData(species_order,
                                    species_positions,
                                    all_indices,
                                    ground_indices,
                                    density_indices,
                                    charge_states,
                                    is_electronic_resolved)
end

function _validate_wall_loss_species(wall_cfg::WallLossConfig,
                                     species_order::Vector{String})
    active_species = Set(species_order)
    invalid_model_species = sort!(String[name
                                         for name in keys(wall_cfg.species_models)
                                         if !(name in active_species)])
    invalid_product_species = sort!(unique(String[product_name
                                                  for model in values(wall_cfg.species_models)
                                                  for product_name in keys(_wall_model_products(model))
                                                  if !(product_name in active_species)]))

    if !isempty(invalid_model_species) || !isempty(invalid_product_species)
        throw(ArgumentError("WallLossConfig species references must match the active non-electron TERRA species ordering. " *
                            "Invalid species_models keys: $(isempty(invalid_model_species) ? "none" : join(invalid_model_species, ", ")); " *
                            "invalid product species: $(isempty(invalid_product_species) ? "none" : join(invalid_product_species, ", "))."))
    end

    return nothing
end

function _build_prepared_wall_reaction(reactant::String,
                                       model::SpeciesWallModel,
                                       reactant_indices::Vector{Int},
                                       species_index_data::WallLossSpeciesIndexData,
                                       molecular_weights::Dict{String, Float64},
                                       charge_state::Int)
    haskey(species_index_data.ground_indices, reactant) ||
        throw(ArgumentError("SpeciesWallModel for `$reactant`: no ground reservoir could be identified in the compact state layout."))
    isempty(reactant_indices) &&
        throw(ArgumentError("SpeciesWallModel for `$reactant`: no compact continuity indices were found for the configured reactant."))

    product_indices = Dict{String, Int}()
    product_molecular_weights = Dict{String, Float64}()
    product_branching = Dict{String, Float64}()

    for (product, coefficient) in pairs(_wall_model_products(model))
        haskey(species_index_data.ground_indices, product) ||
            throw(ArgumentError("SpeciesWallModel for `$reactant`: no ground reservoir could be identified for product `$product`."))
        product_indices[product] = species_index_data.ground_indices[product]
        product_molecular_weights[product] = molecular_weights[product]
        product_branching[product] = coefficient
    end

    reaction = PreparedWallReactionData(reactant,
                                        reactant_indices,
                                        species_index_data.ground_indices[reactant],
                                        product_indices,
                                        product_branching,
                                        molecular_weights[reactant],
                                        product_molecular_weights,
                                        charge_state)
    return PreparedWallReaction(reaction, model)
end

function _prepare_wall_species_model(reactant::String,
                                     model::IonNeutralizationWallModel,
                                     species_index_data::WallLossSpeciesIndexData,
                                     molecular_weights::Dict{String, Float64})
    reactant_charge = get(species_index_data.charge_states, reactant, 0)
    reactant_charge > 0 ||
        throw(ArgumentError("SpeciesWallModel for `$reactant`: `:ion_neutralization` requires a positively charged reactant."))
    reactant_indices = copy(species_index_data.all_indices[reactant])
    return _build_prepared_wall_reaction(reactant,
                                         model,
                                         reactant_indices,
                                         species_index_data,
                                         molecular_weights,
                                         reactant_charge)
end

function _prepare_wall_species_model(reactant::String,
                                     model::BallisticNeutralRecombinationWallModel,
                                     species_index_data::WallLossSpeciesIndexData,
                                     molecular_weights::Dict{String, Float64})
    reactant_charge = get(species_index_data.charge_states, reactant, 0)
    reactant_charge == 0 ||
        throw(ArgumentError("SpeciesWallModel for `$reactant`: `:neutral_recombination` requires a neutral reactant."))
    reactant_indices = [species_index_data.ground_indices[reactant]]
    return _build_prepared_wall_reaction(reactant,
                                         model,
                                         reactant_indices,
                                         species_index_data,
                                         molecular_weights,
                                         reactant_charge)
end

function _prepare_wall_species_model(reactant::String,
                                     model::ConstantNeutralRecombinationWallModel,
                                     species_index_data::WallLossSpeciesIndexData,
                                     molecular_weights::Dict{String, Float64})
    reactant_charge = get(species_index_data.charge_states, reactant, 0)
    reactant_charge == 0 ||
        throw(ArgumentError("SpeciesWallModel for `$reactant`: `:neutral_recombination` requires a neutral reactant."))
    reactant_indices = [species_index_data.ground_indices[reactant]]
    return _build_prepared_wall_reaction(reactant,
                                         model,
                                         reactant_indices,
                                         species_index_data,
                                         molecular_weights,
                                         reactant_charge)
end

@inline function _prepare_wall_losses(layout::ApiLayout,
                                      config::Config,
                                      ::Nothing;
                                      wall_inputs::Union{Nothing, SegmentWallInputs} = nothing)
    return nothing
end

function _prepare_wall_losses(layout::ApiLayout,
                              config::Config,
                              wall_cfg::WallLossConfig;
                              wall_inputs::Union{Nothing, SegmentWallInputs} = nothing)
    wall_cfg.enabled || return nothing
    wall_inputs === nothing &&
        throw(ArgumentError("WallLossConfig is enabled, but no segment wall inputs were provided. " *
                            "Wall losses are profile-driven and require per-segment `wall_profile` data."))

    species_index_data = _build_wall_loss_index_data(layout,
                                                     config.reactor.composition.species)
    _validate_wall_loss_species(wall_cfg, species_index_data.species_order)

    molecular_weights = get_molecular_weights(config.reactor.composition.species)
    molecular_weight_dict = Dict{String, Float64}()
    for isp in eachindex(config.reactor.composition.species)
        isp == layout.esp && continue
        species_name = String(config.reactor.composition.species[isp])
        molecular_weight_dict[species_name] = molecular_weights[isp]
    end

    prepared_models = PreparedWallSpeciesModel[]
    for reactant in sort!(collect(keys(wall_cfg.species_models)))
        push!(prepared_models,
              _prepare_wall_species_model(reactant,
                                          wall_cfg.species_models[reactant],
                                          species_index_data,
                                          molecular_weight_dict))
    end

    return PreparedWallLossData(wall_cfg,
                                wall_inputs,
                                species_index_data,
                                prepared_models,
                                Float64(config.reactor.thermal.Tt),
                                Float64(config.reactor.thermal.Te))
end

@inline _wall_mass_density_to_number_density_cgs(rho_cgs::Real, molecular_weight::Real) =
    Float64(rho_cgs) * AVOGADRO / Float64(molecular_weight)

@inline _wall_number_density_to_mass_density_cgs(n_cgs::Real, molecular_weight::Real) =
    Float64(n_cgs) * Float64(molecular_weight) / AVOGADRO

@inline _molecular_mass_kg(molecular_weight_g_mol::Real) = Float64(molecular_weight_g_mol) *
                                                           1.0e-3 / AVOGADRO

function _wall_rate_1_s(model::PreparedWallReaction{IonNeutralizationWallModel},
                        wall_losses::PreparedWallLossData)
    reaction = _reaction(model)
    h_i = wall_losses.wall_inputs.ion_edge_to_center_ratio === nothing ?
          DEFAULT_ION_EDGE_TO_CENTER_RATIO :
          wall_losses.wall_inputs.ion_edge_to_center_ratio
    u_bohm = model.model.bohm_scale *
             sqrt(reaction.charge_state * BOLTZMANN_J_PER_K * wall_losses.segment_te_K /
                  _molecular_mass_kg(reaction.reactant_molecular_weight))
    return wall_losses.wall_inputs.a_wall_over_v_m_inv * h_i * u_bohm
end

function _wall_rate_1_s(model::PreparedWallReaction{BallisticNeutralRecombinationWallModel},
                        wall_losses::PreparedWallLossData)
    reaction = _reaction(model)
    cbar = sqrt(8.0 * BOLTZMANN_J_PER_K * wall_losses.segment_tt_K /
                (pi * _molecular_mass_kg(reaction.reactant_molecular_weight)))
    return model.model.gamma * cbar * 0.25 * wall_losses.wall_inputs.a_wall_over_v_m_inv
end

function _wall_rate_1_s(model::PreparedWallReaction{ConstantNeutralRecombinationWallModel},
                        wall_losses::PreparedWallLossData)
    return model.model.k_wall_1_s
end

function _accumulate_wall_loss_sources!(du::Union{Nothing, Vector{Float64}},
                                        u::Vector{Float64},
                                        wall_losses::PreparedWallLossData;
                                        species_net_drho::Union{Nothing, Vector{Float64}} = nothing,
                                        k_wall_dict::Union{Nothing, Dict{String, Float64}} = nothing)
    for model in wall_losses.models
        reaction = _reaction(model)
        k_wall = _wall_rate_1_s(model, wall_losses)
        k_wall_dict === nothing || (k_wall_dict[reaction.reactant] = k_wall)

        total_lost_number_rate = 0.0
        total_lost_drho = 0.0
        @inbounds for idx in reaction.reactant_indices
            rho_val = Float64(u[idx])
            rho_val <= 0.0 && continue
            drho_val = -k_wall * rho_val
            du === nothing || (du[idx] += drho_val)
            total_lost_drho -= drho_val
            total_lost_number_rate += _wall_mass_density_to_number_density_cgs(-drho_val,
                                                                               reaction.reactant_molecular_weight)
        end

        if species_net_drho !== nothing && total_lost_drho > 0.0
            species_net_drho[wall_losses.species_index_data.species_positions[reaction.reactant]] -=
                total_lost_drho
        end

        total_lost_number_rate <= 0.0 && continue
        for (product, coefficient) in pairs(reaction.product_branching)
            coefficient == 0.0 && continue
            product_drho = _wall_number_density_to_mass_density_cgs(coefficient *
                                                                    total_lost_number_rate,
                                                                    reaction.product_molecular_weights[product])
            du === nothing || (du[reaction.product_indices[product]] += product_drho)
            if species_net_drho !== nothing
                species_net_drho[wall_losses.species_index_data.species_positions[product]] +=
                    product_drho
            end
        end
    end
    return nothing
end

@inline _apply_wall_losses!(du::Vector{Float64}, u::Vector{Float64}, ::Nothing) = nothing

function _apply_wall_losses!(du::Vector{Float64},
                             u::Vector{Float64},
                             wall_losses::PreparedWallLossData)
    isempty(wall_losses.models) && return nothing
    _accumulate_wall_loss_sources!(du, u, wall_losses)
    return nothing
end

function _wall_loss_frame_snapshot(u::Vector{Float64},
                                   wall_losses::Union{Nothing, PreparedWallLossData},
                                   unit_system::Symbol)
    wall_losses === nothing && return nothing
    isempty(wall_losses.models) && return nothing

    species_drho = zeros(Float64, length(wall_losses.species_index_data.species_order))
    k_wall = Dict{String, Float64}()
    _accumulate_wall_loss_sources!(nothing, u, wall_losses;
                                   species_net_drho = species_drho,
                                   k_wall_dict = k_wall)

    if unit_system == :SI
        species_drho = convert_density_cgs_to_si(species_drho)
    end

    species_rates = Dict{String, Float64}()
    for (i, name) in pairs(wall_losses.species_index_data.species_order)
        species_rates[name] = species_drho[i]
    end

    return (species_mass_density_rates = species_rates,
            k_wall_1_s = k_wall)
end

function _wall_loss_metadata(wall_losses::Union{Nothing, PreparedWallLossData})
    wall_losses === nothing && return nothing

    species_models = Dict{String, Any}()
    for model in wall_losses.models
        reaction = _reaction(model)
        config_model = _config_wall_model(model)
        species_models[reaction.reactant] = Dict{String, Any}("model_type" => _wall_model_type_name(config_model),
                                                              "charge_state" => reaction.charge_state,
                                                              "parameters" => _wall_parameter_dict(config_model),
                                                              "products" => copy(reaction.product_branching),
                                                              "reactant_indices" => copy(reaction.reactant_indices),
                                                              "reactant_ground_index" => reaction.reactant_ground_index,
                                                              "product_indices" => copy(reaction.product_indices))
    end

    return Dict{String, Any}("species_order" => copy(wall_losses.species_index_data.species_order),
                             "segment_inputs" => Dict{String, Any}("a_wall_over_v_m_inv" => wall_losses.wall_inputs.a_wall_over_v_m_inv,
                                                                   "channel_gap_m" => wall_losses.wall_inputs.channel_gap_m,
                                                                   "wall_temperature_K" => wall_losses.wall_inputs.wall_temperature_K,
                                                                   "ion_edge_to_center_ratio" => wall_losses.wall_inputs.ion_edge_to_center_ratio,
                                                                   "tt_K" => wall_losses.segment_tt_K,
                                                                   "te_K" => wall_losses.segment_te_K),
                             "species_models" => species_models)
end
