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
            isfinite(value_f64) && value_f64 > 0.0 || throw(ArgumentError(
                "SegmentWallInputs: $field_name must be finite and strictly positive when provided."
            ))
            return value_f64
        end

        a_wall_over_v_val = Float64(a_wall_over_v_m_inv)
        isfinite(a_wall_over_v_val) && a_wall_over_v_val > 0.0 || throw(ArgumentError(
            "SegmentWallInputs: a_wall_over_v_m_inv must be finite and strictly positive."
        ))

        return new(
            a_wall_over_v_val,
            _coerce_optional_scalar(channel_gap_m, "channel_gap_m"),
            _coerce_optional_scalar(wall_temperature_K, "wall_temperature_K"),
            _coerce_optional_scalar(ion_edge_to_center_ratio, "ion_edge_to_center_ratio"),
        )
    end
end

struct PreparedWallLossData
    wall_config::WallLossConfig
    wall_inputs::SegmentWallInputs
    species_names::Vector{String}
    species_indices::Dict{String, Vector{Int}}
    molecular_weights::Dict{String, Float64}
end

@inline function _wall_losses_enabled(sources::Union{Nothing, SourceTermsConfig})
    return sources !== nothing &&
           sources.wall_losses !== nothing &&
           sources.wall_losses.enabled
end

function _validate_direct_wall_loss_usage(sources::Union{Nothing, SourceTermsConfig})
    _wall_losses_enabled(sources) || return nothing
    throw(ArgumentError(
        "WallLossConfig is currently supported only for profile-driven chain runs. " *
        "Use `solve_terra_chain_steady` with a `terra_chain_profile_v4` profile containing `wall_profile`; " *
        "direct 0D solves do not accept segment wall inputs in Phase 1."
    ))
end

@inline function _build_wall_loss_continuity_indices(layout::ApiLayout,
        species_names::AbstractVector{<:AbstractString})
    return _build_residence_time_continuity_indices(layout, species_names)
end

function _validate_wall_loss_species(wall_cfg::WallLossConfig, species_order::Vector{String})
    active_species = Set(species_order)
    invalid_model_species = sort!(String[
        name for name in keys(wall_cfg.species_models) if !(name in active_species)
    ])
    invalid_product_species = sort!(unique(String[
        product_name
        for model in values(wall_cfg.species_models)
        for product_name in keys(model.products)
        if !(product_name in active_species)
    ]))

    if !isempty(invalid_model_species) || !isempty(invalid_product_species)
        throw(ArgumentError(
            "WallLossConfig species references must match the active non-electron TERRA species ordering. " *
            "Invalid species_models keys: $(isempty(invalid_model_species) ? "none" : join(invalid_model_species, ", ")); " *
            "invalid product species: $(isempty(invalid_product_species) ? "none" : join(invalid_product_species, ", "))."
        ))
    end
end

function _prepare_wall_loss_data(layout::ApiLayout, config::Config, wall_cfg::WallLossConfig;
        wall_inputs::Union{Nothing, SegmentWallInputs} = nothing)
    wall_cfg.enabled || return nothing
    wall_inputs === nothing && throw(ArgumentError(
        "WallLossConfig is enabled, but no segment wall inputs were provided. " *
        "Phase 1 wall losses are profile-driven and require per-segment `wall_profile` data."
    ))

    continuity = _build_wall_loss_continuity_indices(layout, config.reactor.composition.species)
    species_order = continuity.species_order
    _validate_wall_loss_species(wall_cfg, species_order)

    referenced_species = Set{String}()
    for (name, model) in pairs(wall_cfg.species_models)
        push!(referenced_species, name)
        for product_name in keys(model.products)
            push!(referenced_species, product_name)
        end
    end
    referenced_species_vec = sort!(collect(referenced_species))

    species_indices = Dict{String, Vector{Int}}()
    for name in referenced_species_vec
        species_indices[name] = copy(get(continuity.species_indices, name, Int[]))
    end

    molecular_weights = get_molecular_weights(config.reactor.composition.species)
    molecular_weight_dict = Dict{String, Float64}()
    for isp in eachindex(config.reactor.composition.species)
        isp == layout.esp && continue
        species_name = String(config.reactor.composition.species[isp])
        species_name in referenced_species || continue
        molecular_weight_dict[species_name] = molecular_weights[isp]
    end

    return PreparedWallLossData(
        wall_cfg,
        wall_inputs,
        referenced_species_vec,
        species_indices,
        molecular_weight_dict,
    )
end

@inline function _apply_wall_loss_term!(du::Vector{Float64}, u::Vector{Float64}, wall_losses)
    wall_losses === nothing && return nothing
    return nothing
end

function _chain_wall_profile_to_dict(wall_profile::ChainWallProfile)
    payload = Dict{String, Any}(
        "a_wall_over_v_m_inv" => copy(wall_profile.a_wall_over_v_m_inv),
    )
    if wall_profile.channel_gap_m !== nothing
        payload["channel_gap_m"] = copy(wall_profile.channel_gap_m)
    end
    if wall_profile.wall_temperature_K !== nothing
        payload["wall_temperature_K"] = copy(wall_profile.wall_temperature_K)
    end
    if wall_profile.ion_edge_to_center_ratio !== nothing
        payload["ion_edge_to_center_ratio"] = copy(wall_profile.ion_edge_to_center_ratio)
    end
    return payload
end
