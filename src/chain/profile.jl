function _validated_positive_profile_vector(values,
                                            field_name::AbstractString,
                                            context::AbstractString;
                                            expected_length::Union{Nothing, Integer} = nothing,
                                            allow_empty::Bool = false)
    values_vec = Float64.(values)
    isempty(values_vec) && !allow_empty &&
        throw(ArgumentError("$context: $field_name must be non-empty."))
    expected_length !== nothing && length(values_vec) != expected_length &&
        throw(ArgumentError("$context: $field_name length $(length(values_vec)) does not match required length $expected_length."))
    for (i, value) in pairs(values_vec)
        isfinite(value) ||
            throw(ArgumentError("$context: $field_name[$i] must be finite."))
        value > 0.0 ||
            throw(ArgumentError("$context: $field_name[$i] must be strictly positive."))
    end
    return values_vec
end

function _validated_optional_positive_profile_vector(values,
                                                     field_name::AbstractString,
                                                     context::AbstractString,
                                                     expected_length::Integer)
    values === nothing && return nothing
    return _validated_positive_profile_vector(values, field_name, context;
                                              expected_length = expected_length)
end

function _validate_strictly_increasing(values::AbstractVector{<:Real},
                                       field_name::AbstractString,
                                       context::AbstractString)
    for i in 2:length(values)
        values[i] > values[i - 1] ||
            throw(ArgumentError("$context: $field_name must be strictly increasing (failed at indices $(i - 1), $i)."))
    end
    return values
end

function _validated_profile_species_velocity_dict(species_u_m_s::AbstractDict,
                                                  n::Integer)
    species_u_dict = Dict{String, Vector{Float64}}()
    for (name, values) in pairs(species_u_m_s)
        name_str = String(name)
        isempty(strip(name_str)) &&
            throw(ArgumentError("AxialChainProfile: species_u_m_s keys must be non-empty."))
        lowercase(strip(name_str)) in ("e", "e-", "electron") &&
            throw(ArgumentError("AxialChainProfile: species_u_m_s must not include electron species."))
        species_u_dict[name_str] = _validated_positive_profile_vector(values,
                                                                      "species_u_m_s[$name_str]",
                                                                      "AxialChainProfile";
                                                                      expected_length = n)
    end
    isempty(species_u_dict) &&
        throw(ArgumentError("AxialChainProfile: species_u_m_s must contain at least one species."))
    return species_u_dict
end

function _validated_profile_diagnostics_dict(diagnostics::AbstractDict, n::Integer)
    diagnostics_dict = Dict{String, Vector{Float64}}()
    for (key, values) in pairs(diagnostics)
        key_str = String(key)
        values_vec = Float64.(values)
        length(values_vec) == n ||
            throw(ArgumentError("AxialChainProfile: diagnostic `$(key_str)` length $(length(values_vec)) does not match required profile length $n."))
        for (i, value) in pairs(values_vec)
            isfinite(value) ||
                throw(ArgumentError("AxialChainProfile: diagnostic `$(key_str)[$i]` must be finite."))
        end
        diagnostics_dict[key_str] = values_vec
    end
    return diagnostics_dict
end

function _validate_wall_profile_length(wall_profile, n::Integer)
    length(wall_profile.a_wall_over_v_m_inv) == n ||
        throw(ArgumentError("AxialChainProfile: wall_profile.a_wall_over_v_m_inv length $(length(wall_profile.a_wall_over_v_m_inv)) does not match required profile length $n."))

    for (field_name, values) in (("channel_gap_m", wall_profile.channel_gap_m),
                                 ("wall_temperature_K", wall_profile.wall_temperature_K),
                                 ("ion_edge_to_center_ratio",
                                  wall_profile.ion_edge_to_center_ratio))
        values === nothing && continue
        length(values) == n ||
            throw(ArgumentError("AxialChainProfile: wall_profile.$field_name length $(length(values)) does not match required profile length $n."))
    end

    return wall_profile
end

"""
$(SIGNATURES)

Composition payload for the self-contained chain-profile inlet.
"""
struct ChainProfileInletComposition
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density_m3::Float64

    function ChainProfileInletComposition(; species,
                                          mole_fractions,
                                          total_number_density_m3)
        species_vec, mole_frac_vec, n_tot = _validated_species_composition(species,
                                                                           mole_fractions,
                                                                           total_number_density_m3;
                                                                           context = "ChainProfileInletComposition",
                                                                           density_field = "total_number_density_m3")
        return new(species_vec, mole_frac_vec, n_tot)
    end
end

"""
$(SIGNATURES)

Self-contained inlet state used to initialize segment 1 of the chain solver.
"""
struct ChainProfileInlet
    composition::ChainProfileInletComposition
    thermal::ReactorThermalState
    source_compact_index::Int

    function ChainProfileInlet(; composition::ChainProfileInletComposition,
                               thermal::ReactorThermalState,
                               source_compact_index::Integer)
        source_idx = Int(source_compact_index)
        source_idx >= 1 ||
            throw(ArgumentError("ChainProfileInlet: source_compact_index must be >= 1."))
        return new(composition, thermal, source_idx)
    end
end

"""
$(SIGNATURES)

Normalized wall-geometry/profile inputs for chain-driven wall losses.
"""
struct ChainWallProfile
    a_wall_over_v_m_inv::Vector{Float64}
    channel_gap_m::Union{Nothing, Vector{Float64}}
    wall_temperature_K::Union{Nothing, Vector{Float64}}
    ion_edge_to_center_ratio::Union{Nothing, Vector{Float64}}

    function ChainWallProfile(; a_wall_over_v_m_inv,
                              channel_gap_m::Union{Nothing, AbstractVector} = nothing,
                              wall_temperature_K::Union{Nothing, AbstractVector} = nothing,
                              ion_edge_to_center_ratio::Union{Nothing, AbstractVector} = nothing)
        a_wall_over_v_vec = _validated_positive_profile_vector(a_wall_over_v_m_inv,
                                                               "a_wall_over_v_m_inv",
                                                               "ChainWallProfile")
        n = length(a_wall_over_v_vec)

        return new(a_wall_over_v_vec,
                   _validated_optional_positive_profile_vector(channel_gap_m,
                                                              "channel_gap_m",
                                                              "ChainWallProfile",
                                                              n),
                   _validated_optional_positive_profile_vector(wall_temperature_K,
                                                              "wall_temperature_K",
                                                              "ChainWallProfile",
                                                              n),
                   _validated_optional_positive_profile_vector(ion_edge_to_center_ratio,
                                                              "ion_edge_to_center_ratio",
                                                              "ChainWallProfile",
                                                              n))
    end
end

"""
$(SIGNATURES)

Normalized axial chain profile used by the TERRA chain-of-CSTR interface.

# Fields
- `z_m::Vector{Float64}`: Axial coordinate points (m), strictly increasing
- `dx_m::Vector{Float64}`: Point-aligned control-volume lengths (m), strictly positive
- `te_K::Vector{Float64}`: Electron temperature profile (K), strictly positive
- `species_u_m_s::Dict{String, Vector{Float64}}`: Per-species convective velocities (m/s), strictly positive
- `wall_profile::Union{Nothing, ChainWallProfile}`: Optional validated wall-profile inputs for wall-loss source terms
- `inlet::ChainProfileInlet`: Self-contained segment-1 inlet state
- `diagnostics::Dict{String, Vector{Float64}}`: Optional diagnostics arrays
- `generator::Dict{String, Any}`: Generator metadata from interchange artifact
- `selection::Dict{String, Any}`: Selection metadata from interchange artifact
- `schema_version::String`: Interchange schema version
- `source_snapshot::Union{Nothing, Dict{String, Any}}`: Optional source provenance snapshot
"""
struct AxialChainProfile
    z_m::Vector{Float64}
    dx_m::Vector{Float64}
    te_K::Vector{Float64}
    species_u_m_s::Dict{String, Vector{Float64}}
    wall_profile::Union{Nothing, ChainWallProfile}
    inlet::ChainProfileInlet
    diagnostics::Dict{String, Vector{Float64}}
    generator::Dict{String, Any}
    selection::Dict{String, Any}
    schema_version::String
    source_snapshot::Union{Nothing, Dict{String, Any}}

    function AxialChainProfile(; z_m,
                               dx_m,
                               te_K,
                               species_u_m_s,
                               wall_profile::Union{Nothing, ChainWallProfile} = nothing,
                               inlet::ChainProfileInlet,
                               diagnostics::AbstractDict = Dict{String, Vector{Float64}}(),
                               generator::AbstractDict = Dict{String, Any}(),
                               selection::AbstractDict = Dict{String, Any}(),
                               schema_version::AbstractString = "terra_chain_profile_v4",
                               source_snapshot::Union{Nothing, AbstractDict} = nothing)
        z_m_vec = Float64.(z_m)
        isempty(z_m_vec) &&
            throw(ArgumentError("AxialChainProfile: profile arrays must be non-empty."))
        for (i, value) in pairs(z_m_vec)
            isfinite(value) ||
                throw(ArgumentError("AxialChainProfile: z_m[$i] must be finite."))
        end
        _validate_strictly_increasing(z_m_vec, "z_m", "AxialChainProfile")

        n = length(z_m_vec)
        dx_m_vec = _validated_positive_profile_vector(dx_m, "dx_m",
                                                      "AxialChainProfile";
                                                      expected_length = n)
        te_K_vec = _validated_positive_profile_vector(te_K, "te_K",
                                                      "AxialChainProfile";
                                                      expected_length = n)
        species_u_dict = _validated_profile_species_velocity_dict(species_u_m_s, n)
        diagnostics_dict = _validated_profile_diagnostics_dict(diagnostics, n)
        wall_profile !== nothing && _validate_wall_profile_length(wall_profile, n)
        inlet.source_compact_index <= n ||
            throw(ArgumentError("AxialChainProfile: inlet.source_compact_index ($(inlet.source_compact_index)) must be <= retained profile length $(n)."))

        generator_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(generator))
        selection_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(selection))
        source_snapshot_dict = source_snapshot === nothing ? nothing :
                               Dict{String, Any}(String(k) => v
                                                 for (k, v) in pairs(source_snapshot))

        return new(z_m_vec,
                   dx_m_vec,
                   te_K_vec,
                   species_u_dict,
                   wall_profile,
                   inlet,
                   diagnostics_dict,
                   generator_dict,
                   selection_dict,
                   String(schema_version),
                   source_snapshot_dict)
    end
end

"""
$(SIGNATURES)

Validate an `AxialChainProfile` against the chain-profile contract.

# Arguments
- `profile::AxialChainProfile`: Profile to validate

# Returns
- `true` if validation passes

# Throws
- `ArgumentError` if the profile is invalid
"""
function validate_axial_chain_profile(profile::AxialChainProfile)
    isempty(profile.z_m) &&
        throw(ArgumentError("AxialChainProfile must contain at least one axial point."))

    for (i, value) in pairs(profile.z_m)
        isfinite(value) ||
            throw(ArgumentError("AxialChainProfile: z_m[$i] must be finite."))
    end
    _validate_strictly_increasing(profile.z_m, "z_m", "AxialChainProfile")
    _validated_positive_profile_vector(profile.dx_m, "dx_m", "AxialChainProfile";
                                       expected_length = length(profile.z_m))
    _validated_positive_profile_vector(profile.te_K, "te_K", "AxialChainProfile";
                                       expected_length = length(profile.z_m))
    _validated_profile_species_velocity_dict(profile.species_u_m_s, length(profile.z_m))
    _validated_profile_diagnostics_dict(profile.diagnostics, length(profile.z_m))
    profile.wall_profile !== nothing &&
        _validate_wall_profile_length(profile.wall_profile, length(profile.z_m))

    source_idx = profile.inlet.source_compact_index
    source_idx <= length(profile.z_m) ||
        throw(ArgumentError("AxialChainProfile: inlet.source_compact_index ($(source_idx)) must be <= retained profile length $(length(profile.z_m))."))

    return true
end

function _coerce_positive_int(value)::Union{Nothing, Int}
    if value isa Integer
        int_value = Int(value)
        return int_value >= 1 ? int_value : nothing
    end
    if value isa Real && isfinite(value) && value >= 1 && isinteger(value)
        return Int(round(value))
    end
    return nothing
end

function _coerce_nonnegative_int(value)::Union{Nothing, Int}
    if value isa Integer
        int_value = Int(value)
        return int_value >= 0 ? int_value : nothing
    end
    if value isa Real && isfinite(value) && value >= 0 && isinteger(value)
        return Int(round(value))
    end
    return nothing
end

function _selection_index_vector(selection::AbstractDict,
                                 key::AbstractString,
                                 expected_length::Integer)::Union{Nothing, Vector{Int}}
    haskey(selection, key) || return nothing
    values = selection[key]
    values isa AbstractVector || return nothing
    length(values) == expected_length || return nothing

    indices = Vector{Int}(undef, expected_length)
    for i in eachindex(values)
        index = _coerce_positive_int(values[i])
        index === nothing && return nothing
        indices[i] = index
    end
    return indices
end

function _resolve_compact_to_source_index(profile::AxialChainProfile)
    n_segments = length(profile.z_m)
    selection = profile.selection

    for key in ("compact_to_source_index", "source_cell_indices")
        source_indices = _selection_index_vector(selection, key, n_segments)
        source_indices === nothing || return source_indices
    end

    trim_start = haskey(selection, "trim_start_index") ?
                 _coerce_positive_int(selection["trim_start_index"]) : nothing
    if trim_start === nothing && haskey(selection, "trimmed_point_count")
        trimmed_points = _coerce_nonnegative_int(selection["trimmed_point_count"])
        trim_start = trimmed_points === nothing ? nothing : (trimmed_points + 1)
    end
    trim_start !== nothing && return collect(trim_start:(trim_start + n_segments - 1))

    return collect(1:n_segments)
end

function _resolve_original_point_count(profile::AxialChainProfile)::Union{Nothing, Int}
    retained_points = length(profile.z_m)
    selection = profile.selection

    if haskey(selection, "original_point_count")
        original_points = _coerce_positive_int(selection["original_point_count"])
        if original_points !== nothing && original_points >= retained_points
            return original_points
        end
    end

    if haskey(selection, "trimmed_point_count")
        trimmed_points = _coerce_nonnegative_int(selection["trimmed_point_count"])
        trimmed_points === nothing || return retained_points + trimmed_points
    end

    return nothing
end

function _profile_species_velocity_at_segment(profile::AxialChainProfile,
                                              segment_index::Integer)
    species_u = Dict{String, Float64}()
    for (name, values) in pairs(profile.species_u_m_s)
        species_u[name] = values[segment_index]
    end
    return species_u
end

function _segment_wall_profile_values(profile::AxialChainProfile,
                                      segment_index::Integer,
                                      wall_cfg::Union{Nothing, WallLossConfig} = nothing)
    if profile.wall_profile === nothing
        if wall_cfg !== nothing && wall_cfg.enabled
            throw(ArgumentError("WallLossConfig is enabled, but AxialChainProfile.wall_profile is missing. " *
                                "Wall losses require `wall_profile` in a `terra_chain_profile_v4` chain profile."))
        end
        return nothing
    end

    wall_profile = profile.wall_profile
    return (; a_wall_over_v_m_inv = wall_profile.a_wall_over_v_m_inv[segment_index],
            channel_gap_m = wall_profile.channel_gap_m === nothing ? nothing :
                            wall_profile.channel_gap_m[segment_index],
            wall_temperature_K = wall_profile.wall_temperature_K === nothing ? nothing :
                                 wall_profile.wall_temperature_K[segment_index],
            ion_edge_to_center_ratio = wall_profile.ion_edge_to_center_ratio === nothing ?
                                       nothing :
                                       wall_profile.ion_edge_to_center_ratio[segment_index])
end

function _profile_inlet_reactor(profile::AxialChainProfile, unit_system::Symbol)
    inlet = profile.inlet
    n_total = unit_system == :SI ? inlet.composition.total_number_density_m3 :
              convert_number_density_si_to_cgs(inlet.composition.total_number_density_m3)
    composition = ReactorComposition(;
                                     species = inlet.composition.species,
                                     mole_fractions = inlet.composition.mole_fractions,
                                     total_number_density = n_total)
    return ReactorConfig(; composition = composition, thermal = inlet.thermal)
end

function _chain_wall_profile_to_dict(wall_profile::ChainWallProfile)
    payload = Dict{String, Any}("a_wall_over_v_m_inv" => copy(wall_profile.a_wall_over_v_m_inv))
    wall_profile.channel_gap_m !== nothing &&
        (payload["channel_gap_m"] = copy(wall_profile.channel_gap_m))
    wall_profile.wall_temperature_K !== nothing &&
        (payload["wall_temperature_K"] = copy(wall_profile.wall_temperature_K))
    wall_profile.ion_edge_to_center_ratio !== nothing &&
        (payload["ion_edge_to_center_ratio"] = copy(wall_profile.ion_edge_to_center_ratio))
    return payload
end

function _chain_profile_metadata(profile::AxialChainProfile,
                                 diagnostics::AbstractDict,
                                 compact_to_source_index::AbstractVector{<:Integer})
    return ChainMetadata(;
                         schema_version = profile.schema_version,
                         generator = profile.generator,
                         selection = profile.selection,
                         source_snapshot = profile.source_snapshot,
                         diagnostics = diagnostics,
                         compact_to_source_index = compact_to_source_index,
                         original_point_count = _resolve_original_point_count(profile),
                         retained_point_count = length(profile.z_m))
end
