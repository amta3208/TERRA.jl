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
        species_vec = String.(species)
        mole_frac_vec = Float64.(mole_fractions)
        n_tot = Float64(total_number_density_m3)

        if length(species_vec) != length(mole_frac_vec)
            throw(ArgumentError("ChainProfileInletComposition: species and mole_fractions arrays must have same length."))
        end
        isempty(species_vec) &&
            throw(ArgumentError("ChainProfileInletComposition: at least one species must be specified."))
        any(mole_frac_vec .< 0) &&
            throw(ArgumentError("ChainProfileInletComposition: mole fractions must be non-negative."))
        abs(sum(mole_frac_vec) - 1.0) > 1e-10 &&
            throw(ArgumentError("ChainProfileInletComposition: mole fractions must sum to 1.0, got $(sum(mole_frac_vec))."))
        n_tot > 0.0 ||
            throw(ArgumentError("ChainProfileInletComposition: total_number_density_m3 must be positive."))
        length(unique(species_vec)) == length(species_vec) ||
            throw(ArgumentError("ChainProfileInletComposition: duplicate species names found."))
        for name in species_vec
            isempty(strip(name)) &&
                throw(ArgumentError("ChainProfileInletComposition: species names cannot be empty."))
        end

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
        a_wall_over_v_vec = Float64.(a_wall_over_v_m_inv)
        isempty(a_wall_over_v_vec) &&
            throw(ArgumentError("ChainWallProfile: a_wall_over_v_m_inv must be non-empty."))
        for (i, value) in pairs(a_wall_over_v_vec)
            isfinite(value) && value > 0.0 ||
                throw(ArgumentError("ChainWallProfile: a_wall_over_v_m_inv[$i] must be finite and strictly positive."))
        end

        function _coerce_optional_wall_array(values, field_name::AbstractString)
            values === nothing && return nothing
            values_vec = Float64.(values)
            length(values_vec) == length(a_wall_over_v_vec) ||
                throw(ArgumentError("ChainWallProfile: $field_name length $(length(values_vec)) does not match required length $(length(a_wall_over_v_vec))."))
            for (i, value) in pairs(values_vec)
                isfinite(value) && value > 0.0 ||
                    throw(ArgumentError("ChainWallProfile: $field_name[$i] must be finite and strictly positive."))
            end
            return values_vec
        end

        return new(a_wall_over_v_vec,
                   _coerce_optional_wall_array(channel_gap_m, "channel_gap_m"),
                   _coerce_optional_wall_array(wall_temperature_K, "wall_temperature_K"),
                   _coerce_optional_wall_array(ion_edge_to_center_ratio,
                                               "ion_edge_to_center_ratio"))
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
        dx_m_vec = Float64.(dx_m)
        te_K_vec = Float64.(te_K)

        n = length(z_m_vec)
        if n == 0
            throw(ArgumentError("AxialChainProfile: profile arrays must be non-empty."))
        end
        if length(dx_m_vec) != n || length(te_K_vec) != n
            throw(ArgumentError("AxialChainProfile: all required profile arrays must have identical lengths."))
        end

        for (i, value) in pairs(z_m_vec)
            isfinite(value) ||
                throw(ArgumentError("AxialChainProfile: z_m[$i] must be finite."))
        end
        for i in 2:n
            z_m_vec[i] > z_m_vec[i - 1] ||
                throw(ArgumentError("AxialChainProfile: z_m must be strictly increasing (failed at indices $(i - 1), $i)."))
        end

        for (name, values) in (("dx_m", dx_m_vec),
                               ("te_K", te_K_vec))
            for (i, value) in pairs(values)
                isfinite(value) ||
                    throw(ArgumentError("AxialChainProfile: $(name)[$i] must be finite."))
                value > 0.0 ||
                    throw(ArgumentError("AxialChainProfile: $(name)[$i] must be strictly positive."))
            end
        end

        species_u_dict = Dict{String, Vector{Float64}}()
        for (name, values) in pairs(species_u_m_s)
            name_str = String(name)
            isempty(strip(name_str)) &&
                throw(ArgumentError("AxialChainProfile: species_u_m_s keys must be non-empty."))
            values_vec = Float64.(values)
            if length(values_vec) != n
                throw(ArgumentError("AxialChainProfile: species_u_m_s[$name_str] length $(length(values_vec)) does not match required profile length $n."))
            end
            for (i, value) in pairs(values_vec)
                isfinite(value) ||
                    throw(ArgumentError("AxialChainProfile: species_u_m_s[$name_str][$i] must be finite."))
                value > 0.0 ||
                    throw(ArgumentError("AxialChainProfile: species_u_m_s[$name_str][$i] must be strictly positive."))
            end
            species_u_dict[name_str] = values_vec
        end
        isempty(species_u_dict) &&
            throw(ArgumentError("AxialChainProfile: species_u_m_s must contain at least one species."))

        if wall_profile !== nothing
            length(wall_profile.a_wall_over_v_m_inv) == n ||
                throw(ArgumentError("AxialChainProfile: wall_profile.a_wall_over_v_m_inv length $(length(wall_profile.a_wall_over_v_m_inv)) does not match required profile length $n."))
            if wall_profile.channel_gap_m !== nothing
                length(wall_profile.channel_gap_m) == n ||
                    throw(ArgumentError("AxialChainProfile: wall_profile.channel_gap_m length $(length(wall_profile.channel_gap_m)) does not match required profile length $n."))
            end
            if wall_profile.wall_temperature_K !== nothing
                length(wall_profile.wall_temperature_K) == n ||
                    throw(ArgumentError("AxialChainProfile: wall_profile.wall_temperature_K length $(length(wall_profile.wall_temperature_K)) does not match required profile length $n."))
            end
            if wall_profile.ion_edge_to_center_ratio !== nothing
                length(wall_profile.ion_edge_to_center_ratio) == n ||
                    throw(ArgumentError("AxialChainProfile: wall_profile.ion_edge_to_center_ratio length $(length(wall_profile.ion_edge_to_center_ratio)) does not match required profile length $n."))
            end
        end

        diagnostics_dict = Dict{String, Vector{Float64}}()
        for (key, values) in pairs(diagnostics)
            key_str = String(key)
            values_vec = Float64.(values)
            if length(values_vec) != n
                throw(ArgumentError("AxialChainProfile: diagnostic `$(key_str)` length $(length(values_vec)) does not match required profile length $n."))
            end
            for (i, value) in pairs(values_vec)
                isfinite(value) ||
                    throw(ArgumentError("AxialChainProfile: diagnostic `$(key_str)[$i]` must be finite."))
            end
            diagnostics_dict[key_str] = values_vec
        end

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
    n = length(profile.z_m)
    if n == 0
        throw(ArgumentError("AxialChainProfile must contain at least one axial point."))
    end
    if length(profile.dx_m) != n || length(profile.te_K) != n
        throw(ArgumentError("AxialChainProfile required arrays must have identical lengths."))
    end

    for (i, value) in pairs(profile.z_m)
        isfinite(value) ||
            throw(ArgumentError("AxialChainProfile: z_m[$i] must be finite."))
    end
    for i in 2:n
        profile.z_m[i] > profile.z_m[i - 1] ||
            throw(ArgumentError("AxialChainProfile: z_m must be strictly increasing (failed at indices $(i - 1), $i)."))
    end

    for (name, values) in (("dx_m", profile.dx_m),
                           ("te_K", profile.te_K))
        for (i, value) in pairs(values)
            isfinite(value) ||
                throw(ArgumentError("AxialChainProfile: $(name)[$i] must be finite."))
            value > 0.0 ||
                throw(ArgumentError("AxialChainProfile: $(name)[$i] must be strictly positive."))
        end
    end

    isempty(profile.species_u_m_s) &&
        throw(ArgumentError("AxialChainProfile: species_u_m_s must contain at least one species."))
    for (name, values) in pairs(profile.species_u_m_s)
        lowered = lowercase(strip(name))
        lowered in ("e", "e-", "electron") &&
            throw(ArgumentError("AxialChainProfile: species_u_m_s must not include electron species."))
        if length(values) != n
            throw(ArgumentError("AxialChainProfile species_u_m_s[$name] length $(length(values)) does not match required profile length $n."))
        end
        for (i, value) in pairs(values)
            isfinite(value) ||
                throw(ArgumentError("AxialChainProfile: species_u_m_s[$name][$i] must be finite."))
            value > 0.0 ||
                throw(ArgumentError("AxialChainProfile: species_u_m_s[$name][$i] must be strictly positive."))
        end
    end

    for (name, values) in pairs(profile.diagnostics)
        if length(values) != n
            throw(ArgumentError("AxialChainProfile diagnostic `$name` length $(length(values)) does not match required profile length $n."))
        end
        for (i, value) in pairs(values)
            isfinite(value) ||
                throw(ArgumentError("AxialChainProfile: diagnostic `$name[$i]` must be finite."))
        end
    end

    if profile.wall_profile !== nothing
        wall_profile = profile.wall_profile
        length(wall_profile.a_wall_over_v_m_inv) == n ||
            throw(ArgumentError("AxialChainProfile: wall_profile.a_wall_over_v_m_inv length $(length(wall_profile.a_wall_over_v_m_inv)) does not match required profile length $n."))
        for (i, value) in pairs(wall_profile.a_wall_over_v_m_inv)
            isfinite(value) && value > 0.0 ||
                throw(ArgumentError("AxialChainProfile: wall_profile.a_wall_over_v_m_inv[$i] must be finite and strictly positive."))
        end

        for (name, values) in (("channel_gap_m", wall_profile.channel_gap_m),
                               ("wall_temperature_K", wall_profile.wall_temperature_K),
                               ("ion_edge_to_center_ratio",
                                wall_profile.ion_edge_to_center_ratio))
            values === nothing && continue
            length(values) == n ||
                throw(ArgumentError("AxialChainProfile: wall_profile.$name length $(length(values)) does not match required profile length $n."))
            for (i, value) in pairs(values)
                isfinite(value) && value > 0.0 ||
                    throw(ArgumentError("AxialChainProfile: wall_profile.$name[$i] must be finite and strictly positive."))
            end
        end
    end

    source_idx = profile.inlet.source_compact_index
    source_idx <= n ||
        throw(ArgumentError("AxialChainProfile: inlet.source_compact_index ($(source_idx)) must be <= retained profile length $(n)."))

    return true
end
