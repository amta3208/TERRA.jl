const AVOGADRO = 6.02214076e23  # Avogadro's number (1/mol)

# Common molecular weights (g/mol)
const MOLECULAR_WEIGHT_DB = Dict{String, Float64}(
    "N" => 14.007,
    "N2" => 28.014,
    "N+" => 14.007,
    "N2+" => 28.014,
    "E-" => 5.485799e-4,  # electron mass
    "Ar" => 39.948,
    "Ar+" => 39.948,
    "Xe" => 131.293,
    "Xe+" => 131.293,
    "Kr" => 83.798,
    "Kr+" => 83.798,
    "O" => 15.999,
    "O2" => 31.998,
    "O+" => 15.999,
    "O2+" => 31.998
)

"""
$(SIGNATURES)

Get molecular weights for common species (g/mol).

# Arguments
- `species::Vector{String}`: Species names

# Returns
- `Vector{Float64}`: Molecular weights in g/mol
"""
function get_molecular_weights(species::Vector{String})
    weights = Float64[]
    for sp in species
        if haskey(MOLECULAR_WEIGHT_DB, sp)
            push!(weights, MOLECULAR_WEIGHT_DB[sp])
        else
            error("Molecular weight not found for species: $(sp)")
        end
    end
    return weights
end

"""
$(SIGNATURES)

Validate species data consistency between Julia and TERRA.

# Arguments
- `species_names::Vector{String}`: Species names from configuration
- `terra_species::Vector{String}`: Species names from TERRA
- `densities::Vector{Float64}`: Species densities

# Returns
- `true` if validation passes, throws error otherwise
"""
function validate_species_data(species_names::Vector{String},
        terra_species::Vector{String},
        densities::Vector{Float64})

    # Check array lengths match
    if length(species_names) != length(densities)
        error("Species names and densities arrays have different lengths: $(length(species_names)) vs $(length(densities))")
    end

    # Check for negative densities
    if any(densities .< 0)
        negative_indices = findall(densities .< 0)
        error("Negative densities found for species: $(species_names[negative_indices])")
    end

    # Check species names are valid
    for (i, name) in enumerate(species_names)
        if !(name in terra_species)
            error("Species '$(name)' not found in TERRA database. Available species: $(terra_species)")
        end
    end

    # Check for very small densities that might cause numerical issues
    min_density = 1e-30
    small_density_indices = findall(densities .< min_density)
    if !isempty(small_density_indices)
        @warn "Very small densities detected for species: $(species_names[small_density_indices]). This may cause numerical issues."
    end

    return true
end

"""
$(SIGNATURES)

Create species mapping between HallThruster.jl and TERRA conventions.

# Arguments
- `ht_species::Vector{String}`: Species names in HallThruster.jl format
- `terra_species::Vector{String}`: Species names in TERRA format

# Returns
- Dictionary mapping HallThruster.jl species names to TERRA species names
"""
function create_species_mapping(ht_species::Vector{String}, terra_species::Vector{String})
    mapping = Dict{String, String}()

    # Common mappings (this may need to be expanded based on actual usage)
    common_mappings = Dict(
        "N" => "N",
        "N2" => "N2",
        "N+" => "N+",
        "N2+" => "N2+",
        "e-" => "E-",
        "E-" => "E-",
        "Ar" => "Ar",
        "Ar+" => "Ar+",
        "Xe" => "Xe",
        "Xe+" => "Xe+",
        "Kr" => "Kr",
        "Kr+" => "Kr+"
    )

    for ht_name in ht_species
        if haskey(common_mappings, ht_name)
            terra_name = common_mappings[ht_name]
            if terra_name in terra_species
                mapping[ht_name] = terra_name
            else
                error("TERRA species '$(terra_name)' not found in database for HallThruster species '$(ht_name)'")
            end
        else
            # Try direct mapping
            if ht_name in terra_species
                mapping[ht_name] = ht_name
            else
                error("No mapping found for HallThruster species '$(ht_name)' to TERRA database")
            end
        end
    end

    return mapping
end

"""
$(SIGNATURES)

Convert mole fractions to mass densities.

# Arguments
- `mole_fractions::Vector{Float64}`: Species mole fractions
- `molecular_weights::Vector{Float64}`: Species molecular weights (g/mol)
- `total_number_density::Float64`: Total number density (1/cm³)

# Returns
- Vector of mass densities (g/cm³)
"""
function mole_fractions_to_mass_densities(mole_fractions::AbstractVector{<:Real},
        molecular_weights::AbstractVector{<:Real},
        total_number_density::Real)
    if length(mole_fractions) != length(molecular_weights)
        error("Mole fractions and molecular weights arrays must have same length")
    end

    # Per-element validity within numeric tolerance
    if any(mole_fractions .< -eps(Float64)) || any(mole_fractions .> 1 + eps(Float64))
        error("Mole fractions must each be in [0,1] within tolerance")
    end

    if abs(sum(mole_fractions) - 1.0) > 1e-10
        error("Mole fractions must sum to 1.0, got $(sum(mole_fractions))")
    end

    if any(molecular_weights .<= 0)
        error("Molecular weights must be positive")
    end

    # Convert to number densities
    number_densities = mole_fractions .* total_number_density

    # Convert to mass densities (g/cm³)
    # number_density (1/cm³) * molecular_weight (g/mol) / Avogadro (1/mol) = mass_density (g/cm³)
    mass_densities = number_densities .* molecular_weights ./ AVOGADRO

    return mass_densities
end

"""
$(SIGNATURES)

Convert mass densities to mole fractions.

# Arguments
- `mass_densities::Vector{Float64}`: Species mass densities (g/cm³)
- `molecular_weights::Vector{Float64}`: Species molecular weights (g/mol)

# Returns
- Vector of mole fractions
"""
function mass_densities_to_mole_fractions(mass_densities::AbstractVector{<:Real},
        molecular_weights::AbstractVector{<:Real})
    if length(mass_densities) != length(molecular_weights)
        error("Mass densities and molecular weights arrays must have same length")
    end

    if any(mass_densities .< 0)
        error("Mass densities must be non-negative")
    end

    if any(molecular_weights .<= 0)
        error("Molecular weights must be positive")
    end

    # Convert to number densities
    number_densities = mass_densities .* AVOGADRO ./ molecular_weights

    # Convert to mole fractions
    total_number_density = sum(number_densities)

    if isapprox(total_number_density, 0.0; atol = eps(Float64))
        error("Total number density is zero - cannot compute mole fractions")
    end

    mole_fractions = number_densities ./ total_number_density

    return mole_fractions
end
