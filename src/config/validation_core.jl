"""
$(SIGNATURES)

Validate TERRA configuration parameters.

# Arguments
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Initial mole fractions
- `total_number_density::Float64`: Total number density
- `temperatures::TemperatureConfig`: Temperature configuration
- `time_params::TimeIntegrationConfig`: Time integration parameters
- `case_path::String`: Working directory path
- `unit_system::Symbol`: Unit system specification

# Throws
- `ArgumentError` if validation fails
"""
function validate_config(species::Vector{String},
        mole_fractions::Vector{Float64},
        total_number_density::Float64,
        temperatures::TemperatureConfig,
        time_params::TimeIntegrationConfig,
        case_path::String,
        unit_system::Symbol)

    # Check species and mole fractions consistency
    if length(species) != length(mole_fractions)
        throw(ArgumentError("Species and mole_fractions arrays must have same length"))
    end

    if isempty(species)
        throw(ArgumentError("At least one species must be specified"))
    end

    # Check mole fractions sum to 1
    if abs(sum(mole_fractions) - 1.0) > 1e-10
        throw(ArgumentError("Mole fractions must sum to 1.0, got $(sum(mole_fractions))"))
    end

    # Check for negative mole fractions
    if any(mole_fractions .< 0)
        throw(ArgumentError("Mole fractions must be non-negative"))
    end

    # Check total number density
    if total_number_density <= 0
        throw(ArgumentError("Total number density must be positive"))
    end

    # Check for duplicate species
    if length(unique(species)) != length(species)
        throw(ArgumentError("Duplicate species names found"))
    end

    # Validate species names format
    for species_name in species
        if isempty(strip(species_name))
            throw(ArgumentError("Species names cannot be empty"))
        end
    end

    # Validate case path
    if !isdir(case_path)
        throw(ArgumentError("Case path directory does not exist: $case_path"))
    end

    # Validate unit system
    if !(unit_system in [:SI, :CGS])
        throw(ArgumentError("Unit system must be :SI or :CGS, got :$unit_system"))
    end

    return true
end
