"""
$(SIGNATURES)

Reactor composition inputs for a simulation case.

# Fields
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Species mole fractions (must sum to 1)
- `total_number_density::Float64`: Total number density
"""
function _validated_species_composition(species, mole_fractions, total_number_density;
                                        context::AbstractString,
                                        density_field::AbstractString)
    species_vec = String.(species)
    mole_frac_vec = Float64.(mole_fractions)
    n_tot = Float64(total_number_density)

    if length(species_vec) != length(mole_frac_vec)
        throw(ArgumentError("$context: species and mole_fractions arrays must have same length."))
    end
    isempty(species_vec) &&
        throw(ArgumentError("$context: at least one species must be specified."))
    any(mole_frac_vec .< 0) &&
        throw(ArgumentError("$context: mole fractions must be non-negative."))
    abs(sum(mole_frac_vec) - 1.0) > 1e-10 &&
        throw(ArgumentError("$context: mole fractions must sum to 1.0, got $(sum(mole_frac_vec))."))
    n_tot > 0.0 ||
        throw(ArgumentError("$context: $density_field must be positive."))
    length(unique(species_vec)) == length(species_vec) ||
        throw(ArgumentError("$context: duplicate species names found."))
    for name in species_vec
        isempty(strip(name)) &&
            throw(ArgumentError("$context: species names cannot be empty."))
    end

    return species_vec, mole_frac_vec, n_tot
end

struct ReactorComposition
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density::Float64

    function ReactorComposition(species, mole_fractions, total_number_density)
        species_vec, mole_frac_vec, n_tot = _validated_species_composition(species,
                                                                           mole_fractions,
                                                                           total_number_density;
                                                                           context = "ReactorComposition",
                                                                           density_field = "total_number_density")
        return new(species_vec, mole_frac_vec, n_tot)
    end
end

function ReactorComposition(; species, mole_fractions, total_number_density)
    ReactorComposition(species, mole_fractions, total_number_density)
end

"""
$(SIGNATURES)

Thermal state for the reactor.

# Fields
- `Tt::Float64`: Translational temperature (K)
- `Tv::Float64`: Vibrational temperature (K)
- `Tee::Float64`: Electron-electronic temperature (K)
- `Te::Float64`: Electron temperature (K)
"""
struct ReactorThermalState
    Tt::Float64
    Tv::Float64
    Tee::Float64
    Te::Float64

    function ReactorThermalState(Tt, Tv, Tee, Te)
        if any([Tt, Tv, Tee, Te] .<= 0)
            throw(ArgumentError("All temperatures must be positive"))
        end
        return new(Float64(Tt), Float64(Tv), Float64(Tee), Float64(Te))
    end
end

ReactorThermalState(; Tt, Tv, Tee, Te) = ReactorThermalState(Tt, Tv, Tee, Te)

"""
$(SIGNATURES)

Combined reactor-state configuration.
"""
struct ReactorConfig
    composition::ReactorComposition
    thermal::ReactorThermalState
end

function ReactorConfig(; composition::ReactorComposition,
                       thermal::ReactorThermalState)
    return ReactorConfig(composition, thermal)
end
