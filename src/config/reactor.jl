"""
$(SIGNATURES)

Reactor composition inputs for a simulation case.

# Fields
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Species mole fractions (must sum to 1)
- `total_number_density::Float64`: Total number density
"""
struct ReactorComposition
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density::Float64

    function ReactorComposition(species, mole_fractions, total_number_density)
        species_vec = String.(species)
        mole_frac_vec = Float64.(mole_fractions)
        n_tot = Float64(total_number_density)

        if length(species_vec) != length(mole_frac_vec)
            throw(ArgumentError("Species and mole_fractions arrays must have same length"))
        end
        if isempty(species_vec)
            throw(ArgumentError("At least one species must be specified"))
        end
        if any(mole_frac_vec .< 0)
            throw(ArgumentError("Mole fractions must be non-negative"))
        end
        if abs(sum(mole_frac_vec) - 1.0) > 1e-10
            throw(ArgumentError("Mole fractions must sum to 1.0, got $(sum(mole_frac_vec))"))
        end
        if n_tot <= 0
            throw(ArgumentError("Total number density must be positive"))
        end
        if length(unique(species_vec)) != length(species_vec)
            throw(ArgumentError("Duplicate species names found"))
        end
        for name in species_vec
            if isempty(strip(name))
                throw(ArgumentError("Species names cannot be empty"))
            end
        end

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
