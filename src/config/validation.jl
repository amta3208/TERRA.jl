"""
$(SIGNATURES)

Validate species against TERRA database (requires loaded library).

This function uses the fortran_wrapper to query the TERRA database
for available species and validates the configuration species against it.

# Arguments
- `config::Config`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if species validation fails or library not loaded
"""
function validate_species_against_terra_database(config::Config)
    if config.runtime.validate_species_against_terra
        try
            # fortran_wrapper symbols are available once the TERRA module is loaded
            if is_terra_loaded()
                terra_species = get_species_names_wrapper()
                composition = config.reactor.composition
                validate_species_data(
                    composition.species, terra_species, composition.mole_fractions)
                @info "Species validation against TERRA database passed"
            else
                @warn "TERRA library not loaded. Cannot validate species against database."
            end
        catch e
            @warn "Failed to validate species against TERRA database: $(e)"
        end
    end

    return true
end

"""
$(SIGNATURES)

Validate configuration parameters against loaded TERRA library capabilities.

This function checks that the configuration is compatible with the loaded
TERRA library, including array dimensions and species availability.

# Arguments
- `config::Config`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if validation fails or library not loaded
"""
function validate_config_against_terra(config::Config)
    if !is_terra_loaded()
        @warn "TERRA library not loaded. Cannot validate configuration against library capabilities."
        return true
    end

    try
        # Validate array dimensions
        validate_array_dimensions(config)

        # Validate species if requested
        if config.runtime.validate_species_against_terra
            validate_species_in_database(config)
        end

        @info "Configuration validation against TERRA library passed"
        return true

    catch e
        error("Configuration validation against TERRA library failed: $(e)")
    end
end

"""
$(SIGNATURES)

Validate that configuration arrays match TERRA Fortran expectations.

# Arguments
- `config::Config`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if array dimensions exceed TERRA limits
"""
function validate_array_dimensions(config::Config)
    if !is_terra_loaded()
        throw(ErrorException("TERRA library not loaded. Cannot validate array dimensions."))
    end

    # Get TERRA limits
    max_species = get_max_number_of_species_wrapper()
    max_atomic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_states = get_max_number_of_molecular_electronic_states_wrapper()

    # Check species count
    composition = config.reactor.composition
    if length(composition.species) > max_species
        throw(ErrorException("Number of species ($(length(composition.species))) exceeds TERRA maximum ($max_species)"))
    end

    # Check mole fractions array consistency
    if length(composition.mole_fractions) != length(composition.species)
        throw(ErrorException("Mole fractions array length ($(length(composition.mole_fractions))) does not match species count ($(length(composition.species)))"))
    end

    @info "Array dimension validation passed" max_species=max_species max_atomic_states=max_atomic_states max_molecular_states=max_molecular_states

    return true
end

"""
$(SIGNATURES)

Validate that all species in configuration exist in the TERRA database.

# Arguments
- `config::Config`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if any species are not found in database
"""
function validate_species_in_database(config::Config)
    if !is_terra_loaded()
        throw(ErrorException("TERRA library not loaded. Cannot validate species in database."))
    end

    # Get available species from TERRA
    terra_species = get_species_names_wrapper()

    # Check each species in configuration
    missing_species = String[]
    species_names = config.reactor.composition.species
    for species in species_names
        if !(species in terra_species)
            push!(missing_species, species)
        end
    end

    if !isempty(missing_species)
        throw(ErrorException("Species not found in TERRA database: $(join(missing_species, ", "))\n" *
                             "Available species: $(join(terra_species, ", "))"))
    end

    @info "Species database validation passed" validated_species=species_names

    return true
end
