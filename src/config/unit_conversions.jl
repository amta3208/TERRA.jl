"""
$(SIGNATURES)

Convert configuration units if needed.

# Arguments
- `config::TERRAConfig`: Configuration to convert
- `target_unit_system::Symbol`: Target unit system (:SI or :CGS)

# Returns
- `TERRAConfig`: Configuration with converted units
"""
function convert_config_units(config::TERRAConfig, target_unit_system::Symbol)
    if config.unit_system == target_unit_system
        return config  # No conversion needed
    end

    if config.unit_system == :SI && target_unit_system == :CGS
        # Convert SI to CGS
        new_total_number_density = convert_number_density_si_to_cgs(config.total_number_density)

        # Create new config with converted units
        return TERRAConfig(
            species = config.species,
            mole_fractions = config.mole_fractions,
            total_number_density = new_total_number_density,
            temperatures = config.temperatures,  # Temperature units are the same
            time_params = config.time_params,     # Time units are the same
            physics = config.physics,
            processes = config.processes,
            database_path = config.database_path,
            case_path = config.case_path,
            unit_system = target_unit_system,
            validate_species_against_terra = config.validate_species_against_terra,
            print_source_terms = config.print_source_terms,
            write_native_outputs = config.write_native_outputs,
            print_integration_output = config.print_integration_output
        )

    elseif config.unit_system == :CGS && target_unit_system == :SI
        # Convert CGS to SI
        new_total_number_density = convert_number_density_cgs_to_si(config.total_number_density)

        # Create new config with converted units
        return TERRAConfig(
            species = config.species,
            mole_fractions = config.mole_fractions,
            total_number_density = new_total_number_density,
            temperatures = config.temperatures,  # Temperature units are the same
            time_params = config.time_params,     # Time units are the same
            physics = config.physics,
            processes = config.processes,
            database_path = config.database_path,
            case_path = config.case_path,
            unit_system = target_unit_system,
            validate_species_against_terra = config.validate_species_against_terra,
            print_source_terms = config.print_source_terms,
            write_native_outputs = config.write_native_outputs,
            print_integration_output = config.print_integration_output
        )
    else
        error("Unsupported unit conversion: $(config.unit_system) to $target_unit_system")
    end
end

"""
$(SIGNATURES)

Convert nested `Config` units if needed.

# Arguments
- `config::Config`: Configuration to convert
- `target_unit_system::Symbol`: Target unit system (:SI or :CGS)

# Returns
- `Config`: Configuration with converted units
"""
function convert_config_units(config::Config, target_unit_system::Symbol)
    current_unit_system = config.runtime.unit_system
    if current_unit_system == target_unit_system
        return config
    end

    composition = config.reactor.composition
    new_total_number_density = if current_unit_system == :SI && target_unit_system == :CGS
        convert_number_density_si_to_cgs(composition.total_number_density)
    elseif current_unit_system == :CGS && target_unit_system == :SI
        convert_number_density_cgs_to_si(composition.total_number_density)
    else
        error("Unsupported unit conversion: $current_unit_system to $target_unit_system")
    end

    new_reactor = ReactorConfig(;
        composition = ReactorComposition(;
            species = composition.species,
            mole_fractions = composition.mole_fractions,
            total_number_density = new_total_number_density),
        thermal = config.reactor.thermal)

    new_runtime = RuntimeConfig(;
        database_path = config.runtime.database_path,
        case_path = config.runtime.case_path,
        unit_system = target_unit_system,
        validate_species_against_terra = config.runtime.validate_species_against_terra,
        print_source_terms = config.runtime.print_source_terms,
        write_native_outputs = config.runtime.write_native_outputs,
        print_integration_output = config.runtime.print_integration_output)

    return Config(;
        reactor = new_reactor,
        models = config.models,
        numerics = config.numerics,
        runtime = new_runtime)
end
