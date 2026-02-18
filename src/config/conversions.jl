"""
$(SIGNATURES)

Convert legacy `TERRAConfig` to nested `Config`.
"""
function to_config(config::TERRAConfig;
        residence_time::Union{Nothing, ResidenceTimeConfig} = nothing)
    reactor = ReactorConfig(;
        composition = ReactorComposition(;
            species = config.species,
            mole_fractions = config.mole_fractions,
            total_number_density = config.total_number_density),
        thermal = ReactorThermalState(;
            Tt = config.temperatures.Tt,
            Tv = config.temperatures.Tv,
            Tee = config.temperatures.Tee,
            Te = config.temperatures.Te))

    models = ModelConfig(;
        physics = config.physics,
        processes = config.processes)

    numerics = NumericsConfig(;
        time = TimeConfig(;
            dt = config.time_params.dt,
            dt_output = config.time_params.dtm,
            duration = config.time_params.tlim,
            nstep = config.time_params.nstep,
            method = config.time_params.method),
        residence_time = residence_time === nothing ? ResidenceTimeConfig(; enabled = false) :
                         residence_time)

    runtime = RuntimeConfig(;
        database_path = config.database_path,
        case_path = config.case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = config.validate_species_against_terra,
        print_source_terms = config.print_source_terms,
        write_native_outputs = config.write_native_outputs,
        print_integration_output = config.print_integration_output)

    return Config(;
        reactor = reactor,
        models = models,
        numerics = numerics,
        runtime = runtime)
end

"""
$(SIGNATURES)

Convert nested `Config` to legacy `TERRAConfig`.
"""
function to_legacy_config(config::Config)
    return TERRAConfig(
        species = config.reactor.composition.species,
        mole_fractions = config.reactor.composition.mole_fractions,
        total_number_density = config.reactor.composition.total_number_density,
        temperatures = TemperatureConfig(;
            Tt = config.reactor.thermal.Tt,
            Tv = config.reactor.thermal.Tv,
            Tee = config.reactor.thermal.Tee,
            Te = config.reactor.thermal.Te),
        time_params = TimeIntegrationConfig(;
            dt = config.numerics.time.dt,
            dtm = config.numerics.time.dt_output,
            tlim = config.numerics.time.duration,
            nstep = config.numerics.time.nstep,
            method = config.numerics.time.method),
        physics = config.models.physics,
        processes = config.models.processes,
        database_path = config.runtime.database_path,
        case_path = config.runtime.case_path,
        unit_system = config.runtime.unit_system,
        validate_species_against_terra = config.runtime.validate_species_against_terra,
        print_source_terms = config.runtime.print_source_terms,
        write_native_outputs = config.runtime.write_native_outputs,
        print_integration_output = config.runtime.print_integration_output
    )
end

"""
$(SIGNATURES)

Return a copy of `config` with an updated runtime block.
"""
function with_runtime(config::Config;
        database_path::AbstractString = config.runtime.database_path,
        case_path::AbstractString = config.runtime.case_path,
        unit_system::Symbol = config.runtime.unit_system,
        validate_species_against_terra::Bool = config.runtime.validate_species_against_terra,
        print_source_terms::Bool = config.runtime.print_source_terms,
        write_native_outputs::Bool = config.runtime.write_native_outputs,
        print_integration_output::Bool = config.runtime.print_integration_output)
    runtime = RuntimeConfig(;
        database_path = String(database_path),
        case_path = String(case_path),
        unit_system = unit_system,
        validate_species_against_terra = validate_species_against_terra,
        print_source_terms = print_source_terms,
        write_native_outputs = write_native_outputs,
        print_integration_output = print_integration_output)

    return Config(;
        reactor = config.reactor,
        models = config.models,
        numerics = config.numerics,
        runtime = runtime)
end

"""
$(SIGNATURES)

Return a copy of `config` with an updated case path.
"""
function with_case_path(config::Config, case_path::AbstractString)
    return with_runtime(config; case_path = case_path)
end

"""
$(SIGNATURES)

Return a copy of `config` with updated time-integration controls.
"""
function with_time(config::Config;
        dt::Real = config.numerics.time.dt,
        dt_output::Real = config.numerics.time.dt_output,
        duration::Real = config.numerics.time.duration,
        nstep::Integer = config.numerics.time.nstep,
        method::Integer = config.numerics.time.method)
    time = TimeConfig(;
        dt = dt,
        dt_output = dt_output,
        duration = duration,
        nstep = nstep,
        method = method)
    numerics = NumericsConfig(;
        time = time,
        solver = config.numerics.solver,
        space = config.numerics.space,
        residence_time = config.numerics.residence_time)
    return Config(;
        reactor = config.reactor,
        models = config.models,
        numerics = numerics,
        runtime = config.runtime)
end

"""
$(SIGNATURES)

Legacy wrapper for `with_runtime`.
"""
function with_runtime(config::TERRAConfig; kwargs...)
    return to_legacy_config(with_runtime(to_config(config); kwargs...))
end

"""
$(SIGNATURES)

Legacy wrapper for `with_case_path`.
"""
function with_case_path(config::TERRAConfig, case_path::AbstractString)
    return to_legacy_config(with_case_path(to_config(config), case_path))
end

"""
$(SIGNATURES)

Legacy wrapper for `with_time`.
"""
function with_time(config::TERRAConfig;
        dt::Real = config.time_params.dt,
        dtm::Real = config.time_params.dtm,
        tlim::Real = config.time_params.tlim,
        nstep::Integer = config.time_params.nstep,
        method::Integer = config.time_params.method)
    nested = with_time(to_config(config);
        dt = dt,
        dt_output = dtm,
        duration = tlim,
        nstep = nstep,
        method = method)
    return to_legacy_config(nested)
end
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
