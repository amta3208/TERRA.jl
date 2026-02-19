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
