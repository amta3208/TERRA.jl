"""
$(SIGNATURES)

Return a copy of `logging` with updated logging controls.
"""
function with_logging(logging::LoggingConfig;
                      console_mode::Symbol = logging.console_mode,
                      progress_mode::Symbol = logging.progress_mode,
                      native_stream_mode::Symbol = logging.native_stream_mode,
                      integration_detail_mode::Symbol = logging.integration_detail_mode,
                      chain_detail_mode::Symbol = logging.chain_detail_mode,
                      log_dir::Union{Nothing, AbstractString} = logging.log_dir)
    return _with_logging_config(logging;
                                console_mode = console_mode,
                                progress_mode = progress_mode,
                                native_stream_mode = native_stream_mode,
                                integration_detail_mode = integration_detail_mode,
                                chain_detail_mode = chain_detail_mode,
                                log_dir = log_dir)
end

"""
$(SIGNATURES)

Return a copy of `config` with an updated logging block.
"""
function with_logging(config::Config;
                      console_mode::Symbol = config.runtime.logging.console_mode,
                      progress_mode::Symbol = config.runtime.logging.progress_mode,
                      native_stream_mode::Symbol = config.runtime.logging.native_stream_mode,
                      integration_detail_mode::Symbol = config.runtime.logging.integration_detail_mode,
                      chain_detail_mode::Symbol = config.runtime.logging.chain_detail_mode,
                      log_dir::Union{Nothing, AbstractString} = config.runtime.logging.log_dir)
    logging = with_logging(config.runtime.logging;
                           console_mode = console_mode,
                           progress_mode = progress_mode,
                           native_stream_mode = native_stream_mode,
                           integration_detail_mode = integration_detail_mode,
                           chain_detail_mode = chain_detail_mode,
                           log_dir = log_dir)
    return with_runtime(config; logging = logging)
end

"""
$(SIGNATURES)

Return a copy of `runtime` with an updated runtime block.
"""
function with_runtime(runtime::RuntimeConfig;
                      database_path::AbstractString = runtime.database_path,
                      case_path::AbstractString = runtime.case_path,
                      unit_system::Symbol = runtime.unit_system,
                      validate_species_against_terra::Bool = runtime.validate_species_against_terra,
                      print_source_terms::Bool = runtime.print_source_terms,
                      write_native_state_files::Bool = runtime.write_native_state_files,
                      logging::LoggingConfig = runtime.logging,
                      write_native_outputs::Union{Nothing, Bool} = nothing,
                      print_integration_output::Union{Nothing, Bool} = nothing)
    native_state_files = write_native_outputs === nothing ? write_native_state_files :
                         write_native_outputs
    return RuntimeConfig(;
                         database_path = String(database_path),
                         case_path = String(case_path),
                         unit_system = unit_system,
                         validate_species_against_terra = validate_species_against_terra,
                         print_source_terms = print_source_terms,
                         write_native_state_files = native_state_files,
                         logging = logging,
                         print_integration_output = print_integration_output)
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
                      write_native_state_files::Bool = config.runtime.write_native_state_files,
                      logging::LoggingConfig = config.runtime.logging,
                      write_native_outputs::Union{Nothing, Bool} = nothing,
                      print_integration_output::Union{Nothing, Bool} = nothing)
    runtime = with_runtime(config.runtime;
                           database_path = database_path,
                           case_path = case_path,
                           unit_system = unit_system,
                           validate_species_against_terra = validate_species_against_terra,
                           print_source_terms = print_source_terms,
                           write_native_state_files = write_native_state_files,
                           logging = logging,
                           write_native_outputs = write_native_outputs,
                           print_integration_output = print_integration_output)

    return Config(;
                  reactor = config.reactor,
                  models = config.models,
                  sources = config.sources,
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
                              space = config.numerics.space)

    return Config(;
                  reactor = config.reactor,
                  models = config.models,
                  sources = config.sources,
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
                                write_native_state_files = config.runtime.write_native_state_files,
                                logging = config.runtime.logging)

    return Config(;
                  reactor = new_reactor,
                  models = config.models,
                  sources = config.sources,
                  numerics = config.numerics,
                  runtime = new_runtime)
end
