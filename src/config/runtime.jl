"""
$(SIGNATURES)

Runtime and I/O configuration controls.
"""
const _LOGGING_CONSOLE_MODES = (:quiet, :minimal, :verbose)
const _LOGGING_PROGRESS_MODES = (:auto, :off, :summary)
const _LOGGING_STREAM_MODES = (:off, :file, :console, :both)

function _validate_logging_mode(name::AbstractString, value::Symbol, allowed_modes)
    value in allowed_modes ||
        throw(ArgumentError("$(name) must be one of $(collect(allowed_modes)), got :$value"))
    return value
end

function _normalize_log_dir(log_dir::Union{Nothing, AbstractString})
    if log_dir === nothing
        return nothing
    end

    path = String(strip(log_dir))
    isempty(path) &&
        throw(ArgumentError("log_dir must be a non-empty path when provided."))
    return normpath(expanduser(path))
end

_logging_mode_file_only(mode::Symbol) = mode == :off ? :off : :file

"""
$(SIGNATURES)

Logging policy controls for terminal and file-backed runtime diagnostics.
"""
struct LoggingConfig
    console_mode::Symbol
    progress_mode::Symbol
    native_stream_mode::Symbol
    integration_detail_mode::Symbol
    chain_detail_mode::Symbol
    log_dir::Union{Nothing, String}

    function LoggingConfig(; console_mode::Symbol = :minimal,
                           progress_mode::Symbol = :auto,
                           native_stream_mode::Symbol = :file,
                           integration_detail_mode::Symbol = :file,
                           chain_detail_mode::Symbol = :file,
                           log_dir::Union{Nothing, AbstractString} = nothing)
        _validate_logging_mode("console_mode", console_mode, _LOGGING_CONSOLE_MODES)
        _validate_logging_mode("progress_mode", progress_mode, _LOGGING_PROGRESS_MODES)
        _validate_logging_mode("native_stream_mode", native_stream_mode,
                               _LOGGING_STREAM_MODES)
        _validate_logging_mode("integration_detail_mode", integration_detail_mode,
                               _LOGGING_STREAM_MODES)
        _validate_logging_mode("chain_detail_mode", chain_detail_mode,
                               _LOGGING_STREAM_MODES)

        return new(console_mode,
                   progress_mode,
                   native_stream_mode,
                   integration_detail_mode,
                   chain_detail_mode,
                   _normalize_log_dir(log_dir))
    end
end

function _with_logging_config(logging::LoggingConfig;
                              console_mode::Symbol = logging.console_mode,
                              progress_mode::Symbol = logging.progress_mode,
                              native_stream_mode::Symbol = logging.native_stream_mode,
                              integration_detail_mode::Symbol = logging.integration_detail_mode,
                              chain_detail_mode::Symbol = logging.chain_detail_mode,
                              log_dir::Union{Nothing, AbstractString} = logging.log_dir)
    return LoggingConfig(; console_mode = console_mode,
                         progress_mode = progress_mode,
                         native_stream_mode = native_stream_mode,
                         integration_detail_mode = integration_detail_mode,
                         chain_detail_mode = chain_detail_mode,
                         log_dir = log_dir)
end

"""
$(SIGNATURES)

Runtime controls for TERRA wrapper execution.
"""
struct RuntimeConfig
    database_path::String
    case_path::String
    unit_system::Symbol
    validate_species_against_terra::Bool
    print_source_terms::Bool
    write_native_state_files::Bool
    logging::LoggingConfig

    function RuntimeConfig(;
                           database_path::String = "../../databases/n2/elec_sts_expanded_electron_fits",
                           case_path::String = pwd(),
                           unit_system::Symbol = :CGS,
                           validate_species_against_terra::Bool = false,
                           print_source_terms::Bool = true,
                           write_native_state_files::Bool = false,
                           logging::LoggingConfig = LoggingConfig())
        if !isdir(case_path)
            throw(ArgumentError("Case path directory does not exist: $case_path"))
        end
        if !(unit_system in [:SI, :CGS])
            throw(ArgumentError("Unit system must be :SI or :CGS, got :$unit_system"))
        end

        return new(database_path, case_path, unit_system,
                   validate_species_against_terra, print_source_terms,
                   write_native_state_files, logging)
    end
end
