const DEFAULT_PROGRESS_FRACTION_STEP = 0.1
const ACTIVE_RUNTIME_FOR_LOGGING = Ref{Union{Nothing, RuntimeConfig}}(nothing)

_stream_mode_writes_file(mode::Symbol) = mode == :file || mode == :both
_stream_mode_writes_console(mode::Symbol) = mode == :console || mode == :both

function _truncate_log_file(path::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, "")
    end
    return String(path)
end

function _runtime_with_case_path(runtime::RuntimeConfig, case_path::AbstractString)
    runtime.case_path == case_path && return runtime
    return RuntimeConfig(; database_path = runtime.database_path,
                         case_path = String(case_path),
                         unit_system = runtime.unit_system,
                         validate_species_against_terra = runtime.validate_species_against_terra,
                         print_source_terms = runtime.print_source_terms,
                         write_native_state_files = runtime.write_native_state_files,
                         logging = runtime.logging)
end

_active_runtime_for_logging() = ACTIVE_RUNTIME_FOR_LOGGING[]

function _set_active_runtime_for_logging!(runtime::RuntimeConfig)
    ACTIVE_RUNTIME_FOR_LOGGING[] = runtime
    return runtime
end

function _clear_active_runtime_for_logging!()
    ACTIVE_RUNTIME_FOR_LOGGING[] = nothing
    return nothing
end

function _format_log_fields(fields::Pair...)
    isempty(fields) && return ""
    return join(["$(first(field))=$(repr(last(field)))" for field in fields], " ")
end

function _timestamp_string()
    return Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")
end

function _append_log_text(path::AbstractString, text::AbstractString)
    open(path, "a") do io
        write(io, text)
        endswith(text, "\n") || write(io, "\n")
    end
    return nothing
end

function _write_run_log(runtime::RuntimeConfig, text::AbstractString)
    _ensure_log_dir(runtime)
    _append_log_text(_run_log_path(runtime), text)
    return nothing
end

function _native_console_level(runtime::RuntimeConfig)
    logging = runtime.logging
    if !_stream_mode_writes_console(logging.native_stream_mode) ||
       logging.console_mode == :quiet
        return API_NATIVE_LOG_OFF
    elseif logging.console_mode == :verbose
        return API_NATIVE_LOG_VERBOSE
    end
    return API_NATIVE_LOG_MINIMAL
end

function _native_file_level(runtime::RuntimeConfig)
    return _stream_mode_writes_file(runtime.logging.native_stream_mode) ?
           API_NATIVE_LOG_VERBOSE :
           API_NATIVE_LOG_OFF
end

function _prepare_native_logging(runtime::RuntimeConfig)
    console_level = _native_console_level(runtime)
    file_level = _native_file_level(runtime)
    log_path = nothing
    if file_level != API_NATIVE_LOG_OFF
        _ensure_log_dir(runtime)
        log_path = _truncate_log_file(_native_log_path(runtime))
    end
    return (console_level = console_level, file_level = file_level, log_path = log_path)
end

function _format_run_log_line(level::Symbol, message::AbstractString, fields::Pair...)
    field_text = _format_log_fields(fields...)
    prefix = "[$(_timestamp_string())] $(uppercase(String(level)))"
    return isempty(field_text) ? "$prefix $message" : "$prefix $message $field_text"
end

function _console_visibility_enabled(logging::LoggingConfig, visibility::Symbol)
    if visibility == :never
        return false
    elseif visibility == :minimal
        return logging.console_mode != :quiet
    elseif visibility == :verbose
        return logging.console_mode == :verbose
    else
        throw(ArgumentError("Unsupported console visibility mode: :$visibility"))
    end
end

function _emit_console_event(logging::LoggingConfig, visibility::Symbol,
                             message::AbstractString, fields::Pair...)
    _console_visibility_enabled(logging, visibility) || return nothing
    if visibility == :minimal || isempty(fields)
        println(message)
        return nothing
    end

    field_text = _format_log_fields(fields...)
    println(isempty(field_text) ? message : "$message $field_text")
    return nothing
end

function _log_run_event(runtime::RuntimeConfig, level::Symbol, message::AbstractString;
                        console::Symbol = :never, fields...)
    line = _format_run_log_line(level, message, fields...)
    _write_run_log(runtime, line)
    _emit_console_event(runtime.logging, console, message, fields...)
    return nothing
end

function _log_run_exception(runtime::RuntimeConfig, level::Symbol, message::AbstractString,
                            exception; console::Symbol = :never, fields...)
    merged_fields = (fields..., :exception => sprint(showerror, exception))
    _log_run_event(runtime, level, message; console = console, merged_fields...)
    return nothing
end

function _emit_integration_detail(runtime::RuntimeConfig, detail_text::AbstractString)
    mode = runtime.logging.integration_detail_mode
    mode == :off && return nothing

    if _stream_mode_writes_file(mode)
        header = _format_run_log_line(:detail, "0D integration snapshot")
        _write_run_log(runtime,
                       string(header, "\n", detail_text,
                              endswith(detail_text, "\n") ?
                              "" : "\n"))
    end
    if _stream_mode_writes_console(mode)
        print(detail_text)
        endswith(detail_text, "\n") || println()
    end
    return nothing
end

function _progress_mode_enabled(logging::LoggingConfig)
    if logging.progress_mode == :off
        return false
    elseif logging.progress_mode == :summary
        return true
    end
    return logging.console_mode != :quiet && stdout isa Base.TTY
end

mutable struct FractionProgressReporter
    tlim::Float64
    fraction_step::Float64
    next_fraction::Float64
end

function FractionProgressReporter(tlim::Real;
                                  fraction_step::Real = DEFAULT_PROGRESS_FRACTION_STEP)
    step = Float64(fraction_step)
    step > 0.0 || throw(ArgumentError("fraction_step must be positive."))
    return FractionProgressReporter(Float64(tlim), step, step)
end

function _progress_message(progress::FractionProgressReporter, t::Real)
    tlim = progress.tlim
    fraction = tlim > 0.0 ? clamp(Float64(t) / tlim, 0.0, 1.0) : 1.0
    pct = round(Int, 100 * fraction)
    return @sprintf("-> %3d%% (t = %.3e / %.3e s)", pct, Float64(t), tlim)
end

function _report_progress!(runtime::RuntimeConfig, progress::FractionProgressReporter,
                           t::Real)
    message = _progress_message(progress, t)
    _log_run_event(runtime, :info, message; console = :minimal, :t => Float64(t),
                   :tlim => progress.tlim)
    return nothing
end
