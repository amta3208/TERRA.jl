const DEFAULT_PROGRESS_FRACTION_STEP = 0.1

# -----------------------------------------------------------------------------
# Log targets and entries
# -----------------------------------------------------------------------------

abstract type AbstractLogTarget end

struct RunLog <: AbstractLogTarget end
struct NativeLog <: AbstractLogTarget end

const RUN_LOG = RunLog()
const NATIVE_LOG = NativeLog()

abstract type AbstractLogEntry end

struct EventEntry <: AbstractLogEntry
    level::Symbol
    message::String
    console::Symbol
    fields::Vector{Pair{Symbol, Any}}
end

struct ExceptionEntry <: AbstractLogEntry
    level::Symbol
    message::String
    exception::Any
    console::Symbol
    fields::Vector{Pair{Symbol, Any}}
end

struct IntegrationDetailEntry <: AbstractLogEntry
    text::String
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

function EventEntry(level::Symbol, message::AbstractString;
                    console::Symbol = :never, fields...)
    return EventEntry(level,
                      String(message),
                      console,
                      _field_pairs(fields))
end

function ExceptionEntry(level::Symbol, message::AbstractString, exception;
                        console::Symbol = :never,
                        fields...)
    return ExceptionEntry(level,
                          String(message),
                          exception,
                          console,
                          _field_pairs(fields))
end

IntegrationDetailEntry(text::AbstractString) = IntegrationDetailEntry(String(text))

# -----------------------------------------------------------------------------
# File and console routing
# -----------------------------------------------------------------------------

_writes_file(mode::Symbol) = mode == :file || mode == :both
_writes_console(mode::Symbol) = mode == :console || mode == :both

log_path(::RunLog, runtime::RuntimeConfig) = run_log_path(runtime)
log_path(::NativeLog, runtime::RuntimeConfig) = native_log_path(runtime)

stream_mode(::RunLog, ::LoggingConfig) = :file
stream_mode(::NativeLog, logging::LoggingConfig) = logging.native_stream_mode

function _append_log_text(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    open(path, "a") do io
        write(io, text)
        endswith(text, "\n") || write(io, "\n")
    end
    return nothing
end

prepare!(::RunLog, runtime::RuntimeConfig) = (ensure_case_layout!(runtime); nothing)

function _native_console_level(logging::LoggingConfig)
    if !_writes_console(logging.native_stream_mode) || logging.console_mode == :quiet
        return API_NATIVE_LOG_OFF
    elseif logging.console_mode == :verbose
        return API_NATIVE_LOG_VERBOSE
    end
    return API_NATIVE_LOG_MINIMAL
end

function _native_file_level(logging::LoggingConfig)
    _writes_file(logging.native_stream_mode) ?
    API_NATIVE_LOG_VERBOSE :
    API_NATIVE_LOG_OFF
end

function prepare!(::NativeLog, runtime::RuntimeConfig)
    console_level = _native_console_level(runtime.logging)
    file_level = _native_file_level(runtime.logging)
    ensure_case_layout!(runtime)
    path = file_level == API_NATIVE_LOG_OFF ? nothing :
           begin
               native_path = log_path(NATIVE_LOG, runtime)
               mkpath(dirname(native_path))
               open(native_path, "w") do io
                   write(io, "")
               end
               String(native_path)
           end
    return (console_level = console_level, file_level = file_level, log_path = path)
end

function _field_pairs(fields)
    pairs = Pair{Symbol, Any}[]
    for field in fields
        push!(pairs, Symbol(first(field)) => last(field))
    end
    return pairs
end

function _format_log_fields(fields::AbstractVector{<:Pair})
    isempty(fields) && return ""
    return join(("$(first(field))=$(repr(last(field)))" for field in fields), " ")
end

_timestamp_string() = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")

function _format_run_log_line(level::Symbol, message::AbstractString,
                              fields::AbstractVector{<:Pair})
    field_text = _format_log_fields(fields)
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
    end
    throw(ArgumentError("Unsupported console visibility mode: :$visibility"))
end

function _emit_console_text(text::AbstractString)
    print(text)
    endswith(text, "\n") || println()
    return nothing
end

function _emit_console_event(logging::LoggingConfig, entry::EventEntry)
    _console_visibility_enabled(logging, entry.console) || return nothing
    if entry.console == :minimal || isempty(entry.fields)
        println(entry.message)
        return nothing
    end

    field_text = _format_log_fields(entry.fields)
    println(isempty(field_text) ? entry.message : "$(entry.message) $(field_text)")
    return nothing
end

function _emit_console_event(logging::LoggingConfig, entry::ExceptionEntry)
    event = EventEntry(entry.level,
                       entry.message;
                       console = entry.console,
                       _exception_fields(entry)...)
    _emit_console_event(logging, event)
    return nothing
end

_render(entry::AbstractLogEntry) = sprint(show, MIME("text/plain"), entry)

function Base.show(io::IO, ::MIME"text/plain", entry::EventEntry)
    print(io, _format_run_log_line(entry.level, entry.message, entry.fields))
    return nothing
end

function _exception_fields(entry::ExceptionEntry)
    return [entry.fields...,
            :exception => sprint(showerror, entry.exception)]
end

function Base.show(io::IO, ::MIME"text/plain", entry::ExceptionEntry)
    print(io, _format_run_log_line(entry.level, entry.message, _exception_fields(entry)))
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", entry::IntegrationDetailEntry) = print(io, entry.text)

# -----------------------------------------------------------------------------
# Runtime entrypoints and progress reporting
# -----------------------------------------------------------------------------

function emit!(::RunLog, runtime::RuntimeConfig, entry::EventEntry)
    prepare!(RUN_LOG, runtime)
    _append_log_text(log_path(RUN_LOG, runtime), _render(entry))
    _emit_console_event(runtime.logging, entry)
    return nothing
end

function emit!(::RunLog, runtime::RuntimeConfig, entry::ExceptionEntry)
    prepare!(RUN_LOG, runtime)
    _append_log_text(log_path(RUN_LOG, runtime), _render(entry))
    _emit_console_event(runtime.logging, entry)
    return nothing
end

function emit!(::RunLog, runtime::RuntimeConfig, entry::IntegrationDetailEntry)
    mode = runtime.logging.integration_detail_mode
    mode == :off && return nothing

    if _writes_file(mode)
        header = _render(EventEntry(:detail, "reactor integration snapshot"))
        _append_log_text(log_path(RUN_LOG, runtime), string(header, "\n", _render(entry)))
    end
    if _writes_console(mode)
        _emit_console_text(entry.text)
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

function _progress_reporter(runtime::RuntimeConfig, tlim::Real)
    return _progress_mode_enabled(runtime.logging) ? FractionProgressReporter(tlim) :
           nothing
end

function _progress_message(progress::FractionProgressReporter, t::Real)
    fraction = progress.tlim > 0.0 ? clamp(Float64(t) / progress.tlim, 0.0, 1.0) : 1.0
    pct = round(Int, 100 * fraction)
    return @sprintf("-> %3d%% (t = %.3e / %.3e s)", pct, Float64(t), progress.tlim)
end

function _report_progress!(runtime::RuntimeConfig, progress::FractionProgressReporter,
                           t::Real)
    emit!(RUN_LOG, runtime,
          EventEntry(:info, _progress_message(progress, t);
                     console = :minimal,
                     :t => Float64(t),
                     :tlim => progress.tlim))
    return nothing
end
