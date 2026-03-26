const DEFAULT_PROGRESS_FRACTION_STEP = 0.1
const CHAIN_RUN_BANNER = "\n" * "="^9 * " TERRA 1D Chain Simulation " * "="^9
const CHAIN_RUN_FOOTER = "="^(length(CHAIN_RUN_BANNER) - 1)

abstract type AbstractLogTarget end

struct RunLog <: AbstractLogTarget end
struct ChainLog <: AbstractLogTarget end
struct NativeLog <: AbstractLogTarget end

const RUN_LOG = RunLog()
const CHAIN_LOG = ChainLog()
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
    exception
    console::Symbol
    fields::Vector{Pair{Symbol, Any}}
end

struct IntegrationDetailEntry <: AbstractLogEntry
    text::String
end

struct ChainHeaderEntry <: AbstractLogEntry
    config::Config
    profile::AxialChainProfile
    marching::AxialMarchingConfig
    compact_to_source_index::Vector{Int}
end

struct ChainSegmentEntry <: AbstractLogEntry
    base_config::Config
    profile::AxialChainProfile
    compact_to_source_index::Vector{Int}
    segment_index::Int
    segment_case_path::String
    result::ReactorResult
    endpoint_reactor::Union{Nothing, ReactorConfig}
    state_cache_used::Bool
end

struct ChainResultEntry <: AbstractLogEntry
    result::ChainSimulationResult
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

EventEntry(level::Symbol, message::AbstractString;
           console::Symbol = :never, fields...) = EventEntry(level,
                                                             String(message),
                                                             console,
                                                             _field_pairs(fields))

ExceptionEntry(level::Symbol, message::AbstractString, exception;
               console::Symbol = :never,
               fields...) = ExceptionEntry(level,
                                          String(message),
                                          exception,
                                          console,
                                          _field_pairs(fields))

IntegrationDetailEntry(text::AbstractString) = IntegrationDetailEntry(String(text))

function ChainHeaderEntry(config::Config,
                          profile::AxialChainProfile,
                          marching::AxialMarchingConfig,
                          compact_to_source_index::AbstractVector{<:Integer})
    return ChainHeaderEntry(config, profile, marching, Int.(compact_to_source_index))
end

function ChainSegmentEntry(base_config::Config,
                           profile::AxialChainProfile,
                           compact_to_source_index::AbstractVector{<:Integer},
                           segment_index::Integer,
                           segment_case_path::AbstractString,
                           result::ReactorResult;
                           endpoint_reactor::Union{Nothing, ReactorConfig} = nothing,
                           state_cache_used::Bool = false)
    return ChainSegmentEntry(base_config,
                             profile,
                             Int.(compact_to_source_index),
                             Int(segment_index),
                             String(segment_case_path),
                             result,
                             endpoint_reactor,
                             state_cache_used)
end

_writes_file(mode::Symbol) = mode == :file || mode == :both
_writes_console(mode::Symbol) = mode == :console || mode == :both

log_path(::RunLog, runtime::RuntimeConfig) = run_log_path(runtime)
log_path(::NativeLog, runtime::RuntimeConfig) = native_log_path(runtime)
log_path(::ChainLog, runtime::RuntimeConfig) = chain_log_path(runtime)

stream_mode(::RunLog, ::LoggingConfig) = :file
stream_mode(::NativeLog, logging::LoggingConfig) = logging.native_stream_mode
stream_mode(::ChainLog, logging::LoggingConfig) = logging.chain_detail_mode

function _prepare_log_path(path::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, "")
    end
    return String(path)
end

function _append_log_text(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    open(path, "a") do io
        write(io, text)
        endswith(text, "\n") || write(io, "\n")
    end
    return nothing
end

function prepare!(::RunLog, runtime::RuntimeConfig)
    ensure_case_layout!(runtime)
    return nothing
end

function prepare!(::ChainLog, runtime::RuntimeConfig)
    mode = stream_mode(CHAIN_LOG, runtime.logging)
    mode == :off && return nothing
    _writes_file(mode) || return nothing
    ensure_case_layout!(runtime)
    return _prepare_log_path(log_path(CHAIN_LOG, runtime))
end

function _native_console_level(logging::LoggingConfig)
    if !_writes_console(logging.native_stream_mode) || logging.console_mode == :quiet
        return API_NATIVE_LOG_OFF
    elseif logging.console_mode == :verbose
        return API_NATIVE_LOG_VERBOSE
    end
    return API_NATIVE_LOG_MINIMAL
end

_native_file_level(logging::LoggingConfig) = _writes_file(logging.native_stream_mode) ?
                                             API_NATIVE_LOG_VERBOSE :
                                             API_NATIVE_LOG_OFF

function prepare!(::NativeLog, runtime::RuntimeConfig)
    console_level = _native_console_level(runtime.logging)
    file_level = _native_file_level(runtime.logging)
    ensure_case_layout!(runtime)
    path = file_level == API_NATIVE_LOG_OFF ? nothing :
           _prepare_log_path(log_path(NATIVE_LOG, runtime))
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

function _render(entry::AbstractLogEntry)
    return sprint(show, MIME("text/plain"), entry)
end

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

function Base.show(io::IO, ::MIME"text/plain", entry::IntegrationDetailEntry)
    print(io, entry.text)
    return nothing
end

function _chain_scalar(value)
    if value === nothing
        return "null"
    elseif value isa Symbol
        return String(value)
    elseif value isa Bool
        return lowercase(string(value))
    elseif value isa Integer
        return string(value)
    elseif value isa Real
        x = Float64(value)
        if !isfinite(x)
            return string(x)
        elseif x == 0.0
            return "0.0"
        elseif abs(x) >= 1.0e4 || abs(x) < 1.0e-3
            return @sprintf("%.6e", x)
        end
        return @sprintf("%.6f", x)
    elseif value isa AbstractString
        return String(value)
    end
    return repr(value)
end

function _write_chain_line(io::IO, indent::Integer, key::AbstractString, value)
    print(io, repeat(" ", indent), key, ": ", _chain_scalar(value), "\n")
    return nothing
end

function _number_density_to_si(value::Real, unit_system::Symbol)
    number_density = Float64(value)
    return unit_system == :SI ? number_density :
           convert_number_density_cgs_to_si(number_density)
end

function _species_number_density_pairs_m3(composition::ReactorComposition,
                                          unit_system::Symbol)
    total_number_density_m3 = _number_density_to_si(composition.total_number_density,
                                                    unit_system)
    return [name => mole_fraction * total_number_density_m3
            for (name, mole_fraction) in zip(composition.species,
                                             composition.mole_fractions)]
end

function _ordered_chain_velocity_pairs(profile::AxialChainProfile,
                                       segment_index::Integer)
    ordered_species = String[name
                             for name in profile.inlet.composition.species
                             if haskey(profile.species_u_m_s, name)]
    seen = Set(ordered_species)
    for name in sort!(collect(keys(profile.species_u_m_s)))
        name in seen && continue
        push!(ordered_species, name)
    end
    return [name => profile.species_u_m_s[name][segment_index] for name in ordered_species]
end

function Base.show(io::IO, ::MIME"text/plain", entry::ChainHeaderEntry)
    config = entry.config
    profile = entry.profile
    marching = entry.marching

    print(io, "# chain summary generated at ", _timestamp_string(), "\n")
    print(io, "chain:\n")
    _write_chain_line(io, 2, "case_path", config.runtime.case_path)
    _write_chain_line(io, 2, "unit_system", config.runtime.unit_system)
    _write_chain_line(io, 2, "segments", length(profile.z_m))
    _write_chain_line(io, 2, "chain_detail_mode", config.runtime.logging.chain_detail_mode)
    if _writes_file(stream_mode(CHAIN_LOG, config.runtime.logging))
        _write_chain_line(io, 2, "chain_log_path", chain_log_path(config.runtime))
    end
    _write_chain_line(io, 2, "handoff_mode", marching.handoff_mode)
    _write_chain_line(io, 2, "termination_mode", marching.termination_mode)
    _write_chain_line(io, 2, "is_isothermal_teex", marching.is_isothermal_teex)
    _write_chain_line(io, 2, "compact_to_source_index", repr(entry.compact_to_source_index))
    _write_chain_line(io, 2, "species", repr(config.reactor.composition.species))

    print(io, "  inlet:\n")
    _write_chain_line(io, 4, "source_compact_index", profile.inlet.source_compact_index)
    _write_chain_line(io, 4, "total_number_density_m3",
                      profile.inlet.composition.total_number_density_m3)

    print(io, "    temperatures_K:\n")
    _write_chain_line(io, 6, "Tt", profile.inlet.thermal.Tt)
    _write_chain_line(io, 6, "Tv", profile.inlet.thermal.Tv)
    _write_chain_line(io, 6, "Tee", profile.inlet.thermal.Tee)
    _write_chain_line(io, 6, "Te", profile.inlet.thermal.Te)

    print(io, "    mole_fractions:\n")
    for (name, mole_fraction) in zip(profile.inlet.composition.species,
                                     profile.inlet.composition.mole_fractions)
        _write_chain_line(io, 6, name, mole_fraction)
    end

    print(io, "  profile:\n")
    _write_chain_line(io, 4, "z_start_m", first(profile.z_m))
    _write_chain_line(io, 4, "z_end_m", last(profile.z_m))
    _write_chain_line(io, 4, "retained_point_count", length(profile.z_m))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", entry::ChainSegmentEntry)
    print(io, "segment ", entry.segment_index, "/", length(entry.profile.z_m), ":\n")
    _write_chain_line(io, 2, "source_cell_index",
                      entry.compact_to_source_index[entry.segment_index])
    _write_chain_line(io, 2, "case_path", entry.segment_case_path)
    _write_chain_line(io, 2, "z_m", entry.profile.z_m[entry.segment_index])
    _write_chain_line(io, 2, "dx_m", entry.profile.dx_m[entry.segment_index])
    _write_chain_line(io, 2, "te_target_K", entry.profile.te_K[entry.segment_index])
    _write_chain_line(io, 2, "rho_ex_handoff_used", entry.state_cache_used)
    _write_chain_line(io, 2, "status", entry.result.success ? "success" : "failed")
    _write_chain_line(io, 2, "message", entry.result.message)

    print(io, "  species_u_m_s:\n")
    for (name, velocity) in _ordered_chain_velocity_pairs(entry.profile, entry.segment_index)
        _write_chain_line(io, 4, name, velocity)
    end

    if entry.endpoint_reactor !== nothing
        total_number_density_m3 = _number_density_to_si(entry.endpoint_reactor.composition.total_number_density,
                                                        entry.base_config.runtime.unit_system)
        print(io, "  endpoint:\n")
        _write_chain_line(io, 4, "total_number_density_m3", total_number_density_m3)

        print(io, "    temperatures_K:\n")
        _write_chain_line(io, 6, "Tt", entry.endpoint_reactor.thermal.Tt)
        _write_chain_line(io, 6, "Tv", entry.endpoint_reactor.thermal.Tv)
        _write_chain_line(io, 6, "Tee", entry.endpoint_reactor.thermal.Tee)
        _write_chain_line(io, 6, "Te", entry.endpoint_reactor.thermal.Te)

        print(io, "    number_density_m3:\n")
        for (name, number_density) in _species_number_density_pairs_m3(entry.endpoint_reactor.composition,
                                                                       entry.base_config.runtime.unit_system)
            _write_chain_line(io, 6, name, number_density)
        end
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", entry::ChainResultEntry)
    print(io, "result:\n")
    _write_chain_line(io, 2, "success", entry.result.success)
    _write_chain_line(io, 2, "failed_cell", entry.result.failed_cell)
    _write_chain_line(io, 2, "message", entry.result.message)
    return nothing
end

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
        header = _render(EventEntry(:detail, "0D integration snapshot"))
        _append_log_text(log_path(RUN_LOG, runtime), string(header, "\n", _render(entry)))
    end
    if _writes_console(mode)
        _emit_console_text(entry.text)
    end
    return nothing
end

function _chain_summary_entry(entry::ChainHeaderEntry)
    profile = entry.profile
    marching = entry.marching
    lines = [CHAIN_RUN_BANNER,
             @sprintf("segments: %d", length(profile.z_m)),
             "handoff_mode: $(marching.handoff_mode)",
             "termination_mode: $(marching.termination_mode)",
             "is_isothermal_teex: $(marching.is_isothermal_teex)"]
    return EventEntry(:info, join(lines, "\n"); console = :minimal)
end

function _chain_summary_entry(entry::ChainSegmentEntry)
    message = @sprintf("\n===== Segment %d/%d 0D Simulation =====",
                       entry.segment_index, length(entry.profile.z_m))
    return EventEntry(:info, message; console = :minimal)
end

function _chain_summary_entry(entry::ChainResultEntry)
    level = entry.result.success ? :info : :warn
    message = entry.result.success ? "\n$(entry.result.message)" : entry.result.message
    return EventEntry(level, string(message, "\n", CHAIN_RUN_FOOTER); console = :minimal)
end

function emit!(::RunLog, runtime::RuntimeConfig, entry::Union{ChainHeaderEntry,
                                                              ChainSegmentEntry,
                                                              ChainResultEntry})
    emit!(RUN_LOG, runtime, _chain_summary_entry(entry))
    return nothing
end

function emit!(::ChainLog, runtime::RuntimeConfig, entry::Union{ChainHeaderEntry,
                                                                ChainSegmentEntry,
                                                                ChainResultEntry})
    mode = stream_mode(CHAIN_LOG, runtime.logging)
    mode == :off && return nothing

    text = _render(entry)
    if _writes_file(mode)
        _append_log_text(log_path(CHAIN_LOG, runtime), text)
    end
    if _writes_console(mode)
        _emit_console_text(text)
    end
    return nothing
end

function _log_run_event(runtime::RuntimeConfig, level::Symbol, message::AbstractString;
                        console::Symbol = :never, fields...)
    emit!(RUN_LOG, runtime, EventEntry(level, message; console = console, fields...))
    return nothing
end

function _log_run_exception(runtime::RuntimeConfig, level::Symbol, message::AbstractString,
                            exception; console::Symbol = :never, fields...)
    emit!(RUN_LOG, runtime,
          ExceptionEntry(level, message, exception; console = console, fields...))
    return nothing
end

function _emit_integration_detail(runtime::RuntimeConfig, detail_text::AbstractString)
    emit!(RUN_LOG, runtime, IntegrationDetailEntry(detail_text))
    return nothing
end

function _prepare_chain_logging(runtime::RuntimeConfig)
    prepare!(CHAIN_LOG, runtime)
    return nothing
end

function _emit_chain_summary(runtime::RuntimeConfig, entry::Union{ChainHeaderEntry,
                                                                  ChainSegmentEntry,
                                                                  ChainResultEntry})
    emit!(RUN_LOG, runtime, entry)
    return nothing
end

function _emit_chain_detail(runtime::RuntimeConfig, entry::Union{ChainHeaderEntry,
                                                                 ChainSegmentEntry,
                                                                 ChainResultEntry})
    emit!(CHAIN_LOG, runtime, entry)
    return nothing
end

function _prepare_native_logging(runtime::RuntimeConfig)
    return prepare!(NATIVE_LOG, runtime)
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
    return _progress_mode_enabled(runtime.logging) ? FractionProgressReporter(tlim) : nothing
end

function _progress_message(progress::FractionProgressReporter, t::Real)
    fraction = progress.tlim > 0.0 ? clamp(Float64(t) / progress.tlim, 0.0, 1.0) : 1.0
    pct = round(Int, 100 * fraction)
    return @sprintf("-> %3d%% (t = %.3e / %.3e s)", pct, Float64(t), progress.tlim)
end

function _report_progress!(runtime::RuntimeConfig, progress::FractionProgressReporter, t::Real)
    _log_run_event(runtime, :info, _progress_message(progress, t);
                   console = :minimal,
                   :t => Float64(t),
                   :tlim => progress.tlim)
    return nothing
end
