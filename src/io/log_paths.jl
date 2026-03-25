const DEFAULT_OUTPUT_LOG_DIRNAME = "logs"
const DEFAULT_RUN_LOG_FILENAME = "run.log"
const DEFAULT_NATIVE_LOG_FILENAME = "native.log"
const DEFAULT_CHAIN_LOG_FILENAME = "chain.log"

function _resolve_log_dir(case_path::AbstractString, logging::LoggingConfig)::String
    if logging.log_dir === nothing
        return normpath(joinpath(case_path, "output", DEFAULT_OUTPUT_LOG_DIRNAME))
    end
    return isabspath(logging.log_dir) ? normpath(logging.log_dir) :
           normpath(joinpath(case_path, logging.log_dir))
end

_resolve_log_dir(runtime::RuntimeConfig) = _resolve_log_dir(runtime.case_path, runtime.logging)

function _ensure_log_dir(runtime::RuntimeConfig)::String
    log_dir = _resolve_log_dir(runtime)
    mkpath(log_dir)
    return log_dir
end

function _run_log_path(runtime::RuntimeConfig)::String
    return joinpath(_resolve_log_dir(runtime), DEFAULT_RUN_LOG_FILENAME)
end

function _native_log_path(runtime::RuntimeConfig)::String
    return joinpath(_resolve_log_dir(runtime), DEFAULT_NATIVE_LOG_FILENAME)
end

function _chain_log_path(runtime::RuntimeConfig)::String
    return joinpath(_resolve_log_dir(runtime), DEFAULT_CHAIN_LOG_FILENAME)
end

function _segment_logging_log_dir(base_runtime::RuntimeConfig,
                                  segment_case_path::AbstractString)::Union{Nothing, String}
    base_runtime.logging.log_dir === nothing && return nothing

    segment_suffix = relpath(segment_case_path, base_runtime.case_path)
    return normpath(joinpath(_resolve_log_dir(base_runtime), segment_suffix))
end

function _segment_logging_config(base_runtime::RuntimeConfig,
                                 segment_case_path::AbstractString)::LoggingConfig
    return _with_file_only_detail_logging(base_runtime.logging;
                                          log_dir = _segment_logging_log_dir(base_runtime,
                                                                             segment_case_path))
end
