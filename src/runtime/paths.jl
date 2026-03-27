const DEFAULT_OUTPUT_LOG_DIRNAME = "logs"
const DEFAULT_RUN_LOG_FILENAME = "run.log"
const DEFAULT_NATIVE_LOG_FILENAME = "native.log"
const DEFAULT_CHAIN_LOG_FILENAME = "chain.log"

"""
$(SIGNATURES)

Filesystem layout for a prepared TERRA case directory.
"""
struct CaseLayout
    case_path::String
    input_dir::String
    output_dir::String
    sources_dir::String
    states_dir::String
    log_dir::String
    run_log_path::String
    native_log_path::String
    chain_log_path::String
end

function _resolve_log_dir(case_path::AbstractString, logging::LoggingConfig)::String
    if logging.log_dir === nothing
        return normpath(joinpath(case_path, "output", DEFAULT_OUTPUT_LOG_DIRNAME))
    end
    return isabspath(logging.log_dir) ? normpath(logging.log_dir) :
           normpath(joinpath(case_path, logging.log_dir))
end

function CaseLayout(runtime::RuntimeConfig)
    case_path = normpath(runtime.case_path)
    input_dir = normpath(joinpath(case_path, "input"))
    output_dir = normpath(joinpath(case_path, "output"))
    sources_dir = normpath(joinpath(output_dir, "sources"))
    states_dir = normpath(joinpath(output_dir, "states"))
    log_dir_path = _resolve_log_dir(case_path, runtime.logging)
    return CaseLayout(case_path,
                      input_dir,
                      output_dir,
                      sources_dir,
                      states_dir,
                      log_dir_path,
                      joinpath(log_dir_path, DEFAULT_RUN_LOG_FILENAME),
                      joinpath(log_dir_path, DEFAULT_NATIVE_LOG_FILENAME),
                      joinpath(log_dir_path, DEFAULT_CHAIN_LOG_FILENAME))
end

case_layout(runtime::RuntimeConfig) = CaseLayout(runtime)

log_dir(runtime::RuntimeConfig) = case_layout(runtime).log_dir
run_log_path(runtime::RuntimeConfig) = case_layout(runtime).run_log_path
native_log_path(runtime::RuntimeConfig) = case_layout(runtime).native_log_path
chain_log_path(runtime::RuntimeConfig) = case_layout(runtime).chain_log_path

function ensure_case_layout!(runtime::RuntimeConfig)
    layout = case_layout(runtime)
    for dir in (layout.case_path,
                layout.input_dir,
                layout.output_dir,
                layout.sources_dir,
                layout.states_dir,
                layout.log_dir)
        isdir(dir) || mkpath(dir)
    end
    return layout
end
