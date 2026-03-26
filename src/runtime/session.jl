const ACTIVE_RUNTIME = Ref{Union{Nothing, RuntimeConfig}}(nothing)

_active_runtime_for_logging() = ACTIVE_RUNTIME[]

function with_active_runtime_for_logging(f::Function, runtime::RuntimeConfig;
                                         restore::Union{Nothing, RuntimeConfig} = runtime)
    ACTIVE_RUNTIME[] = runtime
    try
        return f()
    finally
        ACTIVE_RUNTIME[] = restore
    end
end

"""
$(SIGNATURES)

Initialize the TERRA system.

This function must be called before any TERRA calculations can be performed.
It sets up the Fortran library, initializes internal data structures, and
prepares the system for simulation. The TERRA shared library is obtained from
the `TERRA_LIB_PATH` environment variable when it is not already loaded.

# Arguments
- `config::Config`: Configuration for initialization
- `case_path::String`: Case directory path (optional, defaults to `config.runtime.case_path`)

# Returns
- `true` if initialization successful

# Throws
- `ErrorException` if initialization fails
"""
function initialize_terra(config::Config, case_path::String = config.runtime.case_path;
                          lifecycle_console::Symbol = :minimal,
                          preserve_active_runtime::Bool = true)
    previous_runtime = _active_runtime_for_logging()
    runtime = case_path == config.runtime.case_path ? config.runtime :
              with_runtime(config.runtime; case_path = case_path)
    restore_runtime = preserve_active_runtime ? runtime : previous_runtime

    return with_active_runtime_for_logging(runtime; restore = restore_runtime) do
        try
            ensure_case_layout!(runtime)
            _log_run_event(runtime, :info, "Initializing TERRA";
                           console = lifecycle_console,
                           :case_path => case_path)

            if !is_terra_loaded()
                load_terra_library!()
                _log_run_event(runtime, :debug,
                               "TERRA library loaded successfully via TERRA_LIB_PATH";
                               console = :never)
            end

            try
                if is_api_initialized_wrapper()
                    _log_run_event(runtime, :debug,
                                   "Finalizing existing TERRA API before re-initialization";
                                   console = :never)
                    finalize_api_wrapper()
                end
            catch e
                _log_run_exception(runtime, :warn,
                                   "Unable to query/finalize existing TERRA API state",
                                   e;
                                   console = :never)
            end

            generate_input_files(config, case_path)
            _log_run_event(runtime, :debug, "TERRA input files generated";
                           console = :never,
                           :case_path => case_path)

            native_logging = _prepare_native_logging(runtime)
            result = initialize_api_wrapper(case_path = case_path,
                                            native_console_level = native_logging.console_level,
                                            native_file_level = native_logging.file_level,
                                            native_log_path = native_logging.log_path)
            num_species = result.num_species
            num_dimensions = result.num_dimensions

            configured_species = config.reactor.composition.species
            if num_species != length(configured_species)
                error("Configured species count ($(length(configured_species))) does not match TERRA setup ($num_species). " *
                      "Ensure configuration and generated input match the TERRA database.")
            end

            _log_run_event(runtime, :info, "TERRA initialized successfully";
                           console = lifecycle_console,
                           :num_species => num_species,
                           :num_dimensions => num_dimensions)

            try
                flags = get_runtime_flags()
                _log_run_event(runtime, :debug, "TERRA runtime flags";
                               console = :never,
                               :ev_relax_set => flags.ev_relax_set,
                               :vib_noneq => flags.vib_noneq,
                               :eex_noneq => flags.eex_noneq,
                               :rot_noneq => flags.rot_noneq,
                               :bfe => flags.consider_elec_bfe,
                               :bbh => flags.consider_elec_bbh,
                               :bfh => flags.consider_elec_bfh,
                               :bbe => flags.consider_elec_bbe)
            catch e
                _log_run_exception(runtime, :warn,
                                   "Unable to read TERRA runtime flags (rebuild library to enable)",
                                   e;
                                   console = :never)
            end
            return true

        catch e
            _log_run_exception(runtime, :error, "Failed to initialize TERRA", e;
                               console = :minimal)
            rethrow(e)
        end
    end
end

"""
$(SIGNATURES)

Finalize the TERRA system and clean up resources.

This function should be called when TERRA is no longer needed to properly
clean up memory and resources.

# Returns
- `Nothing`
"""
function finalize_terra()
    runtime = _active_runtime_for_logging()
    if runtime === nothing
        finalize_api_wrapper()
        close_terra_library()
        return nothing
    end

    return with_active_runtime_for_logging(runtime; restore = nothing) do
        try
            _log_run_event(runtime, :info, "Finalizing TERRA"; console = :minimal)
            finalize_api_wrapper()
            close_terra_library()
            _log_run_event(runtime, :info, "TERRA finalized successfully";
                           console = :minimal)
        catch e
            _log_run_exception(runtime, :error, "Error during TERRA finalization", e;
                               console = :minimal)
            rethrow(e)
        end
    end
end

"""
$(SIGNATURES)

Check if TERRA is properly initialized.

# Returns
- `true` if TERRA is initialized and ready for use
"""
function is_terra_initialized()
    try
        return is_terra_loaded() && is_api_initialized_wrapper()
    catch
        return false
    end
end
