"""
$(SIGNATURES)

Integrate a reactor simulation from config-defined initial conditions.
"""
function _integrate_reactor(config::Config;
                            sources::Union{Nothing, SourceTermsConfig} = config.sources,
                            wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                            state_cache::Union{Nothing, ReactorStateCache} = nothing)
    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        initial_state = build_initial_state(config; state_cache = state_cache)
        return _integrate_reactor(config, initial_state;
                                  sources = sources,
                                  wall_inputs = wall_inputs)
    catch e
        emit!(RUN_LOG, config.runtime,
              ExceptionEntry(:error, "TERRA simulation failed", e;
                             console = :minimal))
        return _failed_reactor_result("Simulation failed: $(string(e))")
    end
end

function integrate_reactor(config::Config;
                           sources::Union{Nothing, SourceTermsConfig} = config.sources)
    _validate_direct_wall_loss_usage(sources)
    return _integrate_reactor(config; sources = sources)
end
