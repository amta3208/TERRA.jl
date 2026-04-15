"""
$(SIGNATURES)

Solve a 0D TERRA simulation.
"""
function _solve_terra_0d_internal(config::Config;
                                  sources::Union{Nothing, SourceTermsConfig} = config.sources,
                                  wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                                  state_cache::Union{Nothing, ReactorStateCache} = nothing)
    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        initial_state = config_to_initial_state(config; state_cache = state_cache)
        return _integrate_0d_system(config, initial_state;
                                    sources = sources,
                                    wall_inputs = wall_inputs)
    catch e
        emit!(RUN_LOG, config.runtime,
              ExceptionEntry(:error, "TERRA simulation failed", e;
                             console = :minimal))
        return _failed_reactor_result("Simulation failed: $(string(e))")
    end
end

function solve_terra_0d(config::Config;
                        sources::Union{Nothing, SourceTermsConfig} = config.sources)
    _validate_direct_wall_loss_usage(sources)
    return _solve_terra_0d_internal(config; sources = sources)
end
