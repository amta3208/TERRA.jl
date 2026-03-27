"""
$(SIGNATURES)

Solve a 0D TERRA simulation.
"""
function _solve_terra_0d_internal(config::Config;
                                  sources::Union{Nothing, SourceTermsConfig} = config.sources,
                                  wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                                  state_cache::Union{Nothing, ReactorStateCache} = nothing,
                                  presentation = STANDALONE_0D_PRESENTATION)
    presentation_obj = _integration_presentation(presentation)

    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        initial_state = config_to_initial_state(config; state_cache = state_cache)
        results, final_state_cache = _integrate_0d_system(config, initial_state;
                                                          sources = sources,
                                                          wall_inputs = wall_inputs,
                                                          inlet_state_cache = state_cache,
                                                          presentation = presentation_obj)
        return results, final_state_cache
    catch e
        emit!(RUN_LOG, config.runtime,
              ExceptionEntry(:error, "TERRA simulation failed", e;
                             console = :minimal))
        return _failed_reactor_result("Simulation failed: $(string(e))"), nothing
    end
end

function solve_terra_0d(config::Config;
                        sources::Union{Nothing, SourceTermsConfig} = config.sources)
    _validate_direct_wall_loss_usage(sources)
    results, _ = _solve_terra_0d_internal(config; sources = sources)
    return results
end
