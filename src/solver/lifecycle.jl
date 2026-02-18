"""
$(SIGNATURES)

Initialize the TERRA system.

This function must be called before any TERRA calculations can be performed.
It sets up the Fortran library, initializes internal data structures, and
prepares the system for simulation. The TERRA shared library is obtained from
the `TERRA_LIB_PATH` environment variable when it is not already loaded.

# Arguments
- `config::TERRAConfig`: Configuration for initialization
- `case_path::String`: Case directory path (optional, defaults to config.case_path)

# Returns
- `true` if initialization successful

# Throws
- `ErrorException` if initialization fails
"""
function initialize_terra(config::TERRAConfig, case_path::String = config.case_path)
    return initialize_terra(to_config(config), case_path)
end

"""
$(SIGNATURES)

Initialize TERRA using nested `Config`.
"""
function initialize_terra(config::Config, case_path::String = config.runtime.case_path)
    try
        # Ensure the shared library is loaded
        if !is_terra_loaded()
            load_terra_library!()
            @debug "TERRA library loaded successfully via TERRA_LIB_PATH"
        end

        # If the Fortran API is already initialized, finalize it so we can
        # reinitialize using the updated configuration and input files.
        try
            if is_api_initialized_wrapper()
                @debug "Finalizing existing TERRA API before re-initialization"
                finalize_api_wrapper()
            end
        catch e
            @warn "Unable to query/finalize existing TERRA API state" exception=e
        end

        # Generate input files for TERRA based on the provided configuration
        generate_input_files(config, case_path)
        @debug "TERRA input files generated" case_path=case_path

        # Initialize the API - get dimensions from Fortran
        result = initialize_api_wrapper(case_path = case_path)
        num_species = result.num_species
        num_dimensions = result.num_dimensions

        # Hard consistency check: configured species count must match TERRA setup
        configured_species = config.reactor.composition.species
        if num_species != length(configured_species)
            error("Configured species count ($(length(configured_species))) does not match TERRA setup ($num_species). " *
                  "Ensure configuration and generated input match the TERRA database.")
        end

        @debug "TERRA initialized successfully" num_species=num_species num_dimensions=num_dimensions
        # Fetch and log runtime flags from TERRA for verification
        try
            flags = get_runtime_flags()
            @debug "TERRA runtime flags" ev_relax_set=flags.ev_relax_set vib_noneq=flags.vib_noneq eex_noneq=flags.eex_noneq rot_noneq=flags.rot_noneq bfe=flags.consider_elec_bfe bbh=flags.consider_elec_bbh bfh=flags.consider_elec_bfh bbe=flags.consider_elec_bbe
        catch e
            @warn "Unable to read TERRA runtime flags (rebuild library to enable)" exception=e
        end
        return true

    catch e
        @error "Failed to initialize TERRA" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Finalize the TERRA system and clean up resources.

This function should be called when TERRA is no longer needed to properly
clean up memory and resources.

# Arguments
- None

# Returns
- `Nothing`
"""
function finalize_terra()
    try
        # Call the wrapper finalization
        finalize_api_wrapper()

        # Close the library
        close_terra_library()

        @info "TERRA finalized successfully"
    catch e
        @error "Error during TERRA finalization" exception=e
        rethrow(e)
    end
end

"""
$(SIGNATURES)

Check if TERRA is properly initialized.

# Arguments
- None

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
