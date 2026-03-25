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

"""
$(SIGNATURES)

Solve a 0D TERRA simulation.

This is the main high-level interface for running TERRA simulations.
It handles all the complexity of data conversion, Fortran interfacing,
and result processing.

# Arguments
- `config::Config`: Configuration for the simulation

# Returns
- `ReactorResult`: Results of the simulation

# Throws
- `ErrorException` if TERRA not initialized or simulation fails
"""
function _solve_terra_0d_internal(config::Config;
                                  sources::Union{Nothing, SourceTermsConfig} = config.sources,
                                  wall_inputs::Union{Nothing, SegmentWallInputs} = nothing,
                                  state_cache::Union{Nothing, ReactorStateCache} = nothing)
    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        @info "Starting TERRA 0D simulation" species=config.reactor.composition.species

        # Convert configuration to initial conditions (SI to CGS)
        initial_state = config_to_initial_state(config; state_cache = state_cache)

        # Run the time integration
        results, final_state_cache = _integrate_0d_system(config, initial_state;
                                                          sources = sources,
                                                          wall_inputs = wall_inputs,
                                                          inlet_state_cache = state_cache)

        @info "TERRA simulation completed successfully"
        return results, final_state_cache

    catch e
        @error "TERRA simulation failed" exception=e
        return ReactorResult(;
                             t = Float64[],
                             frames = ReactorFrame[],
                             source_terms = nothing,
                             success = false,
                             message = "Simulation failed: $(string(e))"), nothing
    end
end

function solve_terra_0d(config::Config;
                        sources::Union{Nothing, SourceTermsConfig} = config.sources)
    _validate_direct_wall_loss_usage(sources)
    results, _ = _solve_terra_0d_internal(config;
                                          sources = sources)
    return results
end

"""
$(SIGNATURES)

Create a default configuration for the 0D Nitrogen Te=10eV example.

# Returns
- `Config`: Configuration matching the example case

# Throws
- `ErrorException` if `$(TERRA_ENV_VAR_NAME)` is unset/invalid or if required database paths do not exist
"""
function nitrogen_10ev_config(; isothermal::Bool = false)
    species = ["N", "N2", "N+", "N2+", "E-"]
    mole_fractions = [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
    total_number_density = 1.0e13  # 1/cm³

    composition = ReactorComposition(;
                                     species = species,
                                     mole_fractions = mole_fractions,
                                     total_number_density = total_number_density)
    thermal = ReactorThermalState(; Tt = 750.0, Tv = 750.0, Tee = 750.0, Te = 115000.0)
    reactor = ReactorConfig(; composition = composition, thermal = thermal)

    physics = PhysicsConfig(; is_isothermal_teex = isothermal)
    models = ModelConfig(; physics = physics)

    # Time parameters are specified in seconds within the wrapper.
    # The TERRA input file expects microseconds; conversion is handled
    # in generate_input_files(). These values correspond to:
    #   dt   = 0.5e-5 microseconds  -> 5e-12 seconds
    #   dtm  = 5.0   microseconds   -> 5e-6  seconds
    #   tlim = 1.0e3 microseconds   -> 1e-3  seconds
    time = TimeConfig(;
                      dt = 5e-12, dt_output = 5e-6, duration = 1e-3, nstep = 500000,
                      method = 2)

    # Resolve database path relative to package root for portability.
    # This file lives in src/solver, so package root is two levels up.
    pkg_root = normpath(joinpath(@__DIR__, "..", ".."))
    database_path = abspath(joinpath(pkg_root, "database", "n2",
                                     "elec_sts_expanded_electron_fits"))

    # Validate that required paths exist
    resolve_terra_library_path()

    if !isdir(database_path)
        error("TERRA database directory not found: $database_path\n" *
              "Please ensure the TERRA database exists and the path is correct.")
    end

    # Additional check for chemistry.dat file in database
    chemistry_file = joinpath(database_path, "chemistry.dat")
    if !isfile(chemistry_file)
        error("Required chemistry.dat file not found in database directory: $chemistry_file\n" *
              "Please ensure the database is complete.")
    end

    numerics = NumericsConfig(; time = time)
    sources = SourceTermsConfig()
    runtime = RuntimeConfig(; database_path = database_path)

    return Config(;
                  reactor = reactor,
                  models = models,
                  sources = sources,
                  numerics = numerics,
                  runtime = runtime)
end

"""
$(SIGNATURES)

Run the 0D Nitrogen Te=10eV example case.

This function provides a convenient way to run the reference test case
that matches the TERRA example in `/terra/examples/0D_Nitrogen_Te_10eV`.
Requires the `TERRA_LIB_PATH` environment variable to point to the TERRA shared library.

# Arguments
- `case_path::String`: Case directory path (optional, creates temp directory if not provided)

# Returns
- `ReactorResult`: Results of the simulation

# Example
```julia
results = nitrogen_10ev_example()
```
"""
function nitrogen_10ev_example(case_path::String = mktempdir();
                               isothermal::Bool = false)
    # Create configuration for the example case
    config = nitrogen_10ev_config(; isothermal = isothermal)

    # Update config with case path
    config_with_path = with_case_path(config, case_path)

    @info "Running 0D Nitrogen Te=10eV example case"
    @info "Configuration" species=config_with_path.reactor.composition.species mole_fractions=config_with_path.reactor.composition.mole_fractions
    @info "Temperatures" Tt=config_with_path.reactor.thermal.Tt Te=config_with_path.reactor.thermal.Te
    @info "Time parameters" dt=config_with_path.numerics.time.dt tlim=config_with_path.numerics.time.duration
    @info "Case path" case_path=case_path

    try
        # Initialize TERRA with config
        initialize_terra(config_with_path, case_path)

        # Run simulation
        results = solve_terra_0d(config_with_path)

        if results.success
            @info "Example simulation completed successfully"
            @info "Final conditions" time=results.t[end]
            final_frame = results.frames[end]

            # Print final species densities
            for (i, species) in enumerate(config_with_path.reactor.composition.species)
                final_density = final_frame.species_densities[i]
                unit_str = config_with_path.runtime.unit_system == :SI ? "kg/m³" : "g/cm³"
                @info "Final density" species=species density=final_density unit=unit_str
            end

            @info "Final temperatures" Tt=final_frame.temperatures.tt Tv=final_frame.temperatures.tv Te=final_frame.temperatures.te
        else
            @error "Example simulation failed" message=results.message
        end

        return results

    finally
        # Clean up TERRA resources
        try
            finalize_terra()
        catch e
            @warn "Error during cleanup" exception=e
        end
    end
end
