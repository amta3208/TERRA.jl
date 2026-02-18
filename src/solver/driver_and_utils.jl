"""
$(SIGNATURES)

Solve a 0D TERRA simulation.

This is the main high-level interface for running TERRA simulations.
It handles all the complexity of data conversion, Fortran interfacing,
and result processing.

# Arguments
- `config::TERRAConfig`: Configuration for the simulation

# Returns
- `TERRAResults`: Results of the simulation

# Throws
- `ErrorException` if TERRA not initialized or simulation fails
"""
function solve_terra_0d(config::TERRAConfig;
        residence_time::Union{Nothing, ResidenceTimeConfig} = nothing,
        use_residence_time::Union{Nothing, Bool} = nothing)
    return solve_terra_0d(
        to_config(config);
        residence_time = residence_time,
        use_residence_time = use_residence_time)
end

"""
$(SIGNATURES)

Solve a 0D TERRA simulation using nested `Config`.
"""
function solve_terra_0d(config::Config;
        residence_time::Union{Nothing, ResidenceTimeConfig} = config.numerics.residence_time,
        use_residence_time::Union{Nothing, Bool} = nothing)
    if !is_terra_initialized()
        error("TERRA not initialized. Call initialize_terra(config) first.")
    end

    try
        @info "Starting TERRA 0D simulation" species=config.reactor.composition.species

        # Convert configuration to initial conditions (SI to CGS)
        initial_state = config_to_initial_state(config)

        # Run the time integration
        results = integrate_0d_system(config, initial_state;
            residence_time = residence_time, use_residence_time = use_residence_time)

        @info "TERRA simulation completed successfully"
        return results

    catch e
        @error "TERRA simulation failed" exception=e
        return TERRAResults(
            Float64[], zeros(0, 0), (;), Float64[], nothing, false,
            "Simulation failed: $(string(e))"
        )
    end
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
- `TERRAResults`: Results of the simulation

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
    config_with_path = TERRAConfig(
        species = config.species,
        mole_fractions = config.mole_fractions,
        total_number_density = config.total_number_density,
        temperatures = config.temperatures,
        time_params = config.time_params,
        physics = config.physics,
        processes = config.processes,
        database_path = config.database_path,
        case_path = case_path,
        unit_system = config.unit_system,
        validate_species_against_terra = config.validate_species_against_terra,
        print_source_terms = config.print_source_terms,
        write_native_outputs = config.write_native_outputs,
        print_integration_output = config.print_integration_output
    )

    @info "Running 0D Nitrogen Te=10eV example case"
    @info "Configuration" species=config_with_path.species mole_fractions=config_with_path.mole_fractions
    @info "Temperatures" Tt=config_with_path.temperatures.Tt Te=config_with_path.temperatures.Te
    @info "Time parameters" dt=config_with_path.time_params.dt tlim=config_with_path.time_params.tlim
    @info "Case path" case_path=case_path

    try
        # Initialize TERRA with config
        initialize_terra(config_with_path, case_path)

        # Run simulation
        results = solve_terra_0d(config_with_path)

        if results.success
            @info "Example simulation completed successfully"
            @info "Final conditions" time=results.time[end]

            # Print final species densities
            for (i, species) in enumerate(config_with_path.species)
                final_density = results.species_densities[i, end]
                unit_str = config_with_path.unit_system == :SI ? "kg/m³" : "g/cm³"
                @info "Final density" species=species density=final_density unit=unit_str
            end

            @info "Final temperatures" Tt=results.temperatures.tt[end] Tv=results.temperatures.tv[end] Te=results.temperatures.te[end]
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

"""
$(SIGNATURES)

Validate simulation results for physical consistency.

# Arguments
- `results::TERRAResults`: Simulation results to validate

# Returns
- `true` if results pass validation, `false` otherwise
"""
function validate_results(results::TERRAResults)
    if !results.success
        @warn "Simulation was not successful"
        return false
    end

    # Check for negative densities
    if any(results.species_densities .< 0)
        @warn "Negative species densities found"
        return false
    end

    # Check for NaN or Inf values
    if any(isnan.(results.species_densities)) || any(isinf.(results.species_densities))
        @warn "NaN or Inf values found in species densities"
        return false
    end

    if any(isnan.(results.temperatures.tt)) || any(isinf.(results.temperatures.tt))
        @warn "NaN or Inf values found in temperatures"
        return false
    end

    # Check temperature ranges (should be positive and reasonable)
    if any(results.temperatures.tt .<= 0) || any(results.temperatures.tt .> 1e6)
        @warn "Unreasonable translational temperatures found"
        return false
    end

    if any(results.temperatures.te .<= 0) || any(results.temperatures.te .> 1e6)
        @warn "Unreasonable electron temperatures found"
        return false
    end

    @info "Results validation passed"
    return true
end

"""
$(SIGNATURES)

Save TERRA results to file.

# Arguments
- `results::TERRAResults`: Results to save
- `filename::String`: Output filename (CSV format)

# Returns
- `true` if save successful
"""
function save_results(results::TERRAResults, filename::String)
    try
        # Prepare data for CSV output
        n_times = length(results.time)
        n_species = size(results.species_densities, 1)

        # Create header
        header = ["time", "total_energy", "T_trans", "T_electron", "T_vib"]
        for i in 1:n_species
            push!(header, "species_$(i)_density")
        end

        # Create data matrix
        data = zeros(n_times, length(header))
        data[:, 1] = results.time
        data[:, 2] = results.total_energy
        data[:, 3] = results.temperatures.tt
        data[:, 4] = results.temperatures.te
        data[:, 5] = results.temperatures.tv

        for i in 1:n_species
            data[:, 5 + i] = results.species_densities[i, :]
        end

        # Write to file
        open(filename, "w") do io
            # Write header
            println(io, join(header, ","))

            # Write data
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end

        @info "Results saved successfully" filename=filename
        return true

    catch e
        @error "Failed to save results" filename=filename exception=e
        return false
    end
end
