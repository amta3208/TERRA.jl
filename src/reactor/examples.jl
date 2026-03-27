"""
$(SIGNATURES)

Create a default configuration for the 0D Nitrogen Te=10eV example.
"""
function nitrogen_10ev_config(; isothermal::Bool = false)
    species = ["N", "N2", "N+", "N2+", "E-"]
    mole_fractions = [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
    total_number_density = 1.0e13

    composition = ReactorComposition(; species = species,
                                     mole_fractions = mole_fractions,
                                     total_number_density = total_number_density)
    thermal = ReactorThermalState(; Tt = 750.0, Tv = 750.0, Tee = 750.0, Te = 115000.0)
    reactor = ReactorConfig(; composition = composition, thermal = thermal)

    physics = PhysicsConfig(; is_isothermal_teex = isothermal)
    models = ModelConfig(; physics = physics)

    time = TimeConfig(; dt = 5e-12,
                      dt_output = 5e-6,
                      duration = 1e-3,
                      nstep = 500000,
                      method = 2)

    pkg_root = normpath(joinpath(@__DIR__, "..", ".."))
    database_path = abspath(joinpath(pkg_root, "database", "n2",
                                     "elec_sts_expanded_electron_fits"))

    resolve_terra_library_path()

    if !isdir(database_path)
        error("TERRA database directory not found: $database_path\n" *
              "Please ensure the TERRA database exists and the path is correct.")
    end

    chemistry_file = joinpath(database_path, "chemistry.dat")
    if !isfile(chemistry_file)
        error("Required chemistry.dat file not found in database directory: $chemistry_file\n" *
              "Please ensure the database is complete.")
    end

    numerics = NumericsConfig(; time = time)
    sources = SourceTermsConfig()
    runtime = RuntimeConfig(; database_path = database_path)

    return Config(; reactor = reactor,
                  models = models,
                  sources = sources,
                  numerics = numerics,
                  runtime = runtime)
end

"""
$(SIGNATURES)

Run the 0D Nitrogen Te=10eV example case.
"""
function nitrogen_10ev_example(case_path::String = mktempdir();
                               isothermal::Bool = false)
    config = nitrogen_10ev_config(; isothermal = isothermal)
    config_with_path = with_case_path(config, case_path)

    @info "Running 0D Nitrogen Te=10eV example case"
    @info "Configuration" species=config_with_path.reactor.composition.species mole_fractions=config_with_path.reactor.composition.mole_fractions
    @info "Temperatures" Tt=config_with_path.reactor.thermal.Tt Te=config_with_path.reactor.thermal.Te
    @info "Time parameters" dt=config_with_path.numerics.time.dt tlim=config_with_path.numerics.time.duration
    @info "Case path" case_path=case_path

    try
        initialize_terra(config_with_path, case_path)
        results = solve_terra_0d(config_with_path)

        if results.success
            @info "Example simulation completed successfully"
            @info "Final conditions" time=results.t[end]
            final_frame = results.frames[end]

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
        try
            finalize_terra()
        catch e
            @warn "Error during cleanup" exception=e
        end
    end
end
