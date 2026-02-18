"""
$(SIGNATURES)

Create a default configuration for the 0D Nitrogen Te=10eV example.

# Returns
- `TERRAConfig`: Configuration matching the example case

# Throws
- `ErrorException` if `$(TERRA_ENV_VAR_NAME)` is unset/invalid or if required database paths do not exist
"""
function nitrogen_10ev_config(; isothermal::Bool = false)
    species = ["N", "N2", "N+", "N2+", "E-"]
    mole_fractions = [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
    total_number_density = 1.0e13  # 1/cm³

    temperatures = TemperatureConfig(; Tt = 750.0, Tv = 750.0, Tee = 750.0, Te = 115000.0)
    physics = PhysicsConfig(; is_isothermal_teex = isothermal)

    # Time parameters are specified in seconds within the wrapper.
    # The TERRA input file expects microseconds; conversion is handled
    # in generate_input_files(). These values correspond to:
    #   dt   = 0.5e-5 microseconds  -> 5e-12 seconds
    #   dtm  = 5.0   microseconds   -> 5e-6  seconds
    #   tlim = 1.0e3 microseconds   -> 1e-3  seconds
    time_params = TimeIntegrationConfig(;
        dt = 5e-12, dtm = 5e-6, tlim = 1e-3, nstep = 500000, method = 2)

    # Resolve database path relative to package root for portability.
    # This file lives in src/config, so package root is two levels up.
    pkg_root = normpath(joinpath(@__DIR__, "..", ".."))
    database_path = abspath(joinpath(
        pkg_root, "database", "n2", "elec_sts_expanded_electron_fits"))

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

    return TERRAConfig(
        species = species,
        mole_fractions = mole_fractions,
        total_number_density = total_number_density,
        temperatures = temperatures,
        time_params = time_params,
        physics = physics,
        database_path = database_path
    )
end

"""
$(SIGNATURES)

Get molecular weights for common species (g/mol).

# Arguments
- `species::Vector{String}`: Species names

# Returns
- `Vector{Float64}`: Molecular weights in g/mol
"""
function get_molecular_weights(species::Vector{String})
    # Common molecular weights (g/mol)
    molecular_weight_db = Dict(
        "N" => 14.007,
        "N2" => 28.014,
        "N+" => 14.007,
        "N2+" => 28.014,
        "E-" => 5.485799e-4,  # electron mass
        "Ar" => 39.948,
        "Ar+" => 39.948,
        "Xe" => 131.293,
        "Xe+" => 131.293,
        "Kr" => 83.798,
        "Kr+" => 83.798,
        "O" => 15.999,
        "O2" => 31.998,
        "O+" => 15.999,
        "O2+" => 31.998
    )

    weights = Float64[]
    for sp in species
        if haskey(molecular_weight_db, sp)
            push!(weights, molecular_weight_db[sp])
        else
            error("Molecular weight not found for species: $(sp)")
        end
    end

    return weights
end
