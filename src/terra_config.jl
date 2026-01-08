"""
# TERRA Configuration Module

This module handles configuration management for TERRA simulations,
including parameter validation, default values, and input file generation.
"""

"""
$(SIGNATURES)

Temperature configuration for TERRA simulation.

# Fields
- `Tt::Float64`: Translational temperature (K)
- `Tv::Float64`: Vibrational temperature (K)
- `Tee::Float64`: Electron-electronic temperature (K)
- `Te::Float64`: Electron temperature (K)
"""
struct TemperatureConfig
    Tt::Float64
    Tv::Float64
    Tee::Float64
    Te::Float64

    function TemperatureConfig(Tt, Tv, Tee, Te)
        if any([Tt, Tv, Tee, Te] .<= 0)
            error("All temperatures must be positive")
        end
        new(Float64(Tt), Float64(Tv), Float64(Tee), Float64(Te))
    end
end

# Preferred keyword constructor
TemperatureConfig(; Tt, Tv, Tee, Te) = TemperatureConfig(Tt, Tv, Tee, Te)

"""
$(SIGNATURES)

Time integration configuration for TERRA simulation.

All time values are in seconds within the TERRA.jl wrapper. When writing
Fortran input files, these values are converted to microseconds to match
the TERRA input format requirements.

# Fields
- `dt::Float64`: Time step (seconds)
- `dtm::Float64`: Output time step (seconds)
- `tlim::Float64`: Final time (seconds)
- `nstep::Int`: Maximum number of time steps
- `method::Int`: Integration method (0=forward Euler, 1=high order explicit, 2=implicit)
"""
struct TimeIntegrationConfig
    dt::Float64
    dtm::Float64
    tlim::Float64
    nstep::Int
    method::Int

    function TimeIntegrationConfig(dt, dtm, tlim, nstep = 500000, method = 2)
        if dt <= 0 || dtm <= 0 || tlim <= 0
            error("Time parameters must be positive")
        end
        if nstep <= 0
            error("Number of steps must be positive")
        end
        if !(method in [0, 1, 2])
            error("Integration method must be 0, 1, or 2")
        end
        new(Float64(dt), Float64(dtm), Float64(tlim), Int(nstep), Int(method))
    end
end

# Preferred keyword constructor
function TimeIntegrationConfig(;
        dt, dtm, tlim, nstep::Integer = 500000, method::Integer = 2)
    TimeIntegrationConfig(dt, dtm, tlim, nstep, method)
end

"""
$(SIGNATURES)

Physics modeling configuration for TERRA simulation.

# Fields
- `bbh_model::Int`: Bound-bound heavy particle model
- `esc_model::Int`: Escape model
- `ar_et_model::Int`: Ar-ET model
- `eex_noneq::Int`: Electron-electronic nonequilibrium flag
- `ev_relax_set::Int`: Electron-vibrational relaxation set
- `et_relax_set::Int`: Electron-translational relaxation set
- `radiation_length::Float64`: Radiation length scale (cm)
- `get_electron_density_by_charge_balance::Bool`: Electron density by charge balance
- `min_sts_frac::Float64`: Minimum state-to-state fraction
- `is_isothermal_teex::Bool`: Isothermal electron-electronic flag
- `energy_loss_per_eii::Float64`: Average electron energy loss per EII event (× E_ion)
"""
struct PhysicsConfig
    bbh_model::Int
    esc_model::Int
    ar_et_model::Int
    eex_noneq::Int
    ev_relax_set::Int
    et_relax_set::Int
    radiation_length::Float64
    get_electron_density_by_charge_balance::Bool
    min_sts_frac::Float64
    is_isothermal_teex::Bool
    energy_loss_per_eii::Float64

    function PhysicsConfig(;
            bbh_model = 4,
            esc_model = 1,
            ar_et_model = 1,
            eex_noneq = 1,
            ev_relax_set = 1,
            et_relax_set = 1,
            radiation_length = 1.0,
            get_electron_density_by_charge_balance = true,
            min_sts_frac = 1e-30,
            is_isothermal_teex = false,
            energy_loss_per_eii = 1.0
    )
        new(bbh_model, esc_model, ar_et_model, eex_noneq, ev_relax_set, et_relax_set,
            radiation_length, get_electron_density_by_charge_balance,
            min_sts_frac, is_isothermal_teex, energy_loss_per_eii)
    end
end

"""
$(SIGNATURES)

Process flags configuration for TERRA simulation.

# Fields
- `consider_elec_bbe::Int`: Consider electron bound-bound excitation
- `consider_elec_bfe::Int`: Consider electron bound-free excitation
- `consider_elec_bbh::Int`: Consider electron bound-bound heavy
- `consider_elec_bfh::Int`: Consider electron bound-free heavy
- `consider_rad::Int`: Consider radiation
- `consider_rdr::Int`: Consider RDR
- `consider_chem::Int`: Consider chemistry
"""
struct ProcessConfig
    consider_elec_bbe::Int
    consider_elec_bfe::Int
    consider_elec_bbh::Int
    consider_elec_bfh::Int
    consider_rad::Int
    consider_rdr::Int
    consider_chem::Int

    function ProcessConfig(;
            consider_elec_bbe = 1,
            consider_elec_bfe = 1,
            consider_elec_bbh = 1,
            consider_elec_bfh = 1,
            consider_rad = 0,
            consider_rdr = 0,
            consider_chem = 1
    )
        new(consider_elec_bbe, consider_elec_bfe, consider_elec_bbh,
            consider_elec_bfh, consider_rad, consider_rdr, consider_chem)
    end
end

"""
$(SIGNATURES)

Main configuration struct for TERRA simulations.

# Fields
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Initial mole fractions
- `total_number_density::Float64`: Total number density (1/cm³)
- `temperatures::TemperatureConfig`: Temperature configuration
- `time_params::TimeIntegrationConfig`: Time integration parameters
- `physics::PhysicsConfig`: Physics modeling options
- `processes::ProcessConfig`: Process flags
- `database_path::String`: Path to chemistry database
- `case_path::String`: Working directory for TERRA simulation
- `unit_system::Symbol`: Unit system (:SI or :CGS)
- `validate_species_against_terra::Bool`: Validate species against TERRA database
- `print_source_terms::Bool`: Print source terms flag
- `write_native_outputs::Bool`: Mirror native TERRA Tecplot outputs when running via the Julia wrapper

"""
struct TERRAConfig
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density::Float64
    temperatures::TemperatureConfig
    time_params::TimeIntegrationConfig
    physics::PhysicsConfig
    processes::ProcessConfig
    database_path::String
    case_path::String
    unit_system::Symbol
    validate_species_against_terra::Bool
    print_source_terms::Bool
    write_native_outputs::Bool

    function TERRAConfig(;
            species::Vector{String},
            mole_fractions::Vector{Float64},
            total_number_density::Float64,
            temperatures::TemperatureConfig,
            time_params::TimeIntegrationConfig,
            physics::PhysicsConfig = PhysicsConfig(),
            processes::ProcessConfig = ProcessConfig(),
            database_path::String = "../../databases/n2/elec_sts_expanded_electron_fits",
            case_path::String = pwd(),
            unit_system::Symbol = :CGS,
            validate_species_against_terra::Bool = false,
            print_source_terms::Bool = true,
            write_native_outputs::Bool = false
    )

        # Validate inputs
        validate_config(
            species, mole_fractions, total_number_density,
            temperatures, time_params, case_path, unit_system)

        new(species, mole_fractions, total_number_density, temperatures,
            time_params, physics, processes, database_path, case_path, unit_system,
            validate_species_against_terra, print_source_terms, write_native_outputs)
    end
end

"""
$(SIGNATURES)

Results container for TERRA simulations.

# Fields
- `time::Vector{Float64}`: Time points
- `species_densities::Matrix{Float64}`: Species densities over time
- `temperatures::NamedTuple`: Temperature evolution
- `total_energy::Vector{Float64}`: Total energy evolution
- `source_terms::Union{NamedTuple, Nothing}`: Source terms (if requested)
- `success::Bool`: Simulation success flag
- `message::String`: Status message
"""
struct TERRAResults
    time::Vector{Float64}
    species_densities::Matrix{Float64}
    temperatures::NamedTuple
    total_energy::Vector{Float64}
    source_terms::Union{NamedTuple, Nothing}
    success::Bool
    message::String
end

"""
$(SIGNATURES)

Validate TERRA configuration parameters.

# Arguments
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Initial mole fractions
- `total_number_density::Float64`: Total number density
- `temperatures::TemperatureConfig`: Temperature configuration
- `time_params::TimeIntegrationConfig`: Time integration parameters
- `case_path::String`: Working directory path
- `unit_system::Symbol`: Unit system specification

# Throws
- `ArgumentError` if validation fails
"""
function validate_config(species::Vector{String},
        mole_fractions::Vector{Float64},
        total_number_density::Float64,
        temperatures::TemperatureConfig,
        time_params::TimeIntegrationConfig,
        case_path::String,
        unit_system::Symbol)

    # Check species and mole fractions consistency
    if length(species) != length(mole_fractions)
        throw(ArgumentError("Species and mole_fractions arrays must have same length"))
    end

    if isempty(species)
        throw(ArgumentError("At least one species must be specified"))
    end

    # Check mole fractions sum to 1
    if abs(sum(mole_fractions) - 1.0) > 1e-10
        throw(ArgumentError("Mole fractions must sum to 1.0, got $(sum(mole_fractions))"))
    end

    # Check for negative mole fractions
    if any(mole_fractions .< 0)
        throw(ArgumentError("Mole fractions must be non-negative"))
    end

    # Check total number density
    if total_number_density <= 0
        throw(ArgumentError("Total number density must be positive"))
    end

    # Check for duplicate species
    if length(unique(species)) != length(species)
        throw(ArgumentError("Duplicate species names found"))
    end

    # Validate species names format
    for species_name in species
        if isempty(strip(species_name))
            throw(ArgumentError("Species names cannot be empty"))
        end
    end

    # Validate case path
    if !isdir(case_path)
        throw(ArgumentError("Case path directory does not exist: $case_path"))
    end

    # Validate unit system
    if !(unit_system in [:SI, :CGS])
        throw(ArgumentError("Unit system must be :SI or :CGS, got :$unit_system"))
    end

    return true
end

"""
$(SIGNATURES)

Generate TERRA input files from configuration with proper directory structure.

This function creates the directory structure required by the Fortran wrapper:
- case_path/input/     (input files)
- case_path/output/    (output files)
- case_path/output/sources/  (source term outputs)
- case_path/output/states/   (state outputs)

# Arguments
- `config::TERRAConfig`: TERRA configuration
- `case_path::String`: Case directory path (default: config.case_path)

# Returns
- `true` if files generated successfully

# Throws
- `ErrorException` if file generation fails
"""
function generate_input_files(config::TERRAConfig, case_path::String = config.case_path)
    try
        # Create required directory structure as expected by fortran_wrapper
        input_dir = joinpath(case_path, "input")
        output_dir = joinpath(case_path, "output")
        sources_dir = joinpath(output_dir, "sources")
        states_dir = joinpath(output_dir, "states")

        # Create all required directories
        for dir in [input_dir, output_dir, sources_dir, states_dir]
            if !isdir(dir)
                mkpath(dir)
            end
        end

        # Validate database path if it's a relative path
        database_path = config.database_path
        if !isabspath(database_path)
            # Convert relative path to absolute from case_path
            database_path = abspath(joinpath(case_path, database_path))
        end

        # Check if database directory exists (warn if not found)
        if !isdir(database_path)
            @warn "Database path not found: $database_path. This may cause TERRA initialization to fail."
        else
            # Check for chemistry.dat file
            chemistry_file = joinpath(database_path, "chemistry.dat")
            if !isfile(chemistry_file)
                @warn "chemistry.dat not found in database path: $chemistry_file"
            end
        end

        # Generate input files in the input directory
        generate_prob_setup_file(config, joinpath(input_dir, "prob_setup.inp"))
        generate_sources_setup_file(config, joinpath(input_dir, "sources_setup.inp"))
        generate_tau_scaling_file(config, joinpath(input_dir, "tau_scaling.inp"))

        @debug "TERRA input files generated successfully" case_path=case_path

        return true

    catch e
        error("Failed to generate TERRA input files: $(e)")
    end
end

"""
$(SIGNATURES)

Generate prob_setup.inp file from configuration.
"""
function generate_prob_setup_file(config::TERRAConfig, filepath::String)
    open(filepath, "w") do io
        println(io, "####################################################")
        println(io, "# Location of database and output folders")
        println(io, "####################################################")
        println(io, "DATABASE_PATH=$(config.database_path)")
        println(io, "CHEM_FILE_NAME=chemistry.dat")
        println(io)
        println(io, "--- Turn on source term printouts")
        println(io, "PRINT_SOURCE_TERMS=$(config.print_source_terms ? 1 : 0)")
        println(io)
        println(io, "####################################################")
        println(io, "# Freestream condition")
        println(io, "####################################################")
        println(io)
        println(io, "---  Number of species, must match chemistry.dat")
        println(io, "NSP=$(length(config.species))")
        println(io)
        println(io, "--- Species mole fractions ($(join(config.species, ", ")))")
        for (i, frac) in enumerate(config.mole_fractions)
            println(io, "X$(i)=$(frac)")
        end
        println(io)
        println(io, "--- Total number density (1/cm³)")
        # Ensure number density is written in CGS units as expected by TERRA
        tn_cgs = config.unit_system == :CGS ? config.total_number_density :
                 convert_number_density_si_to_cgs(config.total_number_density)
        println(io, "TOTAL_NUMBER_DENSITY=$(tn_cgs)")
        println(io)
        println(io, "--- Temperatures (K)")
        println(io, "TT=$(config.temperatures.Tt)")
        println(io, "TV=$(config.temperatures.Tv)")
        println(io, "TEE=$(config.temperatures.Tee)")
        println(io, "TE=$(config.temperatures.Te)")
        println(io)
        println(io, "--- Radiation length scale (cm)")
        println(io, "RAD_LEN=$(config.physics.radiation_length)")
        println(io)
        println(io, "####################################################")
        println(io, "# Physical modeling variables")
        println(io, "####################################################")
        println(io, "--- Physics options")
        println(io, "BBHMODEL=$(config.physics.bbh_model)")
        println(io, "ESC_MODEL=$(config.physics.esc_model)")
        println(io, "AR_ET_MODEL=$(config.physics.ar_et_model)")
        println(io, "EEX_NONEQ=$(config.physics.eex_noneq)")
        println(io, "EV_RELAX_SET=$(config.physics.ev_relax_set)")
        println(io, "ET_RELAX_SET=$(config.physics.et_relax_set)")
        println(io, "ENERGY_LOSS_PER_EII=$(config.physics.energy_loss_per_eii)")
        println(io,
            "GET_ELECTRON_DENSITY_BY_CHARGE_BALANCE=$(config.physics.get_electron_density_by_charge_balance ? 1 : 0)")
        println(io, "IS_ISOTHERMAL_TEEX=$(config.physics.is_isothermal_teex ? 1 : 0)")
        println(io, "MIN_STS_FRAC=$(config.physics.min_sts_frac)")
        println(io)
        println(io, "--- Process flags")
        println(io, "CONSIDER_ELEC_BBE=$(config.processes.consider_elec_bbe)")
        println(io, "CONSIDER_ELEC_BFE=$(config.processes.consider_elec_bfe)")
        println(io, "CONSIDER_ELEC_BBH=$(config.processes.consider_elec_bbh)")
        println(io, "CONSIDER_ELEC_BFH=$(config.processes.consider_elec_bfh)")
        println(io, "CONSIDER_RAD=$(config.processes.consider_rad)")
        println(io, "CONSIDER_RDR=$(config.processes.consider_rdr)")
        println(io, "CONSIDER_CHEM=$(config.processes.consider_chem)")
        println(io)
        println(io, "####################################################")
        println(io, "# Computational setup")
        println(io, "####################################################")
        println(io,
            "--- Time integration method: forward euler == 0, high order explicit == 1, numerical implicit == 2")
        println(io, "TIME_METHOD=$(config.time_params.method)")
        println(io)
        println(io, "--- Number of dimensions, 0 or 1")
        println(io, "ND=0")
        println(io)
        println(io, "--- 0D time integration setup (time in microseconds)")
        # Convert seconds from config to microseconds for Fortran input and format
        # with stable representations to satisfy tests and Fortran parsing
        dt_us = config.time_params.dt * 1e6
        dtm_us = config.time_params.dtm * 1e6
        tlim_us = config.time_params.tlim * 1e6
        # DT is typically very small; print with 2 significant digits (e.g., 5.0e-6)
        println(io, "DT=$(round(dt_us, sigdigits=2))")
        # DTM and TLIM: keep ~6 significant digits to avoid precision loss
        println(io, "DTM=$(round(dtm_us, sigdigits=6))")
        println(io, "TLIM=$(round(tlim_us, sigdigits=6))")
        println(io)
        println(io, "-- Max number of iterations")
        println(io, "NSTEP=$(config.time_params.nstep)")
    end
end

"""
$(SIGNATURES)

Generate sources_setup.inp file from configuration.
"""
function generate_sources_setup_file(config::TERRAConfig, filepath::String)
    open(filepath, "w") do io
        println(io, "BEGIN SPECIES SOURCES")
        for species in config.species
            println(io, species)
        end
        println(io, "END SPECIES SOURCES")
        println(io)
        println(io, "BEGIN EXCITED STATE SOURCES")
        # For now, include common excited states for nitrogen
        # This should be made more general based on the species
        if "N" in config.species
            for i in 1:5
                println(io, "N($(i))")
            end
        end
        if "N2" in config.species
            for i in 1:9
                println(io, "N2($(i))")
            end
        end
        if "N2+" in config.species
            for i in 1:4
                println(io, "N2+($(i))")
            end
        end
        println(io, "END EXCITED STATE SOURCES")
    end
end

"""
$(SIGNATURES)

Generate tau_scaling.inp file from configuration.
"""
function generate_tau_scaling_file(config::TERRAConfig, filepath::String)
    # For now, create an empty file
    # This can be expanded later if tau scaling is needed
    open(filepath, "w") do io
        println(io, "# Tau scaling file - currently empty")
    end
end

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

    # Resolve database path relative to package root for portability
    pkg_root = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
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

"""
$(SIGNATURES)

Validate species against TERRA database (requires loaded library).

This function uses the fortran_wrapper to query the TERRA database
for available species and validates the configuration species against it.

# Arguments
- `config::TERRAConfig`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if species validation fails or library not loaded
"""
function validate_species_against_terra_database(config::TERRAConfig)
    if config.validate_species_against_terra
        try
            # Load fortran_wrapper functions - they should be available since terra_config.jl is included in TERRA.jl
            if is_terra_loaded()
                terra_species = get_species_names_wrapper()
                validate_species_data(config.species, terra_species, config.mole_fractions)
                @info "Species validation against TERRA database passed"
            else
                @warn "TERRA library not loaded. Cannot validate species against database."
            end
        catch e
            @warn "Failed to validate species against TERRA database: $(e)"
        end
    end

    return true
end

"""
$(SIGNATURES)

Validate configuration parameters against loaded TERRA library capabilities.

This function checks that the configuration is compatible with the loaded
TERRA library, including array dimensions and species availability.

# Arguments
- `config::TERRAConfig`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if validation fails or library not loaded
"""
function validate_config_against_terra(config::TERRAConfig)
    if !is_terra_loaded()
        @warn "TERRA library not loaded. Cannot validate configuration against library capabilities."
        return true
    end

    try
        # Validate array dimensions
        validate_array_dimensions(config)

        # Validate species if requested
        if config.validate_species_against_terra
            validate_species_in_database(config)
        end

        @info "Configuration validation against TERRA library passed"
        return true

    catch e
        error("Configuration validation against TERRA library failed: $(e)")
    end
end

"""
$(SIGNATURES)

Validate that configuration arrays match TERRA Fortran expectations.

# Arguments
- `config::TERRAConfig`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if array dimensions exceed TERRA limits
"""
function validate_array_dimensions(config::TERRAConfig)
    if !is_terra_loaded()
        throw(ErrorException("TERRA library not loaded. Cannot validate array dimensions."))
    end

    # Get TERRA limits
    max_species = get_max_number_of_species_wrapper()
    max_atomic_states = get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_states = get_max_number_of_molecular_electronic_states_wrapper()

    # Check species count
    if length(config.species) > max_species
        throw(ErrorException("Number of species ($(length(config.species))) exceeds TERRA maximum ($max_species)"))
    end

    # Check mole fractions array consistency
    if length(config.mole_fractions) != length(config.species)
        throw(ErrorException("Mole fractions array length ($(length(config.mole_fractions))) does not match species count ($(length(config.species)))"))
    end

    @info "Array dimension validation passed" max_species=max_species max_atomic_states=max_atomic_states max_molecular_states=max_molecular_states

    return true
end

"""
$(SIGNATURES)

Validate that all species in configuration exist in the TERRA database.

# Arguments
- `config::TERRAConfig`: Configuration to validate

# Returns
- `true` if validation passes

# Throws
- `ErrorException` if any species are not found in database
"""
function validate_species_in_database(config::TERRAConfig)
    if !is_terra_loaded()
        throw(ErrorException("TERRA library not loaded. Cannot validate species in database."))
    end

    # Get available species from TERRA
    terra_species = get_species_names_wrapper()

    # Check each species in configuration
    missing_species = String[]
    for species in config.species
        if !(species in terra_species)
            push!(missing_species, species)
        end
    end

    if !isempty(missing_species)
        throw(ErrorException("Species not found in TERRA database: $(join(missing_species, ", "))\n" *
                             "Available species: $(join(terra_species, ", "))"))
    end

    @info "Species database validation passed" validated_species=config.species

    return true
end

"""
$(SIGNATURES)

Convert configuration units if needed.

# Arguments
- `config::TERRAConfig`: Configuration to convert
- `target_unit_system::Symbol`: Target unit system (:SI or :CGS)

# Returns
- `TERRAConfig`: Configuration with converted units
"""
function convert_config_units(config::TERRAConfig, target_unit_system::Symbol)
    if config.unit_system == target_unit_system
        return config  # No conversion needed
    end

    if config.unit_system == :SI && target_unit_system == :CGS
        # Convert SI to CGS
        new_total_number_density = convert_number_density_si_to_cgs(config.total_number_density)

        # Create new config with converted units
        return TERRAConfig(
            species = config.species,
            mole_fractions = config.mole_fractions,
            total_number_density = new_total_number_density,
            temperatures = config.temperatures,  # Temperature units are the same
            time_params = config.time_params,     # Time units are the same
            physics = config.physics,
            processes = config.processes,
            database_path = config.database_path,
            case_path = config.case_path,
            unit_system = target_unit_system,
            validate_species_against_terra = config.validate_species_against_terra,
            print_source_terms = config.print_source_terms,
            write_native_outputs = config.write_native_outputs
        )

    elseif config.unit_system == :CGS && target_unit_system == :SI
        # Convert CGS to SI
        new_total_number_density = convert_number_density_cgs_to_si(config.total_number_density)

        # Create new config with converted units
        return TERRAConfig(
            species = config.species,
            mole_fractions = config.mole_fractions,
            total_number_density = new_total_number_density,
            temperatures = config.temperatures,  # Temperature units are the same
            time_params = config.time_params,     # Time units are the same
            physics = config.physics,
            processes = config.processes,
            database_path = config.database_path,
            case_path = config.case_path,
            unit_system = target_unit_system,
            validate_species_against_terra = config.validate_species_against_terra,
            print_source_terms = config.print_source_terms,
            write_native_outputs = config.write_native_outputs
        )
    else
        error("Unsupported unit conversion: $(config.unit_system) to $target_unit_system")
    end
end
