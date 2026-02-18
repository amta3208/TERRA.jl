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
    return generate_input_files(to_config(config), case_path)
end

function generate_input_files(config::Config, case_path::String = config.runtime.case_path)
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
        database_path = config.runtime.database_path
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
    return generate_prob_setup_file(to_config(config), filepath)
end

function generate_prob_setup_file(config::Config, filepath::String)
    runtime = config.runtime
    composition = config.reactor.composition
    thermal = config.reactor.thermal
    physics = config.models.physics
    processes = config.models.processes
    time = config.numerics.time
    space = config.numerics.space

    open(filepath, "w") do io
        println(io, "####################################################")
        println(io, "# Location of database and output folders")
        println(io, "####################################################")
        println(io, "DATABASE_PATH=$(runtime.database_path)")
        println(io, "CHEM_FILE_NAME=chemistry.dat")
        println(io)
        println(io, "--- Turn on source term printouts")
        println(io, "PRINT_SOURCE_TERMS=$(runtime.print_source_terms ? 1 : 0)")
        println(io)
        println(io, "####################################################")
        println(io, "# Freestream condition")
        println(io, "####################################################")
        println(io)
        println(io, "---  Number of species, must match chemistry.dat")
        println(io, "NSP=$(length(composition.species))")
        println(io)
        println(io, "--- Species mole fractions ($(join(composition.species, ", ")))")
        for (i, frac) in enumerate(composition.mole_fractions)
            println(io, "X$(i)=$(frac)")
        end
        println(io)
        println(io, "--- Total number density (1/cm³)")
        # Ensure number density is written in CGS units as expected by TERRA
        tn_cgs = runtime.unit_system == :CGS ? composition.total_number_density :
                 convert_number_density_si_to_cgs(composition.total_number_density)
        println(io, "TOTAL_NUMBER_DENSITY=$(tn_cgs)")
        println(io)
        println(io, "--- Temperatures (K)")
        println(io, "TT=$(thermal.Tt)")
        println(io, "TV=$(thermal.Tv)")
        println(io, "TEE=$(thermal.Tee)")
        println(io, "TE=$(thermal.Te)")
        println(io)
        println(io, "--- Radiation length scale (cm)")
        println(io, "RAD_LEN=$(physics.radiation_length)")
        println(io)
        println(io, "####################################################")
        println(io, "# Physical modeling variables")
        println(io, "####################################################")
        println(io, "--- Physics options")
        println(io, "BBHMODEL=$(physics.bbh_model)")
        println(io, "ESC_MODEL=$(physics.esc_model)")
        println(io, "AR_ET_MODEL=$(physics.ar_et_model)")
        println(io, "EEX_NONEQ=$(physics.eex_noneq)")
        println(io, "EV_RELAX_SET=$(physics.ev_relax_set)")
        println(io, "ET_RELAX_SET=$(physics.et_relax_set)")
        println(io, "ENERGY_LOSS_PER_EII=$(physics.energy_loss_per_eii)")
        println(io,
            "GET_ELECTRON_DENSITY_BY_CHARGE_BALANCE=$(physics.get_electron_density_by_charge_balance ? 1 : 0)")
        println(io, "IS_ISOTHERMAL_TEEX=$(physics.is_isothermal_teex ? 1 : 0)")
        println(io, "MIN_STS_FRAC=$(physics.min_sts_frac)")
        println(io)
        println(io, "--- Process flags")
        println(io, "CONSIDER_ELEC_BBE=$(processes.consider_elec_bbe)")
        println(io, "CONSIDER_ELEC_BFE=$(processes.consider_elec_bfe)")
        println(io, "CONSIDER_ELEC_BBH=$(processes.consider_elec_bbh)")
        println(io, "CONSIDER_ELEC_BFH=$(processes.consider_elec_bfh)")
        println(io, "CONSIDER_RAD=$(processes.consider_rad)")
        println(io, "CONSIDER_RDR=$(processes.consider_rdr)")
        println(io, "CONSIDER_CHEM=$(processes.consider_chem)")
        println(io)
        println(io, "####################################################")
        println(io, "# Computational setup")
        println(io, "####################################################")
        println(io,
            "--- Time integration method: forward euler == 0, high order explicit == 1, numerical implicit == 2")
        println(io, "TIME_METHOD=$(time.method)")
        println(io)
        println(io, "--- Number of dimensions, 0 or 1")
        println(io, "ND=$(space.nd)")
        println(io)
        println(io, "--- 0D time integration setup (time in microseconds)")
        # Convert seconds from config to microseconds for Fortran input and format
        # with stable representations to satisfy tests and Fortran parsing
        dt_us = time.dt * 1e6
        dtm_us = time.dt_output * 1e6
        tlim_us = time.duration * 1e6
        # DT is typically very small; print with 2 significant digits (e.g., 5.0e-6)
        println(io, "DT=$(round(dt_us, sigdigits=2))")
        # DTM and TLIM: keep ~6 significant digits to avoid precision loss
        println(io, "DTM=$(round(dtm_us, sigdigits=6))")
        println(io, "TLIM=$(round(tlim_us, sigdigits=6))")
        println(io)
        println(io, "-- Max number of iterations")
        println(io, "NSTEP=$(time.nstep)")
    end
end

"""
$(SIGNATURES)

Generate sources_setup.inp file from configuration.
"""
function generate_sources_setup_file(config::TERRAConfig, filepath::String)
    return generate_sources_setup_file(to_config(config), filepath)
end

function generate_sources_setup_file(config::Config, filepath::String)
    species_names = config.reactor.composition.species

    open(filepath, "w") do io
        println(io, "BEGIN SPECIES SOURCES")
        for species in species_names
            println(io, species)
        end
        println(io, "END SPECIES SOURCES")
        println(io)
        println(io, "BEGIN EXCITED STATE SOURCES")
        # For now, include common excited states for nitrogen
        # This should be made more general based on the species
        if "N" in species_names
            for i in 1:5
                println(io, "N($(i))")
            end
        end
        if "N2" in species_names
            for i in 1:9
                println(io, "N2($(i))")
            end
        end
        if "N2+" in species_names
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
    return generate_tau_scaling_file(to_config(config), filepath)
end

function generate_tau_scaling_file(config::Config, filepath::String)
    # For now, create an empty file
    # This can be expanded later if tau scaling is needed
    open(filepath, "w") do io
        println(io, "# Tau scaling file - currently empty")
    end
end
