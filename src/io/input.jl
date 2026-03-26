"""
$(SIGNATURES)

Generate TERRA input files from configuration with the prepared case layout.

This function creates the directory structure required by the Fortran wrapper:
- `case_path/input/`
- `case_path/output/`
- `case_path/output/sources/`
- `case_path/output/states/`

# Arguments
- `config::Config`: TERRA configuration
- `case_path::String`: Case directory path (default: `config.runtime.case_path`)

# Returns
- `true` if files generated successfully
"""
function generate_input_files(config::Config, case_path::String = config.runtime.case_path)
    runtime = case_path == config.runtime.case_path ? config.runtime :
              with_runtime(config.runtime; case_path = case_path)
    layout = ensure_case_layout!(runtime)

    try
        database_path = config.runtime.database_path
        if !isabspath(database_path)
            database_path = abspath(joinpath(case_path, database_path))
        end

        if !isdir(database_path)
            @warn "Database path not found: $database_path. This may cause TERRA initialization to fail."
        else
            chemistry_file = joinpath(database_path, "chemistry.dat")
            if !isfile(chemistry_file)
                @warn "chemistry.dat not found in database path: $chemistry_file"
            end
        end

        write_prob_setup!(config, joinpath(layout.input_dir, "prob_setup.inp"))
        write_sources_setup!(config, joinpath(layout.input_dir, "sources_setup.inp"))
        write_tau_scaling!(config, joinpath(layout.input_dir, "tau_scaling.inp"))

        @debug "TERRA input files generated successfully" case_path=case_path
        return true
    catch e
        error("Failed to generate TERRA input files: $(e)")
    end
end

function write_prob_setup!(config::Config, filepath::String)
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
        dt_us = time.dt * 1e6
        dtm_us = time.dt_output * 1e6
        tlim_us = time.duration * 1e6
        println(io, "DT=$(round(dt_us, sigdigits = 2))")
        println(io, "DTM=$(round(dtm_us, sigdigits = 6))")
        println(io, "TLIM=$(round(tlim_us, sigdigits = 6))")
        println(io)
        println(io, "-- Max number of iterations")
        println(io, "NSTEP=$(time.nstep)")
    end
    return nothing
end

function write_sources_setup!(config::Config, filepath::String)
    species_names = config.reactor.composition.species

    open(filepath, "w") do io
        println(io, "BEGIN SPECIES SOURCES")
        for species in species_names
            println(io, species)
        end
        println(io, "END SPECIES SOURCES")
        println(io)
        println(io, "BEGIN EXCITED STATE SOURCES")
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
    return nothing
end

function write_tau_scaling!(config::Config, filepath::String)
    open(filepath, "w") do io
        println(io, "# Tau scaling file - currently empty")
    end
    return nothing
end
