"""
$(SIGNATURES)

Generate TERRA input files from configuration with proper directory structure.

This function creates the directory structure required by the Fortran wrapper:
- case_path/input/     (input files)
- case_path/output/    (output files)
- case_path/output/sources/  (source term outputs)
- case_path/output/states/   (state outputs)

# Arguments
- `config::Config`: TERRA configuration
- `case_path::String`: Case directory path (default: `config.runtime.case_path`)

# Returns
- `true` if files generated successfully

# Throws
- `ErrorException` if file generation fails
"""
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
