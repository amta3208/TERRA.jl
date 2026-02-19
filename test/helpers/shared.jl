const TEST_CASE_PATH = normpath(joinpath(@__DIR__, "..", "test_case"))

"""
Reset the TERRA state (Fortran + Julia) and initialize from scratch.

- `case_path`: Path to a case directory containing `input/prob_setup.inp`
- `config` (keyword, optional; `Config`): If provided, a temporary case directory is created,
  input files are generated from this config there, and initialization is performed
  from that temporary directory (leaving `case_path` untouched).

Returns a NamedTuple with `(num_species, num_dimensions)`.
"""
function reset_and_init!(case_path::AbstractString;
        config::Union{Nothing, terra.Config} = nothing)
    # Best-effort cleanup (ignore errors if library is not loaded yet)
    try
        terra.finalize_api_wrapper()
    catch
        # ignore
    end
    try
        terra.close_terra_library()
    catch
        # ignore
    end

    # Fresh load + init
    terra.load_terra_library!()

    # For test stability within one process, do not have TERRA finalize MPI
    try
        terra.set_api_finalize_mpi_wrapper(false)
    catch
        # ignore if symbol missing (older library)
    end

    if config === nothing
        return terra.initialize_api_wrapper(case_path = case_path)
    else
        # Use a temporary case directory so we don't touch the shared test_case
        mktempdir() do tmp
            # Generate inputs in the temp directory based on the provided config
            terra.generate_input_files(config, tmp)
            # Initialize from the temp directory
            dims = terra.initialize_api_wrapper(case_path = tmp)
            return dims
        end
    end
end
