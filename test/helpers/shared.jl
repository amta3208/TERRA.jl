const TEST_PACKAGE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TEST_CASES_ROOT = normpath(joinpath(@__DIR__, "..", "cases"))
const TEST_TERRA_FORTRAN_REFERENCE_CASE_PATH = normpath(joinpath(TEST_CASES_ROOT,
                                                                 "terra_fortran",
                                                                 "reference_case"))
const TEST_CASE_PATH = TEST_TERRA_FORTRAN_REFERENCE_CASE_PATH
const TEST_HET_CHAIN_INTERFACE_CASE_PATH = normpath(joinpath(TEST_CASES_ROOT,
                                                             "hallthruster_jl",
                                                             "chain_interface_case"))
const TEST_TERRA_CHAIN_INTERFACE_CASE_PATH = normpath(joinpath(TEST_CASES_ROOT, "terra_jl",
                                                               "chain_interface_case"))
const _HALLTHRUSTER_EXPORT_TOOL = let tool = Module(gensym(:HallThrusterExportTool))
    Base.include(tool,
                 joinpath(TEST_PACKAGE_ROOT, "tools", "hallthruster_jl",
                          "export_chain_profile.jl"))
    tool
end

function hallthruster_export_tool()
    return _HALLTHRUSTER_EXPORT_TOOL
end

"""
This avoids calling `dlclose` on an initialized Fortran/MPI runtime by
disabling MPI finalization for tests, finalizing the API if needed, and
only then unloading the library handle.
"""
function unload_terra_library_for_tests!()
    try
        terra.set_api_finalize_mpi_wrapper(false)
    catch
        # ignore if the library is not loaded yet or the symbol is unavailable
    end
    try
        terra.finalize_api_wrapper()
    catch
        # ignore cleanup errors in tests
    end
    try
        terra.close_terra_library()
    catch
        # ignore cleanup errors in tests
    end
    return nothing
end

"""
Reset the native TERRA library to a loaded-but-not-initialized state for tests.
"""
function reload_terra_library_for_tests!()
    unload_terra_library_for_tests!()
    terra.load_terra_library!()
    return nothing
end

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
    unload_terra_library_for_tests!()

    # Fresh load + init
    terra.load_terra_library!()

    # For test stability within one process, do not have TERRA finalize MPI
    try
        terra.set_api_finalize_mpi_wrapper(false)
    catch
        # ignore if symbol missing (older library)
    end

    if config === nothing
        if !isdir(case_path)
            throw(ArgumentError("Test case directory does not exist: $case_path"))
        end
        prob_setup_path = joinpath(case_path, "input", "prob_setup.inp")
        if !isfile(prob_setup_path)
            throw(ArgumentError("Test case is missing required input file: $prob_setup_path"))
        end
        return terra.initialize_api_wrapper(case_path = case_path)
    else
        # Use a temporary case directory so we don't touch the shared reference fixture
        mktempdir() do tmp
            # Generate inputs in the temp directory based on the provided config
            terra.generate_input_files(config, tmp)
            # Initialize from the temp directory
            dims = terra.initialize_api_wrapper(case_path = tmp)
            return dims
        end
    end
end
