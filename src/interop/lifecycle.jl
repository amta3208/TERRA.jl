"""
$(SIGNATURES)

Initialize the TERRA API system.

# Arguments
- `case_path::String`: Path to directory containing input/ subdirectory (default: current directory)

# Returns
- `NamedTuple`: Contains `num_species` and `num_dimensions` as determined by TERRA from input files

# Throws
- `ErrorException`: If case_path doesn't exist, input file is missing, or Fortran call fails
"""
function initialize_api_wrapper(; case_path::String = pwd())
    # Reconcile Julia/Fortran state first to avoid mismatches in tests
    try
        TERRA_INITIALIZED[] = is_api_initialized_wrapper()
    catch
        # ignore if library not loaded yet
        TERRA_INITIALIZED[] = false
    end

    # Always validate inputs and prepare filesystem, even if Fortran is already initialized.
    # This preserves input validation semantics and directory creation guarantees.
    # Validate inputs
    if !isdir(case_path)
        error("Case path does not exist: $case_path")
    end

    input_file = joinpath(case_path, "input", "prob_setup.inp")
    if !isfile(input_file)
        error("Required input file not found: $input_file")
    end

    TERRA_CASE_PATH[] = case_path
    TERRA_OUTPUTS_OPEN[] = false

    # Ensure output directory structure exists
    output_dir = joinpath(case_path, "output")
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Ensure required output subdirectories exist
    sources_dir = joinpath(output_dir, "sources")
    states_dir = joinpath(output_dir, "states")

    if !isdir(sources_dir)
        mkpath(sources_dir)
    end

    if !isdir(states_dir)
        mkpath(states_dir)
    end

    # If Fortran already initialized, skip reinitialization but still return dims
    if TERRA_INITIALIZED[]
        num_species = get_number_of_active_species_wrapper()
        num_dimensions = get_number_of_dimensions_wrapper()
        return (num_species = num_species, num_dimensions = num_dimensions)
    end

    # Store current directory and change to case path
    original_dir = pwd()

    try
        cd(case_path)

        # Create references for the Fortran call (as output parameters)
        # Fortran interface expects (num_species, num_dimensions) order
        num_species_ref = Ref{Int32}(0)
        num_dimensions_ref = Ref{Int32}(0)

        # Call Fortran function with correct parameter order
        ccall((:initialize_api, get_terra_lib_path()), Cvoid,
            (Ref{Int32}, Ref{Int32}),
            num_species_ref, num_dimensions_ref)

        # If we get here, the call succeeded
        TERRA_INITIALIZED[] = true

        # Return the values determined by Fortran
        return (num_species = num_species_ref[], num_dimensions = num_dimensions_ref[])

    catch e
        # Re-throw with more context about the failure
        if isa(e, Base.SystemError) || isa(e, ErrorException)
            error("TERRA initialization failed in directory $case_path: $(e.msg)")
        else
            rethrow(e)
        end
    finally
        # Always restore original directory
        cd(original_dir)
    end
end

"""
$(SIGNATURES)

Finalize the TERRA API system and clean up resources.
"""
function finalize_api_wrapper()
    # If library isn't loaded, there's nothing to do
    if !is_terra_loaded()
        return nothing
    end
    if TERRA_OUTPUTS_OPEN[]
        try
            close_api_output_files_wrapper()
        catch e
            @warn "Failed to close TERRA output files during finalize" exception=e
        end
    end
    # Query Fortran-side state and finalize quietly if needed
    TERRA_INITIALIZED[] = is_api_initialized_wrapper()
    if TERRA_INITIALIZED[]
        ccall((:finalize_api, get_terra_lib_path()), Cvoid, ())
        TERRA_INITIALIZED[] = false
    end
    TERRA_CASE_PATH[] = ""
    TERRA_OUTPUTS_OPEN[] = false
    return nothing
end

"""
$(SIGNATURES)

Control whether finalize_api() will call MPI_Finalize on the Fortran side.
Default is true; tests may disable it to allow reinitialization in one process.
"""
function set_api_finalize_mpi_wrapper(enable::Bool)
    if !is_terra_loaded()
        return nothing
    end
    flag = enable ? Int32(1) : Int32(0)
    ccall((:set_api_finalize_mpi, get_terra_lib_path()), Cvoid, (Int32,), flag)
    return nothing
end

"""
$(SIGNATURES)

Open native TERRA Tecplot-style output files for the current API session.

Requires a successful call to `initialize_api_wrapper` so the case path is tracked.
"""
function open_api_output_files_wrapper()
    if !is_terra_loaded()
        error("TERRA library not loaded. Call load_terra_library! before opening outputs.")
    end
    if !TERRA_INITIALIZED[]
        error("TERRA not initialized. Call initialize_api_wrapper() before opening outputs.")
    end
    if isempty(TERRA_CASE_PATH[])
        error("Case path not set. Initialize the API with a valid case_path before opening outputs.")
    end
    if TERRA_OUTPUTS_OPEN[]
        return nothing
    end

    case_path = TERRA_CASE_PATH[]
    try
        cd(case_path) do
            ccall((:open_api_output_files, get_terra_lib_path()), Cvoid, ())
        end
        TERRA_OUTPUTS_OPEN[] = true
    catch e
        error("Failed to open TERRA output files: $(e)")
    end

    return nothing
end

"""
$(SIGNATURES)

Write a snapshot of the current API state to the native TERRA output files, using
the `y` layout expected by `rhs_api`.
"""
function write_api_outputs_wrapper(istep::Integer, time::Real, dt::Real,
        state::AbstractVector{<:Real}; dist::Real = 0.0, dx::Real = 0.0)
    if !TERRA_OUTPUTS_OPEN[]
        error("TERRA output files are not open. Call open_api_output_files_wrapper() first.")
    end
    if isempty(TERRA_CASE_PATH[])
        error("Case path not set. Cannot write TERRA outputs without an initialized case path.")
    end

    state_vec = state isa Vector{Float64} ? state : Vector{Float64}(state)
    istep32 = Int32(istep)
    neq32 = Int32(length(state_vec))
    time64 = Float64(time)
    dt64 = Float64(dt)
    dist64 = Float64(dist)
    dx64 = Float64(dx)
    case_path = TERRA_CASE_PATH[]

    GC.@preserve state_vec begin
        ptr = Base.unsafe_convert(Ptr{Float64}, state_vec)
        cd(case_path) do
            ccall((:write_api_outputs, get_terra_lib_path()), Cvoid,
                (Int32, Float64, Float64, Float64, Float64, Int32, Ptr{Float64}),
                istep32, time64, dt64, dist64, dx64, neq32, ptr)
        end
    end

    return nothing
end

"""
$(SIGNATURES)

Close the native TERRA output files opened through the API wrappers.
"""
function close_api_output_files_wrapper()
    if !TERRA_OUTPUTS_OPEN[]
        return nothing
    end
    if !is_terra_loaded()
        TERRA_OUTPUTS_OPEN[] = false
        return nothing
    end
    if isempty(TERRA_CASE_PATH[])
        TERRA_OUTPUTS_OPEN[] = false
        return nothing
    end

    case_path = TERRA_CASE_PATH[]
    try
        cd(case_path) do
            ccall((:close_api_output_files, get_terra_lib_path()), Cvoid, ())
        end
    catch e
        error("Failed to close TERRA output files: $(e)")
    finally
        TERRA_OUTPUTS_OPEN[] = false
    end

    return nothing
end
