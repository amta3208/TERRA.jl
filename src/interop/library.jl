"""
$(SIGNATURES)

Resolve the TERRA shared library path from the `$(TERRA_ENV_VAR_NAME)` environment variable.

# Returns
- `String`: Validated absolute path to the TERRA shared library

# Throws
- `ErrorException`: If the environment variable is unset or points to a missing file
"""
function resolve_terra_library_path()
    path_raw = String(strip(get(ENV, TERRA_ENV_VAR_NAME, "")))
    path = expanduser(path_raw)
    if isempty(path)
        error("TERRA library path not provided. Set $(TERRA_ENV_VAR_NAME) in the environment before running TERRA.")
    elseif !isfile(path)
        # If a relative path is provided, interpret it relative to the package root.
        # This makes `TERRA_LIB_PATH=terra/source/libterra.so` work under `Pkg.test`,
        # where the working directory is a temporary environment directory.
        candidate = if isabspath(path) || isempty(PACKAGE_ROOT)
            ""
        else
            abspath(joinpath(PACKAGE_ROOT, path))
        end
        if !isempty(candidate) && isfile(candidate)
            path = candidate
        else
            error("TERRA library file not found: $path\n" *
                  (isempty(candidate) ? "" : "Also tried: $candidate\n") *
                  "Set $(TERRA_ENV_VAR_NAME) to the full path of the TERRA shared library.")
        end
    end
    return path
end

"""
$(SIGNATURES)

Set the path to the TERRA shared library and load it.

# Arguments
- `path::String`: Path to the TERRA shared library file

# Throws
- `ErrorException`: If the library cannot be loaded
"""
function load_terra_library!(path::String)
    # Close existing handle if open
    if TERRA_HANDLE[] != C_NULL
        Libdl.dlclose(TERRA_HANDLE[])
        TERRA_HANDLE[] = C_NULL
    end

    # Validate path exists
    if !isfile(path)
        error("TERRA library file not found: $path\n" *
              "Set environment variable $(TERRA_ENV_VAR_NAME) to the full path of the TERRA shared library.")
    end

    # Open new library
    try
        TERRA_HANDLE[] = Libdl.dlopen(path)
        LOADED_TERRA_LIB_PATH[] = path  # Store the path for ccall usage
        # Reset initialization state for the freshly loaded library
        try
            TERRA_INITIALIZED[] = is_api_initialized_wrapper()
        catch
            TERRA_INITIALIZED[] = false
        end
    catch e
        error("Failed to load TERRA library from $path: $(e.msg)")
    end
end

"""
$(SIGNATURES)

Load the TERRA shared library using the `$(TERRA_ENV_VAR_NAME)` environment variable.

Throws an error if the environment variable is not set or does not point to an existing file.
"""
function load_terra_library!()
    path = resolve_terra_library_path()
    return load_terra_library!(path)
end

"""
$(SIGNATURES)

Enable or disable verbose wrapper debug logging.
"""
function set_wrapper_debug!(flag::Bool)
    DEBUG_WRAPPER[] = flag
end

"""
$(SIGNATURES)

Get the path to the loaded TERRA library.

# Returns
- `String`: Path to the loaded library

# Throws
- `ErrorException`: If no library is loaded
"""
function get_terra_lib_path()
    if TERRA_HANDLE[] == C_NULL
        error("TERRA library not loaded. Set $(TERRA_ENV_VAR_NAME) or call load_terra_library!(path) first.")
    end
    return LOADED_TERRA_LIB_PATH[]
end

"""
$(SIGNATURES)

Check if the TERRA library is currently loaded.

# Returns
- `Bool`: True if library is loaded, false otherwise
"""
function is_terra_loaded()
    return TERRA_HANDLE[] != C_NULL
end

"""
$(SIGNATURES)

Close the TERRA library and free resources.
"""
function close_terra_library()
    if TERRA_HANDLE[] != C_NULL
        Libdl.dlclose(TERRA_HANDLE[])
        TERRA_HANDLE[] = C_NULL
        LOADED_TERRA_LIB_PATH[] = ""  # Clear the path
        TERRA_INITIALIZED[] = false
    end
end
