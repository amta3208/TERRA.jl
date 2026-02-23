# TERRA library state
const TERRA_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const LOADED_TERRA_LIB_PATH = Ref{String}("")
const TERRA_INITIALIZED = Ref{Bool}(false)
const DEBUG_WRAPPER = Ref{Bool}(false)
const TERRA_CASE_PATH = Ref{String}("")
const TERRA_OUTPUTS_OPEN = Ref{Bool}(false)

# Environment variable used to locate the shared library when not provided explicitly
const TERRA_ENV_VAR_NAME = "TERRA_LIB_PATH"
