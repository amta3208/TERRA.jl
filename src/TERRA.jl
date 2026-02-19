module TERRA

using DocStringExtensions
using Libdl
using DifferentialEquations
using DiffEqBase: DiscreteCallback, set_proposed_dt!, u_modified!
using Printf

const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const TEST_DIR = joinpath(PACKAGE_ROOT, "test")

include("data_conversion.jl")
include("fortran_wrapper.jl")

include("config/types.jl")
include("config/validation.jl")
include("config/conversions.jl")

include("io/prob_setup_writer.jl")
include("io/sources_setup_writer.jl")
include("io/tau_scaling_writer.jl")
include("io/input_generation.jl")

include("solver/api_layout.jl")
include("solver/initial_state.jl")
include("solver/state_vector.jl")
include("solver/residence_time.jl")
include("solver/rhs.jl")
include("solver/integrate_0d.jl")
include("solver/driver.jl")

export initialize_terra, finalize_terra
export TERRAResults
export Config, ReactorConfig, ReactorComposition, ReactorThermalState
export ModelConfig, TimeConfig, ODESolverConfig, SpaceConfig, NumericsConfig, RuntimeConfig
export ResidenceTimeConfig
export with_case_path, with_time, with_runtime
export solve_terra_0d, nitrogen_10ev_example

end
