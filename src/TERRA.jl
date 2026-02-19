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
include("terra_config.jl")
include("terra_solver.jl")

export initialize_terra, finalize_terra
export TERRAResults
export Config, ReactorConfig, ReactorComposition, ReactorThermalState
export ModelConfig, TimeConfig, ODESolverConfig, SpaceConfig, NumericsConfig, RuntimeConfig
export ResidenceTimeConfig
export with_case_path, with_time, with_runtime
export solve_terra_0d, nitrogen_10ev_example

end
