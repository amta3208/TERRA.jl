module TERRA

using Dates
using DocStringExtensions
using Libdl
using JSON
using OrdinaryDiffEq
using SciMLBase: ODEProblem, solve, DiscreteCallback, CallbackSet, set_proposed_dt!,
                 u_modified!
using Printf

const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const TEST_DIR = joinpath(PACKAGE_ROOT, "test")

include("conversion/units.jl")
include("conversion/species.jl")

include("interop/state.jl")
include("interop/library.jl")
include("interop/metadata.jl")
include("interop/lifecycle.jl")
include("interop/thermo.jl")
include("interop/rhs_api.jl")

include("config/reactor.jl")
include("config/models.jl")
include("config/numerics.jl")
include("config/runtime.jl")
include("config/residence.jl")
include("config/wall.jl")
include("config/sources.jl")
include("config/config.jl")
include("config/validation.jl")
include("config/conversions.jl")

include("runtime/paths.jl")
include("runtime/session.jl")

include("results/reactor.jl")
include("results/chain.jl")

include("chain/profile.jl")
include("chain/marching.jl")

include("io/logging.jl")
include("io/input.jl")
include("io/profile.jl")
include("io/results.jl")

include("solver/api_layout.jl")
include("solver/initial_state.jl")
include("solver/state_vector.jl")
include("solver/residence_time.jl")
include("solver/wall_losses.jl")
include("solver/source_terms.jl")
include("solver/rhs.jl")
include("solver/integrate_0d.jl")
include("solver/driver.jl")
include("solver/chain_cstr.jl")

export initialize_terra, finalize_terra
export Config, ReactorConfig, ReactorComposition, ReactorThermalState
export ModelConfig, TimeConfig, ODESolverConfig, SpaceConfig
export NumericsConfig, LoggingConfig, RuntimeConfig, ResidenceTimeConfig,
       SourceTermsConfig
export SpeciesWallModel, WallLossConfig, ChainWallProfile
export ChainProfileInletComposition, ChainProfileInlet
export AxialChainProfile, AxialMarchingConfig
export ReactorFrame, ReactorResult, ChainCellResult, ChainMetadata
export ChainSimulationResult, load_chain_profile
export species_density_matrix, temperature_history, total_energy_history
export with_case_path, with_time, with_runtime, with_logging
export solve_terra_0d, nitrogen_10ev_example, save_results, load_results_chain
export solve_terra_chain_steady

end
