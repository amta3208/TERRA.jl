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
include("interop/api_layout.jl")
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

include("reactor/state.jl")
include("reactor/rhs.jl")
include("reactor/initial.jl")

include("io/logging.jl")
include("io/input.jl")
include("io/profile.jl")
include("io/results.jl")

include("sources/residence.jl")
include("sources/wall.jl")
include("sources/core.jl")

include("reactor/integrate.jl")
include("reactor/solve.jl")
include("reactor/examples.jl")
include("solver/chain_cstr.jl")

export initialize_terra, finalize_terra
export Config, ReactorConfig, ReactorComposition, ReactorThermalState
export ModelConfig, TimeConfig, ODESolverConfig, SpaceConfig
export NumericsConfig, LoggingConfig, RuntimeConfig, ResidenceTimeConfig,
       SourceTermsConfig
export SpeciesWallModel, IonNeutralizationWallModel,
       BallisticNeutralRecombinationWallModel,
       ConstantNeutralRecombinationWallModel, WallLossConfig, ChainWallProfile
export ChainProfileInletComposition, ChainProfileInlet
export AxialChainProfile, AxialMarchingConfig
export ReactorFrame, ReactorResult, ChainCellResult, ChainMetadata
export ChainSimulationResult, load_chain_profile
export species_density_matrix, temperature_history, total_energy_history
export with_case_path, with_time, with_runtime, with_logging
export solve_terra_0d, nitrogen_10ev_example, save_results, load_results_chain
export solve_terra_chain_steady

end
