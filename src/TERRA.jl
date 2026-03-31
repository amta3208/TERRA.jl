module TERRA

using Dates
using DocStringExtensions
using Libdl
using JSON
using OrdinaryDiffEq
using PrecompileTools: @compile_workload, @setup_workload
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

include("io/codec.jl")
include("io/logging.jl")
include("io/input.jl")
include("io/profile.jl")
include("io/results.jl")

include("sources/core.jl")
include("sources/residence.jl")
include("sources/wall.jl")

include("reactor/integrate.jl")
include("reactor/solve.jl")
include("reactor/examples.jl")

include("chain/diagnostics.jl")
include("chain/solve.jl")

export initialize_terra, finalize_terra
export solve_terra_0d, solve_terra_chain_steady
export load_chain_profile, save_results, load_results_chain
export species_density_matrix, temperature_history, total_energy_history
export WallLossConfig
export IonNeutralizationWallModel, BallisticNeutralRecombinationWallModel,
       ConstantNeutralRecombinationWallModel
export ChainWallProfile, ChainProfileInletComposition, ChainProfileInlet
export AxialChainProfile, AxialMarchingConfig
export FullStateHandoff, ReinitializeHandoff, FinalTimeTermination
export ReactorResult, ChainSimulationResult

include("precompile.jl")

end
