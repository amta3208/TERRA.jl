using TERRA: TERRA as terra
using Test
using Aqua

include("helpers/shared.jl")

@testset "TERRA.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(terra; ambiguities = false, persistent_tasks = false)
    end

    @testset "Public API" begin
        include("public/exports.jl")
    end

    @testset "Data Conversion" begin
        include("conversion/units.jl")
        include("conversion/species.jl")
    end

    @testset "Fortran Wrapper" begin
        include("interop/library.jl")
        include("interop/lifecycle.jl")
        include("interop/metadata.jl")
        include("interop/thermo.jl")
    end

    @testset "Config" begin
        include("config/reactor.jl")
        include("config/models.jl")
        include("config/numerics.jl")
        include("config/runtime.jl")
        include("config/sources.jl")
        include("config/config.jl")
        include("config/validation.jl")
        include("config/conversions.jl")
    end

    @testset "Runtime" begin
        include("runtime/paths.jl")
        include("runtime/session.jl")
    end

    @testset "Results" begin
        include("results/reactor.jl")
        include("results/chain.jl")
    end

    @testset "Chain" begin
        include("chain/profile.jl")
        include("chain/marching.jl")
    end

    @testset "IO" begin
        include("io/input.jl")
        include("io/profile.jl")
        include("io/export_chain_profile.jl")
        include("io/results.jl")
    end

    @testset "Solver" begin
        include("solver/initial_state.jl")
        include("solver/state_vector.jl")
        include("solver/rhs.jl")
        include("solver/residence_time.jl")
        include("solver/wall_losses.jl")
        include("solver/driver.jl")
        include("solver/integrate_0d.jl")
        include("solver/chain_cstr.jl")
    end
end
