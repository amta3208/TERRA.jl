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
        include("interop/api_layout.jl")
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
        include("chain/solve.jl")
    end

    @testset "Reactor" begin
        include("reactor/initial.jl")
        include("reactor/state.jl")
        include("reactor/rhs.jl")
        include("reactor/solve.jl")
        include("reactor/integrate.jl")
        if "benchmarks" in ARGS
            include("reactor/benchmarks.jl")
        end
    end

    @testset "Sources" begin
        include("sources/residence.jl")
        include("sources/wall.jl")
    end

    @testset "IO" begin
        include("io/input.jl")
        include("io/logging.jl")
        include("io/profile.jl")
        include("io/export_chain_profile.jl")
        include("io/results.jl")
    end

end
