using TERRA: TERRA as terra
using Test
using Aqua

include("helpers/shared.jl")

@progress_testset "TERRA.jl" begin
    @progress_testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(terra; ambiguities = false, persistent_tasks = false)
    end

    @progress_testset "Public API" begin
        include_with_progress("Public API", "public/exports.jl")
    end

    @progress_testset "Data Conversion" begin
        include_with_progress("Data Conversion", "conversion/units.jl")
        include_with_progress("Data Conversion", "conversion/species.jl")
    end

    @progress_testset "Fortran Wrapper" begin
        include_with_progress("Fortran Wrapper", "interop/library.jl")
        include_with_progress("Fortran Wrapper", "interop/lifecycle.jl")
        include_with_progress("Fortran Wrapper", "interop/metadata.jl")
        include_with_progress("Fortran Wrapper", "interop/api_layout.jl")
        include_with_progress("Fortran Wrapper", "interop/thermo.jl")
    end

    @progress_testset "Config" begin
        include_with_progress("Config", "config/reactor.jl")
        include_with_progress("Config", "config/models.jl")
        include_with_progress("Config", "config/numerics.jl")
        include_with_progress("Config", "config/runtime.jl")
        include_with_progress("Config", "config/sources.jl")
        include_with_progress("Config", "config/config.jl")
        include_with_progress("Config", "config/validation.jl")
        include_with_progress("Config", "config/conversions.jl")
    end

    @progress_testset "Runtime" begin
        include_with_progress("Runtime", "runtime/paths.jl")
        include_with_progress("Runtime", "runtime/session.jl")
    end

    @progress_testset "Results" begin
        include_with_progress("Results", "results/reactor.jl")
        include_with_progress("Results", "results/chain.jl")
    end

    @progress_testset "Chain" begin
        include_with_progress("Chain", "chain/profile.jl")
        include_with_progress("Chain", "chain/marching.jl")
        include_with_progress("Chain", "chain/solve.jl")
    end

    @progress_testset "Reactor" begin
        include_with_progress("Reactor", "reactor/initial.jl")
        include_with_progress("Reactor", "reactor/state.jl")
        include_with_progress("Reactor", "reactor/rhs.jl")
        include_with_progress("Reactor", "reactor/solve.jl")
        include_with_progress("Reactor", "reactor/integrate.jl")
    end

    @progress_testset "Sources" begin
        include_with_progress("Sources", "sources/residence.jl")
        include_with_progress("Sources", "sources/wall.jl")
    end

    @progress_testset "IO" begin
        include_with_progress("IO", "io/input.jl")
        include_with_progress("IO", "io/logging.jl")
        include_with_progress("IO", "io/profile.jl")
        include_with_progress("IO", "io/export_chain_profile.jl")
        include_with_progress("IO", "io/results.jl")
    end

end
