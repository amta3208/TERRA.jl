using TERRA: TERRA as terra
using Test
using Aqua

include("helpers/shared.jl")

@testset "TERRA.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(terra; ambiguities = false, persistent_tasks = false)
    end

    @testset "Data Conversion" begin
        include("data_conversion.jl")
    end

    @testset "Fortran Wrapper" begin
        include("fortran_wrapper.jl")
    end

    @testset "Config" begin
        include("config/types.jl")
        include("config/validation.jl")
        include("config/conversions.jl")
    end

    @testset "IO" begin
        include("io/input_generation.jl")
    end

    @testset "Solver" begin
        include("solver/state_vector.jl")
        include("solver/rhs.jl")
        include("solver/residence_time.jl")
        include("solver/driver.jl")
        include("solver/integrate_0d.jl")
    end
end
