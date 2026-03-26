@testset "ReactorComposition" begin
    @testset "Valid Construction" begin
        composition = terra.ReactorComposition(;
                                               species = ["N", "N2", "E-"],
                                               mole_fractions = [0.1, 0.8, 0.1],
                                               total_number_density = 1e13)
        @test composition.species == ["N", "N2", "E-"]
        @test composition.mole_fractions == [0.1, 0.8, 0.1]
        @test composition.total_number_density == 1e13
    end

    @testset "Invalid Construction" begin
        @test_throws ArgumentError terra.ReactorComposition(;
                                                            species = ["N", "N2"],
                                                            mole_fractions = [0.5],
                                                            total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
                                                            species = String[],
                                                            mole_fractions = Float64[],
                                                            total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
                                                            species = ["N", "N2"],
                                                            mole_fractions = [0.3, 0.3],
                                                            total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
                                                            species = ["N", "N"],
                                                            mole_fractions = [0.5, 0.5],
                                                            total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
                                                            species = ["N", "N2"],
                                                            mole_fractions = [0.5, 0.5],
                                                            total_number_density = 0.0)
    end
end

@testset "ReactorThermalState" begin
    @testset "Valid Construction" begin
        thermal = terra.ReactorThermalState(; Tt = 300.0, Tv = 310.0, Tee = 320.0,
                                            Te = 10000.0)
        @test thermal.Tt == 300.0
        @test thermal.Tv == 310.0
        @test thermal.Tee == 320.0
        @test thermal.Te == 10000.0
    end

    @testset "Invalid Construction" begin
        @test_throws ArgumentError terra.ReactorThermalState(; Tt = 0.0, Tv = 310.0,
                                                             Tee = 320.0, Te = 10000.0)
        @test_throws ArgumentError terra.ReactorThermalState(; Tt = 300.0, Tv = -1.0,
                                                             Tee = 320.0, Te = 10000.0)
    end
end
