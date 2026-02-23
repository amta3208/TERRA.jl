@testset "TERRA Fortran Interface Validation Tests" begin
    build_config(; validate_species_against_terra::Bool = false) = terra.Config(;
        reactor = terra.ReactorConfig(;
            composition = terra.ReactorComposition(;
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e13),
            thermal = terra.ReactorThermalState(;
                Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)),
        numerics = terra.NumericsConfig(;
            time = terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4, duration = 1e-3),
            residence_time = nothing),
        runtime = terra.RuntimeConfig(;
            case_path = pwd(),
            unit_system = :CGS,
            validate_species_against_terra = validate_species_against_terra))

    @testset "validate_config_against_terra" begin
        @testset "Basic Functionality" begin
            config = build_config()

            # Should return true (may show warnings if library not loaded)
            @test terra.validate_config_against_terra(config) == true
        end

        @testset "With Species Validation Enabled" begin
            config = build_config(; validate_species_against_terra = true)

            # Should return true (may show warnings if library not loaded or species not found)
            @test terra.validate_config_against_terra(config) == true
        end
    end

    @testset "validate_species_against_terra_database (Enhanced)" begin
        @testset "Validation Disabled" begin
            config = build_config(; validate_species_against_terra = false)

            @test terra.validate_species_against_terra_database(config) == true
        end

        @testset "Validation Enabled" begin
            config = build_config(; validate_species_against_terra = true)

            # Should return true (may show warnings if library not loaded or species not found)
            @test terra.validate_species_against_terra_database(config) == true
        end
    end
end
