@testset "validate_config" begin
    @testset "Valid Inputs" begin
        species = ["N", "N2", "E-"]
        mole_fractions = [0.1, 0.8, 0.1]
        total_number_density = 1e13
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
        time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
        case_path = pwd()
        unit_system = :CGS

        @test terra.validate_config(species, mole_fractions, total_number_density,
            temperatures, time_params, case_path, unit_system) == true
    end

    @testset "Species and Mole Fractions Validation" begin
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
        time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
        case_path = pwd()
        unit_system = :CGS

        # Test mismatched array lengths
        @test_throws ArgumentError terra.validate_config(
            ["N", "N2"], [0.5], 1e13, temperatures, time_params, case_path, unit_system)

        # Test empty species
        @test_throws ArgumentError terra.validate_config(
            String[], Float64[], 1e13, temperatures, time_params, case_path, unit_system)

        # Test mole fractions don't sum to 1
        @test_throws ArgumentError terra.validate_config(
            ["N", "N2"], [0.3, 0.3], 1e13, temperatures,
            time_params, case_path, unit_system)

        # Test negative mole fractions
        @test_throws ArgumentError terra.validate_config(
            ["N", "N2"], [-0.1, 1.1], 1e13, temperatures,
            time_params, case_path, unit_system)

        # Test duplicate species
        @test_throws ArgumentError terra.validate_config(
            ["N", "N"], [0.5, 0.5], 1e13, temperatures,
            time_params, case_path, unit_system)

        # Test empty species name
        @test_throws ArgumentError terra.validate_config(
            ["N", ""], [0.5, 0.5], 1e13, temperatures,
            time_params, case_path, unit_system)
    end

    @testset "Number Density Validation" begin
        species = ["N", "N2"]
        mole_fractions = [0.5, 0.5]
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
        time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
        case_path = pwd()
        unit_system = :CGS

        # Test negative number density
        @test_throws ArgumentError terra.validate_config(
            species, mole_fractions, -1e13, temperatures, time_params, case_path, unit_system)

        # Test zero number density
        @test_throws ArgumentError terra.validate_config(
            species, mole_fractions, 0.0, temperatures, time_params, case_path, unit_system)
    end

    @testset "Path Validation" begin
        species = ["N", "N2"]
        mole_fractions = [0.5, 0.5]
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
        time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
        unit_system = :CGS

        # Test non-existent case path
        @test_throws ArgumentError terra.validate_config(
            species, mole_fractions, 1e13, temperatures,
            time_params, "/nonexistent/path", unit_system)
    end

    @testset "Unit System Validation" begin
        species = ["N", "N2"]
        mole_fractions = [0.5, 0.5]
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
        time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
        case_path = pwd()

        # Test invalid unit system
        @test_throws ArgumentError terra.validate_config(
            species, mole_fractions, 1e13, temperatures,
            time_params, case_path, :INVALID)
    end
end

@testset "TERRA Fortran Interface Validation Tests" begin
    @testset "validate_config_against_terra" begin
        @testset "Basic Functionality" begin
            config = terra.TERRAConfig(
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e13,
                temperatures = terra.TemperatureConfig(;
                    Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                time_params = terra.TimeIntegrationConfig(;
                    dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            )

            # Should return true (may show warnings if library not loaded)
            @test terra.validate_config_against_terra(config) == true
        end

        @testset "With Species Validation Enabled" begin
            config = terra.TERRAConfig(
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e13,
                temperatures = terra.TemperatureConfig(;
                    Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                time_params = terra.TimeIntegrationConfig(;
                    dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                validate_species_against_terra = true
            )

            # Should return true (may show warnings if library not loaded or species not found)
            @test terra.validate_config_against_terra(config) == true
        end
    end

    @testset "validate_species_against_terra_database (Enhanced)" begin
        @testset "Validation Disabled" begin
            config = terra.TERRAConfig(
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e13,
                temperatures = terra.TemperatureConfig(;
                    Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                time_params = terra.TimeIntegrationConfig(;
                    dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                validate_species_against_terra = false
            )

            @test terra.validate_species_against_terra_database(config) == true
        end

        @testset "Validation Enabled" begin
            config = terra.TERRAConfig(
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e13,
                temperatures = terra.TemperatureConfig(;
                    Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
                time_params = terra.TimeIntegrationConfig(;
                    dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
                validate_species_against_terra = true
            )

            # Should return true (may show warnings if library not loaded or species not found)
            @test terra.validate_species_against_terra_database(config) == true
        end
    end
end
