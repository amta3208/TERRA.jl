@testset "convert_config_units" begin
    @testset "No Conversion Needed" begin
        config = terra.nitrogen_10ev_config()  # Default is CGS
        converted = terra.convert_config_units(config, :CGS)
        @test converted === config  # Should return same object
    end

    @testset "SI to CGS Conversion" begin
        # Create SI config
        config_si = terra.TERRAConfig(
            species = ["N", "N2", "E-"],
            mole_fractions = [0.1, 0.8, 0.1],
            total_number_density = 1e19,  # SI: 1/m³
            temperatures = terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0),
            time_params = terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
            unit_system = :SI
        )

        converted = terra.convert_config_units(config_si, :CGS)
        @test converted.unit_system == :CGS
        @test converted.total_number_density ≈ 1e13  # Converted to 1/cm³
        @test converted.species == config_si.species
        @test converted.mole_fractions == config_si.mole_fractions
        @test converted.write_native_outputs == config_si.write_native_outputs
    end

    @testset "CGS to SI Conversion" begin
        config_cgs = terra.nitrogen_10ev_config()  # CGS by default
        converted = terra.convert_config_units(config_cgs, :SI)

        @test converted.unit_system == :SI
        @test converted.total_number_density ≈ 1e19  # Converted to 1/m³
        @test converted.species == config_cgs.species
        @test converted.mole_fractions == config_cgs.mole_fractions
        @test converted.write_native_outputs == config_cgs.write_native_outputs
    end

    @testset "Nested Config Conversion" begin
        legacy = terra.nitrogen_10ev_config()
        nested = terra.to_config(legacy)

        converted_nested = terra.convert_config_units(nested, :SI)
        @test converted_nested.runtime.unit_system == :SI
        @test converted_nested.reactor.composition.total_number_density ≈ 1e19
        @test converted_nested.reactor.composition.species ==
              nested.reactor.composition.species
        @test converted_nested.models.physics == nested.models.physics
        @test converted_nested.models.processes == nested.models.processes
    end

    @testset "Invalid Conversion" begin
        config = terra.nitrogen_10ev_config()
        @test_throws ErrorException terra.convert_config_units(config, :INVALID)
    end
end
