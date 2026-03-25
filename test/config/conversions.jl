@testset "convert_config_units" begin
    @testset "No Conversion Needed" begin
        config = terra.nitrogen_10ev_config()  # Default is CGS
        converted = terra.convert_config_units(config, :CGS)
        @test converted === config  # Should return same object
    end

    @testset "SI to CGS Conversion" begin
        base = terra.nitrogen_10ev_config()

        config_si = terra.Config(;
            reactor = terra.ReactorConfig(;
                composition = terra.ReactorComposition(;
                    species = base.reactor.composition.species,
                    mole_fractions = base.reactor.composition.mole_fractions,
                    total_number_density = 1e19),  # SI: 1/m^3
                thermal = base.reactor.thermal),
            models = base.models,
            sources = base.sources,
            numerics = base.numerics,
            runtime = terra.RuntimeConfig(;
                database_path = base.runtime.database_path,
                case_path = base.runtime.case_path,
                unit_system = :SI,
                validate_species_against_terra = base.runtime.validate_species_against_terra,
                print_source_terms = base.runtime.print_source_terms,
                write_native_state_files = base.runtime.write_native_state_files,
                logging = base.runtime.logging))

        converted = terra.convert_config_units(config_si, :CGS)
        @test converted.runtime.unit_system == :CGS
        @test converted.reactor.composition.total_number_density ≈ 1e13  # 1/cm^3
        @test converted.reactor.composition.species == config_si.reactor.composition.species
        @test converted.reactor.composition.mole_fractions ==
              config_si.reactor.composition.mole_fractions
        @test converted.runtime.write_native_state_files ==
              config_si.runtime.write_native_state_files
        @test converted.runtime.logging == config_si.runtime.logging
        @test converted.sources.residence_time === config_si.sources.residence_time
        @test converted.sources.wall_losses === config_si.sources.wall_losses
    end

    @testset "CGS to SI Conversion" begin
        config_cgs = terra.nitrogen_10ev_config()  # CGS by default
        converted = terra.convert_config_units(config_cgs, :SI)

        @test converted.runtime.unit_system == :SI
        @test converted.reactor.composition.total_number_density ≈ 1e19  # 1/m^3
        @test converted.reactor.composition.species == config_cgs.reactor.composition.species
        @test converted.reactor.composition.mole_fractions ==
              config_cgs.reactor.composition.mole_fractions
        @test converted.runtime.write_native_state_files ==
              config_cgs.runtime.write_native_state_files
        @test converted.runtime.logging == config_cgs.runtime.logging
        @test converted.sources.residence_time === config_cgs.sources.residence_time
        @test converted.sources.wall_losses === config_cgs.sources.wall_losses
    end

    @testset "Invalid Conversion" begin
        config = terra.nitrogen_10ev_config()
        @test_throws ErrorException terra.convert_config_units(config, :INVALID)
    end
end
