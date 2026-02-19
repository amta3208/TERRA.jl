@testset "Data Validation and Species Utilities" begin
    @testset "Species Data Validation" begin
        # Test valid species data
        species_names = ["N2", "N", "N+", "N2+", "E-"]
        terra_species = ["N2", "N", "N+", "N2+", "E-", "O2", "O"]
        densities = [1e-3, 1e-6, 1e-7, 1e-7, 1e-7]  # g/cm³

        @test terra.validate_species_data(species_names, terra_species, densities) == true

        # Test mismatched array lengths
        @test_throws ErrorException terra.validate_species_data(
            species_names[1:3], terra_species, densities)
        @test_throws ErrorException terra.validate_species_data(
            species_names, terra_species, densities[1:3])

        # Test negative densities
        negative_densities = [1e-3, -1e-6, 1e-7, 1e-7, 1e-7]
        @test_throws ErrorException terra.validate_species_data(
            species_names, terra_species, negative_densities)

        # Test invalid species names
        invalid_species = ["N2", "N", "N+", "N2+", "INVALID"]
        @test_throws ErrorException terra.validate_species_data(
            invalid_species, terra_species, densities)

        # Test very small densities (should warn but not error)
        small_densities = [1e-3, 1e-35, 1e-7, 1e-7, 1e-7]
        @test_logs (:warn, r"Very small densities detected") terra.validate_species_data(
            species_names, terra_species, small_densities)

        # Test with zero densities (valid)
        zero_densities = [1e-3, 0.0, 1e-7, 1e-7, 1e-7]
        @test terra.validate_species_data(species_names, terra_species, zero_densities) ==
              true

        # Test with realistic Hall thruster species
        ht_species = ["Xe", "Xe+"]
        ht_terra_species = ["Xe", "Xe+", "Xe2+", "E-"]
        ht_densities = [1e-4, 1e-5]
        @test terra.validate_species_data(ht_species, ht_terra_species, ht_densities) ==
              true
    end

    @testset "Species Mapping" begin
        # Test common species mappings
        ht_species = ["N", "N2", "N+", "N2+", "e-"]
        terra_species = ["N", "N2", "N+", "N2+", "E-"]

        mapping = terra.create_species_mapping(ht_species, terra_species)

        @test mapping["N"] == "N"
        @test mapping["N2"] == "N2"
        @test mapping["N+"] == "N+"
        @test mapping["N2+"] == "N2+"
        @test mapping["e-"] == "E-"

        # Test with noble gas species
        noble_ht = ["Xe", "Xe+", "Ar", "Ar+", "Kr", "Kr+"]
        noble_terra = ["Xe", "Xe+", "Ar", "Ar+", "Kr", "Kr+"]

        noble_mapping = terra.create_species_mapping(noble_ht, noble_terra)

        @test noble_mapping["Xe"] == "Xe"
        @test noble_mapping["Xe+"] == "Xe+"
        @test noble_mapping["Ar"] == "Ar"
        @test noble_mapping["Ar+"] == "Ar+"
        @test noble_mapping["Kr"] == "Kr"
        @test noble_mapping["Kr+"] == "Kr+"

        # Test direct mapping (same names)
        direct_ht = ["N2", "O2"]
        direct_terra = ["N2", "O2", "N", "O"]

        direct_mapping = terra.create_species_mapping(direct_ht, direct_terra)
        @test direct_mapping["N2"] == "N2"
        @test direct_mapping["O2"] == "O2"

        # Test error for missing species in TERRA database
        missing_ht = ["N2", "MISSING_SPECIES"]
        @test_throws ErrorException terra.create_species_mapping(missing_ht, terra_species)

        # Test error for mapped species missing from TERRA database
        mapped_but_missing_ht = ["N"]  # "N" maps to "N" in common_mappings
        terra_without_mapped = ["N2", "O", "E-"]  # Doesn't include "N"
        @test_throws ErrorException terra.create_species_mapping(
            mapped_but_missing_ht, terra_without_mapped)

        # Test error for unmappable species
        unmappable_ht = ["UNKNOWN_SPECIES"]
        unmappable_terra = ["N2", "N"]
        @test_throws ErrorException terra.create_species_mapping(
            unmappable_ht, unmappable_terra)

        # Test E- vs e- mapping
        electron_ht = ["E-"]
        electron_terra = ["E-"]
        electron_mapping = terra.create_species_mapping(electron_ht, electron_terra)
        @test electron_mapping["E-"] == "E-"
    end
end


@testset "get_molecular_weights" begin
    @testset "Known Species" begin
        # Test nitrogen species
        weights = terra.get_molecular_weights(["N", "N2", "N+", "N2+", "E-"])
        @test weights[1] ≈ 14.007  # N
        @test weights[2] ≈ 28.014  # N2
        @test weights[3] ≈ 14.007  # N+
        @test weights[4] ≈ 28.014  # N2+
        @test weights[5] ≈ 5.485799e-4  # E-

        # Test noble gas species
        weights_ar = terra.get_molecular_weights(["Ar", "Ar+"])
        @test weights_ar[1] ≈ 39.948  # Ar
        @test weights_ar[2] ≈ 39.948  # Ar+

        weights_xe = terra.get_molecular_weights(["Xe", "Xe+"])
        @test weights_xe[1] ≈ 131.293  # Xe
        @test weights_xe[2] ≈ 131.293  # Xe+
    end

    @testset "Unknown Species" begin
        @test_throws ErrorException terra.get_molecular_weights(["UNKNOWN"])
        @test_throws ErrorException terra.get_molecular_weights(["N", "UNKNOWN", "E-"])
    end
end
