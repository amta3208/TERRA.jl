@testset "Unit Conversions" begin
    @testset "Density Conversions" begin
        # Test SI to CGS density conversion
        rho_si = [1000.0, 2000.0, 500.0]  # kg/m³
        rho_cgs_expected = [1.0, 2.0, 0.5]  # g/cm³
        rho_cgs = terra.convert_density_si_to_cgs(rho_si)
        @test rho_cgs ≈ rho_cgs_expected

        # Test CGS to SI density conversion
        rho_si_back = terra.convert_density_cgs_to_si(rho_cgs)
        @test rho_si_back ≈ rho_si

        # Test with zero values
        @test terra.convert_density_si_to_cgs([0.0]) ≈ [0.0]
        @test terra.convert_density_cgs_to_si([0.0]) ≈ [0.0]

        # Test with very small values
        small_si = [1e-10]
        small_cgs = terra.convert_density_si_to_cgs(small_si)
        @test terra.convert_density_cgs_to_si(small_cgs) ≈ small_si

        # Test with very large values
        large_si = [1e10]
        large_cgs = terra.convert_density_si_to_cgs(large_si)
        @test terra.convert_density_cgs_to_si(large_cgs) ≈ large_si

        # Accept integer vectors (type-widening)
        int_si = [1, 2, 3]
        @test terra.convert_density_si_to_cgs(int_si) ≈ [1e-3, 2e-3, 3e-3]
    end

    @testset "Energy Density Conversions" begin
        # Test SI to CGS energy density conversion
        energy_si = 1e6  # J/m³
        energy_cgs_expected = 1e7  # erg/cm³
        energy_cgs = terra.convert_energy_density_si_to_cgs(energy_si)
        @test energy_cgs ≈ energy_cgs_expected

        # Test CGS to SI energy density conversion
        energy_si_back = terra.convert_energy_density_cgs_to_si(energy_cgs)
        @test energy_si_back ≈ energy_si

        # Test with zero
        @test terra.convert_energy_density_si_to_cgs(0.0) ≈ 0.0
        @test terra.convert_energy_density_cgs_to_si(0.0) ≈ 0.0

        # Test roundtrip with various values
        test_values = [1.0, 1e3, 1e6, 1e9, 1e-6]
        for val in test_values
            converted = terra.convert_energy_density_si_to_cgs(val)
            back = terra.convert_energy_density_cgs_to_si(converted)
            @test back≈val rtol=1e-12
        end
    end

    @testset "Pressure Conversions" begin
        # Test SI to CGS pressure conversion
        pressure_si = 101325.0  # Pa (1 atm)
        pressure_cgs_expected = 1.01325e6  # dyne/cm²
        pressure_cgs = terra.convert_pressure_si_to_cgs(pressure_si)
        @test pressure_cgs ≈ pressure_cgs_expected

        # Test CGS to SI pressure conversion
        pressure_si_back = terra.convert_pressure_cgs_to_si(pressure_cgs)
        @test pressure_si_back ≈ pressure_si

        # Test with typical Hall thruster pressures
        ht_pressure_si = 1e-2  # Pa (low pressure)
        ht_pressure_cgs = terra.convert_pressure_si_to_cgs(ht_pressure_si)
        @test terra.convert_pressure_cgs_to_si(ht_pressure_cgs) ≈ ht_pressure_si

        # Test roundtrip consistency
        test_pressures = [1.0, 1e3, 1e5, 1e-3, 1e-6]
        for p in test_pressures
            converted = terra.convert_pressure_si_to_cgs(p)
            back = terra.convert_pressure_cgs_to_si(converted)
            @test back≈p rtol=1e-12
        end
    end

    @testset "Number Density Conversions" begin
        # Test SI to CGS number density conversion
        n_si = 1e20  # 1/m³
        n_cgs_expected = 1e14  # 1/cm³
        n_cgs = terra.convert_number_density_si_to_cgs(n_si)
        @test n_cgs ≈ n_cgs_expected

        # Test CGS to SI number density conversion
        n_si_back = terra.convert_number_density_cgs_to_si(n_cgs)
        @test n_si_back ≈ n_si

        # Test with typical plasma densities
        plasma_densities_si = [1e15, 1e17, 1e19]  # 1/m³
        plasma_densities_cgs_expected = [1e9, 1e11, 1e13]  # 1/cm³

        for (si_val, cgs_expected) in zip(
            plasma_densities_si, plasma_densities_cgs_expected)
            cgs_val = terra.convert_number_density_si_to_cgs(si_val)
            @test cgs_val ≈ cgs_expected
            @test terra.convert_number_density_cgs_to_si(cgs_val) ≈ si_val
        end

        # Test roundtrip with edge cases
        edge_cases = [1.0, 1e25, 1e-10]
        for n in edge_cases
            converted = terra.convert_number_density_si_to_cgs(n)
            back = terra.convert_number_density_cgs_to_si(converted)
            @test back≈n rtol=1e-12
        end
    end
end

@testset "State Vector Conversions" begin
    @testset "SI to CGS State Conversion" begin
        # Test basic state conversion
        rho_sp_si = [1000.0, 2000.0, 500.0]  # kg/m³
        rho_etot_si = 1e6  # J/m³
        number_density_si = 1e20  # 1/m³

        result = terra.convert_state_si_to_cgs(rho_sp_si, rho_etot_si, number_density_si)

        @test result.rho_sp ≈ [1.0, 2.0, 0.5]  # g/cm³
        @test result.rho_etot ≈ 1e7  # erg/cm³
        @test result.number_density ≈ 1e14  # 1/cm³
        @test result.rho_erot === nothing
        @test result.rho_eeex === nothing
        @test result.rho_evib === nothing

        # Test with optional energy components
        rho_erot_si = 1e5  # J/m³
        rho_eeex_si = 2e5  # J/m³
        rho_evib_si = 3e5  # J/m³

        result_full = terra.convert_state_si_to_cgs(
            rho_sp_si, rho_etot_si, number_density_si;
            rho_erot_si = rho_erot_si,
            rho_eeex_si = rho_eeex_si,
            rho_evib_si = rho_evib_si)

        @test result_full.rho_erot ≈ 1e6  # erg/cm³
        @test result_full.rho_eeex ≈ 2e6  # erg/cm³
        @test result_full.rho_evib ≈ 3e6  # erg/cm³

        # Test with nitrogen plasma conditions (from TERRA example)
        # N2 dominant plasma at 10 eV electron temperature
        n2_density_si = 1e17 * 28 * 1.66e-27 * 1e3  # Convert from number density to mass density
        n_density_si = 1e13 * 14 * 1.66e-27 * 1e3
        electron_energy_si = 1e17 * 10 * 1.6e-19 * 1e9  # 10 eV per electron

        nitrogen_state = terra.convert_state_si_to_cgs([n2_density_si, n_density_si],
            electron_energy_si, 1e17)

        @test all(nitrogen_state.rho_sp .> 0)
        @test nitrogen_state.rho_etot > 0
        @test nitrogen_state.number_density > 0
    end

    @testset "CGS to SI State Conversion" begin
        # Test basic state conversion back
        rho_sp_cgs = [1.0, 2.0, 0.5]  # g/cm³
        rho_etot_cgs = 1e4  # erg/cm³

        result = terra.convert_state_cgs_to_si(rho_sp_cgs, rho_etot_cgs)

        @test result.rho_sp ≈ [1000.0, 2000.0, 500.0]  # kg/m³
        @test result.rho_etot ≈ 1e3  # J/m³
        @test result.rho_erot === nothing

        # Test with optional energy components
        result_full = terra.convert_state_cgs_to_si(rho_sp_cgs, rho_etot_cgs;
            rho_erot_cgs = 1e3,
            rho_eeex_cgs = 2e3,
            rho_evib_cgs = 3e3)

        @test result_full.rho_erot ≈ 1e2  # J/m³
        @test result_full.rho_eeex ≈ 2e2  # J/m³
        @test result_full.rho_evib ≈ 3e2  # J/m³

        # Test roundtrip consistency
        original_rho_sp = [1500.0, 3000.0, 750.0]
        original_rho_etot = 5e6
        original_number_density = 2e20

        cgs_state = terra.convert_state_si_to_cgs(
            original_rho_sp, original_rho_etot, original_number_density)
        si_state = terra.convert_state_cgs_to_si(cgs_state.rho_sp, cgs_state.rho_etot)

        @test si_state.rho_sp≈original_rho_sp rtol=1e-12
        @test si_state.rho_etot≈original_rho_etot rtol=1e-12
    end

    @testset "Source Terms Conversion" begin
        # Test source terms conversion (same units as state per time)
        drho_sp_cgs = [0.1, 0.2, -0.05]  # g/cm³/s
        drho_etot_cgs = 100.0  # erg/cm³/s

        result = terra.convert_sources_cgs_to_si(drho_sp_cgs, drho_etot_cgs)

        @test result.drho_sp ≈ [100.0, 200.0, -50.0]  # kg/m³/s
        @test result.drho_etot ≈ 1e1  # J/m³/s

        # Test with optional source terms
        result_full = terra.convert_sources_cgs_to_si(drho_sp_cgs, drho_etot_cgs;
            drho_erot_cgs = 10.0,
            drho_eeex_cgs = 20.0,
            drho_evib_cgs = 30.0)

        @test result_full.drho_erot ≈ 1.0  # J/m³/s
        @test result_full.drho_eeex ≈ 2.0  # J/m³/s
        @test result_full.drho_evib ≈ 3.0  # J/m³/s

        # Test with zero source terms
        zero_result = terra.convert_sources_cgs_to_si([0.0, 0.0], 0.0)
        @test all(zero_result.drho_sp .≈ 0.0)
        @test zero_result.drho_etot ≈ 0.0

        # Test with negative source terms (destruction)
        negative_sources = terra.convert_sources_cgs_to_si([-0.1, -0.2], -50.0)
        @test all(negative_sources.drho_sp .< 0)
        @test negative_sources.drho_etot < 0
    end
end

@testset "Data Validation and Preparation" begin
    @testset "Array Preparation for Fortran" begin
        # Test with various array types
        julia_array = [1.0, 2.0, 3.0]
        scalar_val = 5.0
        nothing_val = nothing

        prepared = terra.prepare_arrays_for_fortran(julia_array, scalar_val, nothing_val)

        @test prepared[1] isa Array{Float64}
        @test prepared[1] == [1.0, 2.0, 3.0]
        @test prepared[2] isa Float64
        @test prepared[2] == 5.0
        @test prepared[3] === nothing

        # Test with integer arrays (should convert to Float64)
        int_array = [1, 2, 3]
        prepared_int = terra.prepare_arrays_for_fortran(int_array)
        @test prepared_int[1] isa Array{Float64}
        @test prepared_int[1] == [1.0, 2.0, 3.0]

        # Test with 2D arrays
        matrix = [1.0 2.0; 3.0 4.0]
        prepared_matrix = terra.prepare_arrays_for_fortran(matrix)
        @test prepared_matrix[1] isa Array{Float64}
        @test size(prepared_matrix[1]) == (2, 2)
        # No copy if already Array{Float64}
        @test prepared_matrix[1] === matrix

        # Test with empty arrays
        empty_array = Float64[]
        prepared_empty = terra.prepare_arrays_for_fortran(empty_array)
        @test prepared_empty[1] isa Array{Float64}
        @test length(prepared_empty[1]) == 0
        @test prepared_empty[1] === empty_array

        # Test with multiple arrays
        arr1 = [1.0, 2.0]
        arr2 = [3.0, 4.0, 5.0]
        scalar = 10.0
        prepared_multi = terra.prepare_arrays_for_fortran(arr1, arr2, scalar)

        @test length(prepared_multi) == 3
        @test prepared_multi[1] == [1.0, 2.0]
        @test prepared_multi[2] == [3.0, 4.0, 5.0]
        @test prepared_multi[3] == 10.0
    end

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

@testset "Chemistry Conversions" begin
    @testset "Mole Fractions to Mass Densities" begin
        # Test with nitrogen species (from TERRA example)
        species_names = ["N2", "N", "N+", "N2+", "E-"]
        mole_fractions = [0.9998, 1e-20, 1e-20, 0.0001, 0.0001]
        molecular_weights = [28.0, 14.0, 14.0, 28.0, 0.000549]  # g/mol (electron mass very small)
        total_number_density = 1e13  # 1/cm³

        mass_densities = terra.mole_fractions_to_mass_densities(
            mole_fractions, molecular_weights, total_number_density)

        # Check that results are positive
        @test all(mass_densities .≥ 0)

        # Check that N2 dominates (highest mole fraction and molecular weight)
        @test mass_densities[1] > mass_densities[2]  # N2 > N
        @test mass_densities[1] > mass_densities[3]  # N2 > N+

        # Check approximate values
        expected_n2_density = 0.9998 * 1e13 * 28.0 / 6.02214076e23
        @test mass_densities[1]≈expected_n2_density rtol=1e-10

        # Test with equal mole fractions
        equal_fractions = [0.25, 0.25, 0.25, 0.25]
        equal_weights = [20.0, 20.0, 20.0, 20.0]
        equal_densities = terra.mole_fractions_to_mass_densities(
            equal_fractions, equal_weights, 1e15)

        @test all(equal_densities .≈ equal_densities[1])  # Should all be equal

        # Test error conditions
        @test_throws ErrorException terra.mole_fractions_to_mass_densities(
            [0.5, 0.5], [20.0], 1e15)  # Mismatched lengths
        @test_throws ErrorException terra.mole_fractions_to_mass_densities(
            [0.6, 0.5], [20.0, 20.0], 1e15)  # Don't sum to 1
        # Per-element bounds violations
        @test_throws ErrorException terra.mole_fractions_to_mass_densities(
            [-1e-8, 1.0 + 1e-8], [20.0, 20.0], 1e15)
        # Non-positive molecular weights
        @test_throws ErrorException terra.mole_fractions_to_mass_densities(
            [0.5, 0.5], [0.0, 20.0], 1e15)
    end

    @testset "Mass Densities to Mole Fractions" begin
        # Test conversion back from mass densities
        molecular_weights = [28.0, 14.0, 14.0, 28.0, 0.000549]
        mass_densities = [1e-9, 1e-12, 1e-13, 1e-12, 1e-16]  # g/cm³

        mole_fractions = terra.mass_densities_to_mole_fractions(
            mass_densities, molecular_weights)

        # Check that mole fractions sum to 1
        @test sum(mole_fractions)≈1.0 rtol=1e-10

        # Check that all fractions are positive
        @test all(mole_fractions .≥ 0)

        # Check that N2 has highest mole fraction (highest mass density, same molecular weight as N2+)
        @test mole_fractions[1] > mole_fractions[2]  # N2 > N
        @test mole_fractions[1] > mole_fractions[3]  # N2 > N+

        # Test with equal mass densities but different molecular weights
        equal_masses = [1e-10, 1e-10, 1e-10]
        different_weights = [10.0, 20.0, 30.0]
        fractions = terra.mass_densities_to_mole_fractions(equal_masses, different_weights)

        # Lighter species should have higher mole fraction
        @test fractions[1] > fractions[2] > fractions[3]

        # Test error conditions
        @test_throws ErrorException terra.mass_densities_to_mole_fractions(
            [1e-10, 1e-10], [20.0])  # Mismatched lengths
        @test_throws ErrorException terra.mass_densities_to_mole_fractions(
            [0.0, 0.0], [20.0, 20.0])  # Zero total density
        # Negative mass densities
        @test_throws ErrorException terra.mass_densities_to_mole_fractions(
            [-1e-10, 1e-10], [20.0, 20.0])
        # Non-positive molecular weights
        @test_throws ErrorException terra.mass_densities_to_mole_fractions(
            [1e-10, 1e-10], [20.0, 0.0])
    end

    @testset "Roundtrip Consistency" begin
        # Test that mole_fractions -> mass_densities -> mole_fractions gives original values
        original_fractions = [0.7, 0.2, 0.05, 0.04, 0.01]
        molecular_weights = [28.0, 14.0, 14.0, 28.0, 0.000549]
        total_number_density = 1e15

        # Forward conversion
        mass_densities = terra.mole_fractions_to_mass_densities(
            original_fractions, molecular_weights, total_number_density)

        # Backward conversion
        recovered_fractions = terra.mass_densities_to_mole_fractions(
            mass_densities, molecular_weights)

        # Check roundtrip consistency
        @test recovered_fractions≈original_fractions rtol=1e-12

        # Test with extreme values
        extreme_fractions = [0.99999, 1e-10, 1e-10, 1e-10, 1e-15]
        extreme_fractions = extreme_fractions ./ sum(extreme_fractions)  # Normalize

        extreme_masses = terra.mole_fractions_to_mass_densities(
            extreme_fractions, molecular_weights, 1e20)
        extreme_recovered = terra.mass_densities_to_mole_fractions(
            extreme_masses, molecular_weights)

        @test extreme_recovered≈extreme_fractions rtol=1e-10

        # Test with single species
        single_fraction = [1.0]
        single_weight = [20.0]
        single_mass = terra.mole_fractions_to_mass_densities(
            single_fraction, single_weight, 1e16)
        single_recovered = terra.mass_densities_to_mole_fractions(
            single_mass, single_weight)

        @test single_recovered≈[1.0] rtol=1e-12

        # Test with realistic Hall thruster conditions
        xe_fractions = [0.95, 0.05]  # Xe, Xe+
        xe_weights = [131.3, 131.3]  # g/mol
        xe_density = 1e16  # 1/cm³

        xe_masses = terra.mole_fractions_to_mass_densities(
            xe_fractions, xe_weights, xe_density)
        xe_recovered = terra.mass_densities_to_mole_fractions(xe_masses, xe_weights)

        @test xe_recovered≈xe_fractions rtol=1e-12
    end
end
