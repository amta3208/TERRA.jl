@testset "Utility Functions" begin
    @testset "C API Accessors" begin
        test_case_path = TEST_CASE_PATH
        reset_and_init!(test_case_path)

        max_species = terra.get_max_number_of_species_wrapper()
        active_species = terra.get_number_of_active_species_wrapper()
        species_names = terra.get_species_names_wrapper()

        @test max_species isa Int32 && max_species > 0
        @test active_species isa Int32 && 0 < active_species <= max_species
        @test species_names isa Vector{String}
        @test length(species_names) == active_species
        @test all(!isempty, species_names)
    end
    @testset "Maximum Species Count" begin
        # Ensure library path is set
        terra.load_terra_library!()

        # Test error when library not loaded
        terra.close_terra_library()
        @test_throws ErrorException terra.get_max_number_of_species_wrapper()

        # Test error message content
        try
            terra.get_max_number_of_species_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library not loaded", e.msg)
            @test occursin("load_terra_library!", e.msg)
        end

        # Reload library for actual test
        terra.load_terra_library!()

        max_species = terra.get_max_number_of_species_wrapper()
        println("Max species: ", max_species)

        # Check return type
        @test max_species isa Int32

        # Check that it's positive and reasonable
        @test max_species > 0
        @test max_species <= 100  # Reasonable upper bound
    end

    @testset "Species Names" begin
        # Ensure library is loaded
        terra.load_terra_library!()

        # Test error when library not loaded
        terra.close_terra_library()
        @test_throws ErrorException terra.get_species_names_wrapper()

        # Test error message content
        try
            terra.get_species_names_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library not loaded", e.msg)
            @test occursin("load_terra_library!", e.msg)
        end

        # Reload library for actual test
        test_case_path = TEST_CASE_PATH
        reset_and_init!(test_case_path)

        # Check species composition and counts
        species_names = terra.get_species_names_wrapper()
        active_species = terra.get_number_of_active_species_wrapper()
        allowed_species = ["N", "N2", "N+", "N2+", "E-"]
        for name in species_names
            @test name in allowed_species
        end

        # Check return type
        @test species_names isa Vector{String}

        # Check that all names are non-empty strings
        @test all(length(name) > 0 for name in species_names)

        # Check that number of species is reasonable
        @test length(species_names) == active_species

        # Check that species names contain only valid characters
        isalnum(c) = isletter(c) || isdigit(c)
        for name in species_names
            @test all(c -> isascii(c) && (isalnum(c) || c in ['+', '-', '_']), name)
        end

        # Check that species names are unique
        @test length(species_names) == length(unique(species_names))
    end

    @testset "Electronic States Parameters" begin
        # Ensure library is loaded
        terra.load_terra_library!()

        # Test error when library not loaded
        terra.close_terra_library()
        @test_throws ErrorException terra.get_max_number_of_atomic_electronic_states_wrapper()
        @test_throws ErrorException terra.get_max_number_of_molecular_electronic_states_wrapper()

        # Test error message content
        try
            terra.get_max_number_of_atomic_electronic_states_wrapper()
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library not loaded", e.msg)
            @test occursin("load_terra_library!", e.msg)
        end

        # Reload library for actual test
        terra.load_terra_library!()

        # Test atomic electronic states
        max_atomic_states = terra.get_max_number_of_atomic_electronic_states_wrapper()
        println("Max atomic electronic states: ", max_atomic_states)

        @test max_atomic_states isa Int32
        @test max_atomic_states > 0
        @test max_atomic_states <= 50  # Reasonable upper bound

        # Test molecular electronic states
        max_molecular_states = terra.get_max_number_of_molecular_electronic_states_wrapper()
        println("Max molecular electronic states: ", max_molecular_states)

        @test max_molecular_states isa Int32
        @test max_molecular_states > 0
        @test max_molecular_states <= 50  # Reasonable upper bound

        # Test that values are reasonable relative to each other
        @test max_molecular_states <= max_atomic_states
    end
end
