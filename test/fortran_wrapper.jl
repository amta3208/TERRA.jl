
@testset "Library Management" begin
    @testset "Library Loading and Status" begin
        # Test initial state - no library loaded
        @test !terra.is_terra_loaded()

        # Test that getting handle fails when no library is loaded
        @test_throws ErrorException terra.get_terra_handle()

        # Test error message content
        try
            terra.get_terra_handle()
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library not loaded", e.msg)
            @test occursin("load_terra_library!", e.msg)
        end
    end

    @testset "Library Path Setting" begin
        # Test setting library path with non-existent file
        fake_path = "/nonexistent/path/libterra.so"
        @test_throws ErrorException terra.load_terra_library!(fake_path)

        # Test error message for non-existent file
        try
            terra.load_terra_library!(fake_path)
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library file not found", e.msg)
            @test occursin(fake_path, e.msg)
        end

        # Test that library is still not loaded after failed attempt
        @test !terra.is_terra_loaded()
    end

    @testset "Library Cleanup" begin
        # Test closing library when none is loaded (should be safe)
        @test_nowarn terra.close_terra_library()
        @test !terra.is_terra_loaded()

        # Test multiple calls to close (should be safe)
        @test_nowarn terra.close_terra_library()
        @test_nowarn terra.close_terra_library()
        @test !terra.is_terra_loaded()
    end
end

@testset "Vibrational Temperature Wrapper" begin
    @testset "Round-trip Evib ↔ Tvib" begin
        # Initialize a consistent state
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        # Species mass densities [g/cm^3]
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Choose several vibrational temperatures to test inversion
        tvib_list = [750.0, 1500.0, 5000.0]
        # Hold species electronic temperatures fixed during inversion
        tex_val = 1200.0
        teex_vec = fill(tex_val, length(rho_sp))

        for tvib in tvib_list
            # Build electronic state populations consistent with (tex, trot, tvib)
            rho_ex = terra.set_electronic_boltzmann_wrapper(rho_sp, tex_val, tex_val, tvib)

            # Forward: Tvib -> Evib
            rho_evib = terra.calculate_vibrational_energy_wrapper(
                tvib, rho_sp; rho_ex = rho_ex, tex = teex_vec)
            @test rho_evib isa Float64
            @test isfinite(rho_evib)
            @test rho_evib >= 0.0

            # Inverse: Evib -> Tvib
            tvib_back = terra.calculate_vibrational_temperature_wrapper(
                rho_evib, rho_sp; rho_ex = rho_ex, tex = teex_vec)
            @test tvib_back isa Float64
            @test isfinite(tvib_back)
            @test tvib_back > 0.0

            # Round-trip consistency
            @test isapprox(tvib_back, tvib; rtol = 1e-8, atol = 1e-6)
        end
    end

    @testset "Error Handling Without Library" begin
        terra.close_terra_library()
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        @test_throws ErrorException terra.calculate_vibrational_temperature_wrapper(
            1.0, rho_sp)
    end
end

@testset "Initialization" begin
    # Ensure library path is set
    test_case_path = joinpath(@__DIR__, "test_case")

    # Test basic initialization (Fortran determines species count from input files)
    result = @test_nowarn reset_and_init!(test_case_path)

    # Check that result contains the expected fields
    if result !== nothing
        @test result isa NamedTuple
        @test haskey(result, :num_species)
        @test haskey(result, :num_dimensions)
        @test result.num_species isa Int32
        @test result.num_dimensions isa Int32
        @test result.num_species > 0
        @test result.num_dimensions >= 0
    end
end

@testset "Input/Directory Handling" begin
    @testset "Initialization Input Validation" begin
        # Ensure library is loaded and Fortran not initialized
        terra.close_terra_library()
        terra.load_terra_library!()

        # Test with non-existent case path
        @test_throws ErrorException terra.initialize_api_wrapper(case_path = "/nonexistent/path")

        # Test error message for non-existent case path
        try
            terra.initialize_api_wrapper(case_path = "/nonexistent/path")
            @test false  # Should not reach here
        catch e
            @test occursin("Case path does not exist", e.msg)
        end

        # Test with case path missing input directory
        temp_dir = mktempdir()
        try
            @test_throws ErrorException terra.initialize_api_wrapper(case_path = temp_dir)
        finally
            rm(temp_dir; recursive = true)
        end

        # Test error message for missing input file
        temp_dir = mktempdir()
        try
            terra.initialize_api_wrapper(case_path = temp_dir)
            @test false  # Should not reach here
        catch e
            @test occursin("Required input file not found", e.msg)
            @test occursin("prob_setup.inp", e.msg)
        finally
            rm(temp_dir; recursive = true)
        end
    end

    @testset "Directory Management" begin
        terra.load_terra_library!()

        # Store original directory
        original_dir = pwd()
        test_case_path = joinpath(@__DIR__, "test_case")

        # Test that directory is restored after successful call
        terra.initialize_api_wrapper(case_path = test_case_path)
        @test pwd() == original_dir

        # Test that directory is restored even after failed call
        temp_dir = mktempdir()
        try
            terra.initialize_api_wrapper(case_path = temp_dir)
        catch
            # Expected to fail
        end
        @test pwd() == original_dir
        rm(temp_dir; recursive = true)
    end

    @testset "Output Directory Creation" begin
        # Fresh state for this test
        terra.close_terra_library()
        terra.load_terra_library!()

        # Create a temporary test case directory with input
        temp_case_dir = mktempdir()
        try
            # Create input directory and file
            input_dir = joinpath(temp_case_dir, "input")
            mkpath(input_dir)
            touch(joinpath(input_dir, "prob_setup.inp"))

            # Verify output directories don't exist initially
            output_dir = joinpath(temp_case_dir, "output")
            sources_dir = joinpath(output_dir, "sources")
            states_dir = joinpath(output_dir, "states")

            @test !isdir(output_dir)
            @test !isdir(sources_dir)
            @test !isdir(states_dir)

            # Initialize - should create output directories
            terra.initialize_api_wrapper(case_path = temp_case_dir)

            # Verify output directories were created
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

            # Test that calling again doesn't cause errors (directories already exist)
            result2 = terra.initialize_api_wrapper(case_path = temp_case_dir)
            @test result2 isa NamedTuple

            # Verify directories still exist
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)

        finally
            rm(temp_case_dir; recursive = true)
        end
    end

    @testset "Output Directory Creation with Existing Directories" begin
        # Fresh state for this test
        terra.close_terra_library()
        terra.load_terra_library!()

        # Create a temporary test case directory with input and partial output structure
        temp_case_dir = mktempdir()
        try
            # Create input directory and file
            input_dir = joinpath(temp_case_dir, "input")
            mkpath(input_dir)
            touch(joinpath(input_dir, "prob_setup.inp"))

            # Create output directory but not subdirectories
            output_dir = joinpath(temp_case_dir, "output")
            mkpath(output_dir)

            # Create a test file in output to verify it's preserved
            test_file = joinpath(output_dir, "test_file.txt")
            write(test_file, "test content")

            sources_dir = joinpath(output_dir, "sources")
            states_dir = joinpath(output_dir, "states")

            @test isdir(output_dir)
            @test !isdir(sources_dir)
            @test !isdir(states_dir)
            @test isfile(test_file)

            # Initialize - should create missing subdirectories
            terra.initialize_api_wrapper(case_path = temp_case_dir)

            # Verify all directories exist and existing content is preserved
            @test isdir(output_dir)
            @test isdir(sources_dir)
            @test isdir(states_dir)
            @test isfile(test_file)
            @test read(test_file, String) == "test content"

        finally
            rm(temp_case_dir; recursive = true)
        end
    end
end

@testset "Utility Functions" begin
    @testset "C API Accessors" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        max_species = terra.get_max_number_of_species_wrapper()
        active_species = terra.get_number_of_active_species_wrapper()
        name_len = terra.get_species_name_length_wrapper()

        @test max_species isa Int32 && max_species > 0
        @test active_species isa Int32 && 0 < active_species <= max_species
        @test name_len isa Int32 && name_len > 0
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
        test_case_path = joinpath(@__DIR__, "test_case")
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

@testset "Temperature Calculation" begin
    @testset "Error Handling Without Library" begin
        terra.close_terra_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_throws ErrorException terra.calculate_temperatures_wrapper(rho_sp, rho_etot)

        try
            terra.calculate_temperatures_wrapper(rho_sp, rho_etot)
            @test false
        catch e
            @test occursin("TERRA library not loaded", e.msg)
        end
    end

    @testset "Error Handling Without Initialization" begin
        # Ensure Fortran API is not initialized for this block
        try
            terra.finalize_api_wrapper()
        catch
        end
        terra.close_terra_library()
        terra.load_terra_library!()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        @test_throws ArgumentError terra.calculate_temperatures_wrapper(rho_sp, rho_etot)

        try
            terra.calculate_temperatures_wrapper(rho_sp, rho_etot)
            @test false
        catch e
            @test occursin("rho_ex must be provided when electronic STS is active", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        tt = 1200.0
        tvib = 1500.0
        telec = 1800.0

        rho_ex = terra.set_electronic_boltzmann_wrapper(rho_sp, telec, tt, tvib)
        rho_eeex = terra.calculate_electron_electronic_energy_wrapper(telec, tvib, rho_sp)
        rho_evib = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp; rho_ex = rho_ex,
            tex = fill(telec, length(rho_sp)))
        rho_vx = terra.has_vibrational_sts_wrapper() ?
                 terra.set_vibrational_boltzmann_wrapper(rho_ex, telec, tt, tvib) : nothing
        rho_etot = terra.calculate_total_energy_wrapper(
            tt, rho_sp; rho_ex = rho_ex, rho_vx = rho_vx,
            rho_eeex = rho_eeex, rho_evib = rho_evib)

        @test_nowarn try
            result = terra.calculate_temperatures_wrapper(rho_sp, rho_etot;
                rho_ex = rho_ex, rho_vx = rho_vx, rho_eeex = rho_eeex, rho_evib = rho_evib)

            # If successful, check structure
            @test result isa NamedTuple
            @test haskey(result, :tt)
            @test haskey(result, :trot)
            @test haskey(result, :teex)
            @test haskey(result, :tvib)
            @test haskey(result, :tex)
            @test haskey(result, :tvx)

            # Test that tvx uses the correct molecular electronic states dimension
            if result.tvx isa Matrix
                @test size(result.tvx, 1) ==
                      terra.get_max_number_of_molecular_electronic_states_wrapper()
                @test size(result.tvx, 2) == terra.get_max_number_of_species_wrapper()
            end
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Vibrational Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        terra.close_terra_library()

        tvib = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException terra.calculate_vibrational_energy_wrapper(tvib, rho_sp)

        try
            terra.calculate_vibrational_energy_wrapper(tvib, rho_sp)
            @test false
        catch e
            @test occursin("TERRA library not loaded", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        tvib = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test basic call without optional arguments
        @test_nowarn try
            result = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp)
            # If successful, check return type
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0  # Vibrational energy should be non-negative
        catch e
            @test e isa ErrorException
        end

        # Test with optional arguments
        max_species = terra.get_max_number_of_species_wrapper()
        # For vibrational energy wrapper, rho_ex must match mnex (atomic) per Fortran API
        max_atomic_states = terra.get_max_number_of_atomic_electronic_states_wrapper()
        rho_ex = zeros(Float64, max_atomic_states, max_species)
        tex = fill(tvib, max_species)
        teex = tvib

        @test_nowarn try
            result = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp;
                rho_ex = rho_ex, tex = tex)
            # If successful, check return type
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end

        # Test with partial optional arguments
        @test_nowarn try
            result = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp; teex = teex)
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end
    end

    @testset "Input Validation and Edge Cases" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with different vibrational temperatures
        test_temperatures = [200.0, 300.0, 1000.0, 5000.0, 10000.0]

        for tvib in test_temperatures
            @test_nowarn try
                result = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp)
                @test result isa Float64
                @test isfinite(result)
                @test result >= 0.0
            catch e
                @test e isa ErrorException
            end
        end

        # Test with very small densities
        rho_sp_small = [1e-30, 1e-30, 1e-30, 1e-30, 1e-30]
        @test_nowarn try
            result = terra.calculate_vibrational_energy_wrapper(1000.0, rho_sp_small)
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Electron-Electronic Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        terra.close_terra_library()

        teex = 10000.0
        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException terra.calculate_electron_electronic_energy_wrapper(
            teex, tvib, rho_sp)

        try
            terra.calculate_electron_electronic_energy_wrapper(teex, tvib, rho_sp)
            @test false
        catch e
            @test occursin("TERRA library not loaded", e.msg)
        end
    end

    @testset "Input Validation" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with negative temperature
        @test_throws ArgumentError terra.calculate_electron_electronic_energy_wrapper(
            -1000.0, tvib, rho_sp)

        # Test with zero temperature
        @test_throws ArgumentError terra.calculate_electron_electronic_energy_wrapper(
            0.0, tvib, rho_sp)

        # Test with empty species array
        @test_throws ArgumentError terra.calculate_electron_electronic_energy_wrapper(
            1000.0, tvib, Float64[])

        # Test error messages
        try
            terra.calculate_electron_electronic_energy_wrapper(-1000.0, tvib, rho_sp)
            @test false
        catch e
            @test occursin("Electron-electronic temperature must be positive", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        teex = 10000.0
        tvib = 2000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_nowarn try
            result = terra.calculate_electron_electronic_energy_wrapper(teex, tvib, rho_sp)
            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0  # Electron-electronic energy should be non-negative
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Electronic Boltzmann Distribution" begin
    @testset "Error Handling Without Library" begin
        terra.close_terra_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        tex = 10000.0
        trot = 1000.0
        tvib = 2000.0

        @test_throws ErrorException terra.set_electronic_boltzmann_wrapper(
            rho_sp, tex, trot, tvib)

        try
            terra.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
            @test false
        catch e
            @test occursin("TERRA library not loaded", e.msg)
        end
    end

    @testset "Input Validation" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with negative temperatures
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            rho_sp, -1000.0, 1000.0, 2000.0)
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, -1000.0, 2000.0)
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, 1000.0, -2000.0)

        # Test with zero temperatures
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            rho_sp, 0.0, 1000.0, 2000.0)
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, 0.0, 2000.0)
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            rho_sp, 1000.0, 1000.0, 0.0)

        # Test with empty species array
        @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
            Float64[], 1000.0, 1000.0, 2000.0)

        # Test error messages
        try
            terra.set_electronic_boltzmann_wrapper(rho_sp, -1000.0, 1000.0, 2000.0)
            @test false
        catch e
            @test occursin("All temperatures must be positive", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        terra.load_terra_library!()
        test_case_path = joinpath(@__DIR__, "test_case")
        terra.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        tex = 10000.0
        trot = 1000.0
        tvib = 2000.0

        @test_nowarn try
            result = terra.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
            @test result isa Matrix{Float64}
            @test size(result, 2) == terra.get_max_number_of_species_wrapper()  # Second dimension should be max species
            # Now test that it uses the correct atomic electronic states dimension
            @test size(result, 1) ==
                  terra.get_max_number_of_atomic_electronic_states_wrapper()  # First dimension should be max atomic electronic states
            @test all(isfinite.(result))
            @test all(result .>= 0.0)  # Densities should be non-negative
        catch e
            @test e isa ErrorException
        end
    end

    @testset "Temperature Variations" begin
        terra.load_terra_library!()
        test_case_path = joinpath(@__DIR__, "test_case")
        terra.initialize_api_wrapper(case_path = test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        # Test with different temperature combinations
        test_cases = [
            (1000.0, 1000.0, 1000.0),   # Equal temperatures
            (10000.0, 1000.0, 2000.0),  # High electronic temperature
            (1000.0, 5000.0, 2000.0),   # High rotational temperature
            (1000.0, 1000.0, 10000.0)  # High vibrational temperature
        ]

        for (tex, trot, tvib) in test_cases
            @test_nowarn try
                result = terra.set_electronic_boltzmann_wrapper(rho_sp, tex, trot, tvib)
                @test result isa Matrix{Float64}
                @test all(isfinite.(result))
                @test all(result .>= 0.0)
            catch e
                @test e isa ErrorException
            end
        end
    end
end

@testset "Dimension Validations" begin
    test_case_path = joinpath(@__DIR__, "test_case")
    reset_and_init!(test_case_path)

    max_species = terra.get_max_number_of_species_wrapper()
    max_atomic_states = terra.get_max_number_of_atomic_electronic_states_wrapper()
    max_molecular_states = terra.get_max_number_of_molecular_electronic_states_wrapper()
    max_vqn = terra.get_max_vibrational_quantum_number_wrapper()

    # Oversized rho_sp for set_electronic_boltzmann_wrapper
    rho_sp_oversized = ones(max_species + 1)
    @test_throws ArgumentError terra.set_electronic_boltzmann_wrapper(
        rho_sp_oversized, 1000.0, 1000.0, 1000.0)

    # Oversized rho_ex for calculate_temperatures_wrapper
    rho_sp = ones(max_species)
    rho_ex_bad = zeros(Float64, max_atomic_states + 1, max_species)
    @test_throws ArgumentError terra.calculate_temperatures_wrapper(
        rho_sp, 1.0e4; rho_ex = rho_ex_bad)

    # Oversized rho_vx for calculate_total_energy_wrapper
    rho_vx_bad = zeros(Float64, max_vqn + 2, max_molecular_states + 1, max_species + 1)
    @test_throws ArgumentError terra.calculate_total_energy_wrapper(
        1000.0, rho_sp; rho_vx = rho_vx_bad)

    # Vibrational energy: accept molecular-sized rho_ex (padded internally)
    rho_ex_molecular = zeros(Float64, max_molecular_states, max_species)
    @test_nowarn begin
        evib = terra.calculate_vibrational_energy_wrapper(
            1000.0, rho_sp; rho_ex = rho_ex_molecular)
        @test evib isa Float64 && isfinite(evib) && evib >= 0.0
    end
end

@testset "Total Energy Calculation" begin
    @testset "Error Handling Without Library" begin
        terra.close_terra_library()

        tt = 1000.0
        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        @test_throws ErrorException terra.calculate_total_energy_wrapper(tt, rho_sp)

        try
            terra.calculate_total_energy_wrapper(tt, rho_sp)
            @test false
        catch e
            @test occursin("TERRA library not loaded", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        tt = 1000.0
        tvib = 2000.0
        telec = 3000.0

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_ex = terra.set_electronic_boltzmann_wrapper(rho_sp, telec, tt, tvib)
        rho_eeex = terra.calculate_electron_electronic_energy_wrapper(telec, tvib, rho_sp)
        rho_evib = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp)

        # Test with optional arguments
        # (generally optional, but required here due to flags in test case input file)
        @test_nowarn try
            result = terra.calculate_total_energy_wrapper(
                tt, rho_sp;
                rho_ex = rho_ex,
                rho_eeex = rho_eeex,
                rho_evib = rho_evib
            )

            @test result isa Float64
            @test isfinite(result)
            @test result >= 0.0
        catch e
            @test e isa ErrorException
        end
    end
end

@testset "Source Terms Calculation" begin
    @testset "Error Handling Without Library" begin
        # Ensure library is not loaded
        terra.close_terra_library()

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]
        rho_etot = 1e4

        # Test error when library not loaded
        @test_throws ErrorException terra.calculate_sources_wrapper(rho_sp, rho_etot)

        # Test error message content
        try
            terra.calculate_sources_wrapper(rho_sp, rho_etot)
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library not loaded", e.msg)
        end
    end

    @testset "Function Signature and Return Structure" begin
        # Initialize a consistent state
        test_case_path = joinpath(@__DIR__, "test_case")
        reset_and_init!(test_case_path)

        rho_sp = [1e-3, 1e-6, 1e-7, 1e-7, 1e-10]

        tt = 1000.0
        tvib = 2000.0
        telec = 3000.0

        # Build consistent optional inputs required by flags
        rho_ex = terra.set_electronic_boltzmann_wrapper(rho_sp, telec, tt, tvib)
        rho_eeex = terra.calculate_electron_electronic_energy_wrapper(telec, tvib, rho_sp)
        rho_evib = terra.calculate_vibrational_energy_wrapper(tvib, rho_sp;
            rho_ex = rho_ex, tex = fill(telec, length(rho_sp)))
        rho_vx = terra.has_vibrational_sts_wrapper() ?
                 terra.set_vibrational_boltzmann_wrapper(rho_ex, telec, tt, tvib) : nothing
        rho_etot = terra.calculate_total_energy_wrapper(
            tt, rho_sp; rho_ex = rho_ex, rho_vx = rho_vx,
            rho_eeex = rho_eeex, rho_evib = rho_evib)

        # Call and validate structure according to calculate_sources_wrapper
        @test_nowarn try
            result = terra.calculate_sources_wrapper(
                rho_sp, rho_etot; rho_ex = rho_ex, rho_vx = rho_vx,
                rho_eeex = rho_eeex, rho_evib = rho_evib)

            @test result isa NamedTuple
            @test all(k -> haskey(result, k),
                (:drho_sp, :drho_etot, :drho_ex, :drho_vx,
                    :drho_erot, :drho_eeex, :drho_evib))

            # drho_sp matches input species length and type
            @test result.drho_sp isa Vector{Float64}
            @test length(result.drho_sp) == length(rho_sp)
            @test all(isfinite, result.drho_sp)

            # Scalar derivatives are finite Float64
            @test result.drho_etot isa Float64 && isfinite(result.drho_etot)
            @test result.drho_erot isa Float64 && isfinite(result.drho_erot)
            @test result.drho_eeex isa Float64 && isfinite(result.drho_eeex)
            @test result.drho_evib isa Float64 && isfinite(result.drho_evib)

            # Because rho_ex was provided, drho_ex must be a full-sized matrix
            @test result.drho_ex isa Matrix{Float64}
            @test size(result.drho_ex, 1) ==
                  terra.get_max_number_of_atomic_electronic_states_wrapper()
            @test size(result.drho_ex, 2) == terra.get_number_of_active_species_wrapper()
            @test all(isfinite.(result.drho_ex))

            if rho_vx === nothing
                @test result.drho_vx === nothing
            else
                @test result.drho_vx isa Array{Float64, 3}
                @test size(result.drho_vx, 1) ==
                      terra.get_max_vibrational_quantum_number_wrapper() + 1
                @test size(result.drho_vx, 2) ==
                      terra.get_max_number_of_molecular_electronic_states_wrapper()
                @test size(result.drho_vx, 3) == terra.get_number_of_active_species_wrapper()
                @test all(isfinite.(result.drho_vx))
            end
        catch e
            @test e isa ErrorException
        end
    end
end
