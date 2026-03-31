 @testset "Initialization" begin
    # Ensure library path is set
    test_case_path = TEST_CASE_PATH

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

     @testset "Initialization Creates Missing Output Layout" begin
        temp_case_path = mktempdir()
        try
            input_dir = joinpath(temp_case_path, "input")
            mkpath(input_dir)

            for name in ("prob_setup.inp", "sources_setup.inp", "tau_scaling.inp")
                cp(joinpath(TEST_CASE_PATH, "input", name), joinpath(input_dir, name))
            end

            prob_setup_path = joinpath(input_dir, "prob_setup.inp")
            prob_setup_text = read(prob_setup_path, String)
            database_path = normpath(joinpath(TEST_CASE_PATH,
                                              "..", "..", "..", "..",
                                              "database", "n2",
                                              "elec_sts_expanded_electron_fits"))
            prob_setup_text = replace(prob_setup_text,
                                      r"(?m)^DATABASE_PATH=.*$" =>
                                      "DATABASE_PATH=$(database_path)")
            write(prob_setup_path, prob_setup_text)

            @test !isdir(joinpath(temp_case_path, "output", "states"))

            dims = @test_nowarn reset_and_init!(temp_case_path)
            try
                @test dims isa NamedTuple
                @test isdir(joinpath(temp_case_path, "output"))
                @test isdir(joinpath(temp_case_path, "output", "sources"))
                @test isdir(joinpath(temp_case_path, "output", "states"))
                @test isdir(joinpath(temp_case_path, "output", "logs"))
            finally
                try
                    terra.finalize_api_wrapper()
                catch
                end
            end
        finally
            rm(temp_case_path; recursive = true)
        end
    end
end

 @testset "Input And Lifecycle Handling" begin
     @testset "Initialization Input Validation" begin
        # Ensure library is loaded and Fortran not initialized
        reload_terra_library_for_tests!()

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
        reload_terra_library_for_tests!()

        # Store original directory
        original_dir = pwd()
        test_case_path = TEST_CASE_PATH

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

     @testset "Native Output Open And Close" begin
        dims = @test_nowarn reset_and_init!(TEST_CASE_PATH)
        try
            @test dims isa NamedTuple

            @test_nowarn terra.open_api_output_files_wrapper()
            @test terra.TERRA_OUTPUTS_OPEN[] == true
            @test_nowarn terra.close_api_output_files_wrapper()
            @test terra.TERRA_OUTPUTS_OPEN[] == false
        finally
            try
                terra.finalize_api_wrapper()
            catch
            end
        end
    end
end
