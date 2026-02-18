@testset "generate_input_files" begin
    @testset "Directory Structure Creation" begin
        # Create a temporary directory for testing
        temp_dir = mktempdir()
        try
            config = terra.nitrogen_10ev_config()

            # Test successful file generation
            @test terra.generate_input_files(config, temp_dir) == true

            # Check that required directories were created
            @test isdir(joinpath(temp_dir, "input"))
            @test isdir(joinpath(temp_dir, "output"))
            @test isdir(joinpath(temp_dir, "output", "sources"))
            @test isdir(joinpath(temp_dir, "output", "states"))

            # Check that input files were created
            @test isfile(joinpath(temp_dir, "input", "prob_setup.inp"))
            @test isfile(joinpath(temp_dir, "input", "sources_setup.inp"))
            @test isfile(joinpath(temp_dir, "input", "tau_scaling.inp"))

        finally
            rm(temp_dir; recursive = true)
        end
    end

    @testset "Nested Config Input Generation" begin
        temp_dir = mktempdir()
        try
            legacy = terra.nitrogen_10ev_config()
            nested = terra.to_config(legacy)

            @test terra.generate_input_files(nested, temp_dir) == true
            @test isfile(joinpath(temp_dir, "input", "prob_setup.inp"))
            @test isfile(joinpath(temp_dir, "input", "sources_setup.inp"))
            @test isfile(joinpath(temp_dir, "input", "tau_scaling.inp"))

            prob_setup_content = read(joinpath(temp_dir, "input", "prob_setup.inp"), String)
            @test occursin("NSP=5", prob_setup_content)
            @test occursin("ND=0", prob_setup_content)
            @test occursin("DT=5.0e-6", prob_setup_content)
        finally
            rm(temp_dir; recursive = true)
        end
    end

    @testset "File Content Validation" begin
        temp_dir = mktempdir()
        try
            config = terra.nitrogen_10ev_config()
            terra.generate_input_files(config, temp_dir)

            # Read and check prob_setup.inp content
            prob_setup_content = read(
                joinpath(temp_dir, "input", "prob_setup.inp"), String)
            @test occursin("NSP=5", prob_setup_content)  # 5 species

            @test occursin("X1=1.0e-20", prob_setup_content)
            @test occursin("X2=0.9998", prob_setup_content)
            @test occursin("X3=1.0e-20", prob_setup_content)
            @test occursin("X4=0.0001", prob_setup_content)
            @test occursin("X5=0.0001", prob_setup_content)

            @test occursin("TOTAL_NUMBER_DENSITY=1.0e13", prob_setup_content)

            @test occursin("TT=750.0", prob_setup_content)
            @test occursin("TV=750.0", prob_setup_content)
            @test occursin("TEE=750.0", prob_setup_content)
            @test occursin("TE=115000.0", prob_setup_content)

            @test occursin("RAD_LEN=1.0", prob_setup_content)

            @test occursin("BBHMODEL=4", prob_setup_content)
            @test occursin("ESC_MODEL=1", prob_setup_content)
            @test occursin("AR_ET_MODEL=1", prob_setup_content)
            @test occursin("EEX_NONEQ=1", prob_setup_content)
            @test occursin("EV_RELAX_SET=1", prob_setup_content)
            @test occursin("ET_RELAX_SET=1", prob_setup_content)

            @test occursin("CONSIDER_ELEC_BBE=1", prob_setup_content)
            @test occursin("CONSIDER_ELEC_BFE=1", prob_setup_content)
            @test occursin("CONSIDER_ELEC_BBH=1", prob_setup_content)
            @test occursin("CONSIDER_ELEC_BFH=1", prob_setup_content)
            @test occursin("CONSIDER_RAD=0", prob_setup_content)
            @test occursin("CONSIDER_RDR=0", prob_setup_content)
            @test occursin("CONSIDER_CHEM=1", prob_setup_content)

            @test occursin("TIME_METHOD=2", prob_setup_content)
            @test occursin("IS_ISOTHERMAL_TEEX=0", prob_setup_content)
            @test occursin("ND=0", prob_setup_content)
            @test occursin("DT=5.0e-6", prob_setup_content)
            @test occursin("DTM=5.0", prob_setup_content)
            @test occursin("TLIM=1000.0", prob_setup_content)
            @test occursin("NSTEP=500000", prob_setup_content)

            # Read and check sources_setup.inp content
            sources_content = read(
                joinpath(temp_dir, "input", "sources_setup.inp"), String)
            @test occursin("BEGIN SPECIES SOURCES", sources_content)
            @test occursin("END SPECIES SOURCES", sources_content)
            @test occursin("BEGIN EXCITED STATE SOURCES", sources_content)
            @test occursin("END EXCITED STATE SOURCES", sources_content)
            for species in config.species
                @test occursin(species, sources_content)
            end

        finally
            rm(temp_dir; recursive = true)
        end
    end
end
