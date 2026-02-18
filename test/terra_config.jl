
@testset "TERRA Configuration Tests" begin
    # Mock infrastructure for fortran_wrapper functions during testing
    mock_library_loaded = Ref{Bool}(false)
    mock_max_species = Ref{Int32}(20)
    mock_max_atomic_states = Ref{Int32}(10)
    mock_max_molecular_states = Ref{Int32}(5)
    mock_species_list = Ref{Vector{String}}(["N", "N2", "N+", "N2+", "E-"])

    # Helper functions to control mock behavior
    function set_mock_library_loaded(loaded::Bool)
        mock_library_loaded[] = loaded
    end

    function set_mock_terra_limits(max_species::Int, max_atomic::Int, max_molecular::Int)
        mock_max_species[] = Int32(max_species)
        mock_max_atomic_states[] = Int32(max_atomic)
        mock_max_molecular_states[] = Int32(max_molecular)
    end

    function set_mock_species_list(species::Vector{String})
        mock_species_list[] = species
    end

    @testset "TemperatureConfig" begin
        @testset "Valid Construction" begin
            # Test basic construction (keywords)
            temp_config = terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 10000.0)
            @test temp_config.Tt == 300.0
            @test temp_config.Tv == 300.0
            @test temp_config.Te == 10000.0
            @test temp_config.Tee == 300.0
        end

        @testset "Invalid Construction" begin
            # Test negative temperatures
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = -100.0, Tv = 300.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 300.0, Tv = -100.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = -400.0, Te = 10000.0)
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = 400.0, Te = -10000.0)

            # Test zero temperatures
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 0.0, Tv = 300.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 300.0, Tv = 0.0, Tee = 400.0, Te = 10000.0)
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = 0.0, Te = 10000.0)
            @test_throws ErrorException terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = 400.0, Te = 0.0)
        end
    end

    @testset "TimeIntegrationConfig" begin
        @testset "Valid Construction" begin
            # Test basic construction
            time_config = terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1000, method = 2)
            @test time_config.dt == 1e-6
            @test time_config.dtm == 1e-4
            @test time_config.tlim == 1e-3
            @test time_config.nstep == 1000
            @test time_config.method == 2

            # Test with defaults
            time_config2 = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3)
            @test time_config2.dt == 1e-6
            @test time_config2.dtm == 1e-4
            @test time_config2.tlim == 1e-3
            @test time_config2.nstep == 500000  # Default
            @test time_config2.method == 2      # Default
        end

        @testset "Invalid Construction" begin
            # Test negative time parameters
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = -1e-6, dtm = 1e-4, tlim = 1e-3)
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = -1e-4, tlim = 1e-3)
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = -1e-3)

            # Test zero time parameters
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 0.0, dtm = 1e-4, tlim = 1e-3)
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 0.0, tlim = 1e-3)
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 0.0)

            # Test invalid nstep
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 0)
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = -100)

            # Test invalid method
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1000, method = 3)
            @test_throws ErrorException terra.TimeIntegrationConfig(;
                dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1000, method = -1)
        end
    end

    @testset "PhysicsConfig" begin
        @testset "Default Construction" begin
            physics = terra.PhysicsConfig()
            @test physics.bbh_model == 4
            @test physics.esc_model == 1
            @test physics.ar_et_model == 1
            @test physics.eex_noneq == 1
            @test physics.ev_relax_set == 1
            @test physics.et_relax_set == 1
        end

        @testset "Custom Construction" begin
            physics = terra.PhysicsConfig(
                bbh_model = 2,
                esc_model = 0,
                ar_et_model = 2,
                eex_noneq = 0,
                ev_relax_set = 2,
                et_relax_set = 2
            )
            @test physics.bbh_model == 2
            @test physics.esc_model == 0
            @test physics.ar_et_model == 2
            @test physics.eex_noneq == 0
            @test physics.ev_relax_set == 2
            @test physics.et_relax_set == 2
        end
    end

    @testset "ProcessConfig" begin
        @testset "Default Construction" begin
            processes = terra.ProcessConfig()
            @test processes.consider_elec_bbe == 1
            @test processes.consider_elec_bfe == 1
            @test processes.consider_elec_bbh == 1
            @test processes.consider_elec_bfh == 1
            @test processes.consider_rad == 0
            @test processes.consider_rdr == 0
            @test processes.consider_chem == 1
        end

        @testset "Custom Construction" begin
            processes = terra.ProcessConfig(
                consider_elec_bbe = 0,
                consider_elec_bfe = 0,
                consider_elec_bbh = 0,
                consider_elec_bfh = 0,
                consider_rad = 1,
                consider_rdr = 1,
                consider_chem = 0
            )
            @test processes.consider_elec_bbe == 0
            @test processes.consider_elec_bfe == 0
            @test processes.consider_elec_bbh == 0
            @test processes.consider_elec_bfh == 0
            @test processes.consider_rad == 1
            @test processes.consider_rdr == 1
            @test processes.consider_chem == 0
        end
    end

    @testset "TERRAConfig" begin
        @testset "Valid Construction" begin
            species = ["N", "N2", "N+", "N2+", "E-"]
            mole_fractions = [1e-20, 0.9998, 1e-20, 0.0001, 0.0001]
            total_number_density = 1e13
            temperatures = terra.TemperatureConfig(;
                Tt = 750.0, Tv = 750.0, Tee = 750.0, Te = 115000.0)
            time_params = terra.TimeIntegrationConfig(; dt = 0.5e-5, dtm = 5.0, tlim = 1e3)

            config = terra.TERRAConfig(
                species = species,
                mole_fractions = mole_fractions,
                total_number_density = total_number_density,
                temperatures = temperatures,
                time_params = time_params
            )

            @test config.species == species
            @test config.mole_fractions == mole_fractions
            @test config.total_number_density == total_number_density
            @test config.temperatures == temperatures
            @test config.time_params == time_params
            @test config.unit_system == :CGS  # Default
            @test config.case_path == pwd()   # Default
            @test config.validate_species_against_terra == false  # Default
        end

        @testset "Custom Construction with All Parameters" begin
            species = ["Ar", "Ar+", "E-"]
            mole_fractions = [0.9998, 0.0001, 0.0001]
            total_number_density = 1e12
            temperatures = terra.TemperatureConfig(;
                Tt = 300.0, Tv = 300.0, Tee = 300.0, Te = 50000.0)
            time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-2)
            physics = terra.PhysicsConfig(
                bbh_model = 2,
                radiation_length = 2.0,
                get_electron_density_by_charge_balance = false,
                min_sts_frac = 1e-25,
                is_isothermal_teex = false
            )
            processes = terra.ProcessConfig(consider_rad = 1)
            test_case_path = joinpath(@__DIR__, "test_case")
            test_database_path = joinpath(test_case_path,
                "database/n2/elec_sts_expanded_electron_fits_ground")

            config = terra.TERRAConfig(
                species = species,
                mole_fractions = mole_fractions,
                total_number_density = total_number_density,
                temperatures = temperatures,
                time_params = time_params,
                physics = physics,
                processes = processes,
                database_path = test_database_path,
                case_path = test_case_path,
                unit_system = :SI,
                validate_species_against_terra = true,
                print_source_terms = false,
                write_native_outputs = true
            )

            @test config.species == species
            @test config.mole_fractions == mole_fractions
            @test config.total_number_density == total_number_density
            @test config.temperatures == temperatures
            @test config.time_params == time_params
            @test config.physics == physics
            @test config.processes == processes
            @test config.database_path == test_database_path
            @test config.case_path == test_case_path
            @test config.unit_system == :SI
            @test config.validate_species_against_terra == true
            @test config.physics.radiation_length == 2.0
            @test config.print_source_terms == false
            @test config.write_native_outputs == true
            @test config.physics.get_electron_density_by_charge_balance == false
            @test config.physics.min_sts_frac == 1e-25
            @test config.physics.is_isothermal_teex == false
        end
    end

    @testset "Nested Config + Adapters" begin
        species = ["N", "N2", "E-"]
        mole_fractions = [0.1, 0.8, 0.1]
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 350.0, Tee = 360.0, Te = 10000.0)
        time_params = terra.TimeIntegrationConfig(;
            dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1234, method = 2)

        legacy = terra.TERRAConfig(
            species = species,
            mole_fractions = mole_fractions,
            total_number_density = 1e13,
            temperatures = temperatures,
            time_params = time_params,
            database_path = ".",
            case_path = pwd(),
            unit_system = :CGS,
            validate_species_against_terra = true,
            print_source_terms = false,
            write_native_outputs = true,
            print_integration_output = false
        )

        nested = terra.to_config(legacy)
        @test nested.reactor.composition.species == legacy.species
        @test nested.reactor.composition.mole_fractions == legacy.mole_fractions
        @test nested.reactor.composition.total_number_density == legacy.total_number_density
        @test nested.reactor.thermal.Tt == legacy.temperatures.Tt
        @test nested.reactor.thermal.Tv == legacy.temperatures.Tv
        @test nested.reactor.thermal.Tee == legacy.temperatures.Tee
        @test nested.reactor.thermal.Te == legacy.temperatures.Te
        @test nested.models.physics == legacy.physics
        @test nested.models.processes == legacy.processes
        @test nested.numerics.time.dt == legacy.time_params.dt
        @test nested.numerics.time.dt_output == legacy.time_params.dtm
        @test nested.numerics.time.duration == legacy.time_params.tlim
        @test nested.runtime.database_path == legacy.database_path
        @test nested.runtime.case_path == legacy.case_path
        @test nested.runtime.unit_system == legacy.unit_system
        @test nested.numerics.residence_time !== nothing
        @test nested.numerics.residence_time.enabled == false

        legacy_roundtrip = terra.to_legacy_config(nested)
        @test legacy_roundtrip.species == legacy.species
        @test legacy_roundtrip.mole_fractions == legacy.mole_fractions
        @test legacy_roundtrip.total_number_density == legacy.total_number_density
        @test legacy_roundtrip.temperatures == legacy.temperatures
        @test legacy_roundtrip.time_params == legacy.time_params
        @test legacy_roundtrip.physics == legacy.physics
        @test legacy_roundtrip.processes == legacy.processes
        @test legacy_roundtrip.database_path == legacy.database_path
        @test legacy_roundtrip.case_path == legacy.case_path
        @test legacy_roundtrip.unit_system == legacy.unit_system
        @test legacy_roundtrip.validate_species_against_terra ==
              legacy.validate_species_against_terra
        @test legacy_roundtrip.print_source_terms == legacy.print_source_terms
        @test legacy_roundtrip.write_native_outputs == legacy.write_native_outputs
        @test legacy_roundtrip.print_integration_output == legacy.print_integration_output

        rt = terra.ResidenceTimeConfig(; enabled = true, L = 2.0, U_neutral = 3.0, U_ion = 4.0)
        nested_with_rt = terra.to_config(legacy; residence_time = rt)
        @test nested_with_rt.numerics.residence_time === rt
    end

    @testset "ResidenceTimeConfig (Phase 1)" begin
        rt_default = terra.ResidenceTimeConfig(1.0, 1.0, 1.0)
        @test rt_default.enabled == true
        @test rt_default.L == 1.0
        @test rt_default.U_neutral == 1.0
        @test rt_default.U_ion == 1.0

        rt_disabled = terra.ResidenceTimeConfig(; enabled = false, L = 1.5, U_neutral = 2.0,
            U_ion = 2.5, U_energy = 3.0)
        @test rt_disabled.enabled == false
        @test rt_disabled.U_energy == 3.0

        legacy_inlet = terra.TERRAConfig(
            species = ["N", "N2", "E-"],
            mole_fractions = [0.2, 0.7, 0.1],
            total_number_density = 2e13,
            temperatures = terra.TemperatureConfig(;
                Tt = 800.0, Tv = 900.0, Tee = 1000.0, Te = 11000.0),
            time_params = terra.TimeIntegrationConfig(; dt = 1e-6, dtm = 1e-4, tlim = 1e-3),
            database_path = ".",
            case_path = pwd()
        )

        rt_legacy_inlet = terra.ResidenceTimeConfig(;
            enabled = true,
            L = 1.0,
            U_neutral = 1.0,
            U_ion = 1.0,
            inlet_config = legacy_inlet)
        @test rt_legacy_inlet.inlet_reactor !== nothing
        @test rt_legacy_inlet.inlet_reactor.composition.species == legacy_inlet.species
        @test rt_legacy_inlet.inlet_reactor.thermal.Te == legacy_inlet.temperatures.Te

        nested_inlet = terra.to_config(legacy_inlet)
        rt_nested_inlet = terra.ResidenceTimeConfig(;
            enabled = true,
            L = 1.0,
            U_neutral = 1.0,
            U_ion = 1.0,
            inlet_reactor = nested_inlet)
        @test rt_nested_inlet.inlet_reactor == nested_inlet.reactor

        @test_throws ArgumentError terra.ResidenceTimeConfig(;
            enabled = true,
            L = 1.0,
            U_neutral = 1.0,
            U_ion = 1.0,
            inlet_reactor = nested_inlet.reactor,
            inlet_config = legacy_inlet)
    end

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

    @testset "nitrogen_10ev_config" begin
        @testset "Default Configuration" begin
            config = terra.nitrogen_10ev_config()

            @test config.species == ["N", "N2", "N+", "N2+", "E-"]
            @test config.mole_fractions == [1.0e-20, 0.9998, 1.0e-20, 0.0001, 0.0001]
            @test config.total_number_density == 1.0e13
            @test config.temperatures.Tt == 750.0
            @test config.temperatures.Tv == 750.0
            @test config.temperatures.Te == 115000.0
            @test config.temperatures.Tee == 750.0
            # Stored in seconds (prob_setup writes microseconds)
            @test config.time_params.dt ≈ 5e-12
            @test config.time_params.dtm ≈ 5e-6
            @test config.time_params.tlim ≈ 1e-3
            @test config.time_params.nstep == 500000
            @test config.time_params.method == 2
        end
    end

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

    @testset "TERRAResults" begin
        @testset "Structure Creation" begin
            time = [0.0, 1.0, 2.0]
            species_densities = [1e-3 1e-3 1e-3; 1e-6 1e-6 1e-6]
            temperatures = (tt = [300.0, 310.0, 320.0], te = [10000.0, 11000.0, 12000.0])
            total_energy = [1e4, 1.1e4, 1.2e4]

            results = terra.TERRAResults(
                time, species_densities, temperatures, total_energy,
                nothing, true, "Success"
            )

            @test results.time == time
            @test results.species_densities == species_densities
            @test results.temperatures == temperatures
            @test results.total_energy == total_energy
            @test results.source_terms === nothing
            @test results.success == true
            @test results.message == "Success"
        end
    end
end
