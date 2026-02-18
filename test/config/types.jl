@testset "ReactorComposition" begin
    @testset "Valid Construction" begin
        composition = terra.ReactorComposition(;
            species = ["N", "N2", "E-"],
            mole_fractions = [0.1, 0.8, 0.1],
            total_number_density = 1e13)
        @test composition.species == ["N", "N2", "E-"]
        @test composition.mole_fractions == [0.1, 0.8, 0.1]
        @test composition.total_number_density == 1e13
    end

    @testset "Invalid Construction" begin
        @test_throws ArgumentError terra.ReactorComposition(;
            species = ["N", "N2"], mole_fractions = [0.5], total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
            species = String[], mole_fractions = Float64[], total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
            species = ["N", "N2"], mole_fractions = [0.3, 0.3], total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
            species = ["N", "N"], mole_fractions = [0.5, 0.5], total_number_density = 1e13)
        @test_throws ArgumentError terra.ReactorComposition(;
            species = ["N", "N2"], mole_fractions = [0.5, 0.5], total_number_density = 0.0)
    end
end

@testset "ReactorThermalState" begin
    @testset "Valid Construction" begin
        thermal = terra.ReactorThermalState(; Tt = 300.0, Tv = 310.0, Tee = 320.0, Te = 10000.0)
        @test thermal.Tt == 300.0
        @test thermal.Tv == 310.0
        @test thermal.Tee == 320.0
        @test thermal.Te == 10000.0
    end

    @testset "Invalid Construction" begin
        @test_throws ArgumentError terra.ReactorThermalState(; Tt = 0.0, Tv = 310.0, Tee = 320.0, Te = 10000.0)
        @test_throws ArgumentError terra.ReactorThermalState(; Tt = 300.0, Tv = -1.0, Tee = 320.0, Te = 10000.0)
    end
end

@testset "TimeConfig" begin
    @testset "Valid Construction" begin
        time = terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4, duration = 1e-3, nstep = 1000, method = 2)
        @test time.dt == 1e-6
        @test time.dt_output == 1e-4
        @test time.duration == 1e-3
        @test time.nstep == 1000
        @test time.method == 2
    end

    @testset "Invalid Construction" begin
        @test_throws ArgumentError terra.TimeConfig(; dt = -1e-6, dt_output = 1e-4, duration = 1e-3)
        @test_throws ArgumentError terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4, duration = 0.0)
        @test_throws ArgumentError terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4, duration = 1e-3, method = 9)
    end
end

@testset "Config (Nested)" begin
    reactor = terra.ReactorConfig(;
        composition = terra.ReactorComposition(;
            species = ["N", "N2", "E-"],
            mole_fractions = [0.1, 0.8, 0.1],
            total_number_density = 1e13),
        thermal = terra.ReactorThermalState(; Tt = 300.0, Tv = 350.0, Tee = 360.0, Te = 10000.0))
    models = terra.ModelConfig()
    numerics = terra.NumericsConfig(;
        time = terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4, duration = 1e-3),
        residence_time = nothing)
    runtime = terra.RuntimeConfig(; case_path = pwd(), unit_system = :CGS)

    config = terra.Config(;
        reactor = reactor,
        models = models,
        numerics = numerics,
        runtime = runtime)

    @test config.reactor == reactor
    @test config.models == models
    @test config.numerics == numerics
    @test config.runtime == runtime
end

@testset "Legacy TemperatureConfig" begin
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

@testset "Legacy TimeIntegrationConfig" begin
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

@testset "Legacy TERRAConfig" begin
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
        test_case_path = TEST_CASE_PATH
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

@testset "Config Adapters" begin
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

@testset "Config Modifiers" begin
    legacy = terra.TERRAConfig(
        species = ["N", "N2", "E-"],
        mole_fractions = [0.1, 0.8, 0.1],
        total_number_density = 1e13,
        temperatures = terra.TemperatureConfig(;
            Tt = 300.0, Tv = 350.0, Tee = 360.0, Te = 10000.0),
        time_params = terra.TimeIntegrationConfig(;
            dt = 1e-6, dtm = 1e-4, tlim = 1e-3, nstep = 1234, method = 2),
        database_path = ".",
        case_path = pwd(),
        unit_system = :CGS,
        validate_species_against_terra = false,
        print_source_terms = false,
        write_native_outputs = false,
        print_integration_output = false
    )
    nested = terra.to_config(legacy)

    temp_case = mktempdir()
    nested_case = terra.with_case_path(nested, temp_case)
    @test nested_case.runtime.case_path == temp_case
    @test nested_case.reactor == nested.reactor
    @test nested.runtime.case_path == pwd()

    nested_time = terra.with_time(nested; dt = 2e-6, duration = 2e-3, method = 1)
    @test nested_time.numerics.time.dt == 2e-6
    @test nested_time.numerics.time.duration == 2e-3
    @test nested_time.numerics.time.dt_output == nested.numerics.time.dt_output
    @test nested_time.numerics.time.method == 1
    @test nested.numerics.time.dt == 1e-6

    nested_runtime = terra.with_runtime(nested;
        unit_system = :SI,
        print_source_terms = true,
        write_native_outputs = true)
    @test nested_runtime.runtime.unit_system == :SI
    @test nested_runtime.runtime.print_source_terms == true
    @test nested_runtime.runtime.write_native_outputs == true
    @test nested_runtime.runtime.case_path == nested.runtime.case_path

    legacy_case = terra.with_case_path(legacy, temp_case)
    @test legacy_case.case_path == temp_case
    @test legacy_case.time_params == legacy.time_params

    legacy_time = terra.with_time(legacy; dt = 3e-6, tlim = 4e-3, nstep = 999)
    @test legacy_time.time_params.dt == 3e-6
    @test legacy_time.time_params.tlim == 4e-3
    @test legacy_time.time_params.nstep == 999
    @test legacy_time.case_path == legacy.case_path

    legacy_runtime = terra.with_runtime(legacy;
        unit_system = :SI,
        print_source_terms = true,
        write_native_outputs = true)
    @test legacy_runtime.unit_system == :SI
    @test legacy_runtime.print_source_terms == true
    @test legacy_runtime.write_native_outputs == true

    @test_throws ArgumentError terra.with_case_path(nested, joinpath(temp_case, "missing"))
    @test_throws ArgumentError terra.with_time(nested; method = 9)
end

@testset "ResidenceTimeConfig" begin
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
