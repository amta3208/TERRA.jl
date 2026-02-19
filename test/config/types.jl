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

@testset "ODESolverConfig" begin
    solver = terra.ODESolverConfig(;
        reltol = 1e-8,
        abstol_density = 1e-10,
        saveat_count = 50,
        ramp_understep_ratio = inv(64),
        ramp_history_steps = 4)
    @test solver.reltol == 1e-8
    @test solver.abstol_density == 1e-10
    @test solver.saveat_count == 50
    @test solver.ramp_understep_ratio == inv(64)
    @test solver.ramp_history_steps == 4

    @test_throws ArgumentError terra.ODESolverConfig(; reltol = 0.0)
    @test_throws ArgumentError terra.ODESolverConfig(; saveat_count = 0)
end

@testset "SpaceConfig" begin
    space = terra.SpaceConfig(; nd = 0, dr = nothing)
    @test space.nd == 0
    @test space.dr === nothing

    space2 = terra.SpaceConfig(; nd = 1, dr = 0.1)
    @test space2.nd == 1
    @test space2.dr == 0.1

    @test_throws ArgumentError terra.SpaceConfig(; nd = -1)
    @test_throws ArgumentError terra.SpaceConfig(; nd = 1, dr = 0.0)
end

@testset "PhysicsConfig" begin
    physics = terra.PhysicsConfig()
    @test physics.bbh_model == 4
    @test physics.esc_model == 1
    @test physics.ar_et_model == 1

    custom = terra.PhysicsConfig(
        bbh_model = 2,
        esc_model = 0,
        ar_et_model = 2,
        eex_noneq = 0,
        ev_relax_set = 2,
        et_relax_set = 2
    )
    @test custom.bbh_model == 2
    @test custom.esc_model == 0
    @test custom.ar_et_model == 2
    @test custom.eex_noneq == 0
    @test custom.ev_relax_set == 2
    @test custom.et_relax_set == 2
end

@testset "ProcessConfig" begin
    processes = terra.ProcessConfig()
    @test processes.consider_elec_bbe == 1
    @test processes.consider_rad == 0

    custom = terra.ProcessConfig(
        consider_elec_bbe = 0,
        consider_elec_bfe = 0,
        consider_elec_bbh = 0,
        consider_elec_bfh = 0,
        consider_rad = 1,
        consider_rdr = 1,
        consider_chem = 0
    )
    @test custom.consider_elec_bbe == 0
    @test custom.consider_elec_bfe == 0
    @test custom.consider_elec_bbh == 0
    @test custom.consider_elec_bfh == 0
    @test custom.consider_rad == 1
    @test custom.consider_rdr == 1
    @test custom.consider_chem == 0
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

@testset "Config Modifiers" begin
    config = terra.Config(;
        reactor = terra.ReactorConfig(;
            composition = terra.ReactorComposition(;
                species = ["N", "N2", "E-"],
                mole_fractions = [0.1, 0.8, 0.1],
                total_number_density = 1e13),
            thermal = terra.ReactorThermalState(;
                Tt = 300.0, Tv = 350.0, Tee = 360.0, Te = 10000.0)),
        numerics = terra.NumericsConfig(;
            time = terra.TimeConfig(;
                dt = 1e-6, dt_output = 1e-4, duration = 1e-3, nstep = 1234, method = 2),
            residence_time = nothing),
        runtime = terra.RuntimeConfig(;
            database_path = ".",
            case_path = pwd(),
            unit_system = :CGS,
            validate_species_against_terra = false,
            print_source_terms = false,
            write_native_outputs = false,
            print_integration_output = false))

    temp_case = mktempdir()
    config_case = terra.with_case_path(config, temp_case)
    @test config_case.runtime.case_path == temp_case
    @test config.runtime.case_path == pwd()

    config_time = terra.with_time(config; dt = 2e-6, duration = 2e-3, method = 1)
    @test config_time.numerics.time.dt == 2e-6
    @test config_time.numerics.time.duration == 2e-3
    @test config_time.numerics.time.dt_output == config.numerics.time.dt_output
    @test config_time.numerics.time.method == 1
    @test config.numerics.time.dt == 1e-6

    config_runtime = terra.with_runtime(config;
        unit_system = :SI,
        print_source_terms = true,
        write_native_outputs = true)
    @test config_runtime.runtime.unit_system == :SI
    @test config_runtime.runtime.print_source_terms == true
    @test config_runtime.runtime.write_native_outputs == true

    @test_throws ArgumentError terra.with_case_path(config, joinpath(temp_case, "missing"))
    @test_throws ArgumentError terra.with_time(config; method = 9)
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

    inlet_config = terra.nitrogen_10ev_config(; isothermal = false)

    rt_inlet_reactor = terra.ResidenceTimeConfig(;
        enabled = true,
        L = 1.0,
        U_neutral = 1.0,
        U_ion = 1.0,
        inlet_reactor = inlet_config.reactor)
    @test rt_inlet_reactor.inlet_reactor == inlet_config.reactor

    rt_inlet_config_alias = terra.ResidenceTimeConfig(;
        enabled = true,
        L = 1.0,
        U_neutral = 1.0,
        U_ion = 1.0,
        inlet_config = inlet_config)
    @test rt_inlet_config_alias.inlet_reactor == inlet_config.reactor

    @test_throws ArgumentError terra.ResidenceTimeConfig(;
        enabled = true,
        L = 1.0,
        U_neutral = 1.0,
        U_ion = 1.0,
        inlet_reactor = inlet_config.reactor,
        inlet_config = inlet_config)
end

@testset "SimulationResult" begin
    time = [0.0, 1.0, 2.0]
    species_densities = [1e-3 1e-3 1e-3; 1e-6 1e-6 1e-6]
    temperatures = (tt = [300.0, 310.0, 320.0], te = [10000.0, 11000.0, 12000.0])
    total_energy = [1e4, 1.1e4, 1.2e4]

    results = terra.SimulationResult(
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
