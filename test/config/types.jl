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

@testset "AxialMarchingConfig" begin
    marching = terra.AxialMarchingConfig()
    @test marching.handoff_mode == :reinitialize
    @test marching.termination_mode == :final_time
    @test marching.tee_policy == :match_te
    @test marching.override_tt_K === nothing
    @test marching.override_tv_K === nothing

    custom = terra.AxialMarchingConfig(;
        handoff_mode = :full_state,
        termination_mode = :steady_state,
        tee_policy = :from_inlet,
        override_tt_K = 900.0,
        override_tv_K = 800.0)
    @test custom.handoff_mode == :full_state
    @test custom.termination_mode == :steady_state
    @test custom.tee_policy == :from_inlet
    @test custom.override_tt_K == 900.0
    @test custom.override_tv_K == 800.0

    @test_throws ArgumentError terra.AxialMarchingConfig(; handoff_mode = :bad_mode)
    @test_throws ArgumentError terra.AxialMarchingConfig(; termination_mode = :bad_mode)
    @test_throws ArgumentError terra.AxialMarchingConfig(; tee_policy = :bad_mode)
    @test_throws ArgumentError terra.AxialMarchingConfig(; override_tt_K = 0.0)
    @test_throws ArgumentError terra.AxialMarchingConfig(; override_tv_K = -10.0)
end

@testset "ReactorResult and ReactorFrame" begin
    frame1 = terra.ReactorFrame(;
        t = 0.0,
        species_densities = [1e-3, 1e-6, 1e-8],
        temperatures = (tt = 300.0, te = 10000.0, tv = 310.0),
        total_energy = 1e4,
    )
    frame2 = terra.ReactorFrame(;
        t = 1.0,
        species_densities = [2e-3, 2e-6, 2e-8],
        temperatures = (tt = 320.0, te = 10500.0, tv = 315.0),
        total_energy = 1.1e4,
    )

    reactor = terra.ReactorResult(;
        t = [0.0, 1.0],
        frames = [frame1, frame2],
        success = true,
        message = "ok",
    )
    @test reactor.success == true
    @test length(reactor.frames) == 2
    @test reactor.frames[1].t == 0.0
    @test reactor.frames[2].temperatures.te == 10500.0

    single = reactor[2]
    @test length(single.t) == 1
    @test length(single.frames) == 1
    @test single.frames[1].t == 1.0
end

@testset "ChainSimulationResult" begin
    reactor_cfg = terra.ReactorConfig(;
        composition = terra.ReactorComposition(;
            species = ["N", "N2", "E-"],
            mole_fractions = [0.1, 0.8, 0.1],
            total_number_density = 1e13),
        thermal = terra.ReactorThermalState(; Tt = 300.0, Tv = 310.0, Tee = 320.0, Te = 10000.0))

    reactor = terra.ReactorResult(;
        t = [0.0, 1.0],
        frames = [
            terra.ReactorFrame(;
                t = 0.0,
                species_densities = [1e-3, 1e-6, 1e-8],
                temperatures = (tt = 300.0, te = 10000.0, tv = 310.0),
                total_energy = 1e4,
            ),
            terra.ReactorFrame(;
                t = 1.0,
                species_densities = [1e-3, 1e-6, 1e-8],
                temperatures = (tt = 301.0, te = 10000.0, tv = 311.0),
                total_energy = 1.1e4,
            ),
        ],
        success = true,
        message = "ok",
    )

    cell = terra.ChainCellResult(;
        compact_cell_index = 1,
        source_cell_index = 13,
        z_m = 0.0,
        dx_m = 0.01,
        te_K = 10000.0,
        u_neutral_m_s = 100.0,
        u_ion_m_s = 1000.0,
        reactor = reactor,
        endpoint_reactor = reactor_cfg,
    )
    metadata = terra.ChainMetadata(;
        diagnostics = Dict{String, Any}("note" => "test"),
        compact_to_source_index = [13],
        retained_point_count = 1,
    )

    chain = terra.ChainSimulationResult(;
        cells = [cell],
        metadata = metadata,
        success = true,
        failed_cell = nothing,
        message = "done",
    )

    @test chain.success == true
    @test chain.failed_cell === nothing
    @test length(chain.cells) == 1
    @test chain.cells[1].endpoint_reactor == reactor_cfg
    @test chain.cells[1].reactor == reactor
end

@testset "Nested Chain Cell Indexing" begin
    reactor = terra.ReactorResult(;
        t = [0.0, 1.0],
        frames = [
            terra.ReactorFrame(;
                t = 0.0,
                species_densities = [1e-3, 1e-6, 1e-8],
                temperatures = (tt = 300.0, te = 10000.0, tv = 310.0),
                total_energy = 1e4,
            ),
            terra.ReactorFrame(;
                t = 1.0,
                species_densities = [1e-3, 1e-6, 1e-8],
                temperatures = (tt = 301.0, te = 10000.0, tv = 311.0),
                total_energy = 1.1e4,
            ),
        ],
        success = true,
        message = "ok",
    )

    cell1 = terra.ChainCellResult(;
        compact_cell_index = 1,
        source_cell_index = 13,
        z_m = 0.0,
        dx_m = 0.01,
        te_K = 10000.0,
        u_neutral_m_s = 100.0,
        u_ion_m_s = 1000.0,
        reactor = reactor
    )
    cell2 = terra.ChainCellResult(;
        compact_cell_index = 2,
        source_cell_index = 14,
        z_m = 0.01,
        dx_m = 0.01,
        te_K = 9500.0,
        u_neutral_m_s = 100.0,
        u_ion_m_s = 1200.0,
        reactor = reactor
    )

    metadata = terra.ChainMetadata(;
        compact_to_source_index = [13, 14],
        original_point_count = 102,
        retained_point_count = 2,
        diagnostics = Dict{String, Any}("trimmed_point_count" => 12)
    )

    chain = terra.ChainSimulationResult(;
        cells = [cell1, cell2],
        metadata = metadata,
        success = true,
        failed_cell = nothing,
        message = "ok"
    )

    @test length(chain.cells) == 2
    @test chain.cells[1].source_cell_index == 13
    @test chain.cells[2].source_cell_index == 14
    @test chain.cells[1].reactor.frames[1].t == 0.0

    # Chain slicing is cell-based
    second = chain[2]
    @test length(second.cells) == 1
    @test second.cells[1].compact_cell_index == 2
    @test second.cells[1].source_cell_index == 14
end
