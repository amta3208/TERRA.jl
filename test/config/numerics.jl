 @progress_testset "TimeConfig" begin
     @progress_testset "Valid Construction" begin
        time = terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4, duration = 1e-3,
                                nstep = 1000, method = 2)
        @test time.dt == 1e-6
        @test time.dt_output == 1e-4
        @test time.duration == 1e-3
        @test time.nstep == 1000
        @test time.method == 2
    end

     @progress_testset "Invalid Construction" begin
        @test_throws ArgumentError terra.TimeConfig(; dt = -1e-6, dt_output = 1e-4,
                                                    duration = 1e-3)
        @test_throws ArgumentError terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4,
                                                    duration = 0.0)
        @test_throws ArgumentError terra.TimeConfig(; dt = 1e-6, dt_output = 1e-4,
                                                    duration = 1e-3, method = 9)
    end
end

 @progress_testset "ODESolverConfig" begin
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

 @progress_testset "SpaceConfig" begin
    space = terra.SpaceConfig(; nd = 0, dr = nothing)
    @test space.nd == 0
    @test space.dr === nothing

    space2 = terra.SpaceConfig(; nd = 1, dr = 0.1)
    @test space2.nd == 1
    @test space2.dr == 0.1

    @test_throws ArgumentError terra.SpaceConfig(; nd = -1)
    @test_throws ArgumentError terra.SpaceConfig(; nd = 1, dr = 0.0)
end
