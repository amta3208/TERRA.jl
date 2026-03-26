@testset "AxialMarchingConfig" begin
    marching = terra.AxialMarchingConfig()
    @test marching.handoff_mode == :full_state
    @test marching.termination_mode == :final_time
    @test !(:tee_policy in fieldnames(typeof(marching)))
    @test marching.override_tt_K === nothing
    @test marching.override_tv_K === nothing
    @test marching.is_isothermal_teex == true

    custom = terra.AxialMarchingConfig(;
                                       handoff_mode = :full_state,
                                       termination_mode = :steady_state,
                                       override_tt_K = 900.0,
                                       override_tv_K = 800.0,
                                       is_isothermal_teex = false)
    @test custom.handoff_mode == :full_state
    @test custom.termination_mode == :steady_state
    @test custom.override_tt_K == 900.0
    @test custom.override_tv_K == 800.0
    @test custom.is_isothermal_teex == false

    @test_throws ArgumentError terra.AxialMarchingConfig(; handoff_mode = :bad_mode)
    @test_throws ArgumentError terra.AxialMarchingConfig(; termination_mode = :bad_mode)
    @test_throws MethodError terra.AxialMarchingConfig(; tee_policy = :from_inlet)
    @test_throws ArgumentError terra.AxialMarchingConfig(; override_tt_K = 0.0)
    @test_throws ArgumentError terra.AxialMarchingConfig(; override_tv_K = -10.0)
end
