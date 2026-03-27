@testset "AxialMarchingConfig" begin
    marching = terra.AxialMarchingConfig()
    @test marching.handoff_policy isa terra.FullStateHandoff
    @test marching.termination_policy isa terra.FinalTimeTermination
    @test !(:tee_policy in fieldnames(typeof(marching)))
    @test marching.override_tt_K === nothing
    @test marching.override_tv_K === nothing
    @test marching.is_isothermal_teex == true

    custom = terra.AxialMarchingConfig(;
                                       handoff_policy = terra.ReinitializeHandoff(),
                                       termination_policy = terra.SteadyStateTermination(),
                                       override_tt_K = 900.0,
                                       override_tv_K = 800.0,
                                       is_isothermal_teex = false)
    @test custom.handoff_policy isa terra.ReinitializeHandoff
    @test custom.termination_policy isa terra.SteadyStateTermination
    @test custom.override_tt_K == 900.0
    @test custom.override_tv_K == 800.0
    @test custom.is_isothermal_teex == false

    @test terra.chain_policy_name(marching.handoff_policy) == "full_state"
    @test terra.chain_policy_name(custom.termination_policy) == "steady_state"

    @test_throws MethodError terra.AxialMarchingConfig(; handoff_mode = :bad_mode)
    @test_throws MethodError terra.AxialMarchingConfig(; termination_mode = :bad_mode)
    @test_throws MethodError terra.AxialMarchingConfig(; tee_policy = :from_inlet)
    @test_throws ArgumentError terra.AxialMarchingConfig(; override_tt_K = 0.0)
    @test_throws ArgumentError terra.AxialMarchingConfig(; override_tv_K = -10.0)
end
