 @testset "PhysicsConfig" begin
    physics = terra.PhysicsConfig()
    @test physics.bbh_model == 4
    @test physics.esc_model == 1
    @test physics.ar_et_model == 1

    custom = terra.PhysicsConfig(bbh_model = 2,
                                 esc_model = 0,
                                 ar_et_model = 2,
                                 eex_noneq = 0,
                                 ev_relax_set = 2,
                                 et_relax_set = 2)
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

    custom = terra.ProcessConfig(consider_elec_bbe = 0,
                                 consider_elec_bfe = 0,
                                 consider_elec_bbh = 0,
                                 consider_elec_bfh = 0,
                                 consider_rad = 1,
                                 consider_rdr = 1,
                                 consider_chem = 0)
    @test custom.consider_elec_bbe == 0
    @test custom.consider_elec_bfe == 0
    @test custom.consider_elec_bbh == 0
    @test custom.consider_elec_bfh == 0
    @test custom.consider_rad == 1
    @test custom.consider_rdr == 1
    @test custom.consider_chem == 0
end
