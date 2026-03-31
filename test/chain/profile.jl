 @testset "ChainWallProfile" begin
    wall_profile = terra.ChainWallProfile(;
                                          a_wall_over_v_m_inv = [100.0, 101.0],
                                          channel_gap_m = [0.02, 0.02],)
    @test wall_profile.channel_gap_m == [0.02, 0.02]

    inlet_composition = terra.ChainProfileInletComposition(;
                                                           species = ["N",
                                                               "N2",
                                                               "N+",
                                                               "N2+",
                                                               "E-"],
                                                           mole_fractions = [0.2,
                                                               0.65,
                                                               0.05,
                                                               0.05,
                                                               0.05],
                                                           total_number_density_m3 = 1e19,)
    inlet = terra.ChainProfileInlet(;
                                    composition = inlet_composition,
                                    thermal = terra.ReactorThermalState(; Tt = 500.0,
                                                                        Tv = 500.0,
                                                                        Tee = 500.0,
                                                                        Te = 20000.0),
                                    source_compact_index = 1,)
    profile = terra.AxialChainProfile(;
                                      z_m = [0.0, 0.01],
                                      dx_m = [0.01, 0.01],
                                      te_K = [20000.0, 21000.0],
                                      species_u_m_s = Dict("N" => [200.0, 210.0],
                                                           "N2" => [180.0, 181.0],
                                                           "N+" => [1000.0, 1100.0],
                                                           "N2+" => [900.0, 950.0]),
                                      wall_profile = wall_profile,
                                      inlet = inlet,
                                      selection = Dict("trim_start_index" => 5,
                                                       "original_point_count" => 7),)
    @test profile.schema_version == "terra_chain_profile_v4"
    @test profile.wall_profile !== nothing
    @test terra._resolve_compact_to_source_index(profile) == [5, 6]
    @test terra._resolve_original_point_count(profile) == 7

    inlet_reactor_si = terra._profile_inlet_reactor(profile, :SI)
    @test inlet_reactor_si.composition.total_number_density == 1e19
    inlet_reactor_cgs = terra._profile_inlet_reactor(profile, :CGS)
    @test inlet_reactor_cgs.composition.total_number_density ≈
          terra.convert_number_density_si_to_cgs(1e19)

    wall_values = terra._segment_wall_profile_values(profile, 2)
    @test wall_values.a_wall_over_v_m_inv == 101.0
    @test wall_values.channel_gap_m == 0.02

    @test_throws ArgumentError terra.ChainWallProfile(;
                                                      a_wall_over_v_m_inv = [100.0, -1.0])
    @test_throws ArgumentError terra.ChainWallProfile(;
                                                      a_wall_over_v_m_inv = [100.0, 101.0],
                                                      channel_gap_m = [0.02])
end
