@testset "ReactorResult and ReactorFrame" begin
    frame1 = terra.ReactorFrame(;
                                t = 0.0,
                                species_densities = [1e-3, 1e-6, 1e-8],
                                temperatures = (tt = 300.0, te = 10000.0, tv = 310.0),
                                total_energy = 1e4,)
    frame2 = terra.ReactorFrame(;
                                t = 1.0,
                                species_densities = [2e-3, 2e-6, 2e-8],
                                temperatures = (tt = 320.0, te = 10500.0, tv = 315.0),
                                total_energy = 1.1e4,)

    reactor = terra.ReactorResult(;
                                  t = [0.0, 1.0],
                                  frames = [frame1, frame2],
                                  success = true,
                                  message = "ok",)
    @test reactor.success == true
    @test length(reactor.frames) == 2
    @test reactor.frames[1].t == 0.0
    @test reactor.frames[2].temperatures.te == 10500.0
    @test terra.validate_results(reactor)

    densities = terra.species_density_matrix(reactor)
    temperatures = terra.temperature_history(reactor)
    total_energy = terra.total_energy_history(reactor)
    @test size(densities) == (3, 2)
    @test densities[:, 2] == [2e-3, 2e-6, 2e-8]
    @test temperatures.tt == [300.0, 320.0]
    @test temperatures.te == [10000.0, 10500.0]
    @test temperatures.tv == [310.0, 315.0]
    @test total_energy == [1e4, 1.1e4]

    single = reactor[2]
    @test length(single.t) == 1
    @test length(single.frames) == 1
    @test single.frames[1].t == 1.0
end
