 @testset "ChainSimulationResult" begin
    reactor_cfg = terra.ReactorConfig(;
                                      composition = terra.ReactorComposition(;
                                                                             species = ["N",
                                                                                 "N2",
                                                                                 "E-"],
                                                                             mole_fractions = [0.1,
                                                                                 0.8,
                                                                                 0.1],
                                                                             total_number_density = 1e13),
                                      thermal = terra.ReactorThermalState(; Tt = 300.0,
                                                                          Tv = 310.0,
                                                                          Tee = 320.0,
                                                                          Te = 10000.0))

    reactor = terra.ReactorResult(;
                                  t = [0.0, 1.0],
                                  frames = [terra.ReactorFrame(;
                                                               t = 0.0,
                                                               species_densities = [1e-3,
                                                                   1e-6,
                                                                   1e-8],
                                                               temperatures = (tt = 300.0,
                                                                               te = 10000.0,
                                                                               tv = 310.0),
                                                               total_energy = 1e4,),
                                      terra.ReactorFrame(;
                                                         t = 1.0,
                                                         species_densities = [1e-3,
                                                             1e-6,
                                                             1e-8],
                                                         temperatures = (tt = 301.0,
                                                                         te = 10000.0,
                                                                         tv = 311.0),
                                                         total_energy = 1.1e4,)],
                                  success = true,
                                  message = "ok",)

    cell = terra.ChainCellResult(;
                                 compact_cell_index = 1,
                                 source_cell_index = 13,
                                 z_m = 0.0,
                                 dx_m = 0.01,
                                 te_K = 10000.0,
                                 species_u_m_s = Dict("N" => 100.0, "N2" => 120.0,
                                                      "N+" => 1000.0),
                                 reactor = reactor,
                                 endpoint_reactor = reactor_cfg,)
    metadata = terra.ChainMetadata(;
                                   diagnostics = Dict{String, Any}("note" => "test"),
                                   compact_to_source_index = [13],
                                   retained_point_count = 1,)

    chain = terra.ChainSimulationResult(;
                                        cells = [cell],
                                        metadata = metadata,
                                        success = true,
                                        failed_cell = nothing,
                                        message = "done",)

    @test chain.success == true
    @test chain.failed_cell === nothing
    @test length(chain.cells) == 1
    @test chain.cells[1].endpoint_reactor == reactor_cfg
    @test chain.cells[1].reactor == reactor
end

 @testset "Nested Chain Cell Indexing" begin
    reactor = terra.ReactorResult(;
                                  t = [0.0, 1.0],
                                  frames = [terra.ReactorFrame(;
                                                               t = 0.0,
                                                               species_densities = [1e-3,
                                                                   1e-6,
                                                                   1e-8],
                                                               temperatures = (tt = 300.0,
                                                                               te = 10000.0,
                                                                               tv = 310.0),
                                                               total_energy = 1e4,),
                                      terra.ReactorFrame(;
                                                         t = 1.0,
                                                         species_densities = [1e-3,
                                                             1e-6,
                                                             1e-8],
                                                         temperatures = (tt = 301.0,
                                                                         te = 10000.0,
                                                                         tv = 311.0),
                                                         total_energy = 1.1e4,)],
                                  success = true,
                                  message = "ok",)

    cell1 = terra.ChainCellResult(;
                                  compact_cell_index = 1,
                                  source_cell_index = 13,
                                  z_m = 0.0,
                                  dx_m = 0.01,
                                  te_K = 10000.0,
                                  species_u_m_s = Dict("N" => 100.0, "N2" => 120.0,
                                                       "N+" => 1000.0),
                                  reactor = reactor)
    cell2 = terra.ChainCellResult(;
                                  compact_cell_index = 2,
                                  source_cell_index = 14,
                                  z_m = 0.01,
                                  dx_m = 0.01,
                                  te_K = 9500.0,
                                  species_u_m_s = Dict("N" => 105.0, "N2" => 125.0,
                                                       "N+" => 1200.0),
                                  reactor = reactor)

    metadata = terra.ChainMetadata(;
                                   compact_to_source_index = [13, 14],
                                   original_point_count = 102,
                                   retained_point_count = 2,
                                   diagnostics = Dict{String, Any}("trimmed_point_count" => 12))

    chain = terra.ChainSimulationResult(;
                                        cells = [cell1, cell2],
                                        metadata = metadata,
                                        success = true,
                                        failed_cell = nothing,
                                        message = "ok")

    @test length(chain.cells) == 2
    @test chain.cells[1].source_cell_index == 13
    @test chain.cells[2].source_cell_index == 14
    @test chain.cells[1].reactor.frames[1].t == 0.0

    second = chain[2]
    @test length(second.cells) == 1
    @test second.cells[1].compact_cell_index == 2
    @test second.cells[1].source_cell_index == 14
end
