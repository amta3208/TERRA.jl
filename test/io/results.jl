@testset "Chain Result Save/Load" begin
    frame1 = terra.ReactorFrame(;
        t = 0.0,
        species_densities = [1e-3, 1e-6, 1e-8],
        temperatures = (tt = 300.0, te = 10000.0, tv = 320.0),
        total_energy = 1.0e4,
        source_terms = (production = 1.0, sink = 0.1),
        diagnostics = Dict{String, Any}("stage" => "initial"),
    )
    frame2 = terra.ReactorFrame(;
        t = 2.0e-6,
        species_densities = [1.1e-3, 1.1e-6, 1.0e-8],
        temperatures = (tt = 320.0, te = 9800.0, tv = 330.0),
        total_energy = 1.1e4,
        diagnostics = Dict{String, Any}("stage" => "final"),
    )
    reactor1 = terra.ReactorResult(;
        t = [0.0, 2.0e-6],
        frames = [frame1, frame2],
        success = true,
        message = "cell-1 ok",
        source_terms = (net = [1.0, 2.0],),
        metadata = Dict{String, Any}(
            "solver" => "ode",
            "step_stats" => Dict{String, Any}("accepted_steps" => 24),
        ),
    )

    endpoint = terra.ReactorConfig(;
        composition = terra.ReactorComposition(;
            species = ["N", "N2", "E-"],
            mole_fractions = [0.1, 0.8, 0.1],
            total_number_density = 1e13),
        thermal = terra.ReactorThermalState(;
            Tt = 310.0,
            Tv = 330.0,
            Tee = 9800.0,
            Te = 9800.0),
    )

    cell1 = terra.ChainCellResult(;
        compact_cell_index = 1,
        source_cell_index = 5,
        z_m = 0.0,
        dx_m = 0.01,
        te_K = 10000.0,
        species_u_m_s = Dict("N" => 150.0, "N2" => 180.0, "N+" => 15000.0),
        reactor = reactor1,
        endpoint_reactor = endpoint,
        success = true,
        message = "cell-1 done",
    )

    failed_frame = terra.ReactorFrame(;
        t = 0.0,
        species_densities = [9e-4, 9e-7, 2e-8],
        temperatures = (tt = 305.0, te = 9200.0, tv = 325.0),
        total_energy = 9.8e3,
        diagnostics = Dict{String, Any}("reason" => "mock failure"),
    )
    reactor2 = terra.ReactorResult(;
        t = [0.0],
        frames = [failed_frame],
        success = false,
        message = "cell-2 failed",
    )
    cell2 = terra.ChainCellResult(;
        compact_cell_index = 2,
        source_cell_index = 6,
        z_m = 0.01,
        dx_m = 0.01,
        te_K = 9200.0,
        species_u_m_s = Dict("N" => 140.0, "N2" => 170.0, "N+" => 14000.0),
        reactor = reactor2,
        endpoint_reactor = nothing,
        success = false,
        message = "cell-2 failed",
    )

    metadata = terra.ChainMetadata(;
        schema_version = "terra_chain_profile_v3",
        generator = Dict{String, Any}("tool" => "unit-test"),
        selection = Dict{String, Any}("trim_start_index" => 5),
        diagnostics = Dict{String, Any}(
            "handoff_mode" => "reinitialize",
            "simulated_total_number_density" => [1e13, 9.5e12],
        ),
        compact_to_source_index = [5, 6],
        original_point_count = 20,
        retained_point_count = 2,
    )

    chain = terra.ChainSimulationResult(;
        cells = [cell1, cell2],
        metadata = metadata,
        success = false,
        failed_cell = 2,
        message = "Chain integration failed at cell 2",
    )

    output_path = tempname() * ".json"
    @test terra.save_results(chain, output_path)
    @test isfile(output_path)

    loaded = terra.load_results_chain(output_path)

    @test loaded.success == false
    @test loaded.failed_cell == 2
    @test loaded.message == "Chain integration failed at cell 2"
    @test length(loaded.cells) == 2
    @test loaded.metadata.compact_to_source_index == [5, 6]
    @test loaded.metadata.original_point_count == 20
    @test loaded.cells[1].source_cell_index == 5
    @test loaded.cells[1].endpoint_reactor !== nothing
    @test loaded.cells[1].endpoint_reactor.composition.species == ["N", "N2", "E-"]
    @test loaded.cells[1].reactor.frames[2].temperatures.te ≈ 9800.0
    @test loaded.cells[1].species_u_m_s["N+"] == 15000.0
    @test loaded.cells[1].reactor.frames[1].source_terms !== nothing
    @test loaded.cells[1].reactor.frames[1].source_terms.production ≈ 1.0
    @test loaded.cells[1].reactor.source_terms !== nothing
    @test loaded.cells[1].reactor.source_terms.net == [1.0, 2.0]
    @test loaded.cells[2].endpoint_reactor === nothing
    @test loaded.cells[2].reactor.success == false
    @test loaded.cells[2].reactor.frames[1].diagnostics["reason"] == "mock failure"

    missing_path = joinpath(mktempdir(), "missing_chain_results.json")
    @test_throws ArgumentError terra.load_results_chain(missing_path)
end
