using JSON

function _chain_result_fixture()
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
            "wall_losses" => Dict{String, Any}(
                "species_order" => ["N", "N2", "N+", "N2+"],
                "segment_inputs" => Dict{String, Any}(
                    "a_wall_over_v_m_inv" => 129.03225806451613,
                    "channel_gap_m" => 0.0155,
                    "wall_temperature_K" => nothing,
                    "ion_edge_to_center_ratio" => nothing,
                    "tt_K" => 310.0,
                    "te_K" => 9800.0,
                ),
                "species_models" => Dict{String, Any}(
                    "N+" => Dict{String, Any}(
                        "model_type" => "IonNeutralizationWallModel",
                        "charge_state" => 1,
                        "parameters" => Dict{String, Any}("bohm_scale" => 1.0),
                        "products" => Dict{String, Any}("N" => 1.0),
                        "reactant_indices" => [3, 4],
                        "reactant_ground_index" => 3,
                        "product_indices" => Dict{String, Any}("N" => 1),
                    ),
                ),
            ),
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
        schema_version = "terra_chain_profile_v4",
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

    return (reactor = reactor1, chain = chain)
end

@testset "Reactor Result Save" begin
    fixture = _chain_result_fixture()
    output_path = tempname() * ".csv"

    @test terra.save_results(fixture.reactor, output_path)
    @test isfile(output_path)

    lines = split(chomp(read(output_path, String)), '\n')
    @test length(lines) == 3
    @test startswith(lines[1], "time,total_energy,T_trans,T_electron,T_vib")
    @test occursin("species_1_density", lines[1])
    @test occursin("0.0,10000.0,300.0,10000.0,320.0", lines[2])
end

@testset "Chain Result Save/Load" begin
    fixture = _chain_result_fixture()
    output_path = tempname() * ".json"

    @test terra.save_results(fixture.chain, output_path)
    @test isfile(output_path)

    saved = JSON.parsefile(output_path)
    @test saved["schema_version"] == "terra_chain_results_v3"

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
    @test terra.temperature_history(loaded.cells[1].reactor).te[2] ≈ 9800.0
    @test size(terra.species_density_matrix(loaded.cells[1].reactor)) == (3, 2)
    @test loaded.cells[1].species_u_m_s["N+"] == 15000.0
    @test loaded.cells[1].reactor.frames[1].source_terms !== nothing
    @test loaded.cells[1].reactor.frames[1].source_terms.production ≈ 1.0
    @test loaded.cells[1].reactor.source_terms !== nothing
    @test loaded.cells[1].reactor.source_terms.net == [1.0, 2.0]
    species_model = loaded.cells[1].reactor.metadata["wall_losses"]["species_models"]["N+"]
    @test species_model["model_type"] == "IonNeutralizationWallModel"
    @test sort!(collect(keys(loaded.cells[1].reactor.metadata["wall_losses"]))) ==
          ["segment_inputs", "species_models", "species_order"]
    @test sort!(collect(keys(species_model))) ==
          ["charge_state", "model_type", "parameters", "product_indices", "products",
           "reactant_ground_index", "reactant_indices"]
    @test loaded.cells[2].endpoint_reactor === nothing
    @test loaded.cells[2].reactor.success == false
    @test loaded.cells[2].reactor.frames[1].diagnostics["reason"] == "mock failure"

    missing_path = joinpath(mktempdir(), "missing_chain_results.json")
    @test_throws ArgumentError terra.load_results_chain(missing_path)
end

@testset "Chain Result Validation" begin
    fixture = _chain_result_fixture()
    output_path = tempname() * ".json"
    @test terra.save_results(fixture.chain, output_path)

    bad_schema_path = tempname() * ".json"
    bad_schema = JSON.parsefile(output_path)
    bad_schema["schema_version"] = "terra_chain_results_v1"
    open(bad_schema_path, "w") do io
        JSON.print(io, bad_schema, 2)
        println(io)
    end
    @test_throws ArgumentError terra.load_results_chain(bad_schema_path)

    bad_count_path = tempname() * ".json"
    bad_count = JSON.parsefile(output_path)
    bad_count["metadata"]["retained_point_count"] = 5
    open(bad_count_path, "w") do io
        JSON.print(io, bad_count, 2)
        println(io)
    end
    @test_throws ArgumentError terra.load_results_chain(bad_count_path)

    bad_metadata_path = tempname() * ".json"
    bad_metadata = JSON.parsefile(output_path)
    delete!(bad_metadata["cells"][1]["reactor"]["metadata"]["wall_losses"]["species_models"]["N+"],
            "model_type")
    open(bad_metadata_path, "w") do io
        JSON.print(io, bad_metadata, 2)
        println(io)
    end
    @test_throws ArgumentError terra.load_results_chain(bad_metadata_path)
end
