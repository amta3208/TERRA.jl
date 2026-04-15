function _reactor_result_fixture()
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
    return reactor1
end

 @testset "Reactor Result Save" begin
    fixture = _reactor_result_fixture()
    output_path = tempname() * ".csv"

    @test terra.save_results(fixture, output_path)
    @test isfile(output_path)

    lines = split(chomp(read(output_path, String)), '\n')
    @test length(lines) == 3
    @test startswith(lines[1], "time,total_energy,T_trans,T_electron,T_vib")
    @test occursin("species_1_density", lines[1])
    @test occursin("0.0,10000.0,300.0,10000.0,320.0", lines[2])
end
