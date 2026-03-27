@testset "Wall losses prepare, validate, and apply on the compact state" begin
    config = terra.nitrogen_10ev_config(; isothermal = false)
    @test_nowarn reset_and_init!(tempname(); config = config)

    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()
    u_base = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
                                     rho_ex = state.rho_ex,
                                     rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
                                     rho_erot = layout.rot_noneq ? 0.0 : nothing,
                                     rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    wall_inputs = terra.SegmentWallInputs(;
                                          a_wall_over_v_m_inv = 2.0 / 0.0155,
                                          channel_gap_m = 0.0155)

    ion_wall_cfg = terra.WallLossConfig(;
                                        species_models = Dict("N+" => terra.IonNeutralizationWallModel(;
                                                                                                        bohm_scale = 1.25,
                                                                                                        products = Dict("N" => 1.0)),
                                                              "N2+" => terra.IonNeutralizationWallModel(;
                                                                                                         products = Dict("N2" => 1.0))))

    neutral_wall_cfg = terra.WallLossConfig(;
                                            species_models = Dict("N" => terra.ConstantNeutralRecombinationWallModel(;
                                                                                                                     k_wall_1_s = 4.0e5,
                                                                                                                     products = Dict("N2" => 0.5))))

    combined_wall_cfg = terra.WallLossConfig(;
                                             species_models = Dict("N+" => terra.IonNeutralizationWallModel(;
                                                                                                            bohm_scale = 1.25,
                                                                                                            products = Dict("N" => 1.0)),
                                                                   "N2+" => terra.IonNeutralizationWallModel(;
                                                                                                             products = Dict("N2" => 1.0)),
                                                                   "N" => terra.ConstantNeutralRecombinationWallModel(;
                                                                                                                      k_wall_1_s = 4.0e5,
                                                                                                                      products = Dict("N2" => 0.5))))

    prepared = terra._prepare_sources(layout, config, u_base,
                                      terra.SourceTermsConfig(; wall_losses = combined_wall_cfg);
                                      wall_inputs = wall_inputs)
    prepared_wall = terra.source_operator(prepared,
                                          terra.PreparedWallLossData)::terra.PreparedWallLossData
    @test terra.source_operator(prepared, terra.PreparedResidenceTimeSource) === nothing
    @test prepared_wall !== nothing
    @test prepared_wall.wall_inputs.channel_gap_m ≈ 0.0155
    @test !isempty(prepared_wall.models)
    @test prepared_wall.segment_te_K ≈ config.reactor.thermal.Te
    @test prepared_wall.species_index_data.ground_indices["N"] in prepared_wall.species_index_data.all_indices["N"]
    @test prepared_wall.species_index_data.ground_indices["N2"] in prepared_wall.species_index_data.all_indices["N2"]
    @test length(first(filter(model -> terra._reaction(model).reactant == "N+",
                              prepared_wall.models)).reaction.reactant_indices) ==
          length(prepared_wall.species_index_data.all_indices["N+"])
    @test first(filter(model -> terra._reaction(model).reactant == "N",
                       prepared_wall.models)).reaction.reactant_indices ==
          [prepared_wall.species_index_data.ground_indices["N"]]

    @test_throws ArgumentError terra._prepare_sources(layout, config, u_base,
                                                      terra.SourceTermsConfig(; wall_losses = combined_wall_cfg))

    bad_wall_cfg = terra.WallLossConfig(;
                                        species_models = Dict("Ar+" => terra.IonNeutralizationWallModel(;
                                                                                                         products = Dict("N" => 1.0))))
    @test_throws ArgumentError terra.prepare_source(layout, config, u_base, bad_wall_cfg;
                                                    wall_inputs = wall_inputs)

    @test_throws MethodError terra.SpeciesWallModel(;
                                                    class = :ion_neutralization,
                                                    rate_model = :constant,
                                                    parameters = Dict("k_wall_1_s" => 1.0),
                                                    products = Dict("N" => 1.0))

    molecular_weights = terra.get_molecular_weights(config.reactor.composition.species)
    species_names = config.reactor.composition.species
    atom_counts = Dict("N" => 1.0, "N2" => 2.0, "N+" => 1.0, "N2+" => 2.0)
    function nitrogen_atom_rate(drho_compact)
        drho_sp = zeros(Float64, layout.nsp)
        terra._reconstruct_drho_sp_from_compact!(drho_sp, drho_compact, layout)
        total_rate = 0.0
        for isp in eachindex(species_names)
            name = species_names[isp]
            name == "E-" && continue
            total_rate += atom_counts[name] * drho_sp[isp] * terra.AVOGADRO /
                          molecular_weights[isp]
        end
        return total_rate
    end

    ion_prepared = terra._prepare_sources(layout, config, u_base,
                                          terra.SourceTermsConfig(; wall_losses = ion_wall_cfg);
                                          wall_inputs = wall_inputs)
    ion_prepared_wall = terra.source_operator(ion_prepared,
                                              terra.PreparedWallLossData)::terra.PreparedWallLossData
    n_plus_model = only(filter(model -> terra._reaction(model).reactant == "N+",
                               ion_prepared_wall.models))
    reaction = terra._reaction(n_plus_model)
    expected_default_bohm_rate = wall_inputs.a_wall_over_v_m_inv *
                                 terra.DEFAULT_ION_EDGE_TO_CENTER_RATIO *
                                 n_plus_model.model.bohm_scale *
                                 sqrt(reaction.charge_state * terra.BOLTZMANN_J_PER_K *
                                      ion_prepared_wall.segment_te_K /
                                      terra._molecular_mass_kg(reaction.reactant_molecular_weight))
    @test terra._wall_rate_1_s(n_plus_model, ion_prepared_wall) ≈
          expected_default_bohm_rate

    ion_indices = ion_prepared_wall.species_index_data
    u_ion = zeros(length(u_base))
    for (j, idx) in pairs(ion_indices.all_indices["N+"])
        u_ion[idx] = 1.0e-12 * j
    end
    for (j, idx) in pairs(ion_indices.all_indices["N2+"])
        u_ion[idx] = 2.0e-12 * j
    end
    n_ground_idx = ion_indices.ground_indices["N"]
    n2_ground_idx = ion_indices.ground_indices["N2"]
    n_excited_idx = first(setdiff(ion_indices.all_indices["N"], [n_ground_idx]))
    n2_excited_idx = first(setdiff(ion_indices.all_indices["N2"], [n2_ground_idx]))
    u_ion[n_excited_idx] = 7.0e-13
    u_ion[n2_excited_idx] = 9.0e-13

    du_ion = zeros(length(u_base))
    terra.apply_source!(du_ion, u_ion, ion_prepared_wall)
    @test all(du_ion[idx] < 0.0 for idx in ion_indices.all_indices["N+"])
    @test all(du_ion[idx] < 0.0 for idx in ion_indices.all_indices["N2+"])
    @test du_ion[n_ground_idx] > 0.0
    @test du_ion[n2_ground_idx] > 0.0
    @test all(du_ion[idx] == 0.0
              for idx in setdiff(ion_indices.all_indices["N"], [n_ground_idx]))
    @test all(du_ion[idx] == 0.0
              for idx in setdiff(ion_indices.all_indices["N2"], [n2_ground_idx]))
    @test all(du_ion[idx] == 0.0 for idx in layout.energy_range)
    @test abs(nitrogen_atom_rate(du_ion)) ≤ 1e-6

    neutral_prepared = terra._prepare_sources(layout, config, u_base,
                                              terra.SourceTermsConfig(; wall_losses = neutral_wall_cfg);
                                              wall_inputs = wall_inputs)
    neutral_prepared_wall = terra.source_operator(neutral_prepared,
                                                  terra.PreparedWallLossData)::terra.PreparedWallLossData
    neutral_indices = neutral_prepared_wall.species_index_data
    u_neutral = zeros(length(u_base))
    u_neutral[neutral_indices.ground_indices["N"]] = 4.0e-12
    neutral_excited_idx = first(setdiff(neutral_indices.all_indices["N"],
                                        [neutral_indices.ground_indices["N"]]))
    u_neutral[neutral_excited_idx] = 6.0e-12

    du_neutral = zeros(length(u_base))
    terra.apply_source!(du_neutral, u_neutral, neutral_prepared_wall)
    @test du_neutral[neutral_indices.ground_indices["N"]] < 0.0
    @test du_neutral[neutral_indices.ground_indices["N2"]] > 0.0
    @test all(du_neutral[idx] == 0.0
              for idx in setdiff(neutral_indices.all_indices["N"],
                                 [neutral_indices.ground_indices["N"]]))
    @test all(du_neutral[idx] == 0.0 for idx in layout.energy_range)
    @test abs(nitrogen_atom_rate(du_neutral)) ≤ 1e-6

    du_base = zeros(length(u_base))
    du_wall = zeros(length(u_base))
    du_expected_wall = zeros(length(u_base))

    p_base = (layout = layout,
              config = config,
              teex_const = state.teex_const,
              teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp))
    p_wall = (layout = layout,
              config = config,
              teex_const = state.teex_const,
              teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp),
              sources = terra.PreparedSources(prepared_wall))

    @test_nowarn terra.terra_ode_system!(du_base, u_base, p_base, 0.0)
    @test_nowarn terra.terra_ode_system!(du_wall, u_base, p_wall, 0.0)
    terra.apply_source!(du_expected_wall, u_base, prepared_wall)
    @test any(abs.(du_expected_wall) .> 0.0)
    @test du_wall .- du_base ≈ du_expected_wall rtol = 1e-12 atol = 1e-13
    @test (du_wall .- du_base)[prepared_wall.species_index_data.ground_indices["N2"]] > 0.0
    @test any((du_wall .- du_base)[idx] < 0.0
              for idx in prepared_wall.species_index_data.all_indices["N2+"])
end
