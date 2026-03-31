 @testset "SourceTermsConfig" begin
    cfg = terra.SourceTermsConfig()
    @test cfg.residence_time === nothing
    @test cfg.wall_losses === nothing
end

 @testset "Concrete wall model types" begin
    ion_model = terra.IonNeutralizationWallModel(;
                                                 bohm_scale = 0.5,
                                                 products = Dict("N" => 1.0))
    @test ion_model.bohm_scale == 0.5
    @test ion_model.products["N"] == 1.0

    ballistic_model = terra.BallisticNeutralRecombinationWallModel(;
                                                                   gamma = 0.25,
                                                                   products = Dict("N2" => 0.5))
    @test ballistic_model.gamma == 0.25
    @test ballistic_model.products["N2"] == 0.5

    constant_model = terra.ConstantNeutralRecombinationWallModel(;
                                                                 k_wall_1_s = 2.0e5,
                                                                 products = Dict("N2" => 0.5))
    @test constant_model.k_wall_1_s == 2.0e5
    @test constant_model.products["N2"] == 0.5

    @test_throws ArgumentError terra.IonNeutralizationWallModel(;
                                                                bohm_scale = -1.0)
    @test_throws ArgumentError terra.BallisticNeutralRecombinationWallModel(;
                                                                            gamma = -1.0)
    @test_throws ArgumentError terra.ConstantNeutralRecombinationWallModel(;
                                                                           k_wall_1_s = -1.0)
    @test_throws MethodError terra.SpeciesWallModel(;
                                                    class = :unsupported,
                                                    rate_model = :bohm_gap)
end

 @testset "WallLossConfig" begin
    ion_model = terra.IonNeutralizationWallModel(; products = Dict("N" => 1.0))
    cfg = terra.WallLossConfig(; species_models = Dict("N+" => ion_model))
    @test haskey(cfg.species_models, "N+")
    @test cfg.species_models["N+"] === ion_model

    @test_throws ArgumentError terra.WallLossConfig(;
                                                    species_models = Dict("" => ion_model))
    @test_throws MethodError terra.WallLossConfig(;
                                                  enabled = true,
                                                  species_models = Dict("N+" => ion_model))
    @test_throws MethodError terra.WallLossConfig(;
                                                  use_ion_losses = false,
                                                  species_models = Dict("N+" => ion_model))
    @test_throws MethodError terra.WallLossConfig(;
                                                  use_electronic_quenching = true,
                                                  species_models = Dict("N+" => ion_model))
end

 @testset "ResidenceTimeConfig" begin
    u_species = Dict("N" => 1.0, "N2" => 1.5, "N+" => 2.0, "N2+" => 2.5)

    rt_default = terra.ResidenceTimeConfig(1.0, u_species)
    @test rt_default.L == 1.0
    @test rt_default.U_species == u_species

    rt_custom = terra.ResidenceTimeConfig(; L = 1.5,
                                          U_species = u_species,
                                          U_energy = 3.0)
    @test rt_custom.U_energy == 3.0

    base_config = terra.nitrogen_10ev_config(; isothermal = false)

    rt_inlet_reactor = terra.ResidenceTimeConfig(;
                                                 L = 1.0,
                                                 U_species = u_species,
                                                 inlet_reactor = base_config.reactor)
    @test rt_inlet_reactor.inlet_reactor == base_config.reactor

    @test_throws MethodError terra.ResidenceTimeConfig(;
                                                       enabled = false,
                                                       L = 1.5,
                                                       U_species = u_species)

    @test_throws MethodError terra.ResidenceTimeConfig(;
                                                       enabled = true,
                                                       L = 1.0,
                                                       U_species = u_species,
                                                       inlet_config = base_config)
    @test_throws ArgumentError terra.ResidenceTimeConfig(1.0, u_species, nothing,
                                                         base_config)
end
