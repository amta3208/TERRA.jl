@testset "SourceTermsConfig" begin
    cfg = terra.SourceTermsConfig()
    @test cfg.residence_time === nothing
    @test cfg.wall_losses === nothing
end

@testset "SpeciesWallModel" begin
    model = terra.SpeciesWallModel(;
                                   class = :ion_neutralization,
                                   rate_model = :bohm_gap,
                                   parameters = Dict("alpha" => 0.5),
                                   products = Dict("N" => 1.0),)
    @test model.class == :ion_neutralization
    @test model.rate_model == :bohm_gap
    @test model.parameters["alpha"] == 0.5
    @test model.products["N"] == 1.0

    @test_throws ArgumentError terra.SpeciesWallModel(;
                                                      class = :unsupported,
                                                      rate_model = :bohm_gap)
    @test_throws ArgumentError terra.SpeciesWallModel(;
                                                      class = :ion_neutralization,
                                                      rate_model = :unsupported)
    @test_throws ArgumentError terra.SpeciesWallModel(;
                                                      class = :ion_neutralization,
                                                      rate_model = :bohm_gap,
                                                      parameters = Dict("alpha" => -1.0))
end

@testset "WallLossConfig" begin
    cfg = terra.WallLossConfig(;
                               species_models = Dict("N+" => terra.SpeciesWallModel(;
                                                                                    class = :ion_neutralization,
                                                                                    rate_model = :bohm_gap,
                                                                                    products = Dict("N" => 1.0),)),)
    @test cfg.enabled == true
    @test cfg.use_ion_losses == true
    @test cfg.use_neutral_recombination == false
    @test cfg.use_electronic_quenching == false
    @test haskey(cfg.species_models, "N+")
    @test_throws ArgumentError terra.WallLossConfig(;
                                                    species_models = Dict("" => cfg.species_models["N+"]))
end

@testset "ResidenceTimeConfig" begin
    u_species = Dict("N" => 1.0, "N2" => 1.5, "N+" => 2.0, "N2+" => 2.5)

    rt_default = terra.ResidenceTimeConfig(1.0, u_species)
    @test rt_default.enabled == true
    @test rt_default.L == 1.0
    @test rt_default.U_species == u_species

    rt_disabled = terra.ResidenceTimeConfig(; enabled = false, L = 1.5,
                                            U_species = u_species,
                                            U_energy = 3.0)
    @test rt_disabled.enabled == false
    @test rt_disabled.U_energy == 3.0

    base_config = terra.nitrogen_10ev_config(; isothermal = false)

    rt_inlet_reactor = terra.ResidenceTimeConfig(;
                                                 enabled = true,
                                                 L = 1.0,
                                                 U_species = u_species,
                                                 inlet_reactor = base_config.reactor)
    @test rt_inlet_reactor.inlet_reactor == base_config.reactor

    @test_throws MethodError terra.ResidenceTimeConfig(;
                                                       enabled = true,
                                                       L = 1.0,
                                                       U_species = u_species,
                                                       inlet_config = base_config)
    @test_throws ArgumentError terra.ResidenceTimeConfig(1.0, u_species, nothing,
                                                         base_config)
end
