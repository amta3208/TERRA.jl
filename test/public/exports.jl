 @testset "Export Surface" begin
    expected_exports = (
        :initialize_terra, :finalize_terra,
        :solve_terra_0d, :save_results,
        :species_density_matrix, :temperature_history, :total_energy_history,
        :IonNeutralizationWallModel,
        :BallisticNeutralRecombinationWallModel, :ConstantNeutralRecombinationWallModel,
        :ReactorResult,
    )

    for name in expected_exports
        @test Base.isexported(terra, name)
        @test isdefined(terra, name)
    end

    qualified_only_names = (
        :Config, :ReactorConfig, :ReactorComposition, :ReactorThermalState,
        :ModelConfig, :TimeConfig, :ODESolverConfig, :SpaceConfig,
        :NumericsConfig, :LoggingConfig, :RuntimeConfig, :ResidenceTimeConfig,
        :SourceTermsConfig, :SpeciesWallModel, :WallLossConfig, :ReactorFrame,
        :with_case_path, :with_time,
        :with_runtime, :with_logging, :nitrogen_10ev_example,
    )

    for name in qualified_only_names
        @test !Base.isexported(terra, name)
        @test isdefined(terra, name)
    end
end
