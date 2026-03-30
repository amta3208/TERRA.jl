 @progress_testset "Export Surface" begin
    expected_exports = (
        :initialize_terra, :finalize_terra,
        :solve_terra_0d, :solve_terra_chain_steady,
        :load_chain_profile, :save_results, :load_results_chain,
        :species_density_matrix, :temperature_history, :total_energy_history,
        :WallLossConfig,
        :IonNeutralizationWallModel,
        :BallisticNeutralRecombinationWallModel, :ConstantNeutralRecombinationWallModel,
        :ChainWallProfile,
        :ChainProfileInletComposition, :ChainProfileInlet,
        :AxialChainProfile, :AxialMarchingConfig,
        :FullStateHandoff, :ReinitializeHandoff, :FinalTimeTermination,
        :ReactorResult, :ChainSimulationResult,
    )

    for name in expected_exports
        @test Base.isexported(terra, name)
        @test isdefined(terra, name)
    end

    qualified_only_names = (
        :Config, :ReactorConfig, :ReactorComposition, :ReactorThermalState,
        :ModelConfig, :TimeConfig, :ODESolverConfig, :SpaceConfig,
        :NumericsConfig, :LoggingConfig, :RuntimeConfig, :ResidenceTimeConfig,
        :SourceTermsConfig, :SpeciesWallModel, :AbstractChainHandoffPolicy,
        :AbstractChainTerminationPolicy, :SteadyStateTermination, :ReactorFrame,
        :ChainCellResult, :ChainMetadata, :with_case_path, :with_time,
        :with_runtime, :with_logging, :nitrogen_10ev_example,
    )

    for name in qualified_only_names
        @test !Base.isexported(terra, name)
        @test isdefined(terra, name)
    end
end
