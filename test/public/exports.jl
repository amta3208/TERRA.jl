@testset "Export Surface" begin
    expected_exports = (
        :initialize_terra, :finalize_terra,
        :Config, :ReactorConfig, :ReactorComposition, :ReactorThermalState,
        :ModelConfig, :TimeConfig, :ODESolverConfig, :SpaceConfig,
        :NumericsConfig, :LoggingConfig, :RuntimeConfig, :ResidenceTimeConfig,
        :SourceTermsConfig, :SpeciesWallModel, :WallLossConfig, :ChainWallProfile,
        :ChainProfileInletComposition, :ChainProfileInlet,
        :AxialChainProfile, :AxialMarchingConfig,
        :ReactorFrame, :ReactorResult, :ChainCellResult, :ChainMetadata,
        :ChainSimulationResult, :load_chain_profile,
        :with_case_path, :with_time, :with_runtime, :with_logging,
        :solve_terra_0d, :nitrogen_10ev_example, :save_results, :load_results_chain,
        :solve_terra_chain_steady,
    )

    for name in expected_exports
        @test Base.isexported(terra, name)
        @test isdefined(terra, name)
    end
end
