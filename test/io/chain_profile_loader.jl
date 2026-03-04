@testset "load_chain_profile" begin
    @testset "Loads Fixture Profile" begin
        profile_path = joinpath(
            TEST_TERRA_CHAIN_INTERFACE_CASE_PATH, "input", "chain_profile_v1.json")
        profile = terra.load_chain_profile(profile_path)

        @test profile isa terra.AxialChainProfile
        @test profile.schema_version == "terra_chain_profile_v1"
        @test length(profile.z_m) == length(profile.dx_m) ==
              length(profile.te_K) == length(profile.u_neutral_m_s) ==
              length(profile.u_ion_m_s)
        @test haskey(profile.selection, "neutral_species")
        @test haskey(profile.generator, "tool")
        @test haskey(profile.diagnostics, "n_total_m3")
        @test profile.source_snapshot !== nothing
    end

    @testset "Supports Missing Optional Diagnostics" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_v1.json")
            write(profile_path, """
{
  "schema_version": "terra_chain_profile_v1",
  "generator": {
    "tool": "test",
    "tool_version": "0.0.0",
    "created_utc": "2026-03-04T00:00:00Z"
  },
  "selection": {
    "neutral_species": "N2",
    "ion_species": "N2",
    "ion_charge_state": 1,
    "average_start_time_s": 0.001
  },
  "profile": {
    "z_m": [0.01, 0.02, 0.03],
    "dx_m": [0.01, 0.01, 0.01],
    "te_K": [20000.0, 21000.0, 22000.0],
    "u_neutral_m_s": [200.0, 210.0, 220.0],
    "u_ion_m_s": [1000.0, 1050.0, 1100.0]
  }
}
""")

            profile = terra.load_chain_profile(profile_path)
            @test isempty(profile.diagnostics)
            @test profile.source_snapshot === nothing
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

    @testset "Rejects Unsupported Schema Version" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_bad_version.json")
            write(profile_path, """
{
  "schema_version": "terra_chain_profile_v0",
  "generator": {"tool": "test", "tool_version": "0.0.0", "created_utc": "2026-03-04T00:00:00Z"},
  "selection": {"neutral_species": "N2", "ion_species": "N2", "ion_charge_state": 1, "average_start_time_s": 0.001},
  "profile": {
    "z_m": [0.01, 0.02],
    "dx_m": [0.01, 0.01],
    "te_K": [20000.0, 21000.0],
    "u_neutral_m_s": [200.0, 210.0],
    "u_ion_m_s": [1000.0, 1100.0]
  }
}
""")

            @test_throws ArgumentError terra.load_chain_profile(profile_path)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

    @testset "Rejects Non-monotone Coordinates" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_bad_z.json")
            write(profile_path, """
{
  "schema_version": "terra_chain_profile_v1",
  "generator": {"tool": "test", "tool_version": "0.0.0", "created_utc": "2026-03-04T00:00:00Z"},
  "selection": {"neutral_species": "N2", "ion_species": "N2", "ion_charge_state": 1, "average_start_time_s": 0.001},
  "profile": {
    "z_m": [0.01, 0.02, 0.019],
    "dx_m": [0.01, 0.01, 0.01],
    "te_K": [20000.0, 21000.0, 22000.0],
    "u_neutral_m_s": [200.0, 210.0, 220.0],
    "u_ion_m_s": [1000.0, 1100.0, 1200.0]
  }
}
""")

            @test_throws ArgumentError terra.load_chain_profile(profile_path)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

    @testset "Rejects Non-positive Ion Velocities" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_bad_ui.json")
            write(profile_path, """
{
  "schema_version": "terra_chain_profile_v1",
  "generator": {"tool": "test", "tool_version": "0.0.0", "created_utc": "2026-03-04T00:00:00Z"},
  "selection": {"neutral_species": "N2", "ion_species": "N2", "ion_charge_state": 1, "average_start_time_s": 0.001},
  "profile": {
    "z_m": [0.01, 0.02, 0.03],
    "dx_m": [0.01, 0.01, 0.01],
    "te_K": [20000.0, 21000.0, 22000.0],
    "u_neutral_m_s": [200.0, 210.0, 220.0],
    "u_ion_m_s": [1000.0, -1.0, 1200.0]
  }
}
""")

            @test_throws ArgumentError terra.load_chain_profile(profile_path)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end
end
