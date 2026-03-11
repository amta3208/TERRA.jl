using JSON

@testset "export_chain_profile" begin
    package_root = normpath(joinpath(@__DIR__, "..", ".."))
    script_path = joinpath(package_root, "tools", "hallthruster_jl", "export_chain_profile.jl")
    Base.include(Main, script_path)

    temp_case = mktempdir()
    input_path = joinpath(
        TEST_HET_CHAIN_INTERFACE_CASE_PATH,
        "output",
        "hallthruster_solution_avg.json",
    )

    output_path = Main.export_chain_profile([
        input_path,
        temp_case,
        "--average_start_time=0.0005",
    ])

    @test output_path == joinpath(temp_case, "input", "chain_profile_v2.json")
    @test isfile(output_path)

    raw = JSON.parsefile(output_path)
    @test raw["schema_version"] == "terra_chain_profile_v2"
    @test raw["selection"]["trim_start_index"] == 15
    @test raw["selection"]["trimmed_point_count"] == 14
    @test sort!(String.(raw["selection"]["exported_species"])) == ["N", "N+", "N2", "N2+"]

    profile = raw["profile"]
    @test haskey(profile, "species_u_m_s")
    @test sort!(collect(keys(profile["species_u_m_s"]))) == ["N", "N+", "N2", "N2+"]
    @test !haskey(profile, "u_neutral_m_s")
    @test !haskey(profile, "u_ion_m_s")
end
