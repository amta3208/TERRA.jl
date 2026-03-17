using JSON

function _write_chain_profile_json(path::AbstractString; schema_version::AbstractString,
        z_m = [0.01, 0.02, 0.03],
        dx_m = [0.01, 0.01, 0.01],
        te_K = [20000.0, 21000.0, 22000.0],
        species_u_m_s = Dict(
            "N" => [200.0, 210.0, 220.0],
            "N2" => [180.0, 181.0, 182.0],
            "N+" => [1000.0, 1050.0, 1100.0],
            "N2+" => [900.0, 950.0, 1000.0],
        ),
        diagnostics = nothing)
    payload = Dict{String, Any}(
        "schema_version" => String(schema_version),
        "generator" => Dict{String, Any}(
            "tool" => "test",
            "tool_version" => "0.0.0",
            "created_utc" => "2026-03-04T00:00:00Z",
        ),
        "selection" => Dict{String, Any}(
            "average_start_time_s" => 0.001,
            "exported_species" => ["N", "N2", "N+", "N2+"],
            "ion_velocity_policy" => "trim_to_first_positive",
            "u_ion_floor" => 0.0,
            "min_consecutive_positive" => 3,
            "trim_start_index" => 5,
            "trim_start_z_m" => 0.01,
            "trimmed_point_count" => 4,
            "original_point_count" => 7,
        ),
        "inlet" => Dict{String, Any}(
            "composition" => Dict{String, Any}(
                "species" => ["N", "N2", "N+", "N2+", "E-"],
                "mole_fractions" => [0.2, 0.65, 0.05, 0.05, 0.05],
                "total_number_density_m3" => 1.0e19,
            ),
            "thermal" => Dict{String, Any}(
                "Tt_K" => 500.0,
                "Tv_K" => 500.0,
                "Tee_K" => 20000.0,
                "Te_K" => 20000.0,
            ),
            "source_compact_index" => 1,
        ),
        "profile" => Dict{String, Any}(
            "z_m" => z_m,
            "dx_m" => dx_m,
            "te_K" => te_K,
            "species_u_m_s" => Dict{String, Any}(species_u_m_s),
        ),
    )
    if diagnostics !== nothing
        payload["diagnostics"] = diagnostics
    end

    open(path, "w") do io
        JSON.print(io, payload, 2)
        write(io, '\n')
    end
end

@testset "load_chain_profile" begin
    package_root = normpath(joinpath(@__DIR__, "..", ".."))
    script_path = joinpath(package_root, "tools", "hallthruster_jl", "export_chain_profile.jl")
    Base.include(Main, script_path)

    @testset "Loads Exported Fixture Profile" begin
        temp_case = mktempdir()
        try
            input_path = joinpath(
                TEST_HET_CHAIN_INTERFACE_CASE_PATH,
                "output",
                "hallthruster_solution_avg.json",
            )
            profile_path = Main.export_chain_profile([
                input_path,
                temp_case,
                "--average_start_time=0.0005",
            ])
            profile = terra.load_chain_profile(profile_path)

            @test profile isa terra.AxialChainProfile
            @test profile.schema_version == "terra_chain_profile_v3"
            @test length(profile.z_m) == length(profile.dx_m) == length(profile.te_K)
            @test sort!(collect(keys(profile.species_u_m_s))) == ["N", "N+", "N2", "N2+"]
            @test haskey(profile.selection, "exported_species")
            @test haskey(profile.generator, "tool")
            @test haskey(profile.diagnostics, "n_total_m3")
            @test profile.source_snapshot !== nothing
            @test profile.inlet.composition.species == ["N", "N2", "N+", "N2+", "E-"]
            @test profile.inlet.source_compact_index == 1
            @test profile.inlet.thermal.Te ≈ profile.te_K[1]
            @test profile.inlet.thermal.Tee ≈ profile.inlet.thermal.Tt
            @test profile.inlet.thermal.Tt ≈ profile.inlet.thermal.Tv
            @test profile.inlet.thermal.Te != profile.inlet.thermal.Tee
        finally
            rm(temp_case; recursive = true, force = true)
        end
    end

    @testset "Supports Missing Optional Diagnostics" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_v3.json")
            _write_chain_profile_json(profile_path; schema_version = "terra_chain_profile_v3")

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
            for schema_version in ("terra_chain_profile_v2", "terra_chain_profile_v1")
                profile_path = joinpath(temp_dir, "chain_profile_bad_version.json")
                _write_chain_profile_json(profile_path; schema_version = schema_version)
                @test_throws ArgumentError terra.load_chain_profile(profile_path)
            end
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

    @testset "Rejects Non-monotone Coordinates" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_bad_z.json")
            _write_chain_profile_json(
                profile_path;
                schema_version = "terra_chain_profile_v3",
                z_m = [0.01, 0.02, 0.019],
            )

            @test_throws ArgumentError terra.load_chain_profile(profile_path)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

    @testset "Rejects Per-Species Length Mismatch" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_bad_species_length.json")
            _write_chain_profile_json(
                profile_path;
                schema_version = "terra_chain_profile_v3",
                species_u_m_s = Dict(
                    "N" => [200.0, 210.0, 220.0],
                    "N2" => [180.0, 181.0],
                    "N+" => [1000.0, 1100.0, 1200.0],
                    "N2+" => [900.0, 1000.0, 1100.0],
                ),
            )

            @test_throws ArgumentError terra.load_chain_profile(profile_path)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end

    @testset "Rejects Non-positive Species Velocities" begin
        temp_dir = mktempdir()
        try
            profile_path = joinpath(temp_dir, "chain_profile_bad_species_velocity.json")
            _write_chain_profile_json(
                profile_path;
                schema_version = "terra_chain_profile_v3",
                species_u_m_s = Dict(
                    "N" => [200.0, 210.0, 220.0],
                    "N2" => [180.0, 181.0, 182.0],
                    "N+" => [1000.0, -1.0, 1200.0],
                    "N2+" => [900.0, 1000.0, 1100.0],
                ),
            )

            @test_throws ArgumentError terra.load_chain_profile(profile_path)
        finally
            rm(temp_dir; recursive = true, force = true)
        end
    end
end
