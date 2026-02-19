@testset "Library Management" begin
    @testset "Library Loading and Status" begin
        # Test initial state - no library loaded
        @test !terra.is_terra_loaded()
    end

    @testset "Library Path Setting" begin
        # Test setting library path with non-existent file
        fake_path = "/nonexistent/path/libterra.so"
        @test_throws ErrorException terra.load_terra_library!(fake_path)

        # Test error message for non-existent file
        try
            terra.load_terra_library!(fake_path)
            @test false  # Should not reach here
        catch e
            @test occursin("TERRA library file not found", e.msg)
            @test occursin(fake_path, e.msg)
        end

        # Test that library is still not loaded after failed attempt
        @test !terra.is_terra_loaded()
    end

    @testset "Library Cleanup" begin
        # Test closing library when none is loaded (should be safe)
        @test_nowarn terra.close_terra_library()
        @test !terra.is_terra_loaded()

        # Test multiple calls to close (should be safe)
        @test_nowarn terra.close_terra_library()
        @test_nowarn terra.close_terra_library()
        @test !terra.is_terra_loaded()
    end
end
