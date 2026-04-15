"""
$(SIGNATURES)

Wrapper-managed additive source-term configuration.
"""
struct SourceTermsConfig
    wall_losses::Union{Nothing, WallLossConfig}

    function SourceTermsConfig(; wall_losses::Union{Nothing, WallLossConfig} = nothing)
        return new(wall_losses)
    end
end
