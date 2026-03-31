"""
$(SIGNATURES)

Wrapper-managed additive source-term configuration.
"""
struct SourceTermsConfig
    residence_time::Union{Nothing, ResidenceTimeConfig}
    wall_losses::Union{Nothing, WallLossConfig}

    function SourceTermsConfig(;
                               residence_time::Union{Nothing, ResidenceTimeConfig} = nothing,
                               wall_losses::Union{Nothing, WallLossConfig} = nothing)
        return new(residence_time, wall_losses)
    end
end
