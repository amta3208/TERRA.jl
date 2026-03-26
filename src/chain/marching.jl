"""
$(SIGNATURES)

Controls for the axial-marching chain-of-CSTR solver.

# Fields
- `handoff_mode::Symbol`: Segment handoff mode (`:reinitialize` or `:full_state`)
- `termination_mode::Symbol`: Segment termination mode (`:final_time` or `:steady_state`)
- `override_tt_K::Union{Nothing, Float64}`: Optional translational temperature override (K)
- `override_tv_K::Union{Nothing, Float64}`: Optional vibrational temperature override (K)
- `is_isothermal_teex::Bool`: Whether chain segments enforce the profile `Te(x)` isothermally
"""
struct AxialMarchingConfig
    handoff_mode::Symbol
    termination_mode::Symbol
    override_tt_K::Union{Nothing, Float64}
    override_tv_K::Union{Nothing, Float64}
    is_isothermal_teex::Bool

    function AxialMarchingConfig(; handoff_mode::Symbol = :full_state,
                                 termination_mode::Symbol = :final_time,
                                 override_tt_K::Union{Nothing, Real} = nothing,
                                 override_tv_K::Union{Nothing, Real} = nothing,
                                 is_isothermal_teex::Bool = true)
        handoff_mode in (:reinitialize, :full_state) ||
            throw(ArgumentError("AxialMarchingConfig: handoff_mode must be :reinitialize or :full_state."))
        termination_mode in (:final_time, :steady_state) ||
            throw(ArgumentError("AxialMarchingConfig: termination_mode must be :final_time or :steady_state."))
        if override_tt_K !== nothing && (!isfinite(override_tt_K) || override_tt_K <= 0.0)
            throw(ArgumentError("AxialMarchingConfig: override_tt_K must be finite and positive when provided."))
        end
        if override_tv_K !== nothing && (!isfinite(override_tv_K) || override_tv_K <= 0.0)
            throw(ArgumentError("AxialMarchingConfig: override_tv_K must be finite and positive when provided."))
        end

        return new(handoff_mode,
                   termination_mode,
                   override_tt_K === nothing ? nothing : Float64(override_tt_K),
                   override_tv_K === nothing ? nothing : Float64(override_tv_K),
                   is_isothermal_teex)
    end
end

"""
$(SIGNATURES)

Validate an `AxialMarchingConfig` for the current chain solver capabilities.

# Arguments
- `marching::AxialMarchingConfig`: Marching controls to validate

# Returns
- `true` if validation passes

# Throws
- `ArgumentError` if unsupported modes are requested
"""
function validate_axial_marching_config(marching::AxialMarchingConfig)
    if marching.termination_mode != :final_time
        throw(ArgumentError("Axial chain solver currently supports termination_mode=:final_time only (got $(marching.termination_mode))."))
    end
    return true
end
