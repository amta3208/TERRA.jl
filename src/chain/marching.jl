"""
$(SIGNATURES)

Abstract root for axial chain handoff policies.
"""
abstract type AbstractChainHandoffPolicy end

"""
$(SIGNATURES)

Reuse the prior segment state when initializing downstream chain cells.
"""
struct FullStateHandoff <: AbstractChainHandoffPolicy end

"""
$(SIGNATURES)

Reinitialize each chain cell from the profile/inlet state without state-cache handoff.
"""
struct ReinitializeHandoff <: AbstractChainHandoffPolicy end

"""
$(SIGNATURES)

Abstract root for axial chain termination policies.
"""
abstract type AbstractChainTerminationPolicy end

"""
$(SIGNATURES)

Integrate each chain segment to the configured final integration time.
"""
struct FinalTimeTermination <: AbstractChainTerminationPolicy end

"""
$(SIGNATURES)

Placeholder policy for future steady-state chain termination support.
"""
struct SteadyStateTermination <: AbstractChainTerminationPolicy end

chain_policy_name(::FullStateHandoff) = "full_state"
chain_policy_name(::ReinitializeHandoff) = "reinitialize"
chain_policy_name(::FinalTimeTermination) = "final_time"
chain_policy_name(::SteadyStateTermination) = "steady_state"

Base.show(io::IO, policy::AbstractChainHandoffPolicy) = print(io, chain_policy_name(policy))
Base.show(io::IO, policy::AbstractChainTerminationPolicy) = print(io, chain_policy_name(policy))

requires_state_cache_handoff(::AbstractChainHandoffPolicy) = false
requires_state_cache_handoff(::FullStateHandoff) = true

requested_handoff_state_cache(::AbstractChainHandoffPolicy, ::Integer, upstream_state_cache) = nothing
requested_handoff_state_cache(::FullStateHandoff, segment_index::Integer, upstream_state_cache) =
    segment_index > 1 ? upstream_state_cache : nothing

state_cache_handoff_used(::AbstractChainHandoffPolicy,
                         requested_cache,
                         full_state_handoff_supported::Union{Nothing, Bool}) = false

state_cache_handoff_used(::FullStateHandoff,
                         requested_cache,
                         full_state_handoff_supported::Union{Nothing, Bool}) =
    requested_cache !== nothing &&
    full_state_handoff_supported === true &&
    requested_cache.rho_ex_cgs !== nothing

segment_tee(::FullStateHandoff,
            segment_index::Integer,
            te_profile::Real,
            inlet_reactor::ReactorConfig) =
    segment_index == 1 ? inlet_reactor.thermal.Tee : Float64(te_profile)

segment_tee(::ReinitializeHandoff,
            ::Integer,
            ::Real,
            inlet_reactor::ReactorConfig) = inlet_reactor.thermal.Tee

"""
$(SIGNATURES)

Controls for the axial-marching chain-of-CSTR solver.

# Fields
- `handoff_policy::AbstractChainHandoffPolicy`: Segment handoff behavior
- `termination_policy::AbstractChainTerminationPolicy`: Segment termination behavior
- `override_tt_K::Union{Nothing, Float64}`: Optional translational temperature override (K)
- `override_tv_K::Union{Nothing, Float64}`: Optional vibrational temperature override (K)
- `is_isothermal_teex::Bool`: Whether chain segments enforce the profile `Te(x)` isothermally
"""
struct AxialMarchingConfig{H <: AbstractChainHandoffPolicy, T <: AbstractChainTerminationPolicy}
    handoff_policy::H
    termination_policy::T
    override_tt_K::Union{Nothing, Float64}
    override_tv_K::Union{Nothing, Float64}
    is_isothermal_teex::Bool

    function AxialMarchingConfig(; handoff_policy::AbstractChainHandoffPolicy = FullStateHandoff(),
                                 termination_policy::AbstractChainTerminationPolicy = FinalTimeTermination(),
                                 override_tt_K::Union{Nothing, Real} = nothing,
                                 override_tv_K::Union{Nothing, Real} = nothing,
                                 is_isothermal_teex::Bool = true)
        if override_tt_K !== nothing && (!isfinite(override_tt_K) || override_tt_K <= 0.0)
            throw(ArgumentError("AxialMarchingConfig: override_tt_K must be finite and positive when provided."))
        end
        if override_tv_K !== nothing && (!isfinite(override_tv_K) || override_tv_K <= 0.0)
            throw(ArgumentError("AxialMarchingConfig: override_tv_K must be finite and positive when provided."))
        end

        handoff = handoff_policy
        termination = termination_policy
        return new{typeof(handoff), typeof(termination)}(handoff,
                   termination,
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
validate_axial_marching_config(marching::AxialMarchingConfig) =
    validate_axial_marching_config(marching.termination_policy)

validate_axial_marching_config(::FinalTimeTermination) = true

function validate_axial_marching_config(::SteadyStateTermination)
    throw(ArgumentError("Axial chain solver currently supports `FinalTimeTermination()` only; `SteadyStateTermination()` is not implemented yet."))
end
