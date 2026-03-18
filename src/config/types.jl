"""
# TERRA Configuration Module

This module handles configuration management for TERRA simulations,
including parameter validation, default values, and input file generation.
"""

"""
$(SIGNATURES)

Physics modeling configuration for TERRA simulation.

# Fields
- `bbh_model::Int`: Bound-bound heavy particle model
- `esc_model::Int`: Escape model
- `ar_et_model::Int`: Ar-ET model
- `eex_noneq::Int`: Electron-electronic nonequilibrium flag
- `ev_relax_set::Int`: Electron-vibrational relaxation set
- `et_relax_set::Int`: Electron-translational relaxation set
- `radiation_length::Float64`: Radiation length scale (cm)
- `get_electron_density_by_charge_balance::Bool`: Electron density by charge balance
- `min_sts_frac::Float64`: Minimum state-to-state fraction
- `is_isothermal_teex::Bool`: Isothermal electron-electronic flag
- `energy_loss_per_eii::Float64`: Average electron energy loss per EII event (× E_ion)
"""
struct PhysicsConfig
    bbh_model::Int
    esc_model::Int
    ar_et_model::Int
    eex_noneq::Int
    ev_relax_set::Int
    et_relax_set::Int
    radiation_length::Float64
    get_electron_density_by_charge_balance::Bool
    min_sts_frac::Float64
    is_isothermal_teex::Bool
    energy_loss_per_eii::Float64

    function PhysicsConfig(;
            bbh_model = 4,
            esc_model = 1,
            ar_et_model = 1,
            eex_noneq = 1,
            ev_relax_set = 1,
            et_relax_set = 1,
            radiation_length = 1.0,
            get_electron_density_by_charge_balance = true,
            min_sts_frac = 1e-30,
            is_isothermal_teex = false,
            energy_loss_per_eii = 1.0
    )
        new(bbh_model, esc_model, ar_et_model, eex_noneq, ev_relax_set, et_relax_set,
            radiation_length, get_electron_density_by_charge_balance,
            min_sts_frac, is_isothermal_teex, energy_loss_per_eii)
    end
end

"""
$(SIGNATURES)

Process flags configuration for TERRA simulation.

# Fields
- `consider_elec_bbe::Int`: Consider electron bound-bound excitation
- `consider_elec_bfe::Int`: Consider electron bound-free excitation
- `consider_elec_bbh::Int`: Consider electron bound-bound heavy
- `consider_elec_bfh::Int`: Consider electron bound-free heavy
- `consider_rad::Int`: Consider radiation
- `consider_rdr::Int`: Consider RDR
- `consider_chem::Int`: Consider chemistry
"""
struct ProcessConfig
    consider_elec_bbe::Int
    consider_elec_bfe::Int
    consider_elec_bbh::Int
    consider_elec_bfh::Int
    consider_rad::Int
    consider_rdr::Int
    consider_chem::Int

    function ProcessConfig(;
            consider_elec_bbe = 1,
            consider_elec_bfe = 1,
            consider_elec_bbh = 1,
            consider_elec_bfh = 1,
            consider_rad = 0,
            consider_rdr = 0,
            consider_chem = 1
    )
        new(consider_elec_bbe, consider_elec_bfe, consider_elec_bbh,
            consider_elec_bfh, consider_rad, consider_rdr, consider_chem)
    end
end

"""
$(SIGNATURES)

Reactor composition inputs for a simulation case.

# Fields
- `species::Vector{String}`: Species names
- `mole_fractions::Vector{Float64}`: Species mole fractions (must sum to 1)
- `total_number_density::Float64`: Total number density
"""
struct ReactorComposition
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density::Float64

    function ReactorComposition(species, mole_fractions, total_number_density)
        species_vec = String.(species)
        mole_frac_vec = Float64.(mole_fractions)
        n_tot = Float64(total_number_density)

        if length(species_vec) != length(mole_frac_vec)
            throw(ArgumentError("Species and mole_fractions arrays must have same length"))
        end
        if isempty(species_vec)
            throw(ArgumentError("At least one species must be specified"))
        end
        if any(mole_frac_vec .< 0)
            throw(ArgumentError("Mole fractions must be non-negative"))
        end
        if abs(sum(mole_frac_vec) - 1.0) > 1e-10
            throw(ArgumentError("Mole fractions must sum to 1.0, got $(sum(mole_frac_vec))"))
        end
        if n_tot <= 0
            throw(ArgumentError("Total number density must be positive"))
        end
        if length(unique(species_vec)) != length(species_vec)
            throw(ArgumentError("Duplicate species names found"))
        end
        for name in species_vec
            if isempty(strip(name))
                throw(ArgumentError("Species names cannot be empty"))
            end
        end

        return new(species_vec, mole_frac_vec, n_tot)
    end
end

ReactorComposition(; species, mole_fractions, total_number_density) = ReactorComposition(
    species, mole_fractions, total_number_density)

"""
$(SIGNATURES)

Thermal state for the reactor.

# Fields
- `Tt::Float64`: Translational temperature (K)
- `Tv::Float64`: Vibrational temperature (K)
- `Tee::Float64`: Electron-electronic temperature (K)
- `Te::Float64`: Electron temperature (K)
"""
struct ReactorThermalState
    Tt::Float64
    Tv::Float64
    Tee::Float64
    Te::Float64

    function ReactorThermalState(Tt, Tv, Tee, Te)
        if any([Tt, Tv, Tee, Te] .<= 0)
            throw(ArgumentError("All temperatures must be positive"))
        end
        return new(Float64(Tt), Float64(Tv), Float64(Tee), Float64(Te))
    end
end

ReactorThermalState(; Tt, Tv, Tee, Te) = ReactorThermalState(Tt, Tv, Tee, Te)

"""
$(SIGNATURES)

Combined reactor-state configuration.
"""
struct ReactorConfig
    composition::ReactorComposition
    thermal::ReactorThermalState
end

function ReactorConfig(;
        composition::ReactorComposition,
        thermal::ReactorThermalState)
    return ReactorConfig(composition, thermal)
end

"""
$(SIGNATURES)

Physics/process model configuration.
"""
struct ModelConfig
    physics::PhysicsConfig
    processes::ProcessConfig
end

function ModelConfig(;
        physics::PhysicsConfig = PhysicsConfig(),
        processes::ProcessConfig = ProcessConfig())
    return ModelConfig(physics, processes)
end

"""
$(SIGNATURES)

Numerical time-integration controls.

# Fields
- `dt::Float64`: Time step (seconds)
- `dt_output::Float64`: Native-output time step (seconds)
- `duration::Float64`: Final time (seconds)
- `nstep::Int`: Maximum number of time steps
- `method::Int`: Integration method (0=forward Euler, 1=high order explicit, 2=implicit)
"""
struct TimeConfig
    dt::Float64
    dt_output::Float64
    duration::Float64
    nstep::Int
    method::Int

    function TimeConfig(dt, dt_output, duration, nstep = 500000, method = 2)
        if dt <= 0 || dt_output <= 0 || duration <= 0
            throw(ArgumentError("Time parameters must be positive"))
        end
        if nstep <= 0
            throw(ArgumentError("Number of steps must be positive"))
        end
        if !(method in [0, 1, 2])
            throw(ArgumentError("Integration method must be 0, 1, or 2"))
        end
        return new(Float64(dt), Float64(dt_output), Float64(duration), Int(nstep), Int(method))
    end
end

function TimeConfig(;
        dt, dt_output, duration, nstep::Integer = 500000, method::Integer = 2)
    return TimeConfig(dt, dt_output, duration, nstep, method)
end

"""
$(SIGNATURES)

ODE-solver controls managed by the Julia wrapper.
"""
struct ODESolverConfig
    reltol::Float64
    abstol_density::Float64
    saveat_count::Int
    ramp_understep_ratio::Float64
    ramp_history_steps::Int

    function ODESolverConfig(;
            reltol::Real = 1e-8,
            abstol_density::Real = 1e-10,
            saveat_count::Integer = 100,
            ramp_understep_ratio::Real = inv(128),
            ramp_history_steps::Integer = 5)
        if reltol <= 0
            throw(ArgumentError("reltol must be positive"))
        end
        if abstol_density <= 0
            throw(ArgumentError("abstol_density must be positive"))
        end
        if saveat_count <= 0
            throw(ArgumentError("saveat_count must be positive"))
        end
        if ramp_understep_ratio <= 0
            throw(ArgumentError("ramp_understep_ratio must be positive"))
        end
        if ramp_history_steps <= 0
            throw(ArgumentError("ramp_history_steps must be positive"))
        end

        return new(
            Float64(reltol),
            Float64(abstol_density),
            Int(saveat_count),
            Float64(ramp_understep_ratio),
            Int(ramp_history_steps))
    end
end

"""
$(SIGNATURES)

Spatial-discretization metadata.
"""
struct SpaceConfig
    nd::Int
    dr::Union{Nothing, Float64}

    function SpaceConfig(; nd::Integer = 0, dr::Union{Nothing, Real} = nothing)
        nd_val = Int(nd)
        if nd_val < 0
            throw(ArgumentError("nd must be non-negative"))
        end
        dr_val = dr === nothing ? nothing : Float64(dr)
        if dr_val !== nothing && dr_val <= 0
            throw(ArgumentError("dr must be positive when provided"))
        end
        return new(nd_val, dr_val)
    end
end

"""
$(SIGNATURES)

Residence-time (CSTR / flow-through) model controls.
"""
struct ResidenceTimeConfig
    enabled::Bool
    L::Float64
    U_species::Dict{String, Float64}
    U_energy::Union{Nothing, Float64}
    inlet_reactor::Union{Nothing, ReactorConfig}

    function ResidenceTimeConfig(enabled::Bool, L, U_species,
            U_energy = nothing, inlet_reactor = nothing)
        if !isfinite(L) || L <= 0
            throw(ArgumentError("ResidenceTimeConfig: L must be finite and positive (got $L)."))
        end
        u_species_dict = Dict{String, Float64}()
        for (name, value) in pairs(U_species)
            name_str = String(name)
            isempty(strip(name_str)) && throw(ArgumentError(
                "ResidenceTimeConfig: species keys in U_species must be non-empty."
            ))
            if !isfinite(value) || value <= 0
                throw(ArgumentError(
                    "ResidenceTimeConfig: U_species[$name_str] must be finite and positive (got $value)."
                ))
            end
            u_species_dict[name_str] = Float64(value)
        end
        if U_energy !== nothing && (!isfinite(U_energy) || U_energy <= 0)
            throw(ArgumentError("ResidenceTimeConfig: U_energy must be finite and positive (got $U_energy)."))
        end
        inlet_reactor_val = _coerce_residence_time_inlet_reactor(inlet_reactor)

        return new(
            enabled,
            Float64(L),
            u_species_dict,
            U_energy === nothing ? nothing : Float64(U_energy),
            inlet_reactor_val)
    end
end

ResidenceTimeConfig(L, U_species, U_energy = nothing, inlet_config = nothing) =
    ResidenceTimeConfig(true, L, U_species, U_energy, inlet_config)

function ResidenceTimeConfig(;
        enabled::Bool = true,
        L::Real = 1.0,
        U_species::AbstractDict = Dict{String, Float64}(),
        U_energy::Union{Nothing, Real} = nothing,
        inlet_reactor = nothing,
        inlet_config = nothing)
    if inlet_reactor !== nothing && inlet_config !== nothing
        throw(ArgumentError("ResidenceTimeConfig: provide only one of `inlet_reactor` or legacy `inlet_config`."))
    end
    inlet_value = inlet_reactor === nothing ? inlet_config : inlet_reactor
    return ResidenceTimeConfig(enabled, L, U_species, U_energy, inlet_value)
end

function _coerce_residence_time_inlet_reactor(inlet)
    if inlet === nothing
        return nothing
    elseif inlet isa ReactorConfig
        return inlet
    elseif inlet isa Config
        return inlet.reactor
    else
        throw(ArgumentError("ResidenceTimeConfig: inlet_reactor/inlet_config must be `ReactorConfig`, `Config`, or `nothing`."))
    end
end

"""
$(SIGNATURES)

Wrapper-managed additive source-term configuration.
"""
struct SourceTermsConfig
    residence_time::Union{Nothing, ResidenceTimeConfig}

    function SourceTermsConfig(;
            residence_time::Union{Nothing, ResidenceTimeConfig} = nothing)
        return new(residence_time)
    end
end

"""
$(SIGNATURES)

Solver and discretization controls for the wrapper runtime.
"""
struct NumericsConfig
    time::TimeConfig
    solver::ODESolverConfig
    space::SpaceConfig

    function NumericsConfig(;
            time::TimeConfig,
            solver::ODESolverConfig = ODESolverConfig(),
            space::SpaceConfig = SpaceConfig())
        return new(time, solver, space)
    end
end

"""
$(SIGNATURES)

Runtime and I/O configuration controls.
"""
struct RuntimeConfig
    database_path::String
    case_path::String
    unit_system::Symbol
    validate_species_against_terra::Bool
    print_source_terms::Bool
    write_native_outputs::Bool
    print_integration_output::Bool

    function RuntimeConfig(;
            database_path::String = "../../databases/n2/elec_sts_expanded_electron_fits",
            case_path::String = pwd(),
            unit_system::Symbol = :CGS,
            validate_species_against_terra::Bool = false,
            print_source_terms::Bool = true,
            write_native_outputs::Bool = false,
            print_integration_output::Bool = true)
        if !isdir(case_path)
            throw(ArgumentError("Case path directory does not exist: $case_path"))
        end
        if !(unit_system in [:SI, :CGS])
            throw(ArgumentError("Unit system must be :SI or :CGS, got :$unit_system"))
        end
        return new(database_path, case_path, unit_system,
            validate_species_against_terra, print_source_terms,
            write_native_outputs, print_integration_output)
    end
end

"""
$(SIGNATURES)

Top-level configuration for refactored TERRA workflows.
"""
struct Config
    reactor::ReactorConfig
    models::ModelConfig
    sources::SourceTermsConfig
    numerics::NumericsConfig
    runtime::RuntimeConfig

    function Config(;
            reactor::ReactorConfig,
            numerics::NumericsConfig,
            sources::SourceTermsConfig = SourceTermsConfig(),
            models::ModelConfig = ModelConfig(),
            runtime::RuntimeConfig = RuntimeConfig())
        return new(reactor, models, sources, numerics, runtime)
    end
end

_coerce_residence_time_inlet_reactor(inlet::Config) = inlet.reactor

"""
$(SIGNATURES)

Time-frame result for a single reactor solve.

# Fields
- `t::Float64`: Frame timestamp
- `species_densities::Vector{Float64}`: Species mass densities at this frame
- `temperatures::NamedTuple`: Temperature state at this frame
- `total_energy::Float64`: Total energy density at this frame
- `source_terms::Union{NamedTuple, Nothing}`: Optional source-term snapshot
- `diagnostics::Dict{String, Any}`: Optional frame-level diagnostics
"""
struct ReactorFrame
    t::Float64
    species_densities::Vector{Float64}
    temperatures::NamedTuple
    total_energy::Float64
    source_terms::Union{NamedTuple, Nothing}
    diagnostics::Dict{String, Any}

    function ReactorFrame(;
            t::Real,
            species_densities,
            temperatures::NamedTuple,
            total_energy::Real,
            source_terms::Union{NamedTuple, Nothing} = nothing,
            diagnostics::AbstractDict = Dict{String, Any}())
        species_densities_vec = Float64.(species_densities)
        diagnostics_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(diagnostics))

        return new(
            Float64(t),
            species_densities_vec,
            temperatures,
            Float64(total_energy),
            source_terms,
            diagnostics_dict
        )
    end
end

"""
$(SIGNATURES)

Time-history container for one reactor (one chain cell).

# Fields
- `t::Vector{Float64}`: Saved times
- `frames::Vector{ReactorFrame}`: Saved per-time reactor frames
- `success::Bool`: Reactor success flag
- `message::String`: Reactor status message
- `source_terms::Union{NamedTuple, Nothing}`: Optional source history payload
- `metadata::Dict{String, Any}`: Optional reactor metadata

# Notes
- Supports HallThruster-like slicing: `reactor[i]` returns a one-frame `ReactorResult`.
"""
struct ReactorResult
    t::Vector{Float64}
    frames::Vector{ReactorFrame}
    success::Bool
    message::String
    source_terms::Union{NamedTuple, Nothing}
    metadata::Dict{String, Any}

    function ReactorResult(;
            t,
            frames,
            success::Bool = true,
            message::AbstractString = "",
            source_terms::Union{NamedTuple, Nothing} = nothing,
            metadata::AbstractDict = Dict{String, Any}())
        t_vec = Float64.(t)
        frames_vec = ReactorFrame[frame for frame in frames]
        length(frames_vec) == length(t_vec) || throw(ArgumentError(
            "ReactorResult: `t` and `frames` must have identical lengths."
        ))

        metadata_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(metadata))

        return new(
            t_vec,
            frames_vec,
            success,
            String(message),
            source_terms,
            metadata_dict
        )
    end
end

Base.firstindex(result::ReactorResult) = 1
Base.lastindex(result::ReactorResult) = length(result.frames)

function Base.getindex(result::ReactorResult, frame::Integer)
    return ReactorResult(;
        t = [result.t[frame]],
        frames = [result.frames[frame]],
        success = result.success,
        message = result.message,
        source_terms = result.source_terms,
        metadata = copy(result.metadata)
    )
end

function Base.getindex(result::ReactorResult, frames::AbstractVector{<:Integer})
    return ReactorResult(;
        t = result.t[frames],
        frames = result.frames[frames],
        success = result.success,
        message = result.message,
        source_terms = result.source_terms,
        metadata = copy(result.metadata)
    )
end

"""
$(SIGNATURES)

Metadata for chain-level result organization and segmentation provenance.
"""
struct ChainMetadata
    schema_version::String
    generator::Dict{String, Any}
    selection::Dict{String, Any}
    source_snapshot::Union{Nothing, Dict{String, Any}}
    diagnostics::Dict{String, Any}
    compact_to_source_index::Vector{Int}
    original_point_count::Union{Nothing, Int}
    retained_point_count::Int

    function ChainMetadata(;
            schema_version::AbstractString = "terra_chain_profile_v3",
            generator::AbstractDict = Dict{String, Any}(),
            selection::AbstractDict = Dict{String, Any}(),
            source_snapshot::Union{Nothing, AbstractDict} = nothing,
            diagnostics::AbstractDict = Dict{String, Any}(),
            compact_to_source_index::AbstractVector{<:Integer} = Int[],
            original_point_count::Union{Nothing, Integer} = nothing,
            retained_point_count::Union{Nothing, Integer} = nothing)
        compact_to_source = Int.(compact_to_source_index)
        retained_points = retained_point_count === nothing ? length(compact_to_source) : Int(retained_point_count)
        retained_points >= 0 || throw(ArgumentError("ChainMetadata: retained_point_count must be non-negative."))
        if !isempty(compact_to_source) && length(compact_to_source) != retained_points
            throw(ArgumentError(
                "ChainMetadata: compact_to_source_index length must match retained_point_count."
            ))
        end

        original_points = original_point_count === nothing ? nothing : Int(original_point_count)
        if original_points !== nothing && original_points < retained_points
            throw(ArgumentError(
                "ChainMetadata: original_point_count must be >= retained_point_count."
            ))
        end

        generator_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(generator))
        selection_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(selection))
        source_snapshot_dict = source_snapshot === nothing ? nothing :
                               Dict{String, Any}(String(k) => v for (k, v) in pairs(source_snapshot))
        diagnostics_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(diagnostics))

        return new(
            String(schema_version),
            generator_dict,
            selection_dict,
            source_snapshot_dict,
            diagnostics_dict,
            compact_to_source,
            original_points,
            retained_points
        )
    end
end

"""
$(SIGNATURES)

Composition payload for the self-contained chain-profile inlet.
"""
struct ChainProfileInletComposition
    species::Vector{String}
    mole_fractions::Vector{Float64}
    total_number_density_m3::Float64

    function ChainProfileInletComposition(;
            species,
            mole_fractions,
            total_number_density_m3)
        species_vec = String.(species)
        mole_frac_vec = Float64.(mole_fractions)
        n_tot = Float64(total_number_density_m3)

        if length(species_vec) != length(mole_frac_vec)
            throw(ArgumentError(
                "ChainProfileInletComposition: species and mole_fractions arrays must have same length."
            ))
        end
        isempty(species_vec) && throw(ArgumentError(
            "ChainProfileInletComposition: at least one species must be specified."
        ))
        any(mole_frac_vec .< 0) && throw(ArgumentError(
            "ChainProfileInletComposition: mole fractions must be non-negative."
        ))
        abs(sum(mole_frac_vec) - 1.0) > 1e-10 && throw(ArgumentError(
            "ChainProfileInletComposition: mole fractions must sum to 1.0, got $(sum(mole_frac_vec))."
        ))
        n_tot > 0.0 || throw(ArgumentError(
            "ChainProfileInletComposition: total_number_density_m3 must be positive."
        ))
        length(unique(species_vec)) == length(species_vec) || throw(ArgumentError(
            "ChainProfileInletComposition: duplicate species names found."
        ))
        for name in species_vec
            isempty(strip(name)) && throw(ArgumentError(
                "ChainProfileInletComposition: species names cannot be empty."
            ))
        end

        return new(species_vec, mole_frac_vec, n_tot)
    end
end

"""
$(SIGNATURES)

Self-contained inlet state used to initialize segment 1 of the chain solver.
"""
struct ChainProfileInlet
    composition::ChainProfileInletComposition
    thermal::ReactorThermalState
    source_compact_index::Int

    function ChainProfileInlet(;
            composition::ChainProfileInletComposition,
            thermal::ReactorThermalState,
            source_compact_index::Integer)
        source_idx = Int(source_compact_index)
        source_idx >= 1 || throw(ArgumentError(
            "ChainProfileInlet: source_compact_index must be >= 1."
        ))
        return new(composition, thermal, source_idx)
    end
end

"""
$(SIGNATURES)

Per-cell chain result, including cell-aligned profile inputs and reactor time history.
"""
struct ChainCellResult
    compact_cell_index::Int
    source_cell_index::Int
    z_m::Float64
    dx_m::Float64
    te_K::Float64
    species_u_m_s::Dict{String, Float64}
    reactor::ReactorResult
    endpoint_reactor::Union{Nothing, ReactorConfig}
    success::Bool
    message::String

    function ChainCellResult(;
            compact_cell_index::Integer,
            source_cell_index::Integer = compact_cell_index,
            z_m::Real,
            dx_m::Real,
            te_K::Real,
            species_u_m_s::AbstractDict,
            reactor::ReactorResult,
            endpoint_reactor::Union{Nothing, ReactorConfig} = nothing,
            success::Bool = reactor.success,
            message::AbstractString = reactor.message)
        compact_idx = Int(compact_cell_index)
        source_idx = Int(source_cell_index)
        compact_idx >= 1 || throw(ArgumentError("ChainCellResult: compact_cell_index must be >= 1."))
        source_idx >= 1 || throw(ArgumentError("ChainCellResult: source_cell_index must be >= 1."))

        z_val = Float64(z_m)
        dx_val = Float64(dx_m)
        te_val = Float64(te_K)
        species_u_dict = Dict{String, Float64}()
        for (name, value) in pairs(species_u_m_s)
            name_str = String(name)
            isempty(strip(name_str)) && throw(ArgumentError(
                "ChainCellResult: species_u_m_s keys must be non-empty."
            ))
            value_f64 = Float64(value)
            isfinite(value_f64) && value_f64 > 0 || throw(ArgumentError(
                "ChainCellResult: species_u_m_s[$name_str] must be finite and positive."
            ))
            species_u_dict[name_str] = value_f64
        end

        isfinite(z_val) || throw(ArgumentError("ChainCellResult: z_m must be finite."))
        isfinite(dx_val) && dx_val > 0 || throw(ArgumentError("ChainCellResult: dx_m must be finite and positive."))
        isfinite(te_val) && te_val > 0 || throw(ArgumentError("ChainCellResult: te_K must be finite and positive."))
        isempty(species_u_dict) && throw(ArgumentError(
            "ChainCellResult: species_u_m_s must contain at least one species."
        ))

        return new(
            compact_idx,
            source_idx,
            z_val,
            dx_val,
            te_val,
            species_u_dict,
            reactor,
            endpoint_reactor,
            success,
            String(message)
        )
    end
end

"""
$(SIGNATURES)

Normalized axial chain profile used by the TERRA chain-of-CSTR interface.

# Fields
- `z_m::Vector{Float64}`: Axial coordinate points (m), strictly increasing
- `dx_m::Vector{Float64}`: Point-aligned control-volume lengths (m), strictly positive
- `te_K::Vector{Float64}`: Electron temperature profile (K), strictly positive
- `species_u_m_s::Dict{String, Vector{Float64}}`: Per-species convective velocities (m/s), strictly positive
- `inlet::ChainProfileInlet`: Self-contained segment-1 inlet state
- `diagnostics::Dict{String, Vector{Float64}}`: Optional diagnostics arrays
- `generator::Dict{String, Any}`: Generator metadata from interchange artifact
- `selection::Dict{String, Any}`: Selection metadata from interchange artifact
- `schema_version::String`: Interchange schema version
- `source_snapshot::Union{Nothing, Dict{String, Any}}`: Optional source provenance snapshot
"""
struct AxialChainProfile
    z_m::Vector{Float64}
    dx_m::Vector{Float64}
    te_K::Vector{Float64}
    species_u_m_s::Dict{String, Vector{Float64}}
    inlet::ChainProfileInlet
    diagnostics::Dict{String, Vector{Float64}}
    generator::Dict{String, Any}
    selection::Dict{String, Any}
    schema_version::String
    source_snapshot::Union{Nothing, Dict{String, Any}}

    function AxialChainProfile(;
            z_m,
            dx_m,
            te_K,
            species_u_m_s,
            inlet::ChainProfileInlet,
            diagnostics::AbstractDict = Dict{String, Vector{Float64}}(),
            generator::AbstractDict = Dict{String, Any}(),
            selection::AbstractDict = Dict{String, Any}(),
            schema_version::AbstractString = "terra_chain_profile_v3",
            source_snapshot::Union{Nothing, AbstractDict} = nothing)
        z_m_vec = Float64.(z_m)
        dx_m_vec = Float64.(dx_m)
        te_K_vec = Float64.(te_K)

        n = length(z_m_vec)
        if n == 0
            throw(ArgumentError("AxialChainProfile: profile arrays must be non-empty."))
        end
        if length(dx_m_vec) != n || length(te_K_vec) != n
            throw(ArgumentError("AxialChainProfile: all required profile arrays must have identical lengths."))
        end

        for (i, value) in pairs(z_m_vec)
            isfinite(value) || throw(ArgumentError("AxialChainProfile: z_m[$i] must be finite."))
        end
        for i in 2:n
            z_m_vec[i] > z_m_vec[i - 1] || throw(ArgumentError(
                "AxialChainProfile: z_m must be strictly increasing (failed at indices $(i - 1), $i)."
            ))
        end

        for (name, values) in (
            ("dx_m", dx_m_vec),
            ("te_K", te_K_vec),
        )
            for (i, value) in pairs(values)
                isfinite(value) || throw(ArgumentError("AxialChainProfile: $(name)[$i] must be finite."))
                value > 0.0 || throw(ArgumentError("AxialChainProfile: $(name)[$i] must be strictly positive."))
            end
        end

        species_u_dict = Dict{String, Vector{Float64}}()
        for (name, values) in pairs(species_u_m_s)
            name_str = String(name)
            isempty(strip(name_str)) && throw(ArgumentError(
                "AxialChainProfile: species_u_m_s keys must be non-empty."
            ))
            values_vec = Float64.(values)
            if length(values_vec) != n
                throw(ArgumentError(
                    "AxialChainProfile: species_u_m_s[$name_str] length $(length(values_vec)) does not match required profile length $n."
                ))
            end
            for (i, value) in pairs(values_vec)
                isfinite(value) || throw(ArgumentError(
                    "AxialChainProfile: species_u_m_s[$name_str][$i] must be finite."
                ))
                value > 0.0 || throw(ArgumentError(
                    "AxialChainProfile: species_u_m_s[$name_str][$i] must be strictly positive."
                ))
            end
            species_u_dict[name_str] = values_vec
        end
        isempty(species_u_dict) && throw(ArgumentError(
            "AxialChainProfile: species_u_m_s must contain at least one species."
        ))

        diagnostics_dict = Dict{String, Vector{Float64}}()
        for (key, values) in pairs(diagnostics)
            key_str = String(key)
            values_vec = Float64.(values)
            if length(values_vec) != n
                throw(ArgumentError(
                    "AxialChainProfile: diagnostic `$(key_str)` length $(length(values_vec)) does not match required profile length $n."
                ))
            end
            for (i, value) in pairs(values_vec)
                isfinite(value) || throw(ArgumentError("AxialChainProfile: diagnostic `$(key_str)[$i]` must be finite."))
            end
            diagnostics_dict[key_str] = values_vec
        end

        generator_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(generator))
        selection_dict = Dict{String, Any}(String(k) => v for (k, v) in pairs(selection))
        source_snapshot_dict = source_snapshot === nothing ? nothing :
                               Dict{String, Any}(String(k) => v for (k, v) in pairs(source_snapshot))

        return new(
            z_m_vec,
            dx_m_vec,
            te_K_vec,
            species_u_dict,
            inlet,
            diagnostics_dict,
            generator_dict,
            selection_dict,
            String(schema_version),
            source_snapshot_dict
        )
    end
end

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

    function AxialMarchingConfig(;
            handoff_mode::Symbol = :full_state,
            termination_mode::Symbol = :final_time,
            override_tt_K::Union{Nothing, Real} = nothing,
            override_tv_K::Union{Nothing, Real} = nothing,
            is_isothermal_teex::Bool = true)
        handoff_mode in (:reinitialize, :full_state) || throw(ArgumentError(
            "AxialMarchingConfig: handoff_mode must be :reinitialize or :full_state."
        ))
        termination_mode in (:final_time, :steady_state) || throw(ArgumentError(
            "AxialMarchingConfig: termination_mode must be :final_time or :steady_state."
        ))
        if override_tt_K !== nothing && (!isfinite(override_tt_K) || override_tt_K <= 0.0)
            throw(ArgumentError("AxialMarchingConfig: override_tt_K must be finite and positive when provided."))
        end
        if override_tv_K !== nothing && (!isfinite(override_tv_K) || override_tv_K <= 0.0)
            throw(ArgumentError("AxialMarchingConfig: override_tv_K must be finite and positive when provided."))
        end

        return new(
            handoff_mode,
            termination_mode,
            override_tt_K === nothing ? nothing : Float64(override_tt_K),
            override_tv_K === nothing ? nothing : Float64(override_tv_K),
            is_isothermal_teex
        )
    end
end

"""
$(SIGNATURES)

Results container for axial-marching chain-of-CSTR simulations.

# Fields
- `cells::Vector{ChainCellResult}`: Solved chain cells (compact retained indexing)
- `metadata::ChainMetadata`: Chain-level metadata, diagnostics, and index mapping
- `success::Bool`: Overall chain success flag
- `failed_cell::Union{Nothing, Int}`: First failing retained-cell index, if any
- `message::String`: Overall status message
"""
struct ChainSimulationResult
    cells::Vector{ChainCellResult}
    metadata::ChainMetadata
    success::Bool
    failed_cell::Union{Nothing, Int}
    message::String

    function ChainSimulationResult(
            cells::Vector{ChainCellResult},
            metadata::ChainMetadata,
            success::Bool,
            failed_cell::Union{Nothing, Integer},
            message::AbstractString)
        failed_cell_val = failed_cell === nothing ? nothing : Int(failed_cell)
        if failed_cell_val !== nothing &&
           (failed_cell_val < 1 || failed_cell_val > length(cells))
            throw(ArgumentError(
                "ChainSimulationResult: failed_cell must be `nothing` or in 1:length(cells)."
            ))
        end
        return new(cells, metadata, success, failed_cell_val, String(message))
    end
end

function ChainSimulationResult(;
        cells,
        metadata::ChainMetadata = ChainMetadata(
            compact_to_source_index = [cell.source_cell_index for cell in cells],
            retained_point_count = length(cells)),
        success::Bool = all(cell.success for cell in cells),
        failed_cell::Union{Nothing, Integer} = nothing,
        message::AbstractString = success ?
                                  "Chain integration completed successfully." :
                                  "Chain integration failed.")
    cells_vec = ChainCellResult[cell for cell in cells]
    return ChainSimulationResult(cells_vec, metadata, success, failed_cell, message)
end

Base.firstindex(chain::ChainSimulationResult) = 1
Base.lastindex(chain::ChainSimulationResult) = length(chain.cells)

function _slice_chain_metadata(metadata::ChainMetadata, cells::Vector{ChainCellResult})
    return ChainMetadata(;
        schema_version = metadata.schema_version,
        generator = copy(metadata.generator),
        selection = copy(metadata.selection),
        source_snapshot = metadata.source_snapshot === nothing ? nothing : copy(metadata.source_snapshot),
        diagnostics = copy(metadata.diagnostics),
        compact_to_source_index = [cell.source_cell_index for cell in cells],
        original_point_count = metadata.original_point_count,
        retained_point_count = length(cells)
    )
end

function _slice_chain_result(chain::ChainSimulationResult, indices::AbstractVector{<:Integer})
    cells = chain.cells[indices]
    metadata = _slice_chain_metadata(chain.metadata, cells)
    failed_cell = findfirst(!, [cell.success for cell in cells])
    success = failed_cell === nothing
    message = success ? chain.message : cells[failed_cell].message
    return ChainSimulationResult(cells, metadata, success, failed_cell, message)
end

function Base.getindex(chain::ChainSimulationResult, cell::Integer)
    return _slice_chain_result(chain, [cell])
end

function Base.getindex(chain::ChainSimulationResult, cells::AbstractVector{<:Integer})
    return _slice_chain_result(chain, cells)
end
