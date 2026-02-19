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
    U_neutral::Float64
    U_ion::Float64
    U_energy::Union{Nothing, Float64}
    inlet_reactor::Union{Nothing, ReactorConfig}

    function ResidenceTimeConfig(enabled::Bool, L, U_neutral, U_ion,
            U_energy = nothing, inlet_reactor = nothing)
        if !isfinite(L) || L <= 0
            throw(ArgumentError("ResidenceTimeConfig: L must be finite and positive (got $L)."))
        end
        if !isfinite(U_neutral) || U_neutral <= 0
            throw(ArgumentError("ResidenceTimeConfig: U_neutral must be finite and positive (got $U_neutral)."))
        end
        if !isfinite(U_ion) || U_ion <= 0
            throw(ArgumentError("ResidenceTimeConfig: U_ion must be finite and positive (got $U_ion)."))
        end
        if U_energy !== nothing && (!isfinite(U_energy) || U_energy <= 0)
            throw(ArgumentError("ResidenceTimeConfig: U_energy must be finite and positive (got $U_energy)."))
        end
        inlet_reactor_val = _coerce_residence_time_inlet_reactor(inlet_reactor)

        return new(
            enabled,
            Float64(L),
            Float64(U_neutral),
            Float64(U_ion),
            U_energy === nothing ? nothing : Float64(U_energy),
            inlet_reactor_val)
    end
end

# Backward-compatible positional constructor (enabled defaults to true)
ResidenceTimeConfig(L, U_neutral, U_ion, U_energy = nothing, inlet_config = nothing) =
    ResidenceTimeConfig(true, L, U_neutral, U_ion, U_energy, inlet_config)

function ResidenceTimeConfig(;
        enabled::Bool = true,
        L::Real = 1.0,
        U_neutral::Real = 1.0,
        U_ion::Real = 1.0,
        U_energy::Union{Nothing, Real} = nothing,
        inlet_reactor = nothing,
        inlet_config = nothing)
    if inlet_reactor !== nothing && inlet_config !== nothing
        throw(ArgumentError("ResidenceTimeConfig: provide only one of `inlet_reactor` or legacy `inlet_config`."))
    end
    inlet_value = inlet_reactor === nothing ? inlet_config : inlet_reactor
    return ResidenceTimeConfig(enabled, L, U_neutral, U_ion, U_energy, inlet_value)
end

function _coerce_residence_time_inlet_reactor(inlet)
    if inlet === nothing
        return nothing
    elseif inlet isa ReactorConfig
        return inlet
    else
        throw(ArgumentError("ResidenceTimeConfig: inlet_reactor/inlet_config must be `ReactorConfig`, `Config`, or `nothing`."))
    end
end

"""
$(SIGNATURES)

Numerical controls for the wrapper runtime.
"""
struct NumericsConfig
    time::TimeConfig
    solver::ODESolverConfig
    space::SpaceConfig
    residence_time::Union{Nothing, ResidenceTimeConfig}

    function NumericsConfig(;
            time::TimeConfig,
            solver::ODESolverConfig = ODESolverConfig(),
            space::SpaceConfig = SpaceConfig(),
            residence_time::Union{Nothing, ResidenceTimeConfig} = ResidenceTimeConfig())
        return new(time, solver, space, residence_time)
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
    numerics::NumericsConfig
    runtime::RuntimeConfig

    function Config(;
            reactor::ReactorConfig,
            numerics::NumericsConfig,
            models::ModelConfig = ModelConfig(),
            runtime::RuntimeConfig = RuntimeConfig())
        return new(reactor, models, numerics, runtime)
    end
end

_coerce_residence_time_inlet_reactor(inlet::Config) = inlet.reactor

"""
$(SIGNATURES)

Results container for TERRA simulations.

# Fields
- `time::Vector{Float64}`: Time points
- `species_densities::Matrix{Float64}`: Species densities over time
- `temperatures::NamedTuple`: Temperature evolution
- `total_energy::Vector{Float64}`: Total energy evolution
- `source_terms::Union{NamedTuple, Nothing}`: Source terms (if requested)
- `success::Bool`: Simulation success flag
- `message::String`: Status message
"""
struct TERRAResults
    time::Vector{Float64}
    species_densities::Matrix{Float64}
    temperatures::NamedTuple
    total_energy::Vector{Float64}
    source_terms::Union{NamedTuple, Nothing}
    success::Bool
    message::String
end
