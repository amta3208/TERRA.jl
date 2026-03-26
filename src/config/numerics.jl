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
        return new(Float64(dt), Float64(dt_output), Float64(duration), Int(nstep),
                   Int(method))
    end
end

function TimeConfig(; dt, dt_output, duration, nstep::Integer = 500000, method::Integer = 2)
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

    function ODESolverConfig(; reltol::Real = 1e-8,
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

        return new(Float64(reltol),
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

Solver and discretization controls for the wrapper runtime.
"""
struct NumericsConfig
    time::TimeConfig
    solver::ODESolverConfig
    space::SpaceConfig

    function NumericsConfig(; time::TimeConfig,
                            solver::ODESolverConfig = ODESolverConfig(),
                            space::SpaceConfig = SpaceConfig())
        return new(time, solver, space)
    end
end
