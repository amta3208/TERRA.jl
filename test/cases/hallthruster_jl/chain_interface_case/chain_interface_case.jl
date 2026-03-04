using HallThruster: HallThruster as het
using Unitful

const DEFAULT_OUTPUT_FILE = joinpath(@__DIR__, "output", "hallthruster_solution_avg.json")
const DEFAULT_AVERAGE_START_TIME = 0.5u"ms"

function _seconds(value::Unitful.AbstractQuantity)
    return ustrip(u"s", value)
end

function _seconds(value::Real)
    return Float64(value)
end

function default_output_file()
    return DEFAULT_OUTPUT_FILE
end

function build_case_config(;
    channel_length = 2.5u"cm",
    domain_length = 6.0u"cm",
    temperature_K = 750.0,
    discharge_voltage = 300.0u"V",
    flow_rate = 5u"mg/s",
    max_charge = 1,
)
    geom = het.Geometry1D(
        channel_length = channel_length,
        inner_radius = 34.5u"mm",
        outer_radius = 0.05u"m",
    )

    bfield = het.load_magnetic_field(joinpath(@__DIR__, "bfield_spt100.csv"))

    thruster = het.Thruster(
        name = "SPT-100",
        geometry = geom,
        magnetic_field = bfield,
    )

    return het.Config(
        thruster = thruster,
        domain = (0.0u"cm", domain_length),
        discharge_voltage = discharge_voltage,
        propellants = [
            het.Propellant(
                het.MolecularNitrogen,
                flow_rate_kg_s = flow_rate,
                temperature_K = temperature_K,
                ion_temperature_K = temperature_K,
                max_charge = max_charge,
            ),
        ],
    )
end

function build_simparams(;
    n_cells = 100,
    dt = 5u"ns",
    duration = 1u"ms",
    num_save = 1000,
)
    return het.SimParams(
        grid = het.EvenGrid(n_cells),
        dt = dt,
        duration = duration,
        num_save = num_save,
    )
end

function run_case(;
    config_kwargs = NamedTuple(),
    simparam_kwargs = NamedTuple(),
)
    config = build_case_config(; config_kwargs...)
    simparams = build_simparams(; simparam_kwargs...)
    return het.run_simulation(config, simparams)
end

function export_case_json(
    sol;
    file::AbstractString = default_output_file(),
    average_start_time = DEFAULT_AVERAGE_START_TIME,
    save_time_resolved::Bool = false,
)
    mkpath(dirname(file))
    average_start_time_s = _seconds(average_start_time)
    het.write_to_json(
        String(file),
        sol;
        average_start_time = average_start_time_s,
        save_time_resolved = save_time_resolved,
    )
    return String(file)
end

function main(;
    file::AbstractString = default_output_file(),
    average_start_time = DEFAULT_AVERAGE_START_TIME,
    save_time_resolved::Bool = false,
    write_output::Bool = true,
    return_solution::Bool = false,
    config_kwargs = NamedTuple(),
    simparam_kwargs = NamedTuple(),
)
    sol = run_case(
        config_kwargs = config_kwargs,
        simparam_kwargs = simparam_kwargs,
    )

    output_file = nothing
    if write_output
        output_file = export_case_json(
            sol;
            file = file,
            average_start_time = average_start_time,
            save_time_resolved = save_time_resolved,
        )
    end

    average_start_time_s = _seconds(average_start_time)

    println("HallThruster chain interface case completed.")
    println("retcode: :", sol.retcode)
    println("average_start_time_s: ", average_start_time_s)
    println("save_time_resolved: ", save_time_resolved)
    if output_file !== nothing
        println("output_file: ", output_file)
    end

    return return_solution ? sol : nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
