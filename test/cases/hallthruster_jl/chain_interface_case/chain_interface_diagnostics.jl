using HallThruster: HallThruster as het
using CairoMakie: Makie as mk
using Unitful
using LaTeXStrings

include("chain_interface_case.jl")

const DEFAULT_PLOT_DIR = joinpath(@__DIR__, "plots")

function diff1d_uniform(y::AbstractVector; dx::Real)
    n = length(y)
    n >= 2 || throw(ArgumentError("y must have length >= 2"))
    T = promote_type(eltype(y), Float64)
    dy = similar(y, T)
    dy[1] = (y[2] - y[1]) / dx
    @inbounds for i in 2:(n - 1)
        dy[i] = (y[i + 1] - y[i - 1]) / (2 * dx)
    end
    dy[n] = (y[n] - y[n - 1]) / dx
    return dy
end

function configure_theme!()
    mk.update_theme!(mk.Theme(
        font = "Latin Modern Roman",
        fonts = (
            regular = "Latin Modern Roman",
            bold = "Latin Modern Roman",
            italic = "Latin Modern Roman"
        ),
        Axis = (
            xlabelsize = 20,
            ylabelsize = 20,
            xticklabelsize = 18,
            yticklabelsize = 18,
            titlesize = 28
        ),
        Legend = (labelsize = 20,)
    ))
    return nothing
end

function find_closest_index(positions::AbstractVector, target::Real)
    return argmin(abs.(positions .- target))
end

function average_solution(sol; average_start_time = DEFAULT_AVERAGE_START_TIME)
    return het.time_average(sol, average_start_time)
end

function collect_profile_data(
        avg;
        channel_length = avg.config.thruster.geometry.channel_length
)
    z_m = if hasproperty(avg.grid, :cell_centers)
        avg.grid.cell_centers
    elseif avg.grid isa AbstractVector
        avg.grid
    else
        throw(ArgumentError("Unsupported HallThruster grid format: $(typeof(avg.grid))."))
    end

    z_cm = z_m .* 100.0
    channel_len_cm = channel_length .* 100.0
    frame = only(avg.frames)
    u_neutral = frame.neutrals[:N2].u
    n_neutral = frame.neutrals[:N2].n
    u_ion = frame.ions[:N2][1].u
    dx_m = (z_cm[2] - z_cm[1]) / 100.0
    du_ion = diff1d_uniform(u_ion; dx = dx_m)
    n_ion = frame.ions[:N2][1].n

    return (
        z_cm = z_cm,
        channel_len_cm = channel_len_cm,
        u_neutral = u_neutral,
        n_neutral = n_neutral,
        u_ion = u_ion,
        du_ion = du_ion,
        n_ion = n_ion,
        ne = frame.ne,
        Te = frame.Tev,
        pe = frame.pe,
        potential = frame.potential,
        electric_field = frame.E,
        channel_area = frame.channel_area
    )
end

function selected_state_report(
        data;
        target_axial_position_cm::Real = 1.5,
        translational_temperature_K::Real = 750.0
)
    idx = find_closest_index(data.z_cm, target_axial_position_cm)
    z_to_L = target_axial_position_cm / data.channel_len_cm

    target_neutral_velocity = data.u_neutral[idx]
    target_ion_velocity = data.u_ion[idx]
    target_ion_acceleration = data.du_ion[idx]
    target_electron_temp = data.Te[idx]
    target_electron_pressure = data.pe[idx]

    target_neutral_num_density = data.n_neutral[idx]
    target_ion_num_density = data.n_ion[idx]
    target_electron_num_density = target_electron_pressure / 8.617e-5 /
                                  target_electron_temp / 11604.0
    target_total_num_density = target_neutral_num_density + target_ion_num_density +
                               target_electron_num_density

    return (
        z_to_L = z_to_L,
        target_axial_position_cm = target_axial_position_cm,
        target_neutral_velocity = target_neutral_velocity,
        target_ion_velocity = target_ion_velocity,
        target_ion_acceleration = target_ion_acceleration,
        target_electron_temp_K = target_electron_temp * 11604.0,
        target_neutral_num_density = target_neutral_num_density,
        target_ion_num_density = target_ion_num_density,
        target_electron_num_density = target_electron_num_density,
        target_total_num_density = target_total_num_density,
        neutral_mole_fraction = target_neutral_num_density / target_total_num_density,
        ion_mole_fraction = target_ion_num_density / target_total_num_density,
        electron_mole_fraction = target_electron_num_density / target_total_num_density,
        ion_channel_timescale_us = (data.channel_len_cm / 100.0) / target_ion_velocity *
                                   1e6,
        translational_temperature_K = translational_temperature_K
    )
end

function print_selected_state(report)
    println()
    println("================ Plasma State @ z/L=", report.z_to_L, " ================")
    println("Ion Number Density (1/cm^3)      : ", report.target_ion_num_density / 1e6)
    println("Electron Number Density (1/cm^3) : ", report.target_electron_num_density / 1e6)
    println("Neutral Number Density (1/cm^3)  : ", report.target_neutral_num_density / 1e6)
    println("Total Number Density (1/cm^3)    : ", report.target_total_num_density / 1e6)
    println()
    println("N2  Mole Fraction                : ", report.neutral_mole_fraction)
    println("N2+ Mole Fraction                : ", report.ion_mole_fraction)
    println("E-  Mole Fraction                : ", report.electron_mole_fraction)
    println()
    println("Translational Temperature (K)    : ", report.translational_temperature_K)
    println("Electron Temperature (K)         : ", report.target_electron_temp_K)
    println()
    println("N2  Velocity (m/s)               : ", report.target_neutral_velocity)
    println("N2+ Velocity (m/s)               : ", report.target_ion_velocity)
    println("N2+ Acceleration (m/s^2)         : ", report.target_ion_acceleration)
    println("N2+ Channel Timescale (us)       : ", report.ion_channel_timescale_us)
    println(repeat("=", 57))
    return nothing
end

function build_diagnostic_figures(
        data;
        target_axial_position_cm::Real = 1.5
)
    configure_theme!()
    norm_ticks = collect(0:0.5:(maximum(data.z_cm) / data.channel_len_cm))

    figures = Dict{Symbol, mk.Figure}()

    f1 = mk.Figure()
    ax1_1 = mk.Axis(f1[1, 1],
        xlabel = "Axial coordinate [cm]",
        ylabel = "Ion Velocity [m/s]",
        ylabelcolor = mk.wong_colors()[1],
        xgridvisible = false,
        ygridvisible = false)
    ax1_2 = mk.Axis(f1[1, 1],
        yaxisposition = :right,
        ylabel = "Ion Acceleration [m/s^2]",
        ylabelcolor = mk.wong_colors()[2],
        xlabelvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false)
    ax1_3 = mk.Axis(f1[1, 1],
        xaxisposition = :top,
        xlabel = L"z / L_{ch}",
        ylabelvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticks = (0:0.5:2.5))
    mk.linkxaxes!(ax1_1, ax1_2, ax1_3)
    ax1_3.xticks = (norm_ticks .* data.channel_len_cm, string.(norm_ticks))
    mk.lines!(ax1_1, data.z_cm, data.u_ion, color = mk.wong_colors()[1])
    mk.lines!(ax1_2, data.z_cm, data.du_ion, color = mk.wong_colors()[2])
    figures[:ion_kinematics] = f1

    f2 = mk.Figure()
    ax2_1 = mk.Axis(f2[1, 1],
        xlabel = "Axial coordinate [cm]",
        ylabel = "Neutral Number Density [1/m^3]",
        ylabelcolor = mk.wong_colors()[1],
        xgridvisible = false,
        ygridvisible = false)
    ax2_2 = mk.Axis(f2[1, 1],
        yaxisposition = :right,
        ylabel = "Ion Number Density [1/m^3]",
        ylabelcolor = mk.wong_colors()[2],
        xlabelvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false)
    ln2_1 = mk.lines!(ax2_1, data.z_cm, data.n_neutral, label = "N2")
    ln2_2 = mk.lines!(ax2_2, data.z_cm, data.n_ion, label = "N2+", color = mk.Cycled(2))
    mk.Legend(f2[1, 1], [ln2_1, ln2_2], ["N2", "N2+"],
        tellwidth = false, tellheight = false,
        halign = :right, valign = :top,
        margin = (10, 10, 10, 10))
    figures[:species_density] = f2

    f3 = mk.Figure()
    ax3_1 = mk.Axis(f3[1, 1],
        xlabel = "Axial coordinate [cm]",
        ylabel = "Electron Temperature [eV]",
        xgridvisible = false,
        ygridvisible = false)
    ax3_2 = mk.Axis(f3[1, 1],
        xaxisposition = :top,
        xlabel = L"z / L_{ch}",
        ylabelvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticks = (0:0.5:2.5))
    mk.lines!(ax3_1, data.z_cm, data.Te, label = "Te")
    mk.lines!(ax3_2, data.z_cm / data.channel_len_cm, data.Te, label = "Te")
    figures[:electron_temperature] = f3

    f4 = mk.Figure()
    ax4_1 = mk.Axis(f4[1, 1],
        xlabel = "Axial Coordinate [cm]",
        ylabel = "Potential [V]",
        ylabelcolor = mk.wong_colors()[1],
        xgridvisible = false,
        ygridvisible = false)
    ax4_2 = mk.Axis(f4[1, 1],
        yaxisposition = :right,
        ylabel = "Electric Field [V/m]",
        ylabelcolor = mk.wong_colors()[2],
        xlabelvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false)
    ax4_3 = mk.Axis(f4[1, 1],
        xaxisposition = :top,
        xlabel = L"z\;/L_{ch}",
        ylabelvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xminorticksvisible = true,
        xticks = (0:0.5:2.5))
    mk.lines!(ax4_1, data.z_cm, data.potential, color = mk.wong_colors()[1])
    mk.lines!(ax4_2, data.z_cm, data.electric_field, color = mk.wong_colors()[2])
    mk.linkxaxes!(ax4_1, ax4_2, ax4_3)
    ax4_3.xticks = (norm_ticks .* data.channel_len_cm, string.(norm_ticks))
    figures[:potential_and_field] = f4

    f5 = mk.Figure()
    ax5_1 = mk.Axis(f5[1, 1],
        xlabel = "Axial Coordinate [cm]",
        ylabel = "Potential [V]",
        ylabelcolor = mk.wong_colors()[1],
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1.0,
        ytickalign = 1.0,
        yminorticksvisible = true,
        yminortickalign = 1.0,
        xminorticksvisible = true,
        xminortickalign = 1.0)
    ax5_2 = mk.Axis(f5[1, 1],
        yaxisposition = :right,
        ylabel = "Electron Temperature [eV]",
        ylabelcolor = mk.wong_colors()[2],
        xlabelvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1.0,
        ytickalign = 1.0,
        yminorticksvisible = true,
        yminortickalign = 1.0)
    ax5_3 = mk.Axis(f5[1, 1],
        xaxisposition = :top,
        xlabel = L"z\;/L_{ch}",
        ylabelvisible = false,
        yticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xminorticksvisible = true,
        yticksvisible = false,
        xticks = (0:0.5:2.5),
        xtickalign = 1.0,
        ytickalign = 1.0,
        xminortickalign = 1.0)
    mk.linkxaxes!(ax5_1, ax5_2, ax5_3)
    mk.xlims!(ax5_1, minimum(data.z_cm), maximum(data.z_cm))
    ax5_3.xticks = (norm_ticks .* data.channel_len_cm, string.(norm_ticks))
    mk.lines!(ax5_1, data.z_cm, data.potential, color = mk.wong_colors()[1])
    mk.lines!(ax5_2, data.z_cm, data.Te, color = mk.wong_colors()[2])
    mk.vlines!(
        ax5_1, [target_axial_position_cm], color = :black, linestyle = :dash, linewidth = 2)
    mk.vlines!(
        ax5_2, [target_axial_position_cm], color = :black, linestyle = :dash, linewidth = 2)
    mk.vlines!(
        ax5_3, [target_axial_position_cm], color = :black, linestyle = :dash, linewidth = 2)
    figures[:potential_and_temperature] = f5

    return figures
end

function save_figures(figures; plot_dir::AbstractString = DEFAULT_PLOT_DIR)
    mkpath(plot_dir)
    for (name, fig) in figures
        mk.save(joinpath(plot_dir, string(name) * ".png"), fig, dpi = 300)
    end
    return String(plot_dir)
end

function run_diagnostics(;
        average_start_time = DEFAULT_AVERAGE_START_TIME,
        save_time_resolved::Bool = false,
        write_output::Bool = true,
        display_plots::Bool = true,
        save_plots::Bool = false,
        return_results::Bool = false,
        plot_dir::AbstractString = DEFAULT_PLOT_DIR,
        target_axial_position_cm::Real = 1.5,
        config_kwargs = NamedTuple(),
        simparam_kwargs = NamedTuple()
)
    sol = main(
        average_start_time = average_start_time,
        save_time_resolved = save_time_resolved,
        write_output = write_output,
        return_solution = true,
        config_kwargs = config_kwargs,
        simparam_kwargs = simparam_kwargs
    )

    avg = average_solution(sol; average_start_time = average_start_time)
    data = collect_profile_data(avg)
    report = selected_state_report(
        data;
        target_axial_position_cm = target_axial_position_cm,
        translational_temperature_K = haskey(config_kwargs, :temperature_K) ?
                                      config_kwargs.temperature_K : 750.0
    )
    figures = build_diagnostic_figures(
        data; target_axial_position_cm = target_axial_position_cm)

    print_selected_state(report)

    if display_plots
        foreach(display, values(figures))
    end

    saved_plot_dir = nothing
    if save_plots
        saved_plot_dir = save_figures(figures; plot_dir = plot_dir)
        println("plot_dir: ", saved_plot_dir)
    end

    if return_results
        return (
            solution = sol,
            average = avg,
            data = data,
            report = report,
            figures = figures,
            plot_dir = saved_plot_dir
        )
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_diagnostics()
end
