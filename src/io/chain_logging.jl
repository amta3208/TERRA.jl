const CHAIN_RUN_BANNER = "\n" * "="^9 * " TERRA 1D Chain Simulation " * "="^9
const CHAIN_RUN_FOOTER = "="^(length(CHAIN_RUN_BANNER) - 1)

function _prepare_chain_logging(runtime::RuntimeConfig)
    mode = runtime.logging.chain_detail_mode
    mode == :off && return nothing

    if _stream_mode_writes_file(mode)
        _truncate_log_file(_chain_log_path(runtime))
    end
    return nothing
end

function _write_chain_log(runtime::RuntimeConfig, text::AbstractString)
    _ensure_log_dir(runtime)
    _append_log_text(_chain_log_path(runtime), text)
    return nothing
end

function _emit_chain_detail(runtime::RuntimeConfig, detail_text::AbstractString)
    mode = runtime.logging.chain_detail_mode
    mode == :off && return nothing

    if _stream_mode_writes_file(mode)
        _write_chain_log(runtime, detail_text)
    end
    if _stream_mode_writes_console(mode)
        print(detail_text)
        endswith(detail_text, "\n") || println()
    end
    return nothing
end

function _emit_chain_summary(runtime::RuntimeConfig, level::Symbol, message::AbstractString;
                             console::Symbol = :never, fields...)
    _log_run_event(runtime, level, message; console = console, fields...)
    return nothing
end

function _chain_summary_header_text(profile::AxialChainProfile,
                                    marching::AxialMarchingConfig)
    lines = String[CHAIN_RUN_BANNER,
                   @sprintf("segments: %d", length(profile.z_m)),
                   "handoff_mode: $(marching.handoff_mode)",
                   "termination_mode: $(marching.termination_mode)",
                   "is_isothermal_teex: $(marching.is_isothermal_teex)"]
    return join(lines, "\n")
end

function _chain_summary_segment_text(segment_index::Integer, total_segments::Integer)
    return @sprintf("\n===== Segment %d/%d 0D Simulation =====", segment_index,
                    total_segments)
end

function _chain_summary_result_text(result::ChainSimulationResult)
    return string(result.message, "\n", CHAIN_RUN_FOOTER)
end

function _chain_log_scalar(value)
    if value === nothing
        return "null"
    elseif value isa Symbol
        return String(value)
    elseif value isa Bool
        return lowercase(string(value))
    elseif value isa Integer
        return string(value)
    elseif value isa Real
        x = Float64(value)
        if !isfinite(x)
            return string(x)
        elseif x == 0.0
            return "0.0"
        elseif abs(x) >= 1.0e4 || abs(x) < 1.0e-3
            return @sprintf("%.6e", x)
        end
        return @sprintf("%.6f", x)
    elseif value isa AbstractString
        return String(value)
    end
    return repr(value)
end

function _push_chain_log_line!(lines::Vector{String},
                               indent::Integer,
                               key::AbstractString,
                               value)
    push!(lines, string(repeat(" ", indent), key, ": ", _chain_log_scalar(value)))
    return lines
end

function _number_density_to_si(value::Real, unit_system::Symbol)
    number_density = Float64(value)
    return unit_system == :SI ? number_density :
           convert_number_density_cgs_to_si(number_density)
end

function _species_number_density_pairs_m3(composition::ReactorComposition,
                                          unit_system::Symbol)
    total_number_density_m3 = _number_density_to_si(composition.total_number_density,
                                                    unit_system)
    return [name => mole_fraction * total_number_density_m3
            for (name, mole_fraction) in zip(composition.species,
                                             composition.mole_fractions)]
end

function _ordered_chain_velocity_pairs(profile::AxialChainProfile,
                                       segment_index::Integer)
    ordered_species = String[name
                             for name in profile.inlet.composition.species
                             if haskey(profile.species_u_m_s, name)]

    seen = Set(ordered_species)
    for name in sort!(collect(keys(profile.species_u_m_s)))
        name in seen && continue
        push!(ordered_species, name)
    end

    return [name => profile.species_u_m_s[name][segment_index] for name in ordered_species]
end

function _chain_log_header_text(config::Config,
                                profile::AxialChainProfile,
                                marching::AxialMarchingConfig,
                                compact_to_source_index::AbstractVector{<:Integer})
    lines = String[]
    push!(lines, "# chain summary generated at $(_timestamp_string())")
    push!(lines, "chain:")
    _push_chain_log_line!(lines, 2, "case_path", config.runtime.case_path)
    _push_chain_log_line!(lines, 2, "unit_system", config.runtime.unit_system)
    _push_chain_log_line!(lines, 2, "segments", length(profile.z_m))
    _push_chain_log_line!(lines, 2, "chain_detail_mode",
                          config.runtime.logging.chain_detail_mode)
    if _stream_mode_writes_file(config.runtime.logging.chain_detail_mode)
        _push_chain_log_line!(lines, 2, "chain_log_path", _chain_log_path(config.runtime))
    end
    _push_chain_log_line!(lines, 2, "handoff_mode", marching.handoff_mode)
    _push_chain_log_line!(lines, 2, "termination_mode", marching.termination_mode)
    _push_chain_log_line!(lines, 2, "is_isothermal_teex", marching.is_isothermal_teex)
    _push_chain_log_line!(lines, 2, "compact_to_source_index",
                          repr(Int.(compact_to_source_index)))
    _push_chain_log_line!(lines, 2, "species", repr(config.reactor.composition.species))

    push!(lines, "  inlet:")
    _push_chain_log_line!(lines, 4, "source_compact_index",
                          profile.inlet.source_compact_index)
    _push_chain_log_line!(lines, 4, "total_number_density_m3",
                          profile.inlet.composition.total_number_density_m3)

    push!(lines, "    temperatures_K:")
    _push_chain_log_line!(lines, 6, "Tt", profile.inlet.thermal.Tt)
    _push_chain_log_line!(lines, 6, "Tv", profile.inlet.thermal.Tv)
    _push_chain_log_line!(lines, 6, "Tee", profile.inlet.thermal.Tee)
    _push_chain_log_line!(lines, 6, "Te", profile.inlet.thermal.Te)

    push!(lines, "    mole_fractions:")
    for (name, mole_fraction) in zip(profile.inlet.composition.species,
                                     profile.inlet.composition.mole_fractions)
        _push_chain_log_line!(lines, 6, name, mole_fraction)
    end

    push!(lines, "  profile:")
    _push_chain_log_line!(lines, 4, "z_start_m", first(profile.z_m))
    _push_chain_log_line!(lines, 4, "z_end_m", last(profile.z_m))
    _push_chain_log_line!(lines, 4, "retained_point_count", length(profile.z_m))

    return string(join(lines, "\n"), "\n\n")
end

function _chain_log_segment_text(base_config::Config,
                                 profile::AxialChainProfile,
                                 compact_to_source_index::AbstractVector{<:Integer},
                                 segment_index::Integer,
                                 segment_case_path::AbstractString,
                                 result::ReactorResult;
                                 endpoint_reactor::Union{Nothing, ReactorConfig} = nothing,
                                 state_cache_used::Bool = false)
    lines = String[]
    push!(lines, "segment $(segment_index)/$(length(profile.z_m)):")
    _push_chain_log_line!(lines, 2, "source_cell_index",
                          compact_to_source_index[segment_index])
    _push_chain_log_line!(lines, 2, "case_path", segment_case_path)
    _push_chain_log_line!(lines, 2, "z_m", profile.z_m[segment_index])
    _push_chain_log_line!(lines, 2, "dx_m", profile.dx_m[segment_index])
    _push_chain_log_line!(lines, 2, "te_target_K", profile.te_K[segment_index])
    _push_chain_log_line!(lines, 2, "rho_ex_handoff_used", state_cache_used)
    _push_chain_log_line!(lines, 2, "status", result.success ? "success" : "failed")
    _push_chain_log_line!(lines, 2, "message", result.message)

    push!(lines, "  species_u_m_s:")
    for (name, velocity) in _ordered_chain_velocity_pairs(profile, segment_index)
        _push_chain_log_line!(lines, 4, name, velocity)
    end

    if endpoint_reactor !== nothing
        total_number_density_m3 = _number_density_to_si(endpoint_reactor.composition.total_number_density,
                                                        base_config.runtime.unit_system)
        push!(lines, "  endpoint:")
        _push_chain_log_line!(lines, 4, "total_number_density_m3", total_number_density_m3)

        push!(lines, "    temperatures_K:")
        _push_chain_log_line!(lines, 6, "Tt", endpoint_reactor.thermal.Tt)
        _push_chain_log_line!(lines, 6, "Tv", endpoint_reactor.thermal.Tv)
        _push_chain_log_line!(lines, 6, "Tee", endpoint_reactor.thermal.Tee)
        _push_chain_log_line!(lines, 6, "Te", endpoint_reactor.thermal.Te)

        push!(lines, "    number_density_m3:")
        for (name, number_density) in _species_number_density_pairs_m3(endpoint_reactor.composition,
                                                                       base_config.runtime.unit_system)
            _push_chain_log_line!(lines, 6, name, number_density)
        end
    end

    return string(join(lines, "\n"), "\n\n")
end

function _chain_log_result_text(result::ChainSimulationResult)
    lines = String[]
    push!(lines, "result:")
    _push_chain_log_line!(lines, 2, "success", result.success)
    _push_chain_log_line!(lines, 2, "failed_cell", result.failed_cell)
    _push_chain_log_line!(lines, 2, "message", result.message)
    return string(join(lines, "\n"), "\n")
end
