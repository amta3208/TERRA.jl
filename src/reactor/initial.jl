"""
$(SIGNATURES)

Concrete 0D initial state prepared from `Config` and the active TERRA API.
"""
struct ReactorInitialState
    rho_sp::Vector{Float64}
    rho_etot::Float64
    rho_energy::Float64
    rho_rem::Float64
    rho_ex::Union{Nothing, Matrix{Float64}}
    rho_vx::Union{Nothing, Array{Float64, 3}}
    rho_eeex::Float64
    rho_evib::Float64
    number_density::Float64
    molecular_weights::Vector{Float64}
    teex_const::Float64
    gas_constants::Vector{Float64}
    pressure::Float64
    initial_temperatures::NamedTuple
    state_cache_used::Bool
end

function ReactorInitialState(; rho_sp,
                             rho_etot,
                             rho_energy,
                             rho_rem,
                             rho_ex::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
                             rho_vx::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
                             rho_eeex,
                             rho_evib,
                             number_density,
                             molecular_weights,
                             teex_const,
                             gas_constants,
                             pressure,
                             initial_temperatures::NamedTuple,
                             state_cache_used::Bool = false)
    rho_sp_vec = Float64.(rho_sp)
    molecular_weights_vec = Float64.(molecular_weights)
    gas_constants_vec = Float64.(gas_constants)
    length(rho_sp_vec) == length(molecular_weights_vec) == length(gas_constants_vec) ||
        throw(ArgumentError("ReactorInitialState: density, molecular-weight, and gas-constant vectors must have identical lengths."))

    rho_ex_val = rho_ex === nothing ? nothing : Matrix{Float64}(rho_ex)
    rho_vx_val = rho_vx === nothing ? nothing : Float64.(rho_vx)

    return ReactorInitialState(rho_sp_vec,
                               Float64(rho_etot),
                               Float64(rho_energy),
                               Float64(rho_rem),
                               rho_ex_val,
                               rho_vx_val,
                               Float64(rho_eeex),
                               Float64(rho_evib),
                               Float64(number_density),
                               molecular_weights_vec,
                               Float64(teex_const),
                               gas_constants_vec,
                               Float64(pressure),
                               initial_temperatures,
                               state_cache_used)
end

function _validated_state_cache(config::Config,
                                state_cache::Union{Nothing, ReactorStateCache})
    state_cache === nothing && return nothing

    species = config.reactor.composition.species
    state_cache.species == species ||
        throw(ArgumentError("ReactorStateCache species ordering must match config species ordering."))
    length(state_cache.rho_sp_cgs) == length(species) ||
        throw(ArgumentError("ReactorStateCache rho_sp length must match config species count."))
    if state_cache.rho_ex_cgs !== nothing &&
       size(state_cache.rho_ex_cgs, 2) < length(species)
        throw(ArgumentError("ReactorStateCache rho_ex must provide one column per active species."))
    end

    return state_cache
end

@inline function _use_state_cache(state_cache::Union{Nothing, ReactorStateCache},
                                  has_electronic_sts::Bool)
    return state_cache !== nothing && has_electronic_sts &&
           state_cache.rho_ex_cgs !== nothing
end

"""
$(SIGNATURES)

Convert nested `Config` to a reactor-owned initial state in CGS units.
"""
function config_to_initial_state(config::Config;
                                 state_cache::Union{Nothing, ReactorStateCache} = nothing)
    species = config.reactor.composition.species
    molecular_weights = get_molecular_weights(species)

    config_cgs = config.runtime.unit_system == :CGS ? config :
                 convert_config_units(config, :CGS)
    composition = config_cgs.reactor.composition
    thermal = config_cgs.reactor.thermal

    state_cache = _validated_state_cache(config, state_cache)
    has_elec_sts = has_electronic_sts_wrapper()
    use_state_cache = _use_state_cache(state_cache, has_elec_sts)

    mass_densities_cgs = if use_state_cache
        copy((state_cache::ReactorStateCache).rho_sp_cgs)
    else
        mole_fractions_to_mass_densities(composition.mole_fractions, molecular_weights,
                                         composition.total_number_density)
    end

    gas_constants_full = get_species_gas_constants_wrapper()
    gas_constants = gas_constants_full[1:length(molecular_weights)]

    initial_electronic_states = if use_state_cache
        copy((state_cache::ReactorStateCache).rho_ex_cgs)
    else
        set_electronic_boltzmann_wrapper(mass_densities_cgs,
                                         thermal.Tee,
                                         thermal.Tt,
                                         thermal.Tv)
    end

    initial_electron_electronic_energy = calculate_electron_electronic_energy_wrapper(thermal.Te,
                                                                                      thermal.Tv,
                                                                                      mass_densities_cgs)

    initial_vibrational_states = set_vibrational_boltzmann_wrapper(initial_electronic_states,
                                                                   thermal.Te,
                                                                   thermal.Tt,
                                                                   thermal.Tv)
    initial_vibrational_energy = calculate_vibrational_energy_wrapper(thermal.Tv,
                                                                      mass_densities_cgs;
                                                                      rho_ex = initial_electronic_states,
                                                                      tex = fill(thermal.Te,
                                                                                 length(mass_densities_cgs)))

    initial_total_energy = calculate_total_energy_wrapper(thermal.Tt, mass_densities_cgs;
                                                          rho_ex = initial_electronic_states,
                                                          u = 0.0,
                                                          v = 0.0,
                                                          w = 0.0,
                                                          rho_eeex = initial_electron_electronic_energy,
                                                          rho_evib = initial_vibrational_energy)

    has_vib_sts = has_vibrational_sts_wrapper()
    rho_ex_arg = has_elec_sts ? initial_electronic_states : nothing
    rho_vx_arg = has_vib_sts ? initial_vibrational_states : nothing

    initial_temperatures = calculate_temperatures_wrapper(mass_densities_cgs,
                                                          initial_total_energy;
                                                          rho_ex = rho_ex_arg,
                                                          rho_vx = rho_vx_arg,
                                                          rho_eeex = initial_electron_electronic_energy,
                                                          rho_evib = initial_vibrational_energy)

    initial_enthalpy, initial_pressure = enthalpy_from_energy(initial_total_energy,
                                                              mass_densities_cgs,
                                                              gas_constants,
                                                              species,
                                                              molecular_weights,
                                                              initial_temperatures.tt,
                                                              initial_temperatures.teex)

    electron_index = findfirst(i -> _is_electron_species(species[i], molecular_weights[i]),
                               eachindex(species))
    electron_enthalpy = electron_index === nothing ? 0.0 :
                        mass_densities_cgs[electron_index] * gas_constants[electron_index] *
                        thermal.Te

    initial_remainder_energy = initial_enthalpy - initial_electron_electronic_energy -
                               electron_enthalpy

    number_density_cgs = use_state_cache ?
                         _total_number_density_from_cgs_density(mass_densities_cgs,
                                                                molecular_weights) :
                         composition.total_number_density

    return ReactorInitialState(; rho_sp = mass_densities_cgs,
                               rho_etot = initial_enthalpy,
                               rho_energy = initial_total_energy,
                               rho_rem = initial_remainder_energy,
                               rho_ex = initial_electronic_states,
                               rho_vx = nothing,
                               rho_eeex = initial_electron_electronic_energy,
                               rho_evib = initial_vibrational_energy,
                               number_density = number_density_cgs,
                               molecular_weights = molecular_weights,
                               teex_const = thermal.Te,
                               gas_constants = gas_constants,
                               pressure = initial_pressure,
                               initial_temperatures = initial_temperatures,
                               state_cache_used = use_state_cache)
end
