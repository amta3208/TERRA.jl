"""
$(SIGNATURES)

Convert nested `Config` to initial state vectors for TERRA.
"""
function config_to_initial_state(config::Config)
    species = config.reactor.composition.species

    # Get molecular weights
    molecular_weights = get_molecular_weights(species)

    # Ensure we're working in CGS - convert if needed
    config_cgs = config.runtime.unit_system == :CGS ? config : convert_config_units(config, :CGS)
    composition = config_cgs.reactor.composition
    thermal = config_cgs.reactor.thermal

    # Convert mole fractions to mass densities (CGS units)
    mass_densities_cgs = mole_fractions_to_mass_densities(
        composition.mole_fractions, molecular_weights, composition.total_number_density
    )

    gas_constants_full = get_species_gas_constants_wrapper()
    gas_constants = gas_constants_full[1:length(molecular_weights)]

    # Calculate electronic state populations using TERRA's Boltzmann distribution
    # Use appropriate temperatures for each mode
    initial_electronic_states = set_electronic_boltzmann_wrapper(
        mass_densities_cgs,
        # Initialize electronic-state-resolved populations using the
        # electron-electronic temperature (TEE). Species that are not
        # electronically resolved do not contribute here (their mex = 0),
        # so using TEE only affects resolved species as intended.
        thermal.Tee,  # Electronic-state populations use TEE
        thermal.Tt,  # Rotational temperature (use translational as proxy)
        thermal.Tv   # Vibrational temperature
    )

    # Calculate electron-electronic energy using TERRA's method.
    # IMPORTANT: Species without electronically resolved states contribute to
    # the electron-electronic mode; their energy must be initialized using the
    # electron temperature (TE), not TEE. TERRA internally excludes species with
    # resolved electronic states from this mode, so passing TE here is correct.
    initial_electron_electronic_energy = calculate_electron_electronic_energy_wrapper(
        thermal.Te, thermal.Tv, mass_densities_cgs
    )

    # Initialize vibrational state populations using TERRA's Boltzmann
    # distribution (for potential diagnostics), but do NOT use this to infer
    # STS. The active STS setting comes from the database/setup, which the
    # wrapper cannot reliably query here. To avoid a mismatch with TERRA,
    # always treat vibrational energy in mode form and do not include rho_vx
    # in Etot.
    initial_vibrational_states = set_vibrational_boltzmann_wrapper(
        initial_electronic_states,
        thermal.Te,
        thermal.Tt,
        thermal.Tv
    )

    # Mode-level vibrational energy at Tv. Provide per-species electronic
    # temperature proxy as TE uniformly; species with resolved electronic
    # states are already accounted through rho_ex above.
    initial_vibrational_energy = calculate_vibrational_energy_wrapper(
        thermal.Tv, mass_densities_cgs;
        rho_ex = initial_electronic_states,
        tex = fill(thermal.Te, length(mass_densities_cgs))
    )

    # Total energy: do not include rho_vx to avoid double counting when TERRA
    # is not in vibrational STS mode.
    initial_total_energy = calculate_total_energy_wrapper(
        thermal.Tt, mass_densities_cgs;
        rho_ex = initial_electronic_states,
        u = 0.0, v = 0.0, w = 0.0,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy
    )

    has_elec_sts = has_electronic_sts_wrapper()
    has_vib_sts = has_vibrational_sts_wrapper()
    rho_ex_arg = has_elec_sts ? initial_electronic_states : nothing
    rho_vx_arg = has_vib_sts ? initial_vibrational_states : nothing

    initial_temperatures = calculate_temperatures_wrapper(
        mass_densities_cgs, initial_total_energy;
        rho_ex = rho_ex_arg,
        rho_vx = rho_vx_arg,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy)

    initial_enthalpy, initial_pressure = enthalpy_from_energy(
        initial_total_energy,
        mass_densities_cgs,
        gas_constants,
        species,
        molecular_weights,
        initial_temperatures.tt,
        initial_temperatures.teex)

    electron_index = findfirst(
        i -> _is_electron_species(species[i], molecular_weights[i]),
        eachindex(species))
    electron_enthalpy = electron_index === nothing ? 0.0 :
                        mass_densities_cgs[electron_index] * gas_constants[electron_index] *
                        thermal.Te

    initial_remainder_energy = initial_enthalpy - initial_electron_electronic_energy -
                               electron_enthalpy

    return (
        rho_sp = mass_densities_cgs,
        rho_etot = initial_enthalpy,
        rho_energy = initial_total_energy,
        rho_rem = initial_remainder_energy,
        rho_ex = initial_electronic_states,
        rho_vx = nothing,
        rho_eeex = initial_electron_electronic_energy,
        rho_evib = initial_vibrational_energy,
        number_density = composition.total_number_density,
        molecular_weights = molecular_weights,
        teex_const = thermal.Te,
        gas_constants = gas_constants,
        pressure = initial_pressure,
        initial_temperatures = initial_temperatures
    )
end
