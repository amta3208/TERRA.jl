# Unit conversion constants
const KG_M3_TO_G_CM3 = 1e-3  # kg/m³ to g/cm³
const G_CM3_TO_KG_M3 = 1e3   # g/cm³ to kg/m³
const J_TO_ERG = 1e7         # J to erg
const ERG_TO_J = 1e-7        # erg to J
const PA_TO_DYNE_CM2 = 10.0  # Pa to dyne/cm²
const DYNE_CM2_TO_PA = 0.1   # dyne/cm² to Pa
const J_M3_TO_ERG_CM3 = 10.0  # J/m³ to erg/cm³
const ERG_CM3_TO_J_M3 = 0.1   # erg/cm³ to J/m³

"""
$(SIGNATURES)

Convert species densities from SI (kg/m³) to CGS (g/cm³) units.
"""
function convert_density_si_to_cgs(rho_si::AbstractVector{<:Real})
    return rho_si .* KG_M3_TO_G_CM3
end

"""
$(SIGNATURES)

Convert species densities from CGS (g/cm³) to SI (kg/m³) units.
"""
function convert_density_cgs_to_si(rho_cgs::AbstractVector{<:Real})
    return rho_cgs .* G_CM3_TO_KG_M3
end

"""
$(SIGNATURES)

Convert energy density from SI (J/m³) to CGS (erg/cm³) units.
"""
function convert_energy_density_si_to_cgs(energy_si::Real)
    return energy_si * J_M3_TO_ERG_CM3
end

"""
$(SIGNATURES)

Convert energy density from CGS (erg/cm³) to SI (J/m³) units.
"""
function convert_energy_density_cgs_to_si(energy_cgs::Real)
    return energy_cgs * ERG_CM3_TO_J_M3
end

"""
$(SIGNATURES)

Convert pressure from SI (Pa) to CGS (dyne/cm²) units.
"""
function convert_pressure_si_to_cgs(pressure_si::Real)
    return pressure_si * PA_TO_DYNE_CM2
end

"""
$(SIGNATURES)

Convert pressure from CGS (dyne/cm²) to SI (Pa) units.
"""
function convert_pressure_cgs_to_si(pressure_cgs::Real)
    return pressure_cgs * DYNE_CM2_TO_PA
end

"""
$(SIGNATURES)

Convert number density from SI (1/m³) to CGS (1/cm³) units.
"""
function convert_number_density_si_to_cgs(n_si::Real)
    return n_si * 1e-6  # 1/m³ to 1/cm³
end

"""
$(SIGNATURES)

Convert number density from CGS (1/cm³) to SI (1/m³) units.
"""
function convert_number_density_cgs_to_si(n_cgs::Real)
    return n_cgs * 1e6  # 1/cm³ to 1/m³
end

"""
$(SIGNATURES)

Convert a complete state vector from SI to CGS units for TERRA input.

# Arguments
- `rho_sp_si::Vector{Float64}`: Species densities in SI units (kg/m³)
- `rho_etot_si::Float64`: Total energy density in SI units (J/m³)
- `number_density_si::Float64`: Total number density in SI units (1/m³)
- Additional optional energy components in SI units

# Returns
- Named tuple with all quantities converted to CGS units
"""
function convert_state_si_to_cgs(rho_sp_si::AbstractVector{<:Real},
        rho_etot_si::Real,
        number_density_si::Real;
        rho_erot_si::Union{Real, Nothing} = nothing,
        rho_eeex_si::Union{Real, Nothing} = nothing,
        rho_evib_si::Union{Real, Nothing} = nothing)
    rho_sp_cgs = convert_density_si_to_cgs(rho_sp_si)
    rho_etot_cgs = convert_energy_density_si_to_cgs(rho_etot_si)
    number_density_cgs = convert_number_density_si_to_cgs(number_density_si)

    rho_erot_cgs = rho_erot_si !== nothing ? convert_energy_density_si_to_cgs(rho_erot_si) :
                   nothing
    rho_eeex_cgs = rho_eeex_si !== nothing ? convert_energy_density_si_to_cgs(rho_eeex_si) :
                   nothing
    rho_evib_cgs = rho_evib_si !== nothing ? convert_energy_density_si_to_cgs(rho_evib_si) :
                   nothing

    return (rho_sp = rho_sp_cgs,
        rho_etot = rho_etot_cgs,
        number_density = number_density_cgs,
        rho_erot = rho_erot_cgs,
        rho_eeex = rho_eeex_cgs,
        rho_evib = rho_evib_cgs)
end

"""
$(SIGNATURES)

Convert a complete state vector from CGS to SI units from TERRA output.

# Arguments
- `rho_sp_cgs::Vector{Float64}`: Species densities in CGS units (g/cm³)
- `rho_etot_cgs::Float64`: Total energy density in CGS units (erg/cm³)
- Additional optional energy components in CGS units

# Returns
- Named tuple with all quantities converted to SI units
"""
function convert_state_cgs_to_si(rho_sp_cgs::AbstractVector{<:Real},
        rho_etot_cgs::Real;
        rho_erot_cgs::Union{Real, Nothing} = nothing,
        rho_eeex_cgs::Union{Real, Nothing} = nothing,
        rho_evib_cgs::Union{Real, Nothing} = nothing)
    rho_sp_si = convert_density_cgs_to_si(rho_sp_cgs)
    rho_etot_si = convert_energy_density_cgs_to_si(rho_etot_cgs)

    rho_erot_si = rho_erot_cgs !== nothing ?
                  convert_energy_density_cgs_to_si(rho_erot_cgs) : nothing
    rho_eeex_si = rho_eeex_cgs !== nothing ?
                  convert_energy_density_cgs_to_si(rho_eeex_cgs) : nothing
    rho_evib_si = rho_evib_cgs !== nothing ?
                  convert_energy_density_cgs_to_si(rho_evib_cgs) : nothing

    return (rho_sp = rho_sp_si,
        rho_etot = rho_etot_si,
        rho_erot = rho_erot_si,
        rho_eeex = rho_eeex_si,
        rho_evib = rho_evib_si)
end

"""
$(SIGNATURES)

Convert source terms from CGS to SI units.

# Arguments
- `drho_sp_cgs::Vector{Float64}`: Species source terms in CGS units (g/cm³/s)
- `drho_etot_cgs::Float64`: Energy source term in CGS units (erg/cm³/s)
- Additional optional source terms in CGS units

# Returns
- Named tuple with all source terms converted to SI units
"""
function convert_sources_cgs_to_si(drho_sp_cgs::AbstractVector{<:Real},
        drho_etot_cgs::Real;
        drho_erot_cgs::Union{Real, Nothing} = nothing,
        drho_eeex_cgs::Union{Real, Nothing} = nothing,
        drho_evib_cgs::Union{Real, Nothing} = nothing)

    # Source terms have units of [quantity]/time, so same conversion as state
    drho_sp_si = convert_density_cgs_to_si(drho_sp_cgs)
    drho_etot_si = convert_energy_density_cgs_to_si(drho_etot_cgs)

    drho_erot_si = drho_erot_cgs !== nothing ?
                   convert_energy_density_cgs_to_si(drho_erot_cgs) : nothing
    drho_eeex_si = drho_eeex_cgs !== nothing ?
                   convert_energy_density_cgs_to_si(drho_eeex_cgs) : nothing
    drho_evib_si = drho_evib_cgs !== nothing ?
                   convert_energy_density_cgs_to_si(drho_evib_cgs) : nothing

    return (drho_sp = drho_sp_si,
        drho_etot = drho_etot_si,
        drho_erot = drho_erot_si,
        drho_eeex = drho_eeex_si,
        drho_evib = drho_evib_si)
end
