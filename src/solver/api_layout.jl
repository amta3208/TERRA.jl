"""
# TERRA Solver Module

This module provides the high-level interface for running TERRA simulations
from Julia, hiding the complexity of the Fortran interface and providing
a clean, Julia-native API.
"""

# Debug: RHS call counter
const _TERRA_ODE_DBG_CALLS = Ref(0)

"""
Compact API layout information for the Fortran `rhs_api` / `calculate_rhs_api` state vectors.

This describes the *compact* ordering used internally by TERRA:
`[rho_vib_states, rho_elec_states, rho_species, rho_u..., rho_etot, rho_eeex?, rho_erot?, rho_evib?]`.
"""
struct ApiLayout
    layout_version::Int
    mnsp::Int
    mnex::Int
    mmnex::Int
    mnv::Int
    mneq::Int

    nsp::Int
    nd::Int
    neq::Int
    esp::Int

    get_electron_density_by_charge_balance::Bool
    eex_noneq::Bool
    rot_noneq::Bool
    vib_noneq::Bool
    is_isothermal::Bool
    is_isothermal_teex::Bool
    is_elec_sts::Bool
    is_vib_sts::Bool

    n_eq_vib::Int
    n_eq_elec::Int
    n_eq_sp::Int
    n_eq_mom::Int
    n_eq_energy::Int

    ih::Vector{Int}
    ie::Vector{Int}
    ies::Vector{Int}
    mex::Vector{Int}
    ivs::Matrix{Int}
    mv::Matrix{Int}
    spwt::Vector{Float64}

    vib_range::UnitRange{Int}
    elec_range::UnitRange{Int}
    sp_range::UnitRange{Int}
    mom_range::UnitRange{Int}
    energy_range::UnitRange{Int}

    idx_etot::Int
    idx_eeex::Int
    idx_erot::Int
    idx_evib::Int
end

function ApiLayout(layout::NamedTuple)
    ih = Int.(layout.ih)
    ie = Int.(layout.ie)
    ies = Int.(layout.ies)
    mex = Int.(layout.mex)
    ivs = Int.(layout.ivs)
    mv = Int.(layout.mv)
    spwt = Vector{Float64}(layout.spwt)

    n_eq_vib = Int(layout.n_eq_vib)
    n_eq_elec = Int(layout.n_eq_elec)
    n_eq_sp = Int(layout.n_eq_sp)
    n_eq_mom = Int(layout.n_eq_mom)
    n_eq_energy = Int(layout.n_eq_energy)
    neq = Int(layout.neq)

    vib_start = 1
    vib_stop = n_eq_vib
    elec_start = vib_stop + 1
    elec_stop = vib_stop + n_eq_elec
    sp_start = elec_stop + 1
    sp_stop = elec_stop + n_eq_sp
    mom_start = sp_stop + 1
    mom_stop = sp_stop + n_eq_mom
    energy_start = mom_stop + 1
    energy_stop = neq

    idx_etot = energy_start
    idx_eeex = (layout.eex_noneq == 1) ? (idx_etot + 1) : 0
    idx_erot = (layout.rot_noneq == 1) ? (idx_etot + 1 + (layout.eex_noneq == 1 ? 1 : 0)) :
               0
    idx_evib = (layout.vib_noneq == 1) ?
               (idx_etot + 1 + (layout.eex_noneq == 1 ? 1 : 0) +
                (layout.rot_noneq == 1 ? 1 : 0)) : 0

    return ApiLayout(
        Int(layout.layout_version),
        Int(layout.mnsp),
        Int(layout.mnex),
        Int(layout.mmnex),
        Int(layout.mnv),
        Int(layout.mneq),
        Int(layout.nsp),
        Int(layout.nd),
        neq,
        Int(layout.esp),
        layout.get_electron_density_by_charge_balance == 1,
        layout.eex_noneq == 1,
        layout.rot_noneq == 1,
        layout.vib_noneq == 1,
        layout.is_isothermal == 1,
        layout.is_isothermal_teex == 1,
        layout.is_elec_sts == 1,
        layout.is_vib_sts == 1,
        n_eq_vib,
        n_eq_elec,
        n_eq_sp,
        n_eq_mom,
        n_eq_energy,
        ih,
        ie,
        ies,
        mex,
        ivs,
        mv,
        spwt,
        vib_start:vib_stop,
        elec_start:elec_stop,
        sp_start:sp_stop,
        mom_start:mom_stop,
        energy_start:energy_stop,
        idx_etot,
        idx_eeex,
        idx_erot,
        idx_evib
    )
end

"""
$(SIGNATURES)

Query the Fortran API for the current `y`/`dy` layout.
"""
function get_api_layout()
    return ApiLayout(get_api_layout_wrapper())
end
