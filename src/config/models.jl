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

    function PhysicsConfig(; bbh_model = 4,
                           esc_model = 1,
                           ar_et_model = 1,
                           eex_noneq = 1,
                           ev_relax_set = 1,
                           et_relax_set = 1,
                           radiation_length = 1.0,
                           get_electron_density_by_charge_balance = true,
                           min_sts_frac = 1e-30,
                           is_isothermal_teex = false,
                           energy_loss_per_eii = 1.0)
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

    function ProcessConfig(; consider_elec_bbe = 1,
                           consider_elec_bfe = 1,
                           consider_elec_bbh = 1,
                           consider_elec_bfh = 1,
                           consider_rad = 0,
                           consider_rdr = 0,
                           consider_chem = 1)
        new(consider_elec_bbe, consider_elec_bfe, consider_elec_bbh,
            consider_elec_bfh, consider_rad, consider_rdr, consider_chem)
    end
end

"""
$(SIGNATURES)

Physics/process model configuration.
"""
struct ModelConfig
    physics::PhysicsConfig
    processes::ProcessConfig
end

function ModelConfig(; physics::PhysicsConfig = PhysicsConfig(),
                     processes::ProcessConfig = ProcessConfig())
    return ModelConfig(physics, processes)
end
