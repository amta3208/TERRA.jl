@testset "RHS (rhs_api)" begin
    test_case_path = TEST_CASE_PATH
    @test_nowarn reset_and_init!(test_case_path)

    config = terra.nitrogen_10ev_config(; isothermal = false)
    state = terra.config_to_initial_state(config)
    layout = terra.get_api_layout()

    y0 = terra.pack_state_vector(layout, state.rho_sp, state.rho_energy;
        rho_ex = state.rho_ex,
        rho_eeex = layout.eex_noneq ? state.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? 0.0 : nothing,
        rho_evib = layout.vib_noneq ? state.rho_evib : nothing)

    du = zeros(length(y0))
    p = (
        layout = layout,
        config = config,
        teex_const = state.teex_const,
        teex_const_vec = fill(config.reactor.thermal.Te, layout.nsp),
        work_y = similar(y0),
        work_dy = similar(y0),
        work_rho_sp = zeros(Float64, layout.nsp),
        work_rho_ex = layout.is_elec_sts ? zeros(Float64, layout.mnex, layout.nsp) :
                      zeros(Float64, 0, 0)
    )

    @test_nowarn terra.terra_ode_system!(du, y0, p, 0.0)
    @test all(isfinite, du)
    @test length(du) == length(y0)

    # If the solver probes intermediate states with small negative components,
    # the wrapper should still compute a meaningful RHS by clamping the state
    # passed to Fortran and applying a Shampine-style positivity rule.
    density_indices = Int[]
    append!(density_indices, collect(layout.vib_range))
    append!(density_indices, collect(layout.elec_range))
    append!(density_indices, collect(layout.sp_range))
    unique!(density_indices)
    @test !isempty(density_indices)
    idx = first(density_indices)

    y_neg = copy(y0)
    y_neg[idx] = -abs(y0[idx]) - 1.0e-30

    density_stop = layout.mom_range.start - 1
    @test density_stop >= 1
    rho_target = sum(@view y_neg[1:density_stop])

    y_proj = copy(y_neg)
    @test y_proj[idx] < 0.0
    y_proj[idx] = 0.0
    rho_clamped = sum(@view y_proj[1:density_stop])
    @test rho_target > 0.0
    @test rho_clamped > 0.0
    scale = rho_target / rho_clamped
    @test scale > 0.0
    @inbounds y_proj[1:density_stop] .*= scale

    du_proj = zeros(length(y0))
    @test_nowarn terra.terra_ode_system!(du_proj, y_proj, p, 0.0)
    @test all(isfinite, du_proj)

    du_neg = zeros(length(y0))
    @test_nowarn terra.terra_ode_system!(du_neg, y_neg, p, 0.0)
    @test all(isfinite, du_neg)
    @test du_neg[idx] >= 0.0

    du_expected = copy(du_proj)
    du_expected[idx] = max(du_expected[idx], 0.0)
    @test du_neg == du_expected
end
