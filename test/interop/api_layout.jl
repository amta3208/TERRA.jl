 @testset "ApiLayout" begin
    reset_and_init!(TEST_CASE_PATH)

    raw = terra.get_api_layout_wrapper()
    layout = terra.get_api_layout()

    @test layout isa terra.ApiLayout
    @test layout.layout_version == raw.layout_version
    @test layout.nsp == raw.nsp
    @test layout.nd == raw.nd
    @test layout.neq == raw.neq
    @test layout.mnsp == raw.mnsp
    @test layout.mnv == raw.mnv
    @test layout.get_electron_density_by_charge_balance isa Bool
    @test layout.eex_noneq isa Bool
    @test layout.rot_noneq isa Bool
    @test layout.vib_noneq isa Bool
    @test layout.is_isothermal isa Bool
    @test layout.is_isothermal_teex isa Bool
    @test layout.is_elec_sts isa Bool
    @test layout.is_vib_sts isa Bool

    @test length(layout.ih) == layout.nsp
    @test length(layout.ie) == layout.nsp
    @test length(layout.ies) == layout.nsp
    @test length(layout.mex) == layout.nsp
    @test length(layout.spwt) == layout.nsp
    @test size(layout.ivs, 2) == layout.nsp
    @test size(layout.mv, 1) == layout.nsp

    @test length(layout.vib_range) == layout.n_eq_vib
    @test length(layout.elec_range) == layout.n_eq_elec
    @test length(layout.sp_range) == layout.n_eq_sp
    @test length(layout.mom_range) == layout.n_eq_mom
    @test length(layout.energy_range) == layout.n_eq_energy
    @test first(layout.energy_range) == layout.idx_etot

    if layout.idx_eeex != 0
        @test layout.idx_eeex in layout.energy_range
    end
    if layout.idx_erot != 0
        @test layout.idx_erot in layout.energy_range
    end
    if layout.idx_evib != 0
        @test layout.idx_evib in layout.energy_range
    end
end
