@testset "State Vector Layout (rhs_api)" begin
    test_case_path = TEST_CASE_PATH
    @test_nowarn reset_and_init!(test_case_path)

    layout = terra.get_api_layout()
    @test layout.neq > 0
    @test layout.neq ==
          layout.n_eq_vib + layout.n_eq_elec + layout.n_eq_sp + layout.n_eq_mom +
          layout.n_eq_energy

    nsp = layout.nsp
    rho_sp = collect(range(1.0, length = nsp))

    rho_ex = nothing
    if layout.is_elec_sts
        rho_ex = zeros(Float64, layout.mnex, nsp)
        @inbounds for isp in 1:nsp
            if layout.ies[isp] == 0
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 1
                    continue
                end
                rho_ex[iex, isp] = 1.0e-3 * (100 * isp + iex)
            end
        end
    end

    rho_vx = nothing
    if layout.n_eq_vib > 0
        rho_vx = zeros(Float64, layout.mnv + 1, layout.mmnex, nsp)
        @inbounds for isp in 1:nsp
            if layout.ies[isp] == 0 || layout.ih[isp] != 2
                continue
            end
            for iex in 1:layout.mex[isp]
                if layout.ivs[iex, isp] == 0
                    continue
                end
                mv_isp_iex = layout.mv[isp, iex]
                for ivx in 0:mv_isp_iex
                    rho_vx[ivx + 1, iex, isp] = 1.0e-6 * (10000 * isp + 100 * iex + ivx)
                end
            end
        end
    end

    rho_etot = 123.0
    rho_eeex = layout.eex_noneq ? 7.0 : nothing
    rho_erot = layout.rot_noneq ? 5.0 : nothing
    rho_evib = layout.vib_noneq ? 9.0 : nothing

    y = terra.pack_state_vector(layout, rho_sp, rho_etot;
        rho_ex = rho_ex,
        rho_vx = rho_vx,
        rho_u = layout.nd >= 1 ? 0.0 : nothing,
        rho_v = layout.nd >= 2 ? 0.0 : nothing,
        rho_w = layout.nd >= 3 ? 0.0 : nothing,
        rho_eeex = rho_eeex,
        rho_erot = rho_erot,
        rho_evib = rho_evib)
    @test length(y) == layout.neq

    st = terra.unpack_state_vector(y, layout)
    y2 = terra.pack_state_vector(layout, st.rho_sp, st.rho_etot;
        rho_ex = st.rho_ex,
        rho_vx = st.rho_vx,
        rho_u = layout.nd >= 1 ? st.rho_u : nothing,
        rho_v = layout.nd >= 2 ? st.rho_v : nothing,
        rho_w = layout.nd >= 3 ? st.rho_w : nothing,
        rho_eeex = layout.eex_noneq ? st.rho_eeex : nothing,
        rho_erot = layout.rot_noneq ? st.rho_erot : nothing,
        rho_evib = layout.vib_noneq ? st.rho_evib : nothing)
    @test y2 ≈ y
end
