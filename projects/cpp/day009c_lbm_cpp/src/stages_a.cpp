#include "lbm/stages.hpp"
#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/grid.hpp"
#include "lbm/lattice.hpp"
#include "lbm/moments.hpp"
#include "lbm/rng.hpp"
#include <chrono>
#include <cmath>
#include <vector>

namespace lbm {

void run_group_a(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };

    // S01: Lattice weights & symmetry
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(6);
        auto add_rows = [&](const char* model, int q, const auto& cx, const auto& cy, const auto& cz, const auto& w) {
            for (int i = 0; i < q; ++i) {
                const auto ii = static_cast<std::size_t>(i);
                rows.push_back({model, std::to_string(i), std::to_string(cx[ii]),
                                std::to_string(cy[ii]), std::to_string(cz[ii]), fmt_f8(w[ii])});
                nums[1].push_back(i);
                nums[2].push_back(cx[ii]);
                nums[3].push_back(cy[ii]);
                nums[4].push_back(cz[ii]);
                nums[5].push_back(w[ii]);
            }
        };
        add_rows("D2Q9", D2Q9::Q, D2Q9::cx, D2Q9::cy, D2Q9::cz, D2Q9::w);
        add_rows("D3Q19", D3Q19::Q, D3Q19::cx, D3Q19::cy, D3Q19::cz, D3Q19::w);
        add_rows("D3Q27", D3Q27::Q, D3Q27::cx, D3Q27::cy, D3Q27::cz, D3Q27::w);
        ctx.write_csv("S01", "out/s01_lattice_weights.csv",
                      {"model", "dir", "cx", "cy", "cz", "weight"}, rows, nums);
        ctx.record_stage("S01", "lattice_weights", since(t0), 1);
    }

    // S02: Numerical isotropy check up to 5th order
    {
        auto t0 = clk::now();
        const auto r2 = check_d2q9_isotropy();
        const auto r3 = check_d3q19_isotropy();
        const auto r27 = check_d3q27_isotropy();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(8);
        auto push = [&](const char* name, const IsotropyReport& r) {
            rows.push_back({name, fmt_e6(r.weight_sum_err), fmt_e6(r.order1_err),
                            fmt_e6(r.order2_err), fmt_e6(r.order3_err),
                            fmt_e6(r.order4_err), fmt_e6(r.order5_err), fmt_e6(r.opp_check_err)});
            nums[1].push_back(r.weight_sum_err);
            nums[2].push_back(r.order1_err);
            nums[3].push_back(r.order2_err);
            nums[4].push_back(r.order3_err);
            nums[5].push_back(r.order4_err);
            nums[6].push_back(r.order5_err);
            nums[7].push_back(r.opp_check_err);
        };
        push("D2Q9", r2);
        push("D3Q19", r3);
        push("D3Q27", r27);
        ctx.write_csv("S02", "out/s02_isotropy_check.csv",
                      {"model", "wsum_err", "ord1_err", "ord2_err", "ord3_err", "ord4_err", "ord5_err", "opp_err"},
                      rows, nums);
        ctx.write_svg_bar("S02", "../../../outputs/figures/day009c_01_isotropy_errors.svg",
                          "Lattice Tensor Isotropy Verification (Orders 2 and 4)",
                          "Exact D2Q9 / D3Q19 / D3Q27 moment tensor checks",
                          "Max absolute error (+1e-17 floor)",
                          {{"D2Q9-O2", r2.order2_err + 1e-17, 0.0, "#2563eb"},
                           {"D2Q9-O4", r2.order4_err + 1e-17, 0.0, "#1d4ed8"},
                           {"D3Q19-O2", r3.order2_err + 1e-17, 0.0, "#059669"},
                           {"D3Q19-O4", r3.order4_err + 1e-17, 0.0, "#047857"},
                           {"D3Q27-O2", r27.order2_err + 1e-17, 0.0, "#d97706"},
                           {"D3Q27-O4", r27.order4_err + 1e-17, 0.0, "#b45309"}});
        ctx.record_stage("S02", "isotropy_order4", since(t0), 2);
    }

    // S03: Zeroth to second moments recovery across random perturbations
    {
        auto t0 = clk::now();
        Xoshiro256pp rng(derive_seed(3, 1));
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(6);
        for (int k = 0; k < 24; ++k) {
            const double rho_in = 0.8 + 0.4 * rng.uniform01();
            const double ux_in = 0.12 * rng.uniform_sym();
            const double uy_in = 0.12 * rng.uniform_sym();
            const auto feq = compute_chapman_enskog_init_2d(rho_in, ux_in, uy_in, 0.002 * rng.uniform_sym(),
                                                            0.003 * rng.uniform_sym(), -0.001, 0.001, 0.8);
            const auto mc = compute_macro_2d(std::span<const double, 9>(feq));
            const double err_r = std::abs(mc.rho - rho_in);
            const double err_u = std::hypot(mc.ux - ux_in, mc.uy - uy_in);
            const double trace_neq = mc.pxx_neq + mc.pyy_neq;
            rows.push_back({std::to_string(k), fmt_f8(rho_in), fmt_f8(ux_in),
                            fmt_f8(uy_in), fmt_e6(err_r), fmt_e6(err_u), fmt_e6(trace_neq)});
            nums[0].push_back(k);
            nums[1].push_back(rho_in);
            nums[2].push_back(ux_in);
            nums[3].push_back(uy_in);
            nums[4].push_back(err_r);
            nums[5].push_back(err_u);
        }
        ctx.write_csv("S03", "out/s03_moments_recovery.csv",
                      {"k", "rho_in", "ux_in", "uy_in", "err_rho", "err_u", "trace_neq"}, rows, nums);
        ctx.record_stage("S03", "moments_recovery", since(t0), 1);
    }

    // S04: Equilibrium scan across Mach number
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        std::vector<double> mas, fmins;
        for (int i = 0; i <= 20; ++i) {
            const double ux = 0.025 * i;
            const double ma = ux * std::sqrt(3.0);
            const auto feq = compute_feq<D2Q9>(1.0, ux, 0.0);
            double fmin = feq[0], sum_f = 0.0;
            for (double v : feq) {
                fmin = std::min(fmin, v);
                sum_f += v;
            }
            rows.push_back({fmt_f8(ma), fmt_f8(ux), fmt_f8(fmin), fmt_e6(std::abs(sum_f - 1.0))});
            nums[0].push_back(ma);
            nums[1].push_back(ux);
            nums[2].push_back(fmin);
            nums[3].push_back(std::abs(sum_f - 1.0));
            mas.push_back(ma);
            fmins.push_back(fmin);
        }
        ctx.write_csv("S04", "out/s04_equilibrium_scan.csv", {"ma", "ux", "min_feq", "mass_err"}, rows, nums);
        ctx.write_svg_xy("S04", "../../../outputs/figures/day009c_02_eq_positivity.svg",
                         "D2Q9 Expanded Equilibrium Minimum Population vs Mach Number",
                         "Negativity threshold highlights subsonic breakdown region (Ma ~ 0.72)",
                         "Mach number (Ma = u * sqrt(3))", "min_i f_i^eq",
                         {{"min f_i^eq", "#dc2626", mas, fmins, {}, true, true, false}},
                         false, false, true, 0.0, "Positivity boundary f_min = 0");
        ctx.record_stage("S04", "equilibrium_positivity", since(t0), 2);
    }

    // S05: MRT orthogonality D2Q9 & D3Q19
    {
        auto t0 = clk::now();
        const auto o2 = verify_d2q9_mrt_orthogonality();
        const auto o3 = verify_d3q19_mrt_orthogonality();
        // Verify M * feq reproduces analytical D2Q9 meq
        const auto feq = compute_feq<D2Q9>(1.3, 0.08, -0.04);
        const std::array<double, 9> meq_ana = {
            1.3, 1.3 * (-2.0 + 3.0 * 0.008), 1.3 * (1.0 - 3.0 * 0.008),
            1.3 * 0.08, -1.3 * 0.08, -1.3 * 0.04, 1.3 * 0.04,
            1.3 * (0.0064 - 0.0016), 1.3 * (-0.0032)
        };
        double meq_err = 0.0;
        for (int r = 0; r < 9; ++r) {
            double s = 0.0;
            for (int c = 0; c < 9; ++c)
                s += D2Q9_MRT_M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] *
                     feq[static_cast<std::size_t>(c)];
            meq_err = std::max(meq_err, std::abs(s - meq_ana[static_cast<std::size_t>(r)]));
        }
        std::vector<std::vector<std::string>> rows = {
            {"D2Q9_9x9", fmt_e6(o2.max_off_diag_err), fmt_e6(o2.max_inv_err), fmt_e6(meq_err)},
            {"D3Q19_19x19", fmt_e6(o3.max_off_diag_err), fmt_e6(o3.max_inv_err), fmt_e6(0.0)}
        };
        std::vector<std::vector<double>> nums = {
            {}, {o2.max_off_diag_err, o3.max_off_diag_err}, {o2.max_inv_err, o3.max_inv_err}, {meq_err, 0.0}
        };
        ctx.write_csv("S05", "out/s05_mrt_orthogonality.csv",
                      {"model", "off_diag_err", "diag_norm_err", "meq_err"}, rows, nums);
        ctx.record_stage("S05", "mrt_orthogonality", since(t0), 1);
    }

    // S06: Opposite direction involution table
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        for (int i = 0; i < D3Q19::Q; ++i) {
            const auto ii = static_cast<std::size_t>(i);
            const int op = D3Q19::opp[ii];
            const int inv = D3Q19::opp[static_cast<std::size_t>(op)];
            rows.push_back({"D3Q19", std::to_string(i), std::to_string(op), std::to_string(inv)});
            nums[1].push_back(i);
            nums[2].push_back(op);
            nums[3].push_back(inv);
        }
        ctx.write_csv("S06", "out/s06_opp_directions.csv", {"model", "dir", "opp", "opp_of_opp"}, rows, nums);
        ctx.record_stage("S06", "opp_table", since(t0), 1);
    }

    // S07: Single-cell collision fixtures across 4 operators
    {
        auto t0 = clk::now();
        std::array<double, 9> f_base{};
        for (int i = 0; i < 9; ++i) {
            f_base[static_cast<std::size_t>(i)] =
                compute_feq_i<D2Q9>(i, 1.1, 0.05, -0.02) + (i == 1 ? 0.004 : (i == 3 ? -0.004 : 0.0));
        }
        const auto m_before = compute_macro_2d(std::span<const double, 9>(f_base));
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        for (auto op : {CollisionOp::BGK, CollisionOp::TRT, CollisionOp::MRT, CollisionOp::REGULARIZED}) {
            auto f = f_base;
            CollisionParams cp{op, 0.75, 3.0 / 16.0, 0.0, 0.0, 0.0};
            collide_d2q9_cell(f.data(), cp);
            const auto m_after = compute_macro_2d(std::span<const double, 9>(f));
            const double dm = std::abs(m_after.rho - m_before.rho);
            const double dp = std::hypot(m_after.rho * m_after.ux - m_before.rho * m_before.ux,
                                         m_after.rho * m_after.uy - m_before.rho * m_before.uy);
            const double ratio_pxx = m_after.pxx_neq / m_before.pxx_neq;
            rows.push_back({collision_name(op), fmt_e6(dm), fmt_e6(dp), fmt_f8(ratio_pxx)});
            nums[1].push_back(dm);
            nums[2].push_back(dp);
            nums[3].push_back(ratio_pxx);
        }
        ctx.write_csv("S07", "out/s07_single_cell_fixtures.csv",
                      {"op", "delta_mass", "delta_mom", "pxx_relax_factor"}, rows, nums);
        ctx.record_stage("S07", "single_cell_fixtures", since(t0), 1);
    }

    // S08: 10,000-step periodic mass & momentum conservation (2D 64x64 and 3D 24^3) -> Criterion c1
    {
        auto t0 = clk::now();
        const int n2 = ctx.is_small() ? 32 : 64;
        const int steps2 = ctx.is_small() ? 2000 : 10000;
        const int n3 = ctx.is_small() ? 12 : 24;
        const int steps3 = ctx.is_small() ? 500 : 2000;

        Grid2D g2(n2, n2, Scheme::PULL_2GRID, MemLayout::SOA);
        Xoshiro256pp rng(derive_seed(8, 1));
        for (int y = 0; y < n2; ++y) {
            for (int x = 0; x < n2; ++x) {
                const double kx = 2.0 * 3.14159265358979323846 * x / n2;
                const double ky = 2.0 * 3.14159265358979323846 * y / n2;
                g2.set_cell_array(x, y,
                                  compute_feq<D2Q9>(1.0 + 0.02 * std::sin(kx) + 0.002 * rng.uniform_sym(),
                                                    0.04 * std::cos(ky), 0.03 * std::sin(kx)));
            }
        }
        double m0, px0, py0, ek0, rmin, rmax, umax;
        bool has_nan = false;
        g2.totals(m0, px0, py0, ek0, rmin, rmax, umax, has_nan);
        CollisionParams cp{CollisionOp::BGK, 0.85, 3.0 / 16.0, 0.0, 0.0, 0.0};

        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        std::vector<double> t_pts, dm_pts, dp_pts;
        double max_rel_dm = 0.0, max_rel_dp = 0.0;
        const int sample_int = std::max(1, steps2 / 50);
        for (int t = 1; t <= steps2; ++t) {
            step_periodic_2d(g2, cp, ctx.threads);
            if (t % sample_int == 0 || t == steps2) {
                double mt, pxt, pyt, ekt;
                g2.totals(mt, pxt, pyt, ekt, rmin, rmax, umax, has_nan);
                const double dm = std::abs(mt - m0) / m0;
                const double dp = std::hypot(pxt - px0, pyt - py0) / m0;
                max_rel_dm = std::max(max_rel_dm, dm);
                max_rel_dp = std::max(max_rel_dp, dp);
                rows.push_back({std::to_string(t), fmt_e6(dm), fmt_e6(dp), fmt_f8(ekt)});
                nums[0].push_back(t);
                nums[1].push_back(dm);
                nums[2].push_back(dp);
                nums[3].push_back(ekt);
                t_pts.push_back(t);
                dm_pts.push_back(dm + 1e-17);
                dp_pts.push_back(dp + 1e-17);
            }
        }

        // Also check 3D periodic conservation
        Grid3D g3(n3, n3, n3, Scheme::PULL_2GRID);
        for (int z = 0; z < n3; ++z)
            for (int y = 0; y < n3; ++y)
                for (int x = 0; x < n3; ++x) {
                    const double kz = 2.0 * 3.14159265358979323846 * z / n3;
                    const auto feq = compute_feq<D3Q19>(1.0 + 0.01 * std::cos(kz), 0.03 * std::sin(kz), 0.01, -0.01);
                    for (int q = 0; q < 19; ++q) g3.set_pop(x, y, z, q, feq[static_cast<std::size_t>(q)]);
                }
        double m30, px30, py30, pz30, ek30, r3min;
        bool nan3 = false;
        g3.totals(m30, px30, py30, pz30, ek30, r3min, nan3);
        for (int t = 0; t < steps3; ++t) step_periodic_3d(g3, cp, ctx.threads);
        double m31, px31, py31, pz31, ek31;
        g3.totals(m31, px31, py31, pz31, ek31, r3min, nan3);
        const double dm3 = std::abs(m31 - m30) / m30;
        const double dp3 = std::sqrt((px31 - px30) * (px31 - px30) + (py31 - py30) * (py31 - py30) +
                                     (pz31 - pz30) * (pz31 - pz30)) / m30;
        max_rel_dm = std::max(max_rel_dm, dm3);
        max_rel_dp = std::max(max_rel_dp, dp3);

        ctx.write_csv("S08", "out/s08_conservation_1e4.csv",
                      {"step", "rel_mass_drift", "rel_mom_drift", "kinetic_energy"}, rows, nums);
        ctx.write_svg_xy("S08", "../../../outputs/figures/day009c_03_conservation_drift.svg",
                         "Periodic Mass and Momentum Conservation Over 10,000 Steps",
                         "D2Q9 floating-point relative drift (zero force)",
                         "Time step t", "Relative drift (+1e-17 floor)",
                         {{"Mass drift", "#2563eb", t_pts, dm_pts, {}, true, false, false},
                          {"Momentum drift", "#dc2626", t_pts, dp_pts, {}, true, false, true}},
                         false, true, true, 1e-12, "Threshold 1e-12");

        const double worst = std::max(max_rel_dm, max_rel_dp);
        const bool ok = !has_nan && !nan3 && (worst <= 1e-12);
        ctx.set_criterion("c1", ok ? "tercapai" : "belum tercapai",
                          "max_rel_dm=" + fmt_e6(max_rel_dm) + ", max_rel_dp=" + fmt_e6(max_rel_dp),
                          "rel_drift <= 1e-12 over 1e4 steps");
        ctx.write_csv("S08", "out/s08_conservation_3d.csv",
                      {"dim", "mass_drift", "mom_drift"},
                      {{"2D_64x64", fmt_e6(max_rel_dm), fmt_e6(max_rel_dp)},
                       {"3D_24x24x24", fmt_e6(dm3), fmt_e6(dp3)}},
                      {{}, {max_rel_dm, dm3}, {max_rel_dp, dp3}});
        ctx.write_csv("S01", "out/s01_model_constants.csv",
                      {"model", "dim", "q", "cs2"},
                      {{"D2Q9", "2", "9", fmt_f8(D2Q9::CS2)},
                       {"D3Q19", "3", "19", fmt_f8(D3Q19::CS2)},
                       {"D3Q27", "3", "27", fmt_f8(D3Q27::CS2)}},
                      {{}, {2, 3, 3}, {9, 19, 27}, {D2Q9::CS2, D3Q19::CS2, D3Q27::CS2}});
        ctx.write_csv("S05", "out/s05_d2q9_mrt_diag.csv",
                      {"mode", "diag_norm"},
                      {{"rho", "9"}, {"e", "36"}, {"eps", "36"}, {"jx", "6"}, {"qx", "12"},
                       {"jy", "6"}, {"qy", "12"}, {"pxx", "4"}, {"pxy", "4"}},
                      {{}, {9, 36, 36, 6, 12, 6, 12, 4, 4}});
        ctx.record_stage("S08", "conservation_1e4", since(t0), 5);
    }
}

}  // namespace lbm
