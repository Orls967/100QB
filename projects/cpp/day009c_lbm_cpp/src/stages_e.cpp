#include "lbm/stages.hpp"
#include "lbm/analysis.hpp"
#include "lbm/boundary.hpp"
#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/grid.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

namespace lbm {

void run_group_e(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };
    constexpr double PI = 3.14159265358979323846;

    // S41-S42: 3D D3Q19 shear wave viscosity vs (tau - 0.5) / 3
    double worst_3d_shear = 0.0;
    {
        auto t0 = clk::now();
        const int n = ctx.is_small() ? 16 : 32;
        const double k = 2.0 * PI / static_cast<double>(n);
        const std::vector<double> taus = {0.65, 0.85, 1.15};
        std::vector<std::vector<std::string>> rows, rows_tr;
        std::vector<std::vector<double>> nums(4), nums_tr(3);

        for (double tau : taus) {
            const double nu_th = (tau - 0.5) / 3.0;
            const int t_min = std::max(8, static_cast<int>(std::ceil(0.05 / (nu_th * k * k))));
            const int t_max = std::max(t_min + 25, static_cast<int>(std::floor(1.20 / (nu_th * k * k))));
            const int samp = std::max(1, (t_max - t_min) / 18);

            Grid3D g(n, n, n, Scheme::PULL_2GRID);
            for (int z = 0; z < n; ++z) {
                const double zc = static_cast<double>(z) + 0.5;
                for (int y = 0; y < n; ++y) {
                    for (int x = 0; x < n; ++x) {
                        const double ux = 0.01 * std::sin(k * zc);
                        const auto feq = compute_feq<D3Q19>(1.0, ux, 0.0, 0.0);
                        for (int q = 0; q < 19; ++q) {
                            g.set_pop(x, y, z, q, feq[static_cast<std::size_t>(q)]);
                        }
                    }
                }
            }
            CollisionParams cp{CollisionOp::BGK, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
            std::vector<double> ts, lna;
            for (int t = 1; t <= t_max; ++t) {
                step_periodic_3d(g, cp, ctx.threads);
                if (t >= t_min && (t - t_min) % samp == 0) {
                    double num = 0.0;
                    std::array<double, 19> loc{};
                    for (int z = 0; z < n; ++z) {
                        for (int q = 0; q < 19; ++q) loc[static_cast<std::size_t>(q)] = g.get_pop(n / 2, n / 2, z, q);
                        num += compute_macro_3d<D3Q19>(loc.data()).ux * std::sin(k * (z + 0.5));
                    }
                    const double amp = std::abs(2.0 * num / static_cast<double>(n));
                    ts.push_back(static_cast<double>(t));
                    lna.push_back(std::log(std::max(amp, 1e-30)));
                    if (std::abs(tau - 0.85) < 1e-9) {
                        rows_tr.push_back({std::to_string(t), fmt_e6(amp), fmt_e6(0.01 * std::exp(-nu_th * k * k * t))});
                        nums_tr[0].push_back(t);
                        nums_tr[1].push_back(amp);
                        nums_tr[2].push_back(0.01 * std::exp(-nu_th * k * k * t));
                    }
                }
            }
            const auto fit = ols_fit(ts, lna);
            const double nu_m = -fit.slope / (k * k);
            const double rel = std::abs(nu_m - nu_th) / nu_th;
            worst_3d_shear = std::max(worst_3d_shear, rel);
            rows.push_back({fmt_f8(tau), fmt_f8(nu_th), fmt_f8(nu_m), fmt_e6(rel)});
            nums[0].push_back(tau);
            nums[1].push_back(nu_th);
            nums[2].push_back(nu_m);
            nums[3].push_back(rel);
        }
        ctx.write_csv("S41", "out/s41_d3q19_shear.csv", {"tau", "nu_th", "nu_meas", "rel_err"}, rows, nums);
        ctx.record_stage("S41", "d3q19_shear_visc", since(t0), 1);

        auto t1 = clk::now();
        ctx.write_csv("S42", "out/s42_d3q19_shear_trace.csv", {"t", "amp_lbm", "amp_exact"}, rows_tr, nums_tr);
        ctx.write_svg_xy("S42", "../../../outputs/figures/day009c_31_d3q19_shear.svg",
                         "3D D3Q19 Shear Wave Decay at tau = 0.85",
                         "Max viscosity error across tau = " + fmt_f8(worst_3d_shear * 100.0) + "%",
                         "Time step t", "Amplitude |u_hat(k, t)|",
                         {{"LBM D3Q19", "#2563eb", nums_tr[0], nums_tr[1], {}, true, true, false},
                          {"Exact", "#dc2626", nums_tr[0], nums_tr[2], {}, true, false, true}},
                         false, true);
        ctx.record_stage("S42", "d3q19_shear_trace", since(t1), 2);
    }

    // S43-S44: 3D Square Duct Poiseuille vs Analytical Fourier Series
    double duct_l2_rel = 1.0;
    {
        auto t0 = clk::now();
        const int nx = 4;
        const int ny = ctx.is_small() ? 20 : 32;
        const int nz = ny;
        const double tau = 0.85;
        const double nu = (tau - 0.5) / 3.0;
        const double height = static_cast<double>(ny);
        const double u_scale = 0.02;
        const double fx = (u_scale * 1.0 * nu) / (0.07367 * height * height);
        const int n_steps = ctx.is_small() ? 6000 : 12000;

        Grid3D g(nx, ny, nz, Scheme::PULL_2GRID);
        g.init_uniform(1.0, 0.0, 0.0, 0.0);
        CollisionParams cp{CollisionOp::TRT, tau, 3.0 / 16.0, fx, 0.0, 0.0};
        double prev_u_mid = 0.0;
        std::array<double, 19> loc_chk{};
        for (int t = 1; t <= n_steps; ++t) {
            step_duct_3d(g, cp, ctx.threads);
            if (t % 500 == 0) {
                for (int q = 0; q < 19; ++q) loc_chk[static_cast<std::size_t>(q)] = g.get_pop(nx / 2, ny / 2, nz / 2, q);
                const double uc = compute_macro_3d<D3Q19>(loc_chk.data(), fx, 0.0, 0.0).ux;
                if (std::abs(uc - prev_u_mid) / u_scale < 1e-10) break;
                prev_u_mid = uc;
            }
        }

        double num_sq = 0.0, den_sq = 0.0;
        std::vector<std::vector<std::string>> rows_diag, rows_slice;
        std::vector<std::vector<double>> nums_diag(4), nums_slice(2);
        std::vector<double> map3d(static_cast<std::size_t>(ny) * nz);
        std::array<double, 19> loc{};

        for (int z = 0; z < nz; ++z) {
            const double zp = static_cast<double>(z) + 0.5;
            for (int y = 0; y < ny; ++y) {
                const double yp = static_cast<double>(y) + 0.5;
                for (int q = 0; q < 19; ++q) loc[static_cast<std::size_t>(q)] = g.get_pop(nx / 2, y, z, q);
                const double u_m = compute_macro_3d<D3Q19>(loc.data(), fx, 0.0, 0.0).ux;
                const double u_ex = exact_square_duct_ux(yp, zp, height, fx, 1.0, nu);
                num_sq += (u_m - u_ex) * (u_m - u_ex);
                den_sq += u_ex * u_ex;
                map3d[static_cast<std::size_t>(z) * ny + y] = u_m;
                if (y == z) {
                    rows_diag.push_back({std::to_string(y), fmt_f8(yp / height), fmt_f8(u_m), fmt_f8(u_ex)});
                    nums_diag[0].push_back(y);
                    nums_diag[1].push_back(yp / height);
                    nums_diag[2].push_back(u_m);
                    nums_diag[3].push_back(u_ex);
                }
                if (z == nz / 2) {
                    rows_slice.push_back({std::to_string(y), fmt_f8(u_m), fmt_f8(u_ex)});
                    nums_slice[0].push_back(y);
                    nums_slice[1].push_back(u_m);
                }
            }
        }
        duct_l2_rel = std::sqrt(num_sq / den_sq);
        ctx.write_csv("S43", "out/s43_d3q19_square_duct_diag.csv",
                      {"idx", "diag_norm", "u_lbm", "u_exact"}, rows_diag, nums_diag);
        ctx.write_svg_xy("S43", "../../../outputs/figures/day009c_32_duct3d_profile.svg",
                         "3D D3Q19 Square Duct Poiseuille Diagonal Velocity Profile",
                         "L2 relative error vs analytical Fourier series = " + fmt_f8(duct_l2_rel * 100.0) + "%",
                         "Diagonal position y/H (= z/H)", "Axial velocity u_x",
                         {{"LBM D3Q19 TRT", "#2563eb", nums_diag[1], nums_diag[2], {}, true, true, false},
                          {"Fourier series", "#dc2626", nums_diag[1], nums_diag[3], {}, true, false, true}});
        ctx.record_stage("S43", "d3q19_square_duct", since(t0), 2);

        auto t1 = clk::now();
        ctx.write_csv("S44", "out/s44_d3q19_duct_midslice.csv", {"y", "u_lbm", "u_exact"}, rows_slice, nums_slice);
        ctx.write_svg_heatmap("S44", "../../../outputs/figures/day009c_33_duct3d_map.svg",
                              "3D D3Q19 Square Duct Axial Velocity Cross-Section",
                              "Ny=Nz=" + std::to_string(ny) + ", half-way bounce-back",
                              "y", "z", ny, nz, map3d, "u_x");
        ctx.record_stage("S44", "d3q19_duct_map", since(t1), 2);
    }

    // Set c13 (3D shear viscosity <= 0.5%, square duct <= 2%)
    {
        const bool ok_c13 = (worst_3d_shear <= 0.005) && (duct_l2_rel <= 0.02);
        ctx.set_criterion("c13", ok_c13 ? "tercapai" : "belum tercapai",
                          "3D_shear_rel=" + fmt_f8(worst_3d_shear * 100.0) + "%, duct_L2=" +
                              fmt_f8(duct_l2_rel * 100.0) + "%",
                          "3D shear <= 0.5% & square duct <= 2%");
    }

    // S45-S46: 3D Taylor-Green vortex early kinetic energy decay at low Re
    {
        auto t0 = clk::now();
        const int n = ctx.is_small() ? 16 : 32;
        const double tau = 0.85;
        const double nu = (tau - 0.5) / 3.0;
        const double k = 2.0 * PI / static_cast<double>(n);
        const double u0 = 0.015;
        Grid3D g(n, n, n, Scheme::PULL_2GRID);
        for (int z = 0; z < n; ++z) {
            const double zc = static_cast<double>(z) + 0.5;
            for (int y = 0; y < n; ++y) {
                const double yc = static_cast<double>(y) + 0.5;
                for (int x = 0; x < n; ++x) {
                    const double xc = static_cast<double>(x) + 0.5;
                    const double ux = u0 * std::sin(k * xc) * std::cos(k * yc) * std::cos(k * zc);
                    const double uy = -u0 * std::cos(k * xc) * std::sin(k * yc) * std::cos(k * zc);
                    const double dp = (u0 * u0 / 16.0) * (std::cos(2 * k * xc) + std::cos(2 * k * yc)) *
                                      (std::cos(2 * k * zc) + 2.0);
                    const auto feq = compute_feq<D3Q19>(1.0 + 3.0 * dp, ux, uy, 0.0);
                    for (int q = 0; q < 19; ++q) g.set_pop(x, y, z, q, feq[static_cast<std::size_t>(q)]);
                }
            }
        }
        double m0, px0, py0, pz0, ek0, rmin0;
        bool hn = false;
        g.totals(m0, px0, py0, pz0, ek0, rmin0, hn);
        CollisionParams cp{CollisionOp::MRT, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
        std::vector<std::vector<std::string>> rows_tg3, rows_s46;
        std::vector<std::vector<double>> nums_tg3(4), nums_s46(3);
        const int steps3 = ctx.is_small() ? 60 : 120;
        for (int t = 1; t <= steps3; ++t) {
            step_periodic_3d(g, cp, ctx.threads);
            if (t % 5 == 0) {
                double mt, pxt, pyt, pzt, ekt, rmt;
                g.totals(mt, pxt, pyt, pzt, ekt, rmt, hn);
                const double ek_ex = ek0 * std::exp(-6.0 * nu * k * k * t);
                rows_tg3.push_back({std::to_string(t), fmt_e6(ekt), fmt_e6(ek_ex), fmt_e6(std::abs(ekt - ek_ex) / ek_ex)});
                nums_tg3[0].push_back(t);
                nums_tg3[1].push_back(ekt);
                nums_tg3[2].push_back(ek_ex);
                nums_tg3[3].push_back(std::abs(ekt - ek_ex) / ek_ex);
            }
        }
        ctx.write_csv("S45", "out/s45_tg3d_energy.csv", {"t", "ek_lbm", "ek_analyt", "rel_err"}, rows_tg3, nums_tg3);
        ctx.write_svg_xy("S45", "../../../outputs/figures/day009c_34_tg3d_decay.svg",
                         "3D D3Q19 Taylor-Green Early Kinetic Energy Decay",
                         "Low-Re analytical decay E_0 exp(-6 nu k^2 t)",
                         "Time step t", "Total kinetic energy",
                         {{"LBM D3Q19", "#2563eb", nums_tg3[0], nums_tg3[1], {}, true, true, false},
                          {"Analytical", "#dc2626", nums_tg3[0], nums_tg3[2], {}, true, false, true}},
                         false, true);
        ctx.record_stage("S45", "tg3d_early_decay", since(t0), 2);

        auto t1 = clk::now();
        rows_s46 = {{"rel_err_end", fmt_e6(nums_tg3[3].back())}, {"ek0", fmt_e6(ek0)}};
        nums_s46 = {{}, {nums_tg3[3].back(), ek0}};
        ctx.write_csv("S46", "out/s46_tg3d_summary.csv", {"metric", "value"}, rows_s46, nums_s46);
        ctx.write_svg_xy("S46", "../../../outputs/figures/day009c_35_tg3d_relerr.svg",
                         "3D Taylor-Green Relative Energy Error Trajectory",
                         "D3Q19 on " + std::to_string(n) + "^3 periodic domain",
                         "Time step t", "Relative error |E_k - E_ex| / E_ex",
                         {{"Relative error", "#7c3aed", nums_tg3[0], nums_tg3[3], {}, true, true, false}});
        ctx.record_stage("S46", "tg3d_summary", since(t1), 2);
    }
}

}  // namespace lbm
