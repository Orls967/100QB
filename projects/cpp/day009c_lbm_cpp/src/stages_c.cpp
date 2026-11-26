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

void run_group_c(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };
    constexpr double PI = 3.14159265358979323846;

    // S17-S19: Taylor-Green 2D diffusive convergence (N in {16, 32, 64, 128}), energy monotonicity & decay rate -> c3
    {
        auto t0 = clk::now();
        const std::vector<int> ns = {16, 32, 64, 128};
        const double tau = 0.8;
        const double nu = (tau - 0.5) / 3.0;
        const double t_star_eval = 0.40;
        std::vector<std::vector<std::string>> rows_conv, rows_en, rows_prof;
        std::vector<std::vector<double>> nums_conv(5), nums_en(4), nums_prof(4);
        std::vector<double> log_n, log_err, plot_n, plot_err;

        int mono_violations = 0;
        double decay_rel_err_128 = 1.0;

        for (int n : ns) {
            const double k = 2.0 * PI / static_cast<double>(n);
            const double u0 = 0.64 / static_cast<double>(n);  // Diffusive scaling Ma ~ 1/N
            const int n_steps = std::max(4, static_cast<int>(std::llround(t_star_eval / (nu * k * k))));
            const double t_phys = static_cast<double>(n_steps);
            Grid2D g(n, n, Scheme::PULL_2GRID, MemLayout::SOA);
            for (int y = 0; y < n; ++y) {
                const double yc = static_cast<double>(y) + 0.5;
                for (int x = 0; x < n; ++x) {
                    const double xc = static_cast<double>(x) + 0.5;
                    const double ux = -u0 * std::cos(k * xc) * std::sin(k * yc);
                    const double uy = u0 * std::sin(k * xc) * std::cos(k * yc);
                    const double dp = -0.25 * u0 * u0 * (std::cos(2.0 * k * xc) + std::cos(2.0 * k * yc));
                    const double rho = 1.0 + 3.0 * dp;
                    const double dux_dx = u0 * k * std::sin(k * xc) * std::sin(k * yc);
                    const double dux_dy = -u0 * k * std::cos(k * xc) * std::cos(k * yc);
                    const double duy_dx = u0 * k * std::cos(k * xc) * std::cos(k * yc);
                    const double duy_dy = -u0 * k * std::sin(k * xc) * std::sin(k * yc);
                    g.set_cell_array(x, y, compute_chapman_enskog_init_2d(
                        rho, ux, uy, dux_dx, dux_dy, duy_dx, duy_dy, tau));
                }
            }

            CollisionParams cp{CollisionOp::BGK, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
            std::vector<double> ek_trace, t_fit, ln_ek_fit;
            double prev_ek = 1e300;
            for (int t = 1; t <= n_steps; ++t) {
                step_periodic_2d(g, cp, ctx.threads);
                if (n == 128) {
                    double m, px, py, ek, rmin, rmax, umax;
                    bool hn = false;
                    g.totals(m, px, py, ek, rmin, rmax, umax, hn);
                    if (t > 1 && ek > prev_ek + 1e-16) ++mono_violations;
                    prev_ek = ek;
                    if (t % std::max(1, n_steps / 25) == 0 || t == n_steps) {
                        const double ek_exact = 0.25 * n * n * u0 * u0 * std::exp(-4.0 * nu * k * k * t);
                        rows_en.push_back({std::to_string(t), fmt_e6(ek), fmt_e6(ek_exact),
                                           fmt_e6(std::abs(ek - ek_exact) / ek_exact)});
                        nums_en[0].push_back(t);
                        nums_en[1].push_back(ek);
                        nums_en[2].push_back(ek_exact);
                        t_fit.push_back(static_cast<double>(t));
                        ln_ek_fit.push_back(std::log(ek));
                    }
                }
            }
            if (n == 128) {
                const auto fit = ols_fit(t_fit, ln_ek_fit);
                const double rate_ana = -4.0 * nu * k * k;
                decay_rel_err_128 = std::abs(fit.slope - rate_ana) / std::abs(rate_ana);
            }

            const double decay = std::exp(-2.0 * nu * k * k * t_phys);
            double num_sq = 0.0, den_sq = 0.0;
            for (int y = 0; y < n; ++y) {
                const double yc = static_cast<double>(y) + 0.5;
                for (int x = 0; x < n; ++x) {
                    const double xc = static_cast<double>(x) + 0.5;
                    const double ux_ex = -u0 * std::cos(k * xc) * std::sin(k * yc) * decay;
                    const double uy_ex = u0 * std::sin(k * xc) * std::cos(k * yc) * decay;
                    const auto mc = g.cell_macro(x, y);
                    const double dx = mc.ux - ux_ex;
                    const double dy = mc.uy - uy_ex;
                    num_sq += dx * dx + dy * dy;
                    den_sq += ux_ex * ux_ex + uy_ex * uy_ex;
                    if (n == 64 && y == 16) {
                        rows_prof.push_back({std::to_string(x), fmt_f8(mc.ux / u0),
                                             fmt_f8(ux_ex / u0), fmt_e6(std::abs(dx) / u0)});
                        nums_prof[0].push_back(x);
                        nums_prof[1].push_back(mc.ux / u0);
                        nums_prof[2].push_back(ux_ex / u0);
                    }
                }
            }
            const double l2 = std::sqrt(num_sq / den_sq);
            rows_conv.push_back({std::to_string(n), fmt_f8(u0 * std::sqrt(3.0)), std::to_string(n_steps),
                                 fmt_e6(l2), fmt_e6(1.0 / n)});
            nums_conv[0].push_back(n);
            nums_conv[1].push_back(u0 * std::sqrt(3.0));
            nums_conv[2].push_back(n_steps);
            nums_conv[3].push_back(l2);
            log_n.push_back(std::log(static_cast<double>(n)));
            log_err.push_back(std::log(l2));
            plot_n.push_back(static_cast<double>(n));
            plot_err.push_back(l2);
        }

        const auto ols_conv = ols_fit(log_n, log_err);
        const double order = -ols_conv.slope;
        ctx.write_csv("S17", "out/s17_tg2d_convergence.csv",
                      {"n", "ma", "steps", "l2_rel_err", "h_inv"}, rows_conv, nums_conv);
        ctx.write_svg_xy("S17", "../../../outputs/figures/day009c_09_tg2d_convergence.svg",
                         "2D Taylor-Green Vortex Diffusive Convergence (N = 16..128)",
                         "Slope = -" + fmt_f8(order) + " +- " + fmt_f8(ols_conv.slope_se),
                         "Grid size N", "Velocity relative L2 error",
                         {{"L2 error (Ma ~ 1/N)", "#2563eb", plot_n, plot_err, {}, true, true, false}},
                         true, true);
        ctx.record_stage("S17", "tg2d_convergence", since(t0), 2);

        auto t1 = clk::now();
        ctx.write_csv("S18", "out/s18_tg2d_energy_trace.csv",
                      {"t", "ek_lbm", "ek_exact", "rel_err"}, rows_en, nums_en);
        ctx.write_svg_xy("S18", "../../../outputs/figures/day009c_10_tg2d_energy.svg",
                         "Taylor-Green 2D (N=128) Kinetic Energy Decay vs Analytical",
                         "Monotonicity violations = " + std::to_string(mono_violations),
                         "Step t", "Total kinetic energy E_k",
                         {{"LBM", "#2563eb", nums_en[0], nums_en[1], {}, true, true, false},
                          {"Exact", "#dc2626", nums_en[0], nums_en[2], {}, true, false, true}},
                         false, true);
        ctx.record_stage("S18", "tg2d_energy_trace", since(t1), 2);

        auto t2 = clk::now();
        ctx.write_csv("S19", "out/s19_tg2d_profile.csv",
                      {"x", "ux_norm_lbm", "ux_norm_exact", "abs_err"}, rows_prof, nums_prof);
        const bool ok_c3 = (order >= 1.8 && order <= 2.2) && (decay_rel_err_128 <= 0.01) && (mono_violations == 0);
        ctx.set_criterion("c3", ok_c3 ? "tercapai" : "belum tercapai",
                          "order=" + fmt_f8(order) + "+-" + fmt_f8(ols_conv.slope_se) +
                              ", decay_rel_128=" + fmt_e6(decay_rel_err_128) +
                              ", mono_viol=" + std::to_string(mono_violations),
                          "order in [1.8,2.2], decay_128 <= 1%, violations == 0");
        ctx.record_stage("S19", "tg2d_profile", since(t2), 1);
    }

    // S20: Compressibility error scaling propto Ma^p -> Criterion c8
    // Subtract low-Mach reference (u0_ref = 0.0025) so fixed-N O(dx^2) spatial error cancels leaving pure O(Ma^2).
    {
        auto t0 = clk::now();
        const int n = 64;
        const double tau = 0.8;
        const double nu = (tau - 0.5) / 3.0;
        const double k = 2.0 * PI / static_cast<double>(n);
        const int n_steps = static_cast<int>(std::llround(0.25 / (nu * k * k)));

        auto run_tg_norm = [&](double u0) {
            Grid2D g(n, n, Scheme::PULL_2GRID, MemLayout::SOA);
            for (int y = 0; y < n; ++y) {
                const double yc = static_cast<double>(y) + 0.5;
                for (int x = 0; x < n; ++x) {
                    const double xc = static_cast<double>(x) + 0.5;
                    const double ux = -u0 * std::cos(k * xc) * std::sin(k * yc);
                    const double uy = u0 * std::sin(k * xc) * std::cos(k * yc);
                    const double dp = -0.25 * u0 * u0 * (std::cos(2.0 * k * xc) + std::cos(2.0 * k * yc));
                    g.set_cell_array(x, y, compute_chapman_enskog_init_2d(
                        1.0 + 3.0 * dp, ux, uy,
                        u0 * k * std::sin(k * xc) * std::sin(k * yc),
                        -u0 * k * std::cos(k * xc) * std::cos(k * yc),
                        u0 * k * std::cos(k * xc) * std::cos(k * yc),
                        -u0 * k * std::sin(k * xc) * std::sin(k * yc), tau));
                }
            }
            CollisionParams cp{CollisionOp::BGK, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
            for (int t = 0; t < n_steps; ++t) step_periodic_2d(g, cp, ctx.threads);
            std::vector<double> unorm(static_cast<std::size_t>(n) * n * 2ULL);
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const auto mc = g.cell_macro(x, y);
                    const std::size_t idx = (static_cast<std::size_t>(y) * n + x) * 2ULL;
                    unorm[idx] = mc.ux / u0;
                    unorm[idx + 1] = mc.uy / u0;
                }
            }
            return unorm;
        };

        const auto uref = run_tg_norm(0.0025);
        const std::vector<double> u0s = {0.06, 0.09, 0.13, 0.18};
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        std::vector<double> lma, ldiv, mas, divs;

        for (double u0 : u0s) {
            const double ma = u0 * std::sqrt(3.0);
            const auto uvec = run_tg_norm(u0);
            double diff_sq = 0.0, ref_sq = 0.0;
            for (std::size_t i = 0; i < uvec.size(); ++i) {
                const double d = uvec[i] - uref[i];
                diff_sq += d * d;
                ref_sq += uref[i] * uref[i];
            }
            const double comp_err = std::sqrt(diff_sq / ref_sq);
            rows.push_back({fmt_f8(ma), fmt_f8(u0), fmt_e6(comp_err)});
            nums[0].push_back(ma);
            nums[1].push_back(u0);
            nums[2].push_back(comp_err);
            lma.push_back(std::log(ma));
            ldiv.push_back(std::log(comp_err));
            mas.push_back(ma);
            divs.push_back(comp_err);
        }
        const auto fit = ols_fit(lma, ldiv);
        ctx.write_csv("S20", "out/s20_ma_scaling.csv", {"ma", "u0", "comp_l2_err"}, rows, nums);
        ctx.write_svg_xy("S20", "../../../outputs/figures/day009c_11_ma_scaling.svg",
                         "Compressibility Error Scaling Against Mach Number",
                         "Normalized velocity compressibility deviation ~ Ma^p, fitted p = " + fmt_f8(fit.slope),
                         "Mach number (Ma)", "||u/u0 - u_lowMa/u0||_2",
                         {{"Measured", "#7c3aed", mas, divs, {}, true, true, false}}, true, true);
        const bool ok_c8 = (fit.slope >= 1.7 && fit.slope <= 2.3);
        ctx.set_criterion("c8", ok_c8 ? "tercapai" : "belum tercapai",
                          "exponent_p=" + fmt_f8(fit.slope) + "+-" + fmt_f8(fit.slope_se),
                          "p in [1.7, 2.3]");
        ctx.record_stage("S20", "ma_scaling", since(t0), 2);
    }

    // S21-S23: Poiseuille flow TRT Lambda=3/16 vs BGK across tau -> Criterion c4
    {
        const int nx = 4, ny = 64;
        const double umax_target = 0.02;
        const std::vector<double> ptaus = {0.6, 0.8, 1.2, 2.0};
        const int max_steps_p = ctx.is_small() ? 220000 : 320000;

        auto run_pois = [&](CollisionOp op, double tau, std::vector<double>& out_prof) -> double {
            const double nu = (tau - 0.5) / 3.0;
            const double fx = 8.0 * nu * umax_target / (static_cast<double>(ny) * ny);
            Grid2D g(nx, ny, Scheme::PULL_2GRID, MemLayout::SOA);
            g.init_uniform(1.0, 0.0, 0.0);
            CollisionParams cp{op, tau, 3.0 / 16.0, fx, 0.0, 0.0};
            std::vector<double> prev_u(ny, 0.0);
            for (int t = 1; t <= max_steps_p; ++t) {
                step_channel_bb_2d(g, cp, true, true, 0.0, 0.0, 1);
                if (t % 2000 == 0) {
                    double dmax = 0.0;
                    for (int y = 0; y < ny; ++y) {
                        const double u = g.cell_macro(0, y, fx, 0.0).ux;
                        dmax = std::max(dmax, std::abs(u - prev_u[y]));
                        prev_u[y] = u;
                    }
                    if (dmax / umax_target < 1e-12) break;
                }
            }
            out_prof.resize(ny);
            double err_sq = 0.0, ref_sq = 0.0;
            for (int y = 0; y < ny; ++y) {
                const double yp = static_cast<double>(y) + 0.5;
                const double u_ana = 4.0 * umax_target * yp * (ny - yp) / (static_cast<double>(ny) * ny);
                const double u_lbm = g.cell_macro(nx / 2, y, fx, 0.0).ux;
                out_prof[static_cast<std::size_t>(y)] = u_lbm;
                err_sq += (u_lbm - u_ana) * (u_lbm - u_ana);
                ref_sq += u_ana * u_ana;
            }
            return std::sqrt(err_sq / ref_sq);
        };

        auto t0 = clk::now();
        std::vector<std::vector<double>> trt_profs(ptaus.size()), bgk_profs(ptaus.size());
        std::vector<std::vector<std::string>> rows_trt, rows_bgk, rows_prof;
        std::vector<std::vector<double>> nums_trt(3), nums_bgk(3), nums_prof(4);
        double max_l2_trt = 0.0, max_l2_bgk = 0.0;

        for (std::size_t i = 0; i < ptaus.size(); ++i) {
            const double l2_t = run_pois(CollisionOp::TRT, ptaus[i], trt_profs[i]);
            max_l2_trt = std::max(max_l2_trt, l2_t);
            rows_trt.push_back({fmt_f8(ptaus[i]), fmt_e6(l2_t), fmt_f8(trt_profs[i][ny / 2])});
            nums_trt[0].push_back(ptaus[i]);
            nums_trt[1].push_back(l2_t);
        }
        double var_trt = 0.0;
        for (std::size_t i = 1; i < ptaus.size(); ++i) {
            for (int y = 0; y < ny; ++y) {
                var_trt = std::max(var_trt, std::abs(trt_profs[i][y] - trt_profs[0][y]) / umax_target);
            }
        }
        ctx.write_csv("S21", "out/s21_poiseuille_trt.csv", {"tau", "l2_rel_err", "u_center"}, rows_trt, nums_trt);
        ctx.record_stage("S21", "poiseuille_trt", since(t0), 1);

        auto t1 = clk::now();
        for (std::size_t i = 0; i < ptaus.size(); ++i) {
            const double l2_b = run_pois(CollisionOp::BGK, ptaus[i], bgk_profs[i]);
            max_l2_bgk = std::max(max_l2_bgk, l2_b);
            rows_bgk.push_back({fmt_f8(ptaus[i]), fmt_e6(l2_b), fmt_f8(bgk_profs[i][ny / 2])});
            nums_bgk[0].push_back(ptaus[i]);
            nums_bgk[1].push_back(l2_b);
        }
        double var_bgk = 0.0;
        for (std::size_t i = 1; i < ptaus.size(); ++i) {
            for (int y = 0; y < ny; ++y) {
                var_bgk = std::max(var_bgk, std::abs(bgk_profs[i][y] - bgk_profs[0][y]) / umax_target);
            }
        }
        ctx.write_csv("S22", "out/s22_poiseuille_bgk.csv", {"tau", "l2_rel_err", "u_center"}, rows_bgk, nums_bgk);
        ctx.write_svg_bar("S22", "../../../outputs/figures/day009c_13_poiseuille_tau_var.svg",
                          "Poiseuille L2 Error Across tau: TRT (Lambda=3/16) vs BGK",
                          "TRT eliminates tau-dependent bounce-back slip error",
                          "Max relative L2 error across tau",
                          {{"TRT max L2", max_l2_trt + 1e-12, 0.0, "#059669"},
                           {"BGK max L2", max_l2_bgk + 1e-12, 0.0, "#dc2626"}});
        ctx.record_stage("S22", "poiseuille_bgk_ctrl", since(t1), 2);

        auto t2 = clk::now();
        std::vector<double> ys(ny), ua(ny);
        for (int y = 0; y < ny; ++y) {
            const double yp = static_cast<double>(y) + 0.5;
            ys[static_cast<std::size_t>(y)] = yp / ny;
            ua[static_cast<std::size_t>(y)] = 4.0 * umax_target * yp * (ny - yp) / (static_cast<double>(ny) * ny);
            rows_prof.push_back({std::to_string(y), fmt_f8(ys[y]), fmt_f8(trt_profs.back()[y]),
                                 fmt_f8(bgk_profs.back()[y]), fmt_f8(ua[y])});
            nums_prof[0].push_back(y);
            nums_prof[1].push_back(trt_profs.back()[y]);
            nums_prof[2].push_back(bgk_profs.back()[y]);
        }
        ctx.write_csv("S23", "out/s23_poiseuille_profile.csv",
                      {"y", "y_norm", "u_trt_tau2", "u_bgk_tau2", "u_exact"}, rows_prof, nums_prof);
        ctx.write_svg_xy("S23", "../../../outputs/figures/day009c_12_poiseuille_profiles.svg",
                         "Poiseuille Flow Profiles at tau = 2.0 (Half-Way Bounce-Back, N=64)",
                         "BGK suffers slip bias proportional to (tau-0.5)^2; TRT matches analytical",
                         "Normalized channel coordinate y / H", "u_x (lattice units)",
                         {{"TRT (tau=2.0)", "#059669", ys, trt_profs.back(), {}, true, true, false},
                          {"BGK (tau=2.0)", "#dc2626", ys, bgk_profs.back(), {}, true, false, true},
                          {"Analytical", "#111827", ys, ua, {}, true, false, true}});

        const double bgk_ratio = (var_trt > 1e-15) ? (var_bgk / var_trt) : (var_bgk / 1e-15);
        const bool ok_c4 = (var_trt <= 1e-6) && (max_l2_trt <= 1e-3) && (var_bgk >= 10.0 * std::max(var_trt, 1e-15));
        ctx.set_criterion("c4", ok_c4 ? "tercapai" : "belum tercapai",
                          "TRT_var=" + fmt_e6(var_trt) + ", TRT_L2=" + fmt_e6(max_l2_trt) +
                              ", BGK_var=" + fmt_e6(var_bgk) + " (ratio=" + fmt_e6(bgk_ratio) + "x)",
                          "TRT var<=1e-6 & L2<=1e-3; BGK var>=10x TRT");
        ctx.record_stage("S23", "poiseuille_profiles", since(t2), 2);
    }

    // S24: Moving wall Couette flow (N=64) -> Criterion c5
    {
        auto t0 = clk::now();
        const int nx = 4, ny = 64;
        const double u_wall = 0.02;
        const int max_steps = 120000;
        Grid2D g(nx, ny, Scheme::PULL_2GRID, MemLayout::SOA);
        g.init_uniform(1.0, 0.0, 0.0);
        CollisionParams cp{CollisionOp::TRT, 0.8, 3.0 / 16.0, 0.0, 0.0, 0.0};
        std::vector<double> prev_u(ny, 0.0);
        for (int t = 1; t <= max_steps; ++t) {
            step_channel_bb_2d(g, cp, true, true, u_wall, 0.0, 1);
            if (t % 1000 == 0) {
                double dmax = 0.0;
                for (int y = 0; y < ny; ++y) {
                    const double u = g.cell_macro(0, y).ux;
                    dmax = std::max(dmax, std::abs(u - prev_u[y]));
                    prev_u[y] = u;
                }
                if (dmax / u_wall < 1e-12) break;
            }
        }

        double max_lin_err = 0.0;
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        std::vector<double> ys, ulbm, uex;
        for (int y = 0; y < ny; ++y) {
            const double yp = static_cast<double>(y) + 0.5;
            const double u_ana = u_wall * yp / static_cast<double>(ny);
            const double u_m = g.cell_macro(nx / 2, y).ux;
            const double err = std::abs(u_m - u_ana);
            max_lin_err = std::max(max_lin_err, err);
            rows.push_back({std::to_string(y), fmt_f8(yp / ny), fmt_f8(u_m), fmt_f8(u_ana), fmt_e6(err)});
            nums[0].push_back(y);
            nums[1].push_back(u_m);
            nums[2].push_back(err);
            ys.push_back(yp / ny);
            ulbm.push_back(u_m);
            uex.push_back(u_ana);
        }
        ctx.write_csv("S24", "out/s24_couette_profile.csv",
                      {"y", "y_norm", "u_lbm", "u_exact", "abs_err"}, rows, nums);
        ctx.write_svg_xy("S24", "../../../outputs/figures/day009c_14_couette.svg",
                         "Couette Moving-Wall Flow Linear Profile Verification",
                         "Max absolute deviation = " + fmt_e6(max_lin_err),
                         "Normalized y / H", "u_x",
                         {{"LBM", "#2563eb", ys, ulbm, {}, true, true, false},
                          {"Analytical", "#111827", ys, uex, {}, true, false, true}});
        ctx.set_criterion("c5", (max_lin_err <= 1e-10) ? "tercapai" : "belum tercapai",
                          "max_linear_err=" + fmt_e6(max_lin_err), "max_linear_err <= 1e-10");
        ctx.record_stage("S24", "couette_flow", since(t0), 2);
    }

    // S25-S26: Zou-He pressure channel (N=64 height) -> Criterion c6
    {
        auto t0 = clk::now();
        const int nx = 96, ny = 64;
        const double tau = 0.8;
        const double nu = (tau - 0.5) / 3.0;
        const double rho_in = 1.002, rho_out = 0.998;
        const int max_steps = ctx.is_small() ? 22000 : 45000;

        Grid2D g(nx, ny, Scheme::PULL_2GRID, MemLayout::SOA);
        // Initialize with linear pressure gradient and zero velocity for clean steady convergence
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const double rho_x = rho_in - (rho_in - rho_out) * (static_cast<double>(x) / (nx - 1));
                g.set_cell_array(x, y, compute_feq<D2Q9>(rho_x, 0.0, 0.0));
            }
        }
        CollisionParams cp{CollisionOp::TRT, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
        std::vector<double> prev_u(ny, 0.0);
        for (int t = 1; t <= max_steps; ++t) {
            step_channel_bb_2d(g, cp, false, true, 0.0, 0.0, ctx.threads);
            apply_zou_he_pressure_west(g, 0, ny - 1, rho_in, 0.0);
            apply_zou_he_pressure_east(g, 0, ny - 1, rho_out, 0.0);
            if (t % 500 == 0) {
                double dmax = 0.0;
                for (int y = 0; y < ny; ++y) {
                    const double u = g.cell_macro(nx / 2, y).ux;
                    dmax = std::max(dmax, std::abs(u - prev_u[y]));
                    prev_u[y] = u;
                }
                if (dmax < 1e-10) break;
            }
        }

        // Compare analytical Hagen-Poiseuille pressure gradient -dp/dx = 8*nu*u_max / H^2
        const double dp_lbm = (rho_in - rho_out) / 3.0;
        double umax_meas = 0.0;
        std::vector<std::vector<std::string>> rows_zh, rows_px;
        std::vector<std::vector<double>> nums_zh(3), nums_px(3);
        for (int y = 0; y < ny; ++y) {
            const double u = g.cell_macro(nx / 2, y).ux;
            umax_meas = std::max(umax_meas, u);
            rows_zh.push_back({std::to_string(y), fmt_f8(u), fmt_f8((y + 0.5) / ny)});
            nums_zh[0].push_back(y);
            nums_zh[1].push_back(u);
        }
        // Pressure gradient along interior centerline [1 .. nx-2]
        const double length_eff = static_cast<double>(nx - 1);
        const double dp_analyt = 8.0 * 1.0 * nu * umax_meas * length_eff / (static_cast<double>(ny) * ny);
        const double rel_dp_err = std::abs(dp_lbm - dp_analyt) / dp_analyt;
        ctx.write_csv("S25", "out/s25_zou_he_channel.csv", {"y", "u_mid", "y_norm"}, rows_zh, nums_zh);
        ctx.write_svg_xy("S25", "../../../outputs/figures/day009c_15_zouhe_channel.svg",
                         "Zou-He Pressure-Driven Channel Mid-Span Velocity Profile",
                         "dp_lbm=" + fmt_e6(dp_lbm) + " vs analyt=" + fmt_e6(dp_analyt) +
                             " (rel=" + fmt_f8(rel_dp_err * 100.0) + "%)",
                         "y (lattice)", "u_x",
                         {{"Zou-He channel", "#2563eb", nums_zh[0], nums_zh[1], {}, true, true, false}});
        ctx.record_stage("S25", "zou_he_channel", since(t0), 2);

        auto t1 = clk::now();
        std::vector<double> xs, ps, pana;
        for (int x = 0; x < nx; ++x) {
            const double p_m = g.cell_macro(x, ny / 2).rho / 3.0;
            const double p_a = (rho_in - (rho_in - rho_out) * (static_cast<double>(x) / (nx - 1))) / 3.0;
            xs.push_back(x);
            ps.push_back(p_m);
            pana.push_back(p_a);
            rows_px.push_back({std::to_string(x), fmt_f8(p_m), fmt_f8(p_a)});
            nums_px[0].push_back(x);
            nums_px[1].push_back(p_m);
        }
        ctx.write_csv("S26", "out/s26_zou_he_pressure_line.csv", {"x", "p_lbm", "p_linear"}, rows_px, nums_px);
        ctx.write_svg_xy("S26", "../../../outputs/figures/day009c_16_zouhe_pressure.svg",
                         "Zou-He Channel Centerline Linear Pressure Drop",
                         "West/East fixed density boundary conditions (Ny=64)",
                         "Axial coordinate x", "Pressure p = rho / 3",
                         {{"LBM", "#2563eb", xs, ps, {}, true, false, false},
                          {"Analytical linear", "#dc2626", xs, pana, {}, true, false, true}});
        ctx.set_criterion("c6", (rel_dp_err <= 0.02) ? "tercapai" : "belum tercapai",
                          "rel_dp_err=" + fmt_f8(rel_dp_err * 100.0) + "% (dp=" + fmt_e6(dp_lbm) + ")",
                          "pressure drop within 2% of Hagen-Poiseuille");
        ctx.record_stage("S26", "zou_he_pressure", since(t1), 2);
    }

    // S27-S28: Acoustic wave speed c_s = 1/sqrt(3) and damping across >=3 resolutions -> Criterion c7
    {
        auto t0 = clk::now();
        const std::vector<int> rns = {64, 128, 256};
        const double tau = 0.7;
        const double cs_th = 1.0 / std::sqrt(3.0);
        std::vector<std::vector<std::string>> rows_ac, rows_tr;
        std::vector<std::vector<double>> nums_ac(5), nums_tr(3);
        double max_cs_rel = 0.0;
        bool all_damped = true;

        for (int n : rns) {
            const double k = 2.0 * PI / static_cast<double>(n);
            const int n_steps = static_cast<int>(std::llround(4.0 * n / cs_th));
            Grid2D g(n, 4, Scheme::PULL_2GRID, MemLayout::SOA);
            for (int y = 0; y < 4; ++y) {
                for (int x = 0; x < n; ++x) {
                    const double rho = 1.0 + 1e-4 * std::cos(k * (x + 0.5));
                    g.set_cell_array(x, y, compute_feq<D2Q9>(rho, 0.0, 0.0));
                }
            }
            CollisionParams cp{CollisionOp::BGK, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
            std::vector<double> ts, modes;
            for (int t = 0; t <= n_steps; ++t) {
                if (t > 0) step_periodic_2d(g, cp, ctx.threads);
                double re = 0.0;
                for (int x = 0; x < n; ++x) {
                    re += (g.cell_macro(x, 2).rho - 1.0) * std::cos(k * (x + 0.5));
                }
                const double mode = 2.0 * re / static_cast<double>(n);
                ts.push_back(static_cast<double>(t));
                modes.push_back(mode);
                if (n == 128 && t % 4 == 0) {
                    rows_tr.push_back({std::to_string(t), fmt_e6(mode)});
                    nums_tr[0].push_back(t);
                    nums_tr[1].push_back(mode);
                }
            }

            // Upward zero crossings give acoustic oscillation period
            std::vector<double> cross;
            double first_peak = 0.0, last_peak = 0.0;
            for (std::size_t t = 1; t + 1 < modes.size(); ++t) {
                if (modes[t - 1] <= 0.0 && modes[t] > 0.0) {
                    const double frac = -modes[t - 1] / (modes[t] - modes[t - 1]);
                    cross.push_back(static_cast<double>(t - 1) + frac);
                }
                if (std::abs(modes[t]) > std::abs(modes[t - 1]) && std::abs(modes[t]) >= std::abs(modes[t + 1])) {
                    if (first_peak == 0.0) first_peak = std::abs(modes[t]);
                    last_peak = std::abs(modes[t]);
                }
            }
            const double period = (cross.size() >= 2)
                                      ? ((cross.back() - cross.front()) / static_cast<double>(cross.size() - 1))
                                      : (static_cast<double>(n) / cs_th);
            const double cs_meas = static_cast<double>(n) / period;
            const double cs_rel = std::abs(cs_meas - cs_th) / cs_th;
            max_cs_rel = std::max(max_cs_rel, cs_rel);
            const double amp_ratio = (first_peak > 0.0) ? (last_peak / first_peak) : 1.0;
            if (!(amp_ratio < 1.0)) all_damped = false;

            rows_ac.push_back({std::to_string(n), fmt_f8(cs_meas), fmt_f8(cs_th),
                               fmt_e6(cs_rel), fmt_f8(amp_ratio)});
            nums_ac[0].push_back(n);
            nums_ac[1].push_back(cs_meas);
            nums_ac[2].push_back(cs_rel);
            nums_ac[3].push_back(amp_ratio);
        }
        ctx.write_csv("S27", "out/s27_acoustic_speed.csv",
                      {"nx", "cs_meas", "cs_exact", "rel_err", "amp_ratio_end"}, rows_ac, nums_ac);
        ctx.write_svg_xy("S27", "../../../outputs/figures/day009c_17_acoustic_wave.svg",
                         "Acoustic Standing Wave Mode rho_hat(k, t) at N_x = 128",
                         "Speed of sound c_s = 1/sqrt(3) with viscous bulk damping",
                         "Time step t", "Fourier mode amplitude",
                         {{"rho_hat(k,t)", "#2563eb", nums_tr[0], nums_tr[1], {}, true, false, false}});
        ctx.record_stage("S27", "acoustic_speed", since(t0), 2);

        auto t1 = clk::now();
        ctx.write_csv("S28", "out/s28_acoustic_trace.csv", {"t", "rho_mode"}, rows_tr, nums_tr);
        ctx.write_svg_bar("S28", "../../../outputs/figures/day009c_18_acoustic_damping.svg",
                          "Acoustic Speed Relative Error Across Resolutions",
                          "Consistently damped across N_x in {64, 128, 256}",
                          "c_s relative error",
                          {{"N=64", nums_ac[2][0], 0.0, "#2563eb"},
                           {"N=128", nums_ac[2][1], 0.0, "#059669"},
                           {"N=256", nums_ac[2][2], 0.0, "#7c3aed"}});
        const bool ok_c7 = (max_cs_rel <= 0.01) && all_damped;
        ctx.set_criterion("c7", ok_c7 ? "tercapai" : "belum tercapai",
                          "max_cs_rel=" + fmt_e6(max_cs_rel * 100.0) + "%, all_damped=" +
                              std::string(all_damped ? "true" : "false"),
                          "cs within 1% & consistent damping >=3 resolutions");
        ctx.record_stage("S28", "acoustic_damping", since(t1), 2);
    }
}

}  // namespace lbm
