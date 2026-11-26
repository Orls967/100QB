#include "lbm/stages.hpp"
#include "lbm/analysis.hpp"
#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/grid.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

namespace lbm {

namespace {

struct ShearRun {
    double tau = 0.8;
    double nu_th = 0.1;
    double nu_meas = 0.1;
    double rel_err = 0.0;
    double slope = 0.0;
    double slope_se = 0.0;
    std::vector<double> ts;
    std::vector<double> amps;
};

ShearRun measure_shear_wave_2d(int n, CollisionOp op, double tau, double u0, double u_wind, int threads) {
    constexpr double PI = 3.14159265358979323846;
    const double k = 2.0 * PI / static_cast<double>(n);
    const double nu_th = (tau - 0.5) / 3.0;
    // Pre-declared dimensionless diffusive time window theta = nu * k^2 * t in [0.05, 1.50]
    const int t_min = std::max(10, static_cast<int>(std::ceil(0.05 / (nu_th * k * k))));
    const int t_max = std::max(t_min + 40, static_cast<int>(std::floor(1.50 / (nu_th * k * k))));
    const int sample_step = std::max(1, (t_max - t_min) / 24);

    Grid2D g(n, n, Scheme::PULL_2GRID, MemLayout::SOA);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            const double xc = static_cast<double>(x) + 0.5;
            // Shear wave u_y(x) advected by background wind u_wind along x (exposes Galilean cubic error if u_wind != 0)
            const double uy = u0 * std::sin(k * xc);
            const double duy_dx = u0 * k * std::cos(k * xc);
            g.set_cell_array(x, y, compute_chapman_enskog_init_2d(1.0, u_wind, uy, 0.0, 0.0, duy_dx, 0.0, tau));
        }
    }

    CollisionParams cp{op, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
    ShearRun res{};
    res.tau = tau;
    res.nu_th = nu_th;
    std::vector<double> ln_amp;

    for (int t = 1; t <= t_max; ++t) {
        step_periodic_2d(g, cp, threads);
        if (t >= t_min && (t - t_min) % sample_step == 0) {
            double re = 0.0, im = 0.0;
            for (int x = 0; x < n; ++x) {
                const double xc = static_cast<double>(x) + 0.5;
                const double uy = g.cell_macro(x, n / 2).uy;
                re += uy * std::cos(k * xc);
                im -= uy * std::sin(k * xc);
            }
            const double amp = 2.0 * std::hypot(re, im) / static_cast<double>(n);
            res.ts.push_back(static_cast<double>(t));
            res.amps.push_back(amp);
            ln_amp.push_back(std::log(std::max(amp, 1e-30)));
        }
    }
    const auto fit = ols_fit(res.ts, ln_amp);
    res.slope = fit.slope;
    res.slope_se = fit.slope_se;
    res.nu_meas = -fit.slope / (k * k);
    res.rel_err = std::abs(res.nu_meas - nu_th) / nu_th;
    return res;
}

}  // namespace

void run_group_b(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };

    const int n_sw = ctx.is_small() ? 32 : 64;
    const std::vector<double> taus = {0.6, 0.8, 1.0, 1.5};
    const std::vector<CollisionOp> ops = {
        CollisionOp::BGK, CollisionOp::TRT, CollisionOp::MRT, CollisionOp::REGULARIZED
    };
    const std::vector<std::string> stage_ids = {"S09", "S10", "S11", "S12"};
    const std::vector<std::string> colors = {"#2563eb", "#059669", "#d97706", "#7c3aed"};

    double worst_shear_rel = 0.0;
    std::vector<SvgSeries> rel_err_series;
    std::vector<SvgSeries> decay_series;

    for (std::size_t oi = 0; oi < ops.size(); ++oi) {
        auto t0 = clk::now();
        const CollisionOp op = ops[oi];
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(6);
        std::vector<double> x_tau, y_rel;

        for (double tau : taus) {
            const auto sr = measure_shear_wave_2d(n_sw, op, tau, 0.01, 0.0, ctx.threads);
            worst_shear_rel = std::max(worst_shear_rel, sr.rel_err);
            rows.push_back({collision_name(op), fmt_f8(tau), fmt_f8(sr.nu_th), fmt_f8(sr.nu_meas),
                            fmt_e6(sr.rel_err), fmt_e6(sr.slope_se)});
            nums[1].push_back(tau);
            nums[2].push_back(sr.nu_th);
            nums[3].push_back(sr.nu_meas);
            nums[4].push_back(sr.rel_err);
            nums[5].push_back(sr.slope_se);
            x_tau.push_back(tau);
            y_rel.push_back(sr.rel_err * 100.0);

            if (std::abs(tau - 0.8) < 1e-9) {
                decay_series.push_back({collision_name(op), colors[oi], sr.ts, sr.amps, {}, true, true, oi % 2 == 1});
            }
        }
        rel_err_series.push_back({collision_name(op), colors[oi], x_tau, y_rel, {}, true, true, false});

        const std::string csv_name = std::string("out/") +
            (oi == 0 ? "s09_shear_bgk.csv" : oi == 1 ? "s10_shear_trt.csv" : oi == 2 ? "s11_shear_mrt.csv" : "s12_shear_reg.csv");
        ctx.write_csv(stage_ids[oi], csv_name,
                      {"op", "tau", "nu_th", "nu_meas", "rel_err", "slope_se"}, rows, nums);

        int n_files = 1;
        if (oi == 1) {
            ctx.write_svg_xy("S10", "../../../outputs/figures/day009c_05_shear_decay.svg",
                             "Shear Wave Amplitude Decay at tau = 0.8 Across All Collision Operators",
                             "N=64 periodic, initialized with equilibrium + non-equilibrium stress",
                             "Step t", "Fundamental Fourier amplitude |u_hat(k, t)|", decay_series, false, true);
            ++n_files;
        } else if (oi == 3) {
            ctx.write_svg_xy("S12", "../../../outputs/figures/day009c_04_shear_viscosity.svg",
                             "Measured Shear Viscosity Error (%) vs Relaxation Time tau",
                             "All 4 collision operators across tau in {0.6, 0.8, 1.0, 1.5}",
                             "Relaxation time tau", "Relative viscosity error (%)", rel_err_series,
                             false, false, true, 0.5, "0.5% threshold");
            ++n_files;
        }
        ctx.record_stage(stage_ids[oi], std::string("shear_") + collision_name(op), since(t0), n_files);
    }

    ctx.set_criterion("c2", (worst_shear_rel <= 0.005) ? "tercapai" : "belum tercapai",
                      "max_rel_err=" + fmt_f8(worst_shear_rel * 100.0) + "%",
                      "rel_err <= 0.5% across 4 ops x 4 tau");

    // S13-S14: Stability Map across (tau, u0) on Double Shear Layer (BGK vs TRT vs MRT) -> Criterion c12
    {
        auto t0 = clk::now();
        const int n_stab = ctx.is_small() ? 32 : 48;
        const int steps_stab = ctx.is_small() ? 400 : 1000;
        const std::vector<double> st_taus = {0.502, 0.505, 0.51, 0.52, 0.55};
        const std::vector<double> st_u0s = {0.06, 0.10, 0.14, 0.18, 0.22};

        std::vector<std::vector<std::string>> rows_map, rows_sum;
        std::vector<std::vector<double>> nums_map(4), nums_sum(4);
        std::vector<double> hm_bgk, hm_trt;
        int bgk_ok = 0, trt_ok = 0, mrt_ok = 0, total_pts = 0;

        auto test_stability = [&](CollisionOp op, double tau, double u0) -> bool {
            Grid2D g(n_stab, n_stab, Scheme::PULL_2GRID, MemLayout::SOA);
            constexpr double PI = 3.14159265358979323846;
            for (int y = 0; y < n_stab; ++y) {
                const double yf = (y + 0.5) / n_stab;
                const double ux = (yf <= 0.5) ? u0 * std::tanh(40.0 * (yf - 0.25))
                                              : u0 * std::tanh(40.0 * (0.75 - yf));
                for (int x = 0; x < n_stab; ++x) {
                    const double xf = (x + 0.5) / n_stab;
                    const double uy = 0.05 * u0 * std::sin(2.0 * PI * (xf + 0.25));
                    g.set_cell_array(x, y, compute_feq<D2Q9>(1.0, ux, uy));
                }
            }
            CollisionParams cp{op, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
            for (int t = 0; t < steps_stab; ++t) {
                step_periodic_2d(g, cp, ctx.threads);
                if (t % 50 == 49) {
                    double m, px, py, ek, rmin, rmax, umax;
                    bool has_nan = false;
                    g.totals(m, px, py, ek, rmin, rmax, umax, has_nan);
                    if (has_nan || rmin <= 0.1 || rmax > 5.0 || umax > 0.9) return false;
                }
            }
            return true;
        };

        for (double tau : st_taus) {
            for (double u0 : st_u0s) {
                ++total_pts;
                const bool s_bgk = test_stability(CollisionOp::BGK, tau, u0);
                const bool s_trt = test_stability(CollisionOp::TRT, tau, u0);
                const bool s_mrt = test_stability(CollisionOp::MRT, tau, u0);
                bgk_ok += s_bgk ? 1 : 0;
                trt_ok += s_trt ? 1 : 0;
                mrt_ok += s_mrt ? 1 : 0;
                hm_bgk.push_back(s_bgk ? 1.0 : 0.0);
                hm_trt.push_back(s_mrt ? 2.0 : (s_trt ? 1.0 : 0.0));
                rows_map.push_back({fmt_f8(tau), fmt_f8(u0 * std::sqrt(3.0)),
                                    std::to_string(s_bgk ? 1 : 0),
                                    std::to_string(s_trt ? 1 : 0),
                                    std::to_string(s_mrt ? 1 : 0)});
                nums_map[0].push_back(tau);
                nums_map[1].push_back(u0 * std::sqrt(3.0));
            }
        }
        ctx.write_csv("S13", "out/s13_stability_map.csv",
                      {"tau", "ma", "bgk_ok", "trt_ok", "mrt_ok"}, rows_map, nums_map);
        ctx.write_svg_heatmap("S13", "../../../outputs/figures/day009c_07_stability_map.svg",
                              "Double Shear Layer Stability Map (TRT/MRT Score)",
                              "0=Blowup, 1=TRT survives, 2=TRT+MRT survive",
                              "Mach number index", "tau index",
                              static_cast<int>(st_u0s.size()), static_cast<int>(st_taus.size()),
                              hm_trt, "Survival");
        ctx.record_stage("S13", "stability_map", since(t0), 2);

        auto t1 = clk::now();
        const double r_bgk = static_cast<double>(bgk_ok) / total_pts;
        const double r_trt = static_cast<double>(trt_ok) / total_pts;
        const double r_mrt = static_cast<double>(mrt_ok) / total_pts;
        rows_sum = {
            {"BGK", std::to_string(bgk_ok), std::to_string(total_pts), fmt_f8(r_bgk)},
            {"TRT", std::to_string(trt_ok), std::to_string(total_pts), fmt_f8(r_trt)},
            {"MRT", std::to_string(mrt_ok), std::to_string(total_pts), fmt_f8(r_mrt)}
        };
        nums_sum[1] = {static_cast<double>(bgk_ok), static_cast<double>(trt_ok), static_cast<double>(mrt_ok)};
        nums_sum[3] = {r_bgk, r_trt, r_mrt};
        ctx.write_csv("S14", "out/s14_stability_summary.csv",
                      {"op", "survived", "total_pts", "survival_rate"}, rows_sum, nums_sum);
        ctx.write_svg_bar("S14", "../../../outputs/figures/day009c_08_stability_rates.svg",
                          "Double Shear Layer Stability Rate by Collision Operator",
                          "Fraction of (tau, Ma) grid surviving 1000 steps",
                          "Survival fraction",
                          {{"BGK", r_bgk, 0.0, "#dc2626"},
                           {"TRT", r_trt, 0.0, "#059669"},
                           {"MRT", r_mrt, 0.0, "#2563eb"}});
        const bool ok_c12 = (bgk_ok > 0 && bgk_ok < total_pts) && (mrt_ok >= bgk_ok || trt_ok >= bgk_ok);
        ctx.set_criterion("c12", ok_c12 ? "tercapai" : "belum tercapai",
                          "survived BGK=" + std::to_string(bgk_ok) + "/" + std::to_string(total_pts) +
                              ", TRT=" + std::to_string(trt_ok) + ", MRT=" + std::to_string(mrt_ok),
                          "BGK ratio in (0, 1) and TRT/MRT survive >= BGK");
        ctx.record_stage("S14", "stability_summary", since(t1), 2);
    }

    // S15-S16: Galilean invariance with uniform wind u_wind in {0.0, 0.05, 0.10, 0.15} -> Criterion c9
    {
        auto t0 = clk::now();
        const std::vector<double> winds = {0.0, 0.05, 0.10, 0.15};
        std::vector<std::vector<std::string>> rows_g;
        std::vector<std::vector<double>> nums_g(5);
        std::vector<double> ws, e_bgk, e_trt;

        for (double uw : winds) {
            const auto b = measure_shear_wave_2d(n_sw, CollisionOp::BGK, 0.65, 0.005, uw, ctx.threads);
            const auto tr = measure_shear_wave_2d(n_sw, CollisionOp::TRT, 0.65, 0.005, uw, ctx.threads);
            const double ratio = (tr.rel_err > 1e-12) ? (b.rel_err / tr.rel_err) : 1.0;
            rows_g.push_back({fmt_f8(uw), fmt_e6(b.rel_err), fmt_e6(tr.rel_err), fmt_f8(ratio)});
            nums_g[0].push_back(uw);
            nums_g[1].push_back(b.rel_err);
            nums_g[2].push_back(tr.rel_err);
            nums_g[3].push_back(ratio);
            ws.push_back(uw);
            e_bgk.push_back(b.rel_err * 100.0);
            e_trt.push_back(tr.rel_err * 100.0);
        }
        ctx.write_csv("S15", "out/s15_galilean_invariance.csv",
                      {"u_wind", "bgk_rel_err", "trt_rel_err", "ratio_bgk_over_trt"}, rows_g, nums_g);
        ctx.write_svg_xy("S15", "../../../outputs/figures/day009c_06_tau_error.svg",
                         "Galilean Invariance Breakdown vs Background Wind Speed",
                         "Standard D2Q9 lacks O(u^3) Hermite terms; error grows ~ u_wind^2",
                         "Background wind speed u_0", "Viscosity relative error (%)",
                         {{"BGK", "#dc2626", ws, e_bgk, {}, true, true, false},
                          {"TRT (Lambda=3/16)", "#059669", ws, e_trt, {}, true, true, true}});
        ctx.record_stage("S15", "galilean_invariance", since(t0), 2);

        auto t1 = clk::now();
        const double bgk_0 = nums_g[1][0];
        const double bgk_hi = nums_g[1].back();
        const double trt_hi = nums_g[2].back();
        const double ratio_hi = (trt_hi > 1e-12) ? (bgk_hi / trt_hi) : 1.0;
        const bool ok_c9 = (bgk_hi > 2.0 * std::max(bgk_0, 1e-4)) && (ratio_hi >= 2.0);
        std::vector<std::vector<std::string>> rows_ctrl = {
            {"neg_ctrl_u0", fmt_e6(bgk_0), fmt_e6(nums_g[2][0])},
            {"wind_0.15", fmt_e6(bgk_hi), fmt_e6(trt_hi)}
        };
        std::vector<std::vector<double>> nums_ctrl = {{}, {bgk_0, bgk_hi}, {nums_g[2][0], trt_hi}};
        ctx.write_csv("S16", "out/s16_galilean_summary.csv", {"case", "bgk_err", "trt_err"}, rows_ctrl, nums_ctrl);
        ctx.set_criterion("c9", ok_c9 ? "tercapai" : "belum tercapai",
                          "u0=0_err=" + fmt_e6(bgk_0) + ", u0=0.15_BGK=" + fmt_e6(bgk_hi) +
                              ", TRT=" + fmt_e6(trt_hi) + ", ratio=" + fmt_f8(ratio_hi),
                          "error grows with wind (ctrl u0=0) & BGK/TRT >= 2");
        ctx.record_stage("S16", "galilean_summary", since(t1), 1);
    }
}

}  // namespace lbm
