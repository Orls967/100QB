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

namespace {

struct CavityRun {
    int re = 100;
    int n = 128;
    double u_lid = 0.05;
    double tau = 0.7;
    int steps_run = 0;
    double final_res = 0.0;
    double rho_min = 1.0;
    VortexCenter vc{};
    std::vector<double> hist_t;
    std::vector<double> hist_res;
    std::vector<double> u_centerline;
    std::vector<double> v_centerline;
    std::vector<double> speed_map;
};

CavityRun run_cavity_case(int re, int n, double u_lid, CollisionOp op, int max_steps,
                          double tol_res, int threads) {
    const double nu = u_lid * static_cast<double>(n) / static_cast<double>(re);
    const double tau = 0.5 + 3.0 * nu;
    Grid2D g(n, n, Scheme::PULL_2GRID, MemLayout::SOA);
    g.init_uniform(1.0, 0.0, 0.0);
    CollisionParams cp{op, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};

    std::vector<double> prev_ux(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> prev_uy(static_cast<std::size_t>(n) * n, 0.0);

    CavityRun out{};
    out.re = re;
    out.n = n;
    out.u_lid = u_lid;
    out.tau = tau;

    const int check_interval = 500;
    double last_res = 1.0;
    int t_done = 0;
    for (int t = 1; t <= max_steps; ++t) {
        step_cavity_2d(g, cp, u_lid, threads);
        t_done = t;
        if (t % check_interval == 0) {
            double num = 0.0, den = 0.0;
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const std::size_t c = static_cast<std::size_t>(y) * n + x;
                    const auto mc = g.cell_macro(x, y);
                    const double dx = mc.ux - prev_ux[c];
                    const double dy = mc.uy - prev_uy[c];
                    num += dx * dx + dy * dy;
                    den += mc.ux * mc.ux + mc.uy * mc.uy;
                    prev_ux[c] = mc.ux;
                    prev_uy[c] = mc.uy;
                }
            }
            last_res = std::sqrt(num / std::max(den, 1e-24));
            out.hist_t.push_back(static_cast<double>(t));
            out.hist_res.push_back(last_res);
            if (last_res <= tol_res) break;
        }
    }

    out.steps_run = t_done;
    out.final_res = last_res;
    double m, px, py, ek, rmin, rmax, umax;
    bool hn = false;
    g.totals(m, px, py, ek, rmin, rmax, umax, hn);
    out.rho_min = hn ? -1.0 : rmin;
    out.vc = analyze_cavity_vortices(g);

    out.u_centerline.resize(n);
    out.v_centerline.resize(n);
    for (int i = 0; i < n; ++i) {
        out.u_centerline[static_cast<std::size_t>(i)] = g.cell_macro(n / 2, i).ux / u_lid;
        out.v_centerline[static_cast<std::size_t>(i)] = g.cell_macro(i, n / 2).uy / u_lid;
    }
    out.speed_map.resize(static_cast<std::size_t>(n) * n);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            const auto mc = g.cell_macro(x, y);
            out.speed_map[static_cast<std::size_t>(y) * n + x] = std::hypot(mc.ux, mc.uy) / u_lid;
        }
    }
    return out;
}

}  // namespace

void run_group_d(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };

    // ---------- S29-S34: Lid-Driven Cavity Re=100, 400, 1000 -> Criterion c10 ----------
    const int n100 = ctx.is_small() ? 64 : 128;
    const int n400 = ctx.is_small() ? 64 : 128;
    const int n1000 = ctx.is_small() ? 80 : 160;
    const auto cav100 = run_cavity_case(100, n100, 0.05, CollisionOp::TRT,
                                        ctx.is_small() ? 8000 : 50000, 1e-9, ctx.threads);
    const auto cav400 = run_cavity_case(400, n400, 0.07, CollisionOp::MRT,
                                        ctx.is_small() ? 8000 : 35000, 1e-8, ctx.threads);
    const auto cav1000 = run_cavity_case(1000, n1000, 0.08, CollisionOp::MRT,
                                         ctx.is_small() ? 8000 : 45000, 1e-7, ctx.threads);

    // S29: Re=100 residual convergence
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(2);
        for (std::size_t i = 0; i < cav100.hist_t.size(); ++i) {
            rows.push_back({std::to_string(static_cast<int>(cav100.hist_t[i])), fmt_e6(cav100.hist_res[i])});
            nums[0].push_back(cav100.hist_t[i]);
            nums[1].push_back(cav100.hist_res[i]);
        }
        ctx.write_csv("S29", "out/s29_cavity_re100_residual.csv", {"step", "rel_residual"}, rows, nums);
        ctx.write_svg_xy("S29", "../../../outputs/figures/day009c_19_cavity_residual.svg",
                         "Lid-Driven Cavity (Re = 100) Steady-State Residual Convergence",
                         "Final residual = " + fmt_e6(cav100.final_res),
                         "Time step t", "Relative L2 velocity change",
                         {{"Re=100 TRT", "#2563eb", cav100.hist_t, cav100.hist_res, {}, true, true, false}},
                         false, true, true, 1e-9, "Steady threshold 1e-9");
        ctx.record_stage("S29", "cavity_re100_conv", since(t0), 2);
    }

    // S30: Centerline u(y) and v(x) profiles for Re=100
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(3);
        std::vector<double> coord(cav100.n);
        for (int i = 0; i < cav100.n; ++i) {
            coord[static_cast<std::size_t>(i)] = (i + 0.5) / cav100.n;
            rows.push_back({fmt_f8(coord[i]), fmt_f8(cav100.u_centerline[i]), fmt_f8(cav100.v_centerline[i])});
            nums[0].push_back(coord[i]);
            nums[1].push_back(cav100.u_centerline[i]);
            nums[2].push_back(cav100.v_centerline[i]);
        }
        ctx.write_csv("S30", "out/s30_cavity_re100_centerline.csv", {"coord", "u_mid_x", "v_mid_y"}, rows, nums);
        ctx.write_svg_xy("S30", "../../../outputs/figures/day009c_20_cavity_centerlines.svg",
                         "Cavity Re = 100 Centerline Velocity Profiles u(0.5, y) and v(x, 0.5)",
                         "Normalized by lid speed U_lid = 0.05",
                         "Normalized coordinate", "Normalized velocity",
                         {{"u(0.5, y) / U", "#2563eb", coord, cav100.u_centerline, {}, true, false, false},
                          {"v(x, 0.5) / U", "#dc2626", coord, cav100.v_centerline, {}, true, false, true}});
        ctx.record_stage("S30", "cavity_re100_centerlines", since(t0), 2);
    }

    // S31: Re=100 speed heatmap
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(3);
        for (int i = 0; i < cav100.n; i += std::max(1, cav100.n / 16)) {
            rows.push_back({std::to_string(i), fmt_f8(cav100.u_centerline[i]), fmt_f8(cav100.speed_map[(cav100.n / 2) * cav100.n + i])});
            nums[0].push_back(i);
            nums[1].push_back(cav100.u_centerline[i]);
            nums[2].push_back(cav100.speed_map[(cav100.n / 2) * cav100.n + i]);
        }
        ctx.write_csv("S31", "out/s31_cavity_re100_slice.csv", {"idx", "u_mid", "speed_mid"}, rows, nums);
        ctx.write_svg_heatmap("S31", "../../../outputs/figures/day009c_21_cavity_re100_map.svg",
                              "Lid-Driven Cavity (Re = 100) Speed Magnitude |u| / U_lid",
                              "Primary vortex center=(" + fmt_f8(cav100.vc.x) + ", " + fmt_f8(cav100.vc.y) + ")",
                              "x", "y", cav100.n, cav100.n, cav100.speed_map, "|u|/U");
        ctx.record_stage("S31", "cavity_re100_map", since(t0), 2);
    }

    // S32: Re=400 cavity
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(3);
        for (int i = 0; i < cav400.n; ++i) {
            rows.push_back({fmt_f8((i + 0.5) / cav400.n), fmt_f8(cav400.u_centerline[i]), fmt_f8(cav400.v_centerline[i])});
            nums[0].push_back((i + 0.5) / cav400.n);
            nums[1].push_back(cav400.u_centerline[i]);
            nums[2].push_back(cav400.v_centerline[i]);
        }
        ctx.write_csv("S32", "out/s32_cavity_re400.csv", {"coord", "u_mid", "v_mid"}, rows, nums);
        ctx.write_svg_heatmap("S32", "../../../outputs/figures/day009c_22_cavity_re400_map.svg",
                              "Lid-Driven Cavity (Re = 400, MRT) Speed Magnitude |u| / U_lid",
                              "Primary vortex center=(" + fmt_f8(cav400.vc.x) + ", " + fmt_f8(cav400.vc.y) + ")",
                              "x", "y", cav400.n, cav400.n, cav400.speed_map, "|u|/U");
        ctx.record_stage("S32", "cavity_re400", since(t0), 2);
    }

    // S33: Re=1000 cavity
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(3);
        for (int i = 0; i < cav1000.n; ++i) {
            rows.push_back({fmt_f8((i + 0.5) / cav1000.n), fmt_f8(cav1000.u_centerline[i]), fmt_f8(cav1000.v_centerline[i])});
            nums[0].push_back((i + 0.5) / cav1000.n);
            nums[1].push_back(cav1000.u_centerline[i]);
            nums[2].push_back(cav1000.v_centerline[i]);
        }
        ctx.write_csv("S33", "out/s33_cavity_re1000.csv", {"coord", "u_mid", "v_mid"}, rows, nums);
        ctx.write_svg_heatmap("S33", "../../../outputs/figures/day009c_23_cavity_re1000_map.svg",
                              "Lid-Driven Cavity (Re = 1000, MRT) Speed Magnitude |u| / U_lid",
                              "Primary vortex center=(" + fmt_f8(cav1000.vc.x) + ", " + fmt_f8(cav1000.vc.y) + ")",
                              "x", "y", cav1000.n, cav1000.n, cav1000.speed_map, "|u|/U");
        ctx.record_stage("S33", "cavity_re1000", since(t0), 2);
    }

    // S34: Cavity summary & Ghia benchmark -> c10
    {
        auto t0 = clk::now();
        const double dist100 = std::hypot(cav100.vc.x - 0.6172, cav100.vc.y - 0.7344);
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(7);
        for (const auto& c : {cav100, cav400, cav1000}) {
            rows.push_back({std::to_string(c.re), std::to_string(c.n), fmt_f8(c.tau),
                            fmt_f8(c.vc.x), fmt_f8(c.vc.y), std::to_string(c.vc.corner_eddies),
                            fmt_f8(c.rho_min), fmt_e6(c.final_res)});
            nums[0].push_back(c.re);
            nums[3].push_back(c.vc.x);
            nums[4].push_back(c.vc.y);
            nums[5].push_back(c.vc.corner_eddies);
            nums[6].push_back(c.rho_min);
        }
        ctx.write_csv("S34", "out/s34_cavity_summary.csv",
                      {"re", "n", "tau", "vc_x", "vc_y", "eddies", "rho_min", "residual"}, rows, nums);
        ctx.write_svg_xy("S34", "../../../outputs/figures/day009c_24_cavity_vortex_centers.svg",
                         "Cavity Primary Vortex Center Migration Re = 100 -> 400 -> 1000",
                         "Ghia et al. (1982) literature reference at Re=100: (0.6172, 0.7344)",
                         "Vortex center x_c", "Vortex center y_c",
                         {{"LBM (100, 400, 1000)", "#2563eb",
                           {cav100.vc.x, cav400.vc.x, cav1000.vc.x},
                           {cav100.vc.y, cav400.vc.y, cav1000.vc.y}, {}, true, true, false},
                          {"Ghia Re=100 lit.", "#dc2626", {0.6172}, {0.7344}, {}, false, true, false}});
        const bool ok_c10 = (dist100 <= 0.03) && (cav100.final_res <= 1e-9) &&
                            (cav400.rho_min > 0.0) && (cav1000.rho_min > 0.0);
        ctx.set_criterion("c10", ok_c10 ? "tercapai" : "belum tercapai",
                          "Re100=(" + fmt_f8(cav100.vc.x) + "," + fmt_f8(cav100.vc.y) +
                              ") dist=" + fmt_f8(dist100) + " res=" + fmt_e6(cav100.final_res) +
                              "; Re400=(" + fmt_f8(cav400.vc.x) + "," + fmt_f8(cav400.vc.y) + ") eddies=" +
                              std::to_string(cav400.vc.corner_eddies) +
                              "; Re1000=(" + fmt_f8(cav1000.vc.x) + "," + fmt_f8(cav1000.vc.y) + ") eddies=" +
                              std::to_string(cav1000.vc.corner_eddies),
                          "Re100 dist<=0.03 of (0.6172,0.7344) & res<=1e-9; Re400/1000 rho_min>0");
        ctx.record_stage("S34", "cavity_summary", since(t0), 2);
    }

    // ---------- S35-S40: Cylinder Re=20, Re=100, and Bouzidi vs Staircase -> Criterion c11 ----------
    auto run_cyl_sim = [&](double re, int d_diam, bool use_bouzidi, int n_steps, double u_in,
                           bool trigger_shedding, std::vector<double>* out_ts = nullptr,
                           std::vector<double>* out_cd = nullptr, std::vector<double>* out_cl = nullptr,
                           std::vector<double>* out_speed = nullptr, int* out_nx = nullptr,
                           int* out_ny = nullptr, double* out_lr = nullptr) -> ForceResult {
        const int nx = 16 * d_diam;
        const int ny = 10 * d_diam;
        const double cx = 4.0 * d_diam + 0.23;
        const double cy = 0.5 * ny + (trigger_shedding ? 0.31 : 0.17);
        const double rad = 0.5 * d_diam;
        const double nu = u_in * d_diam / re;
        const double tau = 0.5 + 3.0 * nu;

        Grid2D g(nx, ny, Scheme::PULL_2GRID, MemLayout::SOA);
        // Start smoothly from rest and ramp inlet velocity over first 600 steps to avoid acoustic shock!
        g.init_uniform(1.0, 0.0, 0.0);
        const auto dom = build_cylinder_domain(g, cx, cy, rad, use_bouzidi);
        CollisionParams cp{CollisionOp::MRT, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};

        ForceResult last{};
        double cd_avg_sum = 0.0;
        int cd_avg_cnt = 0;
        const int avg_win = std::max(400, n_steps / 5);
        for (int t = 1; t <= n_steps; ++t) {
            const double u_cur = u_in * std::min(1.0, static_cast<double>(t) / 600.0);
            last = step_cylinder_2d(g, dom, cp, u_cur, 1.0, ctx.threads);
            if (trigger_shedding && t == 700) {
                for (int y = ny / 2 - d_diam; y <= ny / 2 + d_diam; ++y) {
                    if (y <= 0 || y >= ny - 1) continue;
                    for (int x = static_cast<int>(cx + d_diam); x < static_cast<int>(cx + 2 * d_diam); ++x) {
                        const double dy = (y - cy) / d_diam;
                        g.set_cell_array(x, y, compute_feq<D2Q9>(1.0, u_in, 0.12 * u_in * std::exp(-dy * dy)));
                    }
                }
            }
            if (out_ts && t > n_steps / 2) {
                out_ts->push_back(static_cast<double>(t));
                out_cd->push_back(last.cd);
                out_cl->push_back(last.cl);
            }
            if (t > n_steps - avg_win) {
                cd_avg_sum += last.cd;
                ++cd_avg_cnt;
            }
        }
        last.cd = cd_avg_sum / std::max(1, cd_avg_cnt);
        if (out_lr) *out_lr = measure_recirculation_length(g, cx, cy, static_cast<double>(d_diam));
        if (out_speed) {
            *out_nx = nx;
            *out_ny = ny;
            out_speed->resize(static_cast<std::size_t>(nx) * ny);
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    const std::size_t c = static_cast<std::size_t>(y) * nx + x;
                    if (g.solid[c]) {
                        (*out_speed)[c] = 0.0;
                    } else {
                        const auto mc = g.cell_macro(x, y);
                        (*out_speed)[c] = std::hypot(mc.ux, mc.uy) / u_in;
                    }
                }
            }
        }
        return last;
    };

    // S35-S36: Cylinder Re=20 steady drag & recirculation
    double cd20 = 0.0, lr20 = 0.0;
    {
        auto t0 = clk::now();
        std::vector<double> ts20, cd_tr, cl_tr, spd20;
        int nx20 = 0, ny20 = 0;
        const int steps20 = ctx.is_small() ? 3000 : 10000;
        const auto f20 = run_cyl_sim(20.0, 14, true, steps20, 0.04, false,
                                     &ts20, &cd_tr, &cl_tr, &spd20, &nx20, &ny20, &lr20);
        cd20 = f20.cd;
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(3);
        const std::size_t stride = std::max<std::size_t>(1, ts20.size() / 40);
        for (std::size_t i = 0; i < ts20.size(); i += stride) {
            rows.push_back({std::to_string(static_cast<int>(ts20[i])), fmt_f8(cd_tr[i]), fmt_e6(cl_tr[i])});
            nums[0].push_back(ts20[i]);
            nums[1].push_back(cd_tr[i]);
            nums[2].push_back(cl_tr[i]);
        }
        ctx.write_csv("S35", "out/s35_cylinder_re20_trace.csv", {"step", "cd", "cl"}, rows, nums);
        ctx.write_svg_xy("S35", "../../../outputs/figures/day009c_25_cylinder_re20_drag.svg",
                         "Cylinder Cross-Flow (Re = 20, D=14, Bouzidi) Drag Convergence",
                         "Steady C_d = " + fmt_f8(cd20) + ", wake L_r/D = " + fmt_f8(lr20),
                         "Time step t", "Drag coefficient C_d",
                         {{"C_d(t)", "#2563eb", nums[0], nums[1], {}, true, false, false}});
        ctx.record_stage("S35", "cylinder_re20_drag", since(t0), 2);

        auto t1 = clk::now();
        std::vector<std::vector<std::string>> rows_wake;
        std::vector<std::vector<double>> nums_wake(2);
        for (int x = 0; x < nx20; x += 4) {
            rows_wake.push_back({std::to_string(x), fmt_f8(spd20[static_cast<std::size_t>(ny20 / 2) * nx20 + x])});
            nums_wake[0].push_back(x);
            nums_wake[1].push_back(spd20[static_cast<std::size_t>(ny20 / 2) * nx20 + x]);
        }
        ctx.write_csv("S36", "out/s36_cylinder_re20_centerwake.csv", {"x", "speed_mid"}, rows_wake, nums_wake);
        ctx.write_svg_heatmap("S36", "../../../outputs/figures/day009c_26_cylinder_re20_map.svg",
                              "Cylinder Wake Speed Field at Re = 20 (Steady Twin Eddies)",
                              "Bouzidi curved boundary (H/D = 10)",
                              "x", "y", nx20, ny20, spd20, "|u|/U");
        ctx.record_stage("S36", "cylinder_re20_wake", since(t1), 2);
    }

    // S37 & S39: Cylinder Re=100 Karman vortex shedding & Strouhal number
    double st100 = 0.0, cd100 = 0.0, cl_amp100 = 0.0;
    {
        auto t0 = clk::now();
        std::vector<double> ts100, cd_tr, cl_tr, spd100;
        int nx100 = 0, ny100 = 0;
        const int d100 = 16;
        const double u100 = 0.055;
        const int steps100 = ctx.is_small() ? 5000 : 18000;
        const auto f100 = run_cyl_sim(100.0, d100, true, steps100, u100, true,
                                      &ts100, &cd_tr, &cl_tr, &spd100, &nx100, &ny100, nullptr);
        // Analyze tail 50% of recorded shedding window
        const std::size_t tail_beg = ts100.size() / 2;
        std::vector<double> cl_tail(cl_tr.begin() + static_cast<std::ptrdiff_t>(tail_beg), cl_tr.end());
        const auto st = analyze_strouhal(cl_tail, static_cast<double>(d100), u100);
        st100 = st.strouhal;
        cl_amp100 = st.cl_amp;
        double cd_sum = 0.0;
        for (std::size_t i = tail_beg; i < cd_tr.size(); ++i) cd_sum += cd_tr[i];
        cd100 = cd_sum / std::max<std::size_t>(1, cd_tr.size() - tail_beg);
        (void)f100;

        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(3);
        const std::size_t stride = std::max<std::size_t>(1, ts100.size() / 80);
        for (std::size_t i = 0; i < ts100.size(); i += stride) {
            rows.push_back({std::to_string(static_cast<int>(ts100[i])), fmt_f8(cd_tr[i]), fmt_f8(cl_tr[i])});
            nums[0].push_back(ts100[i]);
            nums[1].push_back(cd_tr[i]);
            nums[2].push_back(cl_tr[i]);
        }
        ctx.write_csv("S37", "out/s37_cylinder_re100_shedding.csv", {"step", "cd", "cl"}, rows, nums);
        ctx.write_svg_xy("S37", "../../../outputs/figures/day009c_27_cylinder_re100_cl.svg",
                         "Cylinder Karman Vortex Street at Re = 100 (Lift Oscillations)",
                         "Strouhal St = " + fmt_f8(st100) + ", C_d_mean = " + fmt_f8(cd100),
                         "Time step t", "Lift coefficient C_l",
                         {{"C_l(t)", "#dc2626", nums[0], nums[2], {}, true, false, false}});
        ctx.record_stage("S37", "cylinder_re100_shedding", since(t0), 2);

        auto t1 = clk::now();
        std::vector<std::vector<std::string>> rows_st = {
            {fmt_f8(st.freq_zc), fmt_f8(st.freq_dft), fmt_f8(st100), fmt_f8(cl_amp100), fmt_f8(cd100)}
        };
        std::vector<std::vector<double>> nums_st = {{st.freq_zc}, {st.freq_dft}, {st100}, {cl_amp100}, {cd100}};
        ctx.write_csv("S39", "out/s39_cylinder_re100_strouhal.csv",
                      {"freq_zc", "freq_dft", "strouhal", "cl_amp", "cd_mean"}, rows_st, nums_st);
        ctx.write_svg_heatmap("S39", "../../../outputs/figures/day009c_29_cylinder_re100_map.svg",
                              "Cylinder Karman Vortex Street Speed Field (Re = 100, D=16)",
                              "Strouhal St = " + fmt_f8(st100),
                              "x", "y", nx100, ny100, spd100, "|u|/U");
        ctx.record_stage("S39", "cylinder_re100_map", since(t1), 2);
    }

    // S38 & S40: Bouzidi vs Staircase across >=3 resolutions (D in {10, 14, 18}) vs fine Bouzidi D=24
    {
        auto t0 = clk::now();
        const int d_ref = ctx.is_small() ? 18 : 24;
        const int st_ref = ctx.is_small() ? 3500 : 8000;
        const double cd_ref = run_cyl_sim(20.0, d_ref, true, st_ref, 0.04, false).cd;
        const std::vector<int> ds = {10, 14, 18};
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(5);
        std::vector<double> xd, eb, es;
        int better_count = 0;

        for (int d : ds) {
            const int st = ctx.is_small() ? 3500 : 7500;
            const double cb = run_cyl_sim(20.0, d, true, st, 0.04, false).cd;
            const double cs = run_cyl_sim(20.0, d, false, st, 0.04, false).cd;
            const double err_b = std::abs(cb - cd_ref);
            const double err_s = std::abs(cs - cd_ref);
            if (err_b <= 0.5 * err_s) ++better_count;
            rows.push_back({std::to_string(d), fmt_f8(cb), fmt_f8(cs), fmt_e6(err_b), fmt_e6(err_s)});
            nums[0].push_back(d);
            nums[1].push_back(cb);
            nums[2].push_back(cs);
            nums[3].push_back(err_b);
            nums[4].push_back(err_s);
            xd.push_back(d);
            eb.push_back(err_b);
            es.push_back(err_s);
        }
        ctx.write_csv("S38", "out/s38_bouzidi_vs_staircase.csv",
                      {"diameter", "cd_bouzidi", "cd_staircase", "err_bouzidi", "err_staircase"}, rows, nums);
        ctx.write_svg_xy("S38", "../../../outputs/figures/day009c_28_bouzidi_convergence.svg",
                         "Cylinder Drag Error Convergence: Bouzidi Curved Wall vs Staircase",
                         "Re = 20 vs fine Bouzidi reference (D=" + std::to_string(d_ref) + ")",
                         "Cylinder diameter D (lattice cells)", "|C_d - C_d_ref|",
                         {{"Bouzidi curved", "#059669", xd, eb, {}, true, true, false},
                          {"Staircase BB", "#dc2626", xd, es, {}, true, true, true}},
                         true, true);
        ctx.record_stage("S38", "bouzidi_vs_staircase", since(t0), 2);

        auto t1 = clk::now();
        std::vector<std::vector<std::string>> r40 = {
            {"Re20_Bouzidi", fmt_f8(cd20), fmt_f8(lr20), "0.0"},
            {"Re100_Bouzidi", fmt_f8(cd100), "shedding", fmt_f8(st100)}
        };
        std::vector<std::vector<double>> n40 = {{}, {cd20, cd100}, {}, {0.0, st100}};
        ctx.write_csv("S40", "out/s40_cylinder_summary.csv", {"case", "cd", "wake_lr_d", "strouhal"}, r40, n40);
        ctx.write_svg_bar("S40", "../../../outputs/figures/day009c_30_bouzidi_ratio.svg",
                          "Bouzidi vs Staircase Drag Error at Three Resolutions",
                          "Bouzidi achieves <= 0.5x staircase error at " + std::to_string(better_count) + "/3 resolutions",
                          "Relative error |C_d - C_d_ref| / C_d_ref",
                          {{"D8-Bouz", eb[0] / cd_ref, 0.0, "#059669"},
                           {"D8-Stair", es[0] / cd_ref, 0.0, "#dc2626"},
                           {"D12-Bouz", eb[1] / cd_ref, 0.0, "#059669"},
                           {"D12-Stair", es[1] / cd_ref, 0.0, "#dc2626"},
                           {"D16-Bouz", eb[2] / cd_ref, 0.0, "#059669"},
                           {"D16-Stair", es[2] / cd_ref, 0.0, "#dc2626"}});
        const bool ok_c11 = (cd20 >= 1.8 && cd20 <= 2.6) && (lr20 >= 0.6 && lr20 <= 1.2) &&
                            (st100 >= 0.14 && st100 <= 0.20) && (better_count == 3);
        ctx.set_criterion("c11", ok_c11 ? "tercapai" : "belum tercapai",
                          "Re20: Cd=" + fmt_f8(cd20) + ", Lr/D=" + fmt_f8(lr20) +
                              "; Re100: St=" + fmt_f8(st100) + "; Bouzidi<=0.5x=" +
                              std::to_string(better_count) + "/3",
                          "Cd in [1.8,2.6], Lr/D in [0.6,1.2], St in [0.14,0.20], Bouzidi<=0.5x in >=3 res");
        ctx.record_stage("S40", "cylinder_summary", since(t1), 2);
    }
}

}  // namespace lbm
