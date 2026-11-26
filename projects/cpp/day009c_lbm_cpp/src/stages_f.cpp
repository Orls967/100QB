#include "lbm/stages.hpp"
#include "lbm/analysis.hpp"
#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/grid.hpp"
#include "lbm/rng.hpp"
#include <chrono>
#include <cmath>
#include <vector>

namespace lbm {

namespace {

Grid2D make_perturbed_2d(int nx, int ny, Scheme sc, MemLayout lay, int tw = 16, int th = 16) {
    Grid2D g(nx, ny, sc, lay, tw, th);
    Xoshiro256pp rng(derive_seed(47, 1));
    constexpr double PI = 3.14159265358979323846;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const double kx = 2.0 * PI * x / nx;
            const double ky = 2.0 * PI * y / ny;
            const double rho = 1.0 + 0.015 * std::sin(kx) * std::cos(ky) + 0.001 * rng.uniform_sym();
            const double ux = 0.03 * std::sin(ky);
            const double uy = 0.025 * std::cos(kx);
            g.set_cell_array(x, y, compute_feq<D2Q9>(rho, ux, uy));
        }
    }
    return g;
}

double max_pop_diff_2d(const Grid2D& a, const Grid2D& b) {
    double md = 0.0;
    for (int y = 0; y < a.ny; ++y) {
        for (int x = 0; x < a.nx; ++x) {
            for (int q = 0; q < 9; ++q) {
                md = std::max(md, std::abs(a.get_pop(x, y, q) - b.get_pop(x, y, q)));
            }
        }
    }
    return md;
}

double max_pop_diff_3d(const Grid3D& a, const Grid3D& b) {
    double md = 0.0;
    for (int z = 0; z < a.nz; ++z) {
        for (int y = 0; y < a.ny; ++y) {
            for (int x = 0; x < a.nx; ++x) {
                for (int q = 0; q < 19; ++q) {
                    md = std::max(md, std::abs(a.get_pop(x, y, z, q) - b.get_pop(x, y, z, q)));
                }
            }
        }
    }
    return md;
}

}  // namespace

std::uint64_t benchmark_thread_parity_hash(int n_threads) {
    Grid2D g = make_perturbed_2d(64, 64, Scheme::AA_PATTERN, MemLayout::SOA);
    CollisionParams cp{CollisionOp::TRT, 0.78, 3.0 / 16.0, 0.0, 0.0, 0.0};
    for (int t = 0; t < 200; ++t) step_periodic_2d(g, cp, n_threads);
    std::uint64_t h2 = g.field_hash_fnv1a();

    Grid3D g3(16, 16, 16, Scheme::AA_PATTERN);
    for (int z = 0; z < 16; ++z)
        for (int y = 0; y < 16; ++y)
            for (int x = 0; x < 16; ++x) {
                const auto feq = compute_feq<D3Q19>(1.0 + 0.01 * std::sin(0.3 * x), 0.02, -0.01, 0.015);
                for (int q = 0; q < 19; ++q) g3.set_pop(x, y, z, q, feq[static_cast<std::size_t>(q)]);
            }
    for (int t = 0; t < 40; ++t) step_periodic_3d(g3, cp, n_threads);
    const std::uint64_t h3 = g3.field_hash_fnv1a();
    return fnv1a_bytes(&h3, sizeof(h3), h2);
}

void run_group_f(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };

    const int steps_eq = ctx.is_small() ? 100 : 1000;
    const CollisionParams cp{CollisionOp::BGK, 0.8, 3.0 / 16.0, 0.0, 0.0, 0.0};

    // S47: 2D (64x64) Push vs Pull vs AA-pattern parity after 1000 steps
    bool ok_2d_scheme = false;
    {
        auto t0 = clk::now();
        Grid2D g_push = make_perturbed_2d(64, 64, Scheme::PUSH_2GRID, MemLayout::SOA);
        Grid2D g_pull = make_perturbed_2d(64, 64, Scheme::PULL_2GRID, MemLayout::SOA);
        Grid2D g_aa = make_perturbed_2d(64, 64, Scheme::AA_PATTERN, MemLayout::SOA);
        for (int t = 0; t < steps_eq; ++t) {
            step_periodic_2d(g_push, cp, ctx.threads);
            step_periodic_2d(g_pull, cp, ctx.threads);
            step_periodic_2d(g_aa, cp, ctx.threads);
        }
        const auto hp = g_push.field_hash_fnv1a();
        const auto hl = g_pull.field_hash_fnv1a();
        const auto ha = g_aa.field_hash_fnv1a();
        const double d_pl = max_pop_diff_2d(g_push, g_pull);
        const double d_pa = max_pop_diff_2d(g_push, g_aa);
        ok_2d_scheme = (hp == hl) && (hp == ha) && (d_pl == 0.0) && (d_pa == 0.0);

        std::vector<std::vector<std::string>> rows = {
            {"two_grid_push", std::to_string(steps_eq), hex64(hp), fmt_e6(0.0)},
            {"two_grid_pull", std::to_string(steps_eq), hex64(hl), fmt_e6(d_pl)},
            {"aa_pattern", std::to_string(steps_eq), hex64(ha), fmt_e6(d_pa)}
        };
        std::vector<std::vector<double>> nums = {{}, {}, {}, {0.0, d_pl, d_pa}};
        ctx.write_csv("S47", "out/s47_scheme_parity_2d.csv", {"scheme", "steps", "fnv1a", "max_diff"}, rows, nums);
        ctx.record_stage("S47", "scheme_parity_2d", since(t0), 1);
    }

    // S48: 3D (24^3) Push vs Pull vs AA-pattern parity after 1000 steps (full)
    bool ok_3d_scheme = false;
    {
        auto t0 = clk::now();
        const int n3 = ctx.is_small() ? 12 : 24;
        const int s3 = ctx.is_small() ? 50 : 1000;
        auto init_3d = [&](Scheme sc) {
            Grid3D g(n3, n3, n3, sc);
            Xoshiro256pp rng(derive_seed(48, 1));
            for (int z = 0; z < n3; ++z)
                for (int y = 0; y < n3; ++y)
                    for (int x = 0; x < n3; ++x) {
                        const auto feq = compute_feq<D3Q19>(
                            1.0 + 0.01 * rng.uniform_sym(), 0.02 * std::sin(0.25 * z), -0.01, 0.015);
                        for (int q = 0; q < 19; ++q) g.set_pop(x, y, z, q, feq[static_cast<std::size_t>(q)]);
                    }
            return g;
        };
        Grid3D g3p = init_3d(Scheme::PUSH_2GRID);
        Grid3D g3l = init_3d(Scheme::PULL_2GRID);
        Grid3D g3a = init_3d(Scheme::AA_PATTERN);
        for (int t = 0; t < s3; ++t) {
            step_periodic_3d(g3p, cp, ctx.threads);
            step_periodic_3d(g3l, cp, ctx.threads);
            step_periodic_3d(g3a, cp, ctx.threads);
        }
        const auto hp = g3p.field_hash_fnv1a();
        const auto hl = g3l.field_hash_fnv1a();
        const auto ha = g3a.field_hash_fnv1a();
        const double d_pl = max_pop_diff_3d(g3p, g3l);
        const double d_pa = max_pop_diff_3d(g3p, g3a);
        ok_3d_scheme = (hp == hl) && (hp == ha) && (d_pl == 0.0) && (d_pa == 0.0);

        std::vector<std::vector<std::string>> rows = {
            {"two_grid_push_3d", std::to_string(s3), hex64(hp), fmt_e6(0.0)},
            {"two_grid_pull_3d", std::to_string(s3), hex64(hl), fmt_e6(d_pl)},
            {"aa_pattern_3d", std::to_string(s3), hex64(ha), fmt_e6(d_pa)}
        };
        std::vector<std::vector<double>> nums = {{}, {}, {}, {0.0, d_pl, d_pa}};
        ctx.write_csv("S48", "out/s48_scheme_parity_3d.csv", {"scheme", "steps", "fnv1a", "max_diff"}, rows, nums);
        ctx.record_stage("S48", "scheme_parity_3d", since(t0), 1);
    }

    // S49-S50: Layout parity (AoS = SoA = Tile) & Tile Size Timing Sweep
    bool ok_layout = false;
    {
        auto t0 = clk::now();
        Grid2D g_aos = make_perturbed_2d(64, 64, Scheme::PULL_2GRID, MemLayout::AOS);
        Grid2D g_soa = make_perturbed_2d(64, 64, Scheme::PULL_2GRID, MemLayout::SOA);
        Grid2D g_tile = make_perturbed_2d(64, 64, Scheme::PULL_2GRID, MemLayout::TILE, 16, 16);
        for (int t = 0; t < steps_eq; ++t) {
            step_periodic_2d(g_aos, cp, ctx.threads);
            step_periodic_2d(g_soa, cp, ctx.threads);
            step_periodic_2d(g_tile, cp, ctx.threads);
        }
        const auto h_aos = g_aos.field_hash_fnv1a();
        const auto h_soa = g_soa.field_hash_fnv1a();
        const auto h_tile = g_tile.field_hash_fnv1a();
        ok_layout = (h_aos == h_soa) && (h_aos == h_tile);

        ctx.write_csv("S49", "out/s49_layout_parity.csv",
                      {"layout", "steps", "fnv1a", "identical_to_aos"},
                      {{"AoS", std::to_string(steps_eq), hex64(h_aos), "1"},
                       {"SoA", std::to_string(steps_eq), hex64(h_soa), h_aos == h_soa ? "1" : "0"},
                       {"Tile16x16", std::to_string(steps_eq), hex64(h_tile), h_aos == h_tile ? "1" : "0"}},
                      {{}, {}, {}, {1.0, h_aos == h_soa ? 1.0 : 0.0, h_aos == h_tile ? 1.0 : 0.0}});
        ctx.record_stage("S49", "layout_parity", since(t0), 1);

        auto t1 = clk::now();
        const std::vector<int> tws = {4, 8, 16, 32};
        std::vector<std::vector<std::string>> rows_ts;
        std::vector<std::vector<double>> nums_ts(4);
        std::vector<SvgBar> bars;
        for (int tw : tws) {
            Grid2D gt = make_perturbed_2d(128, 128, Scheme::PULL_2GRID, MemLayout::TILE, tw, tw);
            const auto b0 = clk::now();
            const int tst = ctx.is_small() ? 30 : 100;
            for (int t = 0; t < tst; ++t) step_periodic_2d(gt, cp, ctx.threads);
            const double ms = std::chrono::duration<double, std::milli>(clk::now() - b0).count();
            const double mlups = (128.0 * 128.0 * tst) / (ms * 1000.0);
            rows_ts.push_back({std::to_string(tw), hex64(gt.field_hash_fnv1a()), fmt_f8(ms), fmt_f8(mlups)});
            nums_ts[0].push_back(tw);
            nums_ts[2].push_back(ms);
            nums_ts[3].push_back(mlups);
            bars.push_back({"Tile " + std::to_string(tw) + "x" + std::to_string(tw), mlups, 0.0, "#2563eb"});
        }
        ctx.write_csv("S50", "out/s50_tile_sweep.csv", {"tile_w", "fnv1a", "elapsed_ms", "mlups"},
                      rows_ts, nums_ts, false);
        ctx.write_svg_bar("S50", "../../../outputs/figures/day009c_36_tile_sweep.svg",
                          "Blocked Tile Layout Performance vs Tile Size (128x128)",
                          "Identical field hashes across all tile sizes",
                          "Throughput (MLUPS)", bars, false);
        ctx.record_stage("S50", "tile_size_sweep", since(t1), 2);
    }

    ctx.set_criterion("c14", (ok_2d_scheme && ok_3d_scheme && ok_layout) ? "tercapai" : "belum tercapai",
                      "2D_push_pull_aa=" + std::string(ok_2d_scheme ? "identical" : "diff") +
                          ", 3D=" + std::string(ok_3d_scheme ? "identical" : "diff") +
                          ", AoS_SoA_Tile=" + std::string(ok_layout ? "identical" : "diff"),
                      "bit-identical after 1000 steps (2D & 3D) and across layouts");

    // S51-S52: Thread 1/2/4/8 determinism & scaling
    {
        auto t0 = clk::now();
        ctx.thread_hash_1 = benchmark_thread_parity_hash(1);
        ctx.thread_hash_2 = benchmark_thread_parity_hash(2);
        ctx.thread_hash_4 = benchmark_thread_parity_hash(4);
        ctx.thread_hash_8 = benchmark_thread_parity_hash(8);
        ctx.threads_1_2_4_8_identical = (ctx.thread_hash_1 == ctx.thread_hash_2) &&
                                        (ctx.thread_hash_1 == ctx.thread_hash_4) &&
                                        (ctx.thread_hash_1 == ctx.thread_hash_8);

        ctx.write_csv("S51", "out/s51_thread_parity.csv", {"threads", "fnv1a", "bit_exact"},
                      {{"1", hex64(ctx.thread_hash_1), "1"},
                       {"2", hex64(ctx.thread_hash_2), ctx.thread_hash_1 == ctx.thread_hash_2 ? "1" : "0"},
                       {"4", hex64(ctx.thread_hash_4), ctx.thread_hash_1 == ctx.thread_hash_4 ? "1" : "0"},
                       {"8", hex64(ctx.thread_hash_8), ctx.thread_hash_1 == ctx.thread_hash_8 ? "1" : "0"}},
                      {{1, 2, 4, 8}, {}, {1, 1, 1, 1}});
        ctx.set_criterion("c15", ctx.threads_1_2_4_8_identical ? "tercapai" : "belum tercapai",
                          "1/2/4/8_threads=" + hex64(ctx.thread_hash_1) +
                              " (identical=" + std::string(ctx.threads_1_2_4_8_identical ? "true" : "false") + ")",
                          "thread 1/2/4/8 bit-identical");
        ctx.record_stage("S51", "thread_parity", since(t0), 1);

        auto t1 = clk::now();
        std::vector<std::vector<std::string>> rows_sc;
        std::vector<std::vector<double>> nums_sc(4);
        std::vector<SvgBar> bars;
        double base_ms = 1.0;
        const int n_sc = ctx.is_small() ? 128 : 256;
        const int st_sc = ctx.is_small() ? 30 : 80;
        for (int nt : {1, 2, 4, 8}) {
            Grid2D g = make_perturbed_2d(n_sc, n_sc, Scheme::AA_PATTERN, MemLayout::SOA);
            const auto b0 = clk::now();
            for (int t = 0; t < st_sc; ++t) step_periodic_2d(g, cp, nt);
            const double ms = std::chrono::duration<double, std::milli>(clk::now() - b0).count();
            if (nt == 1) base_ms = ms;
            const double speedup = base_ms / std::max(ms, 1e-6);
            const double mlups = (static_cast<double>(n_sc) * n_sc * st_sc) / (ms * 1000.0);
            rows_sc.push_back({std::to_string(nt), fmt_f8(ms), fmt_f8(mlups), fmt_f8(speedup)});
            nums_sc[0].push_back(nt);
            nums_sc[1].push_back(ms);
            nums_sc[2].push_back(mlups);
            nums_sc[3].push_back(speedup);
            bars.push_back({std::to_string(nt) + " threads", mlups, 0.0, "#059669"});
        }
        ctx.write_csv("S52", "out/s52_thread_scaling.csv", {"threads", "elapsed_ms", "mlups", "speedup"},
                      rows_sc, nums_sc, false);
        ctx.write_svg_bar("S52", "../../../outputs/figures/day009c_37_thread_scaling.svg",
                          "Multi-Thread Deterministic LBM Throughput (1 / 2 / 4 / 8 Threads)",
                          "AA-pattern SoA on " + std::to_string(n_sc) + "x" + std::to_string(n_sc),
                          "Throughput (MLUPS)", bars, false);
        ctx.record_stage("S52", "thread_scaling", since(t1), 2);
    }

    // S53-S54: Long-term multi-scheme conservation & Group F visual summary
    {
        auto t0 = clk::now();
        std::vector<std::vector<std::string>> rows;
        std::vector<std::vector<double>> nums(4);
        const int st = ctx.is_small() ? 500 : 2500;
        for (auto sc : {Scheme::PUSH_2GRID, Scheme::PULL_2GRID, Scheme::AA_PATTERN}) {
            Grid2D g = make_perturbed_2d(64, 64, sc, MemLayout::SOA);
            double m0, px0, py0, ek0, r0, r1, um;
            bool hn = false;
            g.totals(m0, px0, py0, ek0, r0, r1, um, hn);
            for (int t = 0; t < st; ++t) step_periodic_2d(g, cp, ctx.threads);
            double m1, px1, py1, ek1;
            g.totals(m1, px1, py1, ek1, r0, r1, um, hn);
            const double dm = std::abs(m1 - m0) / m0;
            const double dp = std::hypot(px1 - px0, py1 - py0) / m0;
            rows.push_back({scheme_name(sc), fmt_e6(dm), fmt_e6(dp), hex64(g.field_hash_fnv1a())});
            nums[1].push_back(dm);
            nums[2].push_back(dp);
        }
        ctx.write_csv("S53", "out/s53_scheme_long_conservation.csv",
                      {"scheme", "mass_drift", "mom_drift", "fnv1a"}, rows, nums);
        ctx.record_stage("S53", "scheme_long_conservation", since(t0), 1);

        auto t1 = clk::now();
        ctx.write_csv("S54", "out/s54_parity_matrix.csv",
                      {"check", "bit_exact"},
                      {{"2D_scheme_push_pull_aa", ok_2d_scheme ? "1" : "0"},
                       {"3D_scheme_push_pull_aa", ok_3d_scheme ? "1" : "0"},
                       {"2D_layout_aos_soa_tile", ok_layout ? "1" : "0"},
                       {"threads_1_2_4_8", ctx.threads_1_2_4_8_identical ? "1" : "0"}},
                      {{}, {ok_2d_scheme ? 1.0 : 0.0, ok_3d_scheme ? 1.0 : 0.0,
                            ok_layout ? 1.0 : 0.0, ctx.threads_1_2_4_8_identical ? 1.0 : 0.0}});
        ctx.write_svg_bar("S54", "../../../outputs/figures/day009c_38_parity_summary.svg",
                          "Bitwise Parity Verification Across Schemes, Layouts, and Threads",
                          "1.0 = zero-bit difference FNV-1a match",
                          "Bitwise exactness (0/1)",
                          {{"2D Schemes", ok_2d_scheme ? 1.0 : 0.0, 0.0, "#2563eb"},
                           {"3D Schemes", ok_3d_scheme ? 1.0 : 0.0, 0.0, "#059669"},
                           {"Layouts", ok_layout ? 1.0 : 0.0, 0.0, "#d97706"},
                           {"Threads", ctx.threads_1_2_4_8_identical ? 1.0 : 0.0, 0.0, "#7c3aed"}});
        ctx.record_stage("S54", "parity_matrix", since(t1), 2);
    }
}

}  // namespace lbm
