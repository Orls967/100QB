#include "lbm/stages.hpp"
#include "lbm/analysis.hpp"
#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/grid.hpp"
#include "lbm/parallel.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <vector>

namespace lbm {

namespace {

// Fast branchless interior pull kernel on SoA (V2 / V5) without per-cell modulo wrap on interior [1..nx-2]x[1..ny-2]
void step_fast_soa_pull_bgk(double* __restrict dst, const double* __restrict src,
                            int nx, int ny, double omega, int n_threads) {
    const std::size_t ncells = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    const auto stream_and_collide = [&](int ix, int iy, std::size_t c, bool interior) {
        double f[9];
        if (interior) {
            const std::size_t w = static_cast<std::size_t>(nx);
            f[0] = src[c];
            f[1] = src[ncells + c - 1];
            f[2] = src[2 * ncells + c - w];
            f[3] = src[3 * ncells + c + 1];
            f[4] = src[4 * ncells + c + w];
            f[5] = src[5 * ncells + c - w - 1];
            f[6] = src[6 * ncells + c - w + 1];
            f[7] = src[7 * ncells + c + w + 1];
            f[8] = src[8 * ncells + c + w - 1];
        } else {
            for (int q = 0; q < 9; ++q) {
                const int sx = (ix - D2Q9::cx[static_cast<std::size_t>(q)] + nx) % nx;
                const int sy = (iy - D2Q9::cy[static_cast<std::size_t>(q)] + ny) % ny;
                f[q] = src[static_cast<std::size_t>(q) * ncells +
                           static_cast<std::size_t>(sy) * static_cast<std::size_t>(nx) +
                           static_cast<std::size_t>(sx)];
            }
        }
        const double rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
        const double ux = ((f[1] + f[5] + f[8]) - (f[3] + f[6] + f[7])) / rho;
        const double uy = ((f[2] + f[5] + f[6]) - (f[4] + f[7] + f[8])) / rho;
        const double u2 = 1.5 * (ux * ux + uy * uy);
        const double one_m_om = 1.0 - omega;

        dst[c] = one_m_om * f[0] + omega * (4.0 / 9.0) * rho * (1.0 - u2);
        const double w1r = omega * (1.0 / 9.0) * rho;
        const double w2r = omega * (1.0 / 36.0) * rho;
        dst[ncells + c] = one_m_om * f[1] + w1r * (1.0 + 3.0 * ux + 4.5 * ux * ux - u2);
        dst[2 * ncells + c] = one_m_om * f[2] + w1r * (1.0 + 3.0 * uy + 4.5 * uy * uy - u2);
        dst[3 * ncells + c] = one_m_om * f[3] + w1r * (1.0 - 3.0 * ux + 4.5 * ux * ux - u2);
        dst[4 * ncells + c] = one_m_om * f[4] + w1r * (1.0 - 3.0 * uy + 4.5 * uy * uy - u2);
        const double c5 = ux + uy;
        const double c6 = -ux + uy;
        dst[5 * ncells + c] = one_m_om * f[5] + w2r * (1.0 + 3.0 * c5 + 4.5 * c5 * c5 - u2);
        dst[6 * ncells + c] = one_m_om * f[6] + w2r * (1.0 + 3.0 * c6 + 4.5 * c6 * c6 - u2);
        dst[7 * ncells + c] = one_m_om * f[7] + w2r * (1.0 - 3.0 * c5 + 4.5 * c5 * c5 - u2);
        dst[8 * ncells + c] = one_m_om * f[8] + w2r * (1.0 - 3.0 * c6 + 4.5 * c6 * c6 - u2);
    };

    parallel_for_range(static_cast<std::size_t>(ny), n_threads, [&](std::size_t y0, std::size_t y1, int) {
        for (std::size_t y = y0; y < y1; ++y) {
            const int iy = static_cast<int>(y);
            const std::size_t row = y * static_cast<std::size_t>(nx);
            if (iy > 0 && iy < ny - 1) {
                stream_and_collide(0, iy, row, false);
                for (int ix = 1; ix < nx - 1; ++ix) {
                    stream_and_collide(ix, iy, row + static_cast<std::size_t>(ix), true);
                }
                stream_and_collide(nx - 1, iy, row + static_cast<std::size_t>(nx - 1), false);
            } else {
                for (int ix = 0; ix < nx; ++ix) {
                    stream_and_collide(ix, iy, row + static_cast<std::size_t>(ix), false);
                }
            }
        }
    });
}

// Fast branchless in-place AA-pattern BGK SoA (V4 & V6)
void step_fast_aa_bgk(double* __restrict buf, int nx, int ny, double omega, std::size_t step_idx,
                      int n_threads) {
    const std::size_t ncells = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    const double one_m_om = 1.0 - omega;
    if ((step_idx & 1ULL) == 0ULL) {
        parallel_for_range(static_cast<std::size_t>(ny), n_threads, [&](std::size_t y0, std::size_t y1, int) {
            for (std::size_t y = y0; y < y1; ++y) {
                const std::size_t row = y * static_cast<std::size_t>(nx);
                for (int ix = 0; ix < nx; ++ix) {
                    const std::size_t c = row + static_cast<std::size_t>(ix);
                    const double f0 = buf[c];
                    const double f1 = buf[ncells + c];
                    const double f2 = buf[2 * ncells + c];
                    const double f3 = buf[3 * ncells + c];
                    const double f4 = buf[4 * ncells + c];
                    const double f5 = buf[5 * ncells + c];
                    const double f6 = buf[6 * ncells + c];
                    const double f7 = buf[7 * ncells + c];
                    const double f8 = buf[8 * ncells + c];
                    const double rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
                    const double ux = ((f1 + f5 + f8) - (f3 + f6 + f7)) / rho;
                    const double uy = ((f2 + f5 + f6) - (f4 + f7 + f8)) / rho;
                    const double u2 = 1.5 * (ux * ux + uy * uy);
                    const double w1r = omega * (1.0 / 9.0) * rho;
                    const double w2r = omega * (1.0 / 36.0) * rho;
                    const double c5 = ux + uy, c6 = -ux + uy;
                    buf[c] = one_m_om * f0 + omega * (4.0 / 9.0) * rho * (1.0 - u2);
                    buf[3 * ncells + c] = one_m_om * f1 + w1r * (1.0 + 3.0 * ux + 4.5 * ux * ux - u2);
                    buf[4 * ncells + c] = one_m_om * f2 + w1r * (1.0 + 3.0 * uy + 4.5 * uy * uy - u2);
                    buf[1 * ncells + c] = one_m_om * f3 + w1r * (1.0 - 3.0 * ux + 4.5 * ux * ux - u2);
                    buf[2 * ncells + c] = one_m_om * f4 + w1r * (1.0 - 3.0 * uy + 4.5 * uy * uy - u2);
                    buf[7 * ncells + c] = one_m_om * f5 + w2r * (1.0 + 3.0 * c5 + 4.5 * c5 * c5 - u2);
                    buf[8 * ncells + c] = one_m_om * f6 + w2r * (1.0 + 3.0 * c6 + 4.5 * c6 * c6 - u2);
                    buf[5 * ncells + c] = one_m_om * f7 + w2r * (1.0 - 3.0 * c5 + 4.5 * c5 * c5 - u2);
                    buf[6 * ncells + c] = one_m_om * f8 + w2r * (1.0 - 3.0 * c6 + 4.5 * c6 * c6 - u2);
                }
            }
        });
    } else {
        const std::size_t w = static_cast<std::size_t>(nx);
        parallel_for_range(static_cast<std::size_t>(ny), n_threads, [&](std::size_t y0, std::size_t y1, int) {
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                const std::size_t row = y * w;
                const bool mid_y = (iy > 0 && iy < ny - 1);
                for (int ix = 0; ix < nx; ++ix) {
                    const std::size_t c = row + static_cast<std::size_t>(ix);
                    std::size_t a[9];
                    if (mid_y && ix > 0 && ix < nx - 1) {
                        a[0] = c;
                        a[1] = ncells + c + 1;
                        a[2] = 2 * ncells + c + w;
                        a[3] = 3 * ncells + c - 1;
                        a[4] = 4 * ncells + c - w;
                        a[5] = 5 * ncells + c + w + 1;
                        a[6] = 6 * ncells + c + w - 1;
                        a[7] = 7 * ncells + c - w - 1;
                        a[8] = 8 * ncells + c - w + 1;
                    } else {
                        for (int m = 0; m < 9; ++m) {
                            const int mx = (ix + D2Q9::cx[static_cast<std::size_t>(m)] + nx) % nx;
                            const int my = (iy + D2Q9::cy[static_cast<std::size_t>(m)] + ny) % ny;
                            a[m] = static_cast<std::size_t>(m) * ncells +
                                   static_cast<std::size_t>(my) * w + static_cast<std::size_t>(mx);
                        }
                    }
                    // Load swapped slot opp[m] -> logical direction q
                    const double f0 = buf[a[0]];
                    const double f1 = buf[a[3]];
                    const double f2 = buf[a[4]];
                    const double f3 = buf[a[1]];
                    const double f4 = buf[a[2]];
                    const double f5 = buf[a[7]];
                    const double f6 = buf[a[8]];
                    const double f7 = buf[a[5]];
                    const double f8 = buf[a[6]];
                    const double rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
                    const double ux = ((f1 + f5 + f8) - (f3 + f6 + f7)) / rho;
                    const double uy = ((f2 + f5 + f6) - (f4 + f7 + f8)) / rho;
                    const double u2 = 1.5 * (ux * ux + uy * uy);
                    const double w1r = omega * (1.0 / 9.0) * rho;
                    const double w2r = omega * (1.0 / 36.0) * rho;
                    const double c5 = ux + uy, c6 = -ux + uy;
                    buf[a[0]] = one_m_om * f0 + omega * (4.0 / 9.0) * rho * (1.0 - u2);
                    buf[a[1]] = one_m_om * f1 + w1r * (1.0 + 3.0 * ux + 4.5 * ux * ux - u2);
                    buf[a[2]] = one_m_om * f2 + w1r * (1.0 + 3.0 * uy + 4.5 * uy * uy - u2);
                    buf[a[3]] = one_m_om * f3 + w1r * (1.0 - 3.0 * ux + 4.5 * ux * ux - u2);
                    buf[a[4]] = one_m_om * f4 + w1r * (1.0 - 3.0 * uy + 4.5 * uy * uy - u2);
                    buf[a[5]] = one_m_om * f5 + w2r * (1.0 + 3.0 * c5 + 4.5 * c5 * c5 - u2);
                    buf[a[6]] = one_m_om * f6 + w2r * (1.0 + 3.0 * c6 + 4.5 * c6 * c6 - u2);
                    buf[a[7]] = one_m_om * f7 + w2r * (1.0 - 3.0 * c5 + 4.5 * c5 * c5 - u2);
                    buf[a[8]] = one_m_om * f8 + w2r * (1.0 - 3.0 * c6 + 4.5 * c6 * c6 - u2);
                }
            }
        });
    }
}

}  // namespace

void run_group_g(RunContext& ctx) {
    using clk = std::chrono::steady_clock;
    auto since = [](clk::time_point t0) {
        return std::chrono::duration<double>(clk::now() - t0).count();
    };

    const int nx = ctx.is_small() ? 128 : 256;
    const int ny = nx;
    const int n_steps = ctx.is_small() ? 20 : 40;
    const int repeats = 7;
    const double omega = 1.0 / 0.8;
    const CollisionParams cp{CollisionOp::BGK, 0.8, 3.0 / 16.0, 0.0, 0.0, 0.0};

    // Host reference copy bandwidth (STREAM-like copy kernel)
    const double copy_bw_gbps = measure_copy_bandwidth_gbps(static_cast<std::size_t>(nx) * ny * 9ULL, 9);

    struct BenchVer {
        std::string name;
        double bytes_per_upd;
        std::vector<double> mlups_runs;
        double med_mlups = 0.0;
        double mad_mlups = 0.0;
        double bw_gbps = 0.0;
        double roof_frac = 0.0;
    };
    std::vector<BenchVer> vers = {
        {"V0_AoS_Push_1t", 144.0, {}},
        {"V1_SoA_Push_1t", 144.0, {}},
        {"V2_SoA_BranchlessPull_1t", 144.0, {}},
        {"V3_Tile16x16_Pull_1t", 144.0, {}},
        {"V4_AA_InPlace_1t", 72.0, {}},
        {"V5_SoA_BranchlessPull_mt", 144.0, {}},
        {"V6_AA_CombinedBest_mt", 72.0, {}}
    };

    // S55-S58: Execute optimization ladder V0..V6 (>=5 repeats each)
    std::vector<std::vector<std::string>> rows_rep, rows_sum;
    std::vector<std::vector<double>> nums_rep(3), nums_sum(6);

    auto run_ladder_version = [&](std::size_t vi) {
        auto& bv = vers[vi];
        bv.mlups_runs.reserve(repeats);
        for (int r = 0; r < repeats; ++r) {
            if (vi == 0 || vi == 1 || vi == 3) {
                const Scheme sc = Scheme::PUSH_2GRID;
                const MemLayout lay = (vi == 0) ? MemLayout::AOS : (vi == 1 ? MemLayout::SOA : MemLayout::TILE);
                Grid2D g(nx, ny, (vi == 3) ? Scheme::PULL_2GRID : sc, lay, 16, 16);
                g.init_uniform(1.0, 0.02, 0.01);
                step_periodic_2d(g, cp, 1);  // warmup
                const auto b0 = clk::now();
                for (int t = 0; t < n_steps; ++t) step_periodic_2d(g, cp, 1);
                const double sec = std::chrono::duration<double>(clk::now() - b0).count();
                const double ml = (static_cast<double>(nx) * ny * n_steps) / (sec * 1.0e6);
                bv.mlups_runs.push_back(ml);
            } else if (vi == 2 || vi == 5) {
                const int nt = (vi == 5) ? ctx.threads : 1;
                Grid2D g(nx, ny, Scheme::PULL_2GRID, MemLayout::SOA);
                g.init_uniform(1.0, 0.02, 0.01);
                step_fast_soa_pull_bgk(g.f1.data(), g.f0.data(), nx, ny, omega, nt);
                g.f0.swap(g.f1);
                const auto b0 = clk::now();
                for (int t = 0; t < n_steps; ++t) {
                    step_fast_soa_pull_bgk(g.f1.data(), g.f0.data(), nx, ny, omega, nt);
                    g.f0.swap(g.f1);
                }
                const double sec = std::chrono::duration<double>(clk::now() - b0).count();
                const double ml = (static_cast<double>(nx) * ny * n_steps) / (sec * 1.0e6);
                bv.mlups_runs.push_back(ml);
            } else {
                const int nt = (vi == 6) ? ctx.threads : 1;
                Grid2D g(nx, ny, Scheme::AA_PATTERN, MemLayout::SOA);
                g.init_uniform(1.0, 0.02, 0.01);
                step_fast_aa_bgk(g.f0.data(), nx, ny, omega, 0, nt);
                step_fast_aa_bgk(g.f0.data(), nx, ny, omega, 1, nt);
                const auto b0 = clk::now();
                for (int t = 0; t < n_steps; ++t) {
                    step_fast_aa_bgk(g.f0.data(), nx, ny, omega, static_cast<std::size_t>(t), nt);
                }
                const double sec = std::chrono::duration<double>(clk::now() - b0).count();
                const double ml = (static_cast<double>(nx) * ny * n_steps) / (sec * 1.0e6);
                bv.mlups_runs.push_back(ml);
            }
            rows_rep.push_back({bv.name, std::to_string(r), fmt_f8(bv.mlups_runs.back())});
            nums_rep[1].push_back(r);
            nums_rep[2].push_back(bv.mlups_runs.back());
        }
        const auto mm = median_mad(bv.mlups_runs);
        bv.med_mlups = mm.median;
        bv.mad_mlups = mm.mad;
        bv.bw_gbps = bv.med_mlups * bv.bytes_per_upd / 1000.0;
        bv.roof_frac = bv.bw_gbps / std::max(copy_bw_gbps, 1e-6);
        rows_sum.push_back({bv.name, fmt_f8(bv.med_mlups), fmt_f8(bv.mad_mlups),
                            fmt_f8(bv.bytes_per_upd), fmt_f8(bv.bw_gbps), fmt_f8(bv.roof_frac)});
        nums_sum[1].push_back(bv.med_mlups);
        nums_sum[2].push_back(bv.mad_mlups);
        nums_sum[3].push_back(bv.bytes_per_upd);
        nums_sum[4].push_back(bv.bw_gbps);
        nums_sum[5].push_back(bv.roof_frac);
    };

    {
        auto t0 = clk::now();
        run_ladder_version(0);
        run_ladder_version(1);
        ctx.write_csv("S55", "out/s55_ladder_v0_v1.csv", {"version", "repeat", "mlups"},
                      rows_rep, nums_rep, false);
        ctx.record_stage("S55", "ladder_v0_v1", since(t0), 1);
    }
    {
        auto t0 = clk::now();
        run_ladder_version(2);
        run_ladder_version(3);
        ctx.write_csv("S56", "out/s56_ladder_v2_v3.csv", {"version", "repeat", "mlups"},
                      rows_rep, nums_rep, false);
        ctx.record_stage("S56", "ladder_v2_v3", since(t0), 1);
    }
    {
        auto t0 = clk::now();
        run_ladder_version(4);
        run_ladder_version(5);
        run_ladder_version(6);
        ctx.write_csv("S57", "out/s57_ladder_summary.csv",
                      {"version", "median_mlups", "mad_mlups", "bytes_per_upd", "bw_gbps", "roofline_frac"},
                      rows_sum, nums_sum, false);
        std::vector<SvgBar> m_bars, r_bars;
        const std::vector<std::string> cols = {"#6b7280", "#3b82f6", "#2563eb", "#d97706", "#059669", "#7c3aed", "#047857"};
        for (std::size_t i = 0; i < vers.size(); ++i) {
            m_bars.push_back({vers[i].name.substr(0, 8), vers[i].med_mlups, vers[i].mad_mlups, cols[i]});
            r_bars.push_back({vers[i].name.substr(0, 8), vers[i].roof_frac, 0.0, cols[i]});
        }
        ctx.write_svg_bar("S57", "../../../outputs/figures/day009c_39_opt_ladder_mlups.svg",
                          "Optimization Ladder V0..V6 Throughput (Median +- MAD, 7 Repeats)",
                          std::to_string(nx) + "x" + std::to_string(ny) + " D2Q9 BGK (" + std::to_string(ctx.threads) + " threads for V5/V6)",
                          "Throughput (MLUPS)", m_bars, false);
        ctx.write_svg_bar("S58", "../../../outputs/figures/day009c_40_roofline_bw.svg",
                          "Effective Bandwidth Relative to Copy Bandwidth (" + fmt_f8(copy_bw_gbps) + " GB/s)",
                          "2-grid = 144 B/upd vs AA-pattern footprint = 72 B/upd (bandwidth-fraction view)",
                          "Fraction of copy bandwidth", r_bars, false);
        ctx.record_stage("S57", "ladder_v4_v6", since(t0), 2);
    }
    {
        auto t0 = clk::now();
        ctx.write_csv("S58", "out/s58_copy_bandwidth.csv",
                      {"metric", "value"},
                      {{"host_copy_bw_gbps", fmt_f8(copy_bw_gbps)},
                       {"v0_mlups", fmt_f8(vers[0].med_mlups)},
                       {"v6_mlups", fmt_f8(vers[6].med_mlups)},
                       {"v6_roof_frac", fmt_f8(vers[6].roof_frac)}},
                      {{}, {copy_bw_gbps, vers[0].med_mlups, vers[6].med_mlups, vers[6].roof_frac}}, false);
        std::string meas = "MLUPS=[" + fmt_f8(vers[0].med_mlups) + "," + fmt_f8(vers[1].med_mlups) + "," +
                           fmt_f8(vers[2].med_mlups) + "," + fmt_f8(vers[3].med_mlups) + "," +
                           fmt_f8(vers[4].med_mlups) + "," + fmt_f8(vers[5].med_mlups) + "," +
                           fmt_f8(vers[6].med_mlups) + "], copy_bw=" + fmt_f8(copy_bw_gbps) +
                           "GB/s, best_frac=" + fmt_f8(vers[6].roof_frac);
        ctx.set_criterion("c17", "tercapai", meas, "V0..V6 median+-MAD >=5 repeats + roofline bw fraction");
        ctx.record_stage("S58", "roofline_summary", since(t0), 2);
    }

    // S59-S60: Cell-count time scaling N >= 256^2 -> Criterion c18
    {
        auto t0 = clk::now();
        const std::vector<int> ns = ctx.is_small() ? std::vector<int>{128, 192, 256, 320}
                                                   : std::vector<int>{256, 384, 512, 640};
        const int steps_sc = ctx.is_small() ? 20 : 40;
        std::vector<std::vector<std::string>> rows_sc;
        std::vector<std::vector<double>> nums_sc(4);
        std::vector<double> log_c, log_t, xc, yt, yml;

        for (int n : ns) {
            Grid2D g(n, n, Scheme::AA_PATTERN, MemLayout::SOA);
            g.init_uniform(1.0, 0.02, 0.01);
            step_fast_aa_bgk(g.f0.data(), n, n, omega, 0, ctx.threads);
            step_fast_aa_bgk(g.f0.data(), n, n, omega, 1, ctx.threads);
            std::vector<double> dts;
            for (int r = 0; r < 5; ++r) {
                const auto b0 = clk::now();
                for (int t = 0; t < steps_sc; ++t) {
                    step_fast_aa_bgk(g.f0.data(), n, n, omega, static_cast<std::size_t>(t), ctx.threads);
                }
                dts.push_back(std::chrono::duration<double, std::milli>(clk::now() - b0).count());
            }
            const double med_ms = median_mad(dts).median;
            const double ncells = static_cast<double>(n) * n;
            const double mlups = (ncells * steps_sc) / (med_ms * 1000.0);
            rows_sc.push_back({std::to_string(n), std::to_string(static_cast<std::size_t>(ncells)),
                               fmt_f8(med_ms), fmt_f8(mlups)});
            nums_sc[0].push_back(n);
            nums_sc[1].push_back(ncells);
            nums_sc[2].push_back(med_ms);
            nums_sc[3].push_back(mlups);
            log_c.push_back(std::log(ncells));
            log_t.push_back(std::log(med_ms));
            xc.push_back(ncells);
            yt.push_back(med_ms);
            yml.push_back(mlups);
        }
        const auto fit = ols_fit(log_c, log_t);
        ctx.write_csv("S59", "out/s59_cell_time_scaling.csv",
                      {"n", "n_cells", "median_ms", "mlups"}, rows_sc, nums_sc, false);
        ctx.write_svg_xy("S59", "../../../outputs/figures/day009c_41_cell_scaling.svg",
                         "Grid Size Wall-Clock Time Scaling (N >= 256^2)",
                         "Log-log slope = " + fmt_f8(fit.slope) + " +- " + fmt_f8(fit.slope_se),
                         "Number of cells (N^2)", "Median step-block time (ms)",
                         {{"Measured time", "#2563eb", xc, yt, {}, true, true, false}}, true, true);
        ctx.record_stage("S59", "cell_time_scaling", since(t0), 2);

        auto t1 = clk::now();
        ctx.write_csv("S60", "out/s60_scaling_fit.csv",
                      {"slope", "slope_se", "r2"},
                      {{fmt_f8(fit.slope), fmt_f8(fit.slope_se), fmt_f8(fit.r2)}},
                      {{fit.slope}, {fit.slope_se}, {fit.r2}}, false);
        ctx.write_svg_xy("S60", "../../../outputs/figures/day009c_42_throughput_scale.svg",
                         "Throughput (MLUPS) Across Large Grids",
                         "AA-pattern multi-threaded (" + std::to_string(ctx.threads) + " threads)",
                         "Number of cells (N^2)", "Throughput (MLUPS)",
                         {{"MLUPS", "#059669", xc, yml, {}, true, true, false}}, true, false);
        const bool ok_c18 = (fit.slope >= 0.9 && fit.slope <= 1.15);
        ctx.set_criterion("c18", ok_c18 ? "tercapai" : "belum tercapai",
                          "slope=" + fmt_f8(fit.slope) + "+-" + fmt_f8(fit.slope_se),
                          "log-log slope in [0.9, 1.15] for N >= 256^2");
        ctx.record_stage("S60", "scaling_fit", since(t1), 2);
    }
}

}  // namespace lbm
