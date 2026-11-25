#include "wator.hpp"
#include <cstdio>
#include <cstring>
#include <string>

static std::pair<int, double> acf_second_peak(const std::vector<int> &x, int t0) {
    int n = static_cast<int>(x.size()) - t0;
    double mean = 0.0;
    for (int i = 0; i < n; ++i) mean += x[static_cast<size_t>(t0 + i)];
    mean /= static_cast<double>(n);
    double var = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[static_cast<size_t>(t0 + i)] - mean;
        var += d * d;
    }
    var /= static_cast<double>(n);
    bool passed_zero = false;
    int best_lag = -1;
    double best_acf = -1.0;
    for (int lag = 1; lag < n / 2; ++lag) {
        double c = 0.0;
        for (int i = 0; i < n - lag; ++i)
            c += (x[static_cast<size_t>(t0 + i)] - mean) * (x[static_cast<size_t>(t0 + i + lag)] - mean);
        double rho = (var > 1e-12) ? ((c / static_cast<double>(n - lag)) / var) : 0.0;
        if (rho < 0.0) passed_zero = true;
        if (passed_zero && rho > best_acf) { best_acf = rho; best_lag = lag; }
        if (passed_zero && lag > best_lag + 20 && best_acf > 0.0) break;
    }
    return {best_lag, best_acf};
}

static void write_svg_timeseries(const std::string &path, const std::vector<int> &F, const std::vector<int> &S) {
    FILE *fp = std::fopen(path.c_str(), "w");
    int max_v = 1;
    for (size_t i = 0; i < F.size(); ++i) max_v = std::max(max_v, std::max(F[i], S[i]));
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"420\" viewBox=\"0 0 800 420\">\n");
    std::fprintf(fp, "<rect width=\"800\" height=\"420\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"400\" y=\"24\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Wa-Tor 200x200: Deret Waktu Populasi (5000 Langkah)</text>\n");
    std::fprintf(fp, "<rect x=\"65\" y=\"40\" width=\"710\" height=\"330\" fill=\"none\" stroke=\"#cccccc\"/>\n");
    auto poly = [&](const std::vector<int> &v, const char *col) {
        std::fprintf(fp, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.4\" points=\"", col);
        for (size_t i = 0; i < v.size(); i += 5) {
            double px = 65.0 + 710.0 * static_cast<double>(i) / static_cast<double>(v.size() - 1);
            double py = 370.0 - 320.0 * static_cast<double>(v[i]) / static_cast<double>(max_v);
            std::fprintf(fp, "%.1f,%.1f ", px, py);
        }
        std::fprintf(fp, "\"/>\n");
    };
    poly(F, "#1f77b4");
    poly(S, "#d62728");
    std::fprintf(fp, "<text x=\"650\" y=\"62\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1f77b4\">Ikan (Mangsa)</text>\n");
    std::fprintf(fp, "<text x=\"650\" y=\"80\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#d62728\">Hiu (Pemangsa)</text>\n");
    std::fprintf(fp, "<text x=\"420\" y=\"405\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">Langkah (Kronon)</text>\n</svg>\n");
    std::fclose(fp);
}

static void write_svg_phase(const std::string &path, const std::vector<int> &F, const std::vector<int> &S,
                            const std::vector<double> &lvF, const std::vector<double> &lvS) {
    FILE *fp = std::fopen(path.c_str(), "w");
    double maxF = 1.0, maxS = 1.0;
    for (size_t i = 200; i < F.size(); ++i) { maxF = std::max(maxF, static_cast<double>(F[i])); maxS = std::max(maxS, static_cast<double>(S[i])); }
    for (size_t i = 0; i < lvF.size(); ++i) { maxF = std::max(maxF, lvF[i]); maxS = std::max(maxS, lvS[i]); }
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"680\" height=\"460\" viewBox=\"0 0 680 460\">\n");
    std::fprintf(fp, "<rect width=\"680\" height=\"460\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"340\" y=\"24\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Potret Fase Wa-Tor vs Orbit Lotka-Volterra (RK4 Fit)</text>\n");
    std::fprintf(fp, "<rect x=\"65\" y=\"40\" width=\"585\" height=\"370\" fill=\"none\" stroke=\"#cccccc\"/>\n");
    std::fprintf(fp, "<polyline fill=\"none\" stroke=\"#1f77b4\" stroke-opacity=\"0.7\" stroke-width=\"1.2\" points=\"");
    for (size_t i = 200; i < F.size(); i += 3) {
        double px = 65.0 + 575.0 * static_cast<double>(F[i]) / maxF;
        double py = 410.0 - 360.0 * static_cast<double>(S[i]) / maxS;
        std::fprintf(fp, "%.1f,%.1f ", px, py);
    }
    std::fprintf(fp, "\"/>\n<polyline fill=\"none\" stroke=\"#d62728\" stroke-width=\"2.0\" stroke-dasharray=\"5,3\" points=\"");
    for (size_t i = 0; i < lvF.size(); i += 2) {
        double px = 65.0 + 575.0 * lvF[i] / maxF;
        double py = 410.0 - 360.0 * lvS[i] / maxS;
        std::fprintf(fp, "%.1f,%.1f ", px, py);
    }
    std::fprintf(fp, "\"/>\n<text x=\"470\" y=\"62\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1f77b4\">Simulasi Wa-Tor (t>=200)</text>\n");
    std::fprintf(fp, "<text x=\"470\" y=\"80\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#d62728\">Orbit LV RK4 (Fit OLS)</text>\n");
    std::fprintf(fp, "<text x=\"355\" y=\"442\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">Populasi Ikan</text>\n</svg>\n");
    std::fclose(fp);
}

static void write_svg_snapshots(const std::string &path, const std::vector<std::vector<uint8_t>> &snaps,
                                const int steps[4], int W, int H) {
    FILE *fp = std::fopen(path.c_str(), "w");
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"840\" height=\"240\" viewBox=\"0 0 840 240\">\n");
    std::fprintf(fp, "<rect width=\"840\" height=\"240\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"420\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\" font-weight=\"bold\">Snapshot Grid Wa-Tor 200x200 (Biru=Ikan, Merah=Hiu)</text>\n");
    for (int k = 0; k < 4; ++k) {
        int ox = 20 + k * 205, oy = 42;
        std::fprintf(fp, "<text x=\"%d\" y=\"36\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\">Langkah %d</text>\n", ox + 90, steps[k]);
        std::fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"180\" height=\"180\" fill=\"#f7f7f7\" stroke=\"#999999\"/>\n", ox, oy);
        for (int y = 0; y < H; y += 2) {
            for (int x = 0; x < W; x += 2) {
                uint8_t sp = snaps[static_cast<size_t>(k)][static_cast<size_t>(y * W + x)];
                if (sp == 0) continue;
                const char *col = (sp == 1) ? "#1f77b4" : "#d62728";
                std::fprintf(fp, "<rect x=\"%.1f\" y=\"%.1f\" width=\"1.8\" height=\"1.8\" fill=\"%s\"/>\n",
                             ox + x * 0.9, oy + y * 0.9, col);
            }
        }
    }
    std::fprintf(fp, "</svg>\n");
    std::fclose(fp);
}

static void write_svg_extinction(const std::string &path, const int sizes[4], const double p_ext[4],
                                 const double lo[4], const double hi[4]) {
    FILE *fp = std::fopen(path.c_str(), "w");
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"400\" viewBox=\"0 0 640 400\">\n");
    std::fprintf(fp, "<rect width=\"640\" height=\"400\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"320\" y=\"24\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Proporsi Kepunahan vs Ukuran Grid (40 Replika, Wilson 95%%)</text>\n");
    std::fprintf(fp, "<rect x=\"65\" y=\"40\" width=\"535\" height=\"300\" fill=\"none\" stroke=\"#cccccc\"/>\n");
    std::fprintf(fp, "<polyline fill=\"none\" stroke=\"#2ca02c\" stroke-width=\"2\" points=\"");
    double xs[4];
    for (int i = 0; i < 4; ++i) {
        xs[i] = 110.0 + i * 150.0;
        double py = 340.0 - 280.0 * p_ext[i];
        std::fprintf(fp, "%.1f,%.1f ", xs[i], py);
    }
    std::fprintf(fp, "\"/>\n");
    for (int i = 0; i < 4; ++i) {
        double py = 340.0 - 280.0 * p_ext[i], plo = 340.0 - 280.0 * lo[i], phi = 340.0 - 280.0 * hi[i];
        std::fprintf(fp, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#2ca02c\" stroke-width=\"1.8\"/>\n", xs[i], plo, xs[i], phi);
        std::fprintf(fp, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4.5\" fill=\"#2ca02c\"/>\n", xs[i], py);
        std::fprintf(fp, "<text x=\"%.1f\" y=\"360\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">W=%d</text>\n", xs[i], sizes[i]);
        std::fprintf(fp, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\">%.2f</text>\n", xs[i] + 18, py - 6, p_ext[i]);
    }
    std::fprintf(fp, "</svg>\n");
    std::fclose(fp);
}

int main(int argc, char **argv) {
    std::string out_dir = "data", fig_dir = "../../../outputs/figures";
    int n_threads = 8;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) out_dir = argv[++i];
        else if (std::strcmp(argv[i], "--fig-dir") == 0 && i + 1 < argc) fig_dir = argv[++i];
        else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) n_threads = std::atoi(argv[++i]);
    }

    wator::Params def_p;
    wator::World w(def_p, wator::MASTER_SEED);
    std::vector<int> F(5001), S(5001);
    std::vector<std::vector<uint8_t>> snaps;
    const int snap_steps[4] = {0, 200, 1000, 4000};
    bool inv_ok = w.check_invariants();
    for (int t = 0; t <= 5000; ++t) {
        F[static_cast<size_t>(t)] = w.n_fish;
        S[static_cast<size_t>(t)] = w.n_sharks;
        for (int st : snap_steps) {
            if (t == st) {
                std::vector<uint8_t> buf(w.grid.size());
                for (size_t j = 0; j < w.grid.size(); ++j) buf[j] = static_cast<uint8_t>(w.grid[j].kind);
                snaps.push_back(std::move(buf));
            }
        }
        if (t < 5000) { w.step(); inv_ok = inv_ok && w.check_invariants(); }
    }
    if (!inv_ok) { std::fprintf(stderr, "Invariant failed!\n"); return 1; }

    FILE *fts = std::fopen((out_dir + "/timeseries.csv").c_str(), "w");
    std::fprintf(fts, "step,fish,sharks\n");
    for (int t = 0; t <= 5000; ++t) std::fprintf(fts, "%d,%d,%d\n", t, F[static_cast<size_t>(t)], S[static_cast<size_t>(t)]);
    std::fclose(fts);

    auto lv = wator::fit_lv_ols(F, S, 200);
    auto [sim_period, acf_peak] = acf_second_peak(F, 200);
    double lv_period = (lv.a > 0 && lv.c > 0) ? (2.0 * 3.14159265358979323846 / std::sqrt(lv.a * lv.c)) : 0.0;
    double lf = static_cast<double>(F[200]), ls = static_cast<double>(S[200]);
    std::vector<double> lvF, lvS;
    for (int k = 0; k < 2000; ++k) {
        lvF.push_back(lf); lvS.push_back(ls);
        wator::rk4_lv_step(lf, ls, 0.1, lv);
    }

    FILE *flv = std::fopen((out_dir + "/lv_fit.csv").c_str(), "w");
    std::fprintf(flv, "a,b,c,d,sim_period,acf_peak,lv_period,invariants_ok\n");
    std::fprintf(flv, "%.10g,%.10g,%.10g,%.10g,%d,%.6f,%.6f,%d\n", lv.a, lv.b, lv.c, lv.d, sim_period, acf_peak, lv_period, inv_ok ? 1 : 0);
    std::fclose(flv);

    const int sizes[4] = {20, 40, 80, 160};
    double p_ext[4]{}, lo[4]{}, hi[4]{};
    FILE *fex = std::fopen((out_dir + "/extinction_size.csv").c_str(), "w");
    std::fprintf(fex, "W,n_reps,n_extinct,prop,wilson_lo,wilson_hi\n");
    for (int si = 0; si < 4; ++si) {
        int W = sizes[si], n_reps = 40;
        std::vector<int> ext_flag(static_cast<size_t>(n_reps), 0), rep_inv(static_cast<size_t>(n_reps), 1);
        std::atomic<int> next_r{0};
        std::vector<std::thread> pool;
        for (int th = 0; th < n_threads; ++th) {
            pool.emplace_back([&]() {
                while (true) {
                    int r = next_r.fetch_add(1);
                    if (r >= n_reps) break;
                    wator::Params rp = def_p; rp.W = W; rp.H = W;
                    uint64_t seed = wator::splitmix64_mix(wator::MASTER_SEED ^ static_cast<uint64_t>(W * 1000 + r));
                    wator::World rw(rp, seed);
                    for (int t = 0; t < 5000; ++t) {
                        rw.step();
                        if (!rw.check_invariants()) rep_inv[static_cast<size_t>(r)] = 0;
                        if (rw.n_fish == 0 || rw.n_sharks == 0) { ext_flag[static_cast<size_t>(r)] = 1; break; }
                    }
                }
            });
        }
        for (auto &th : pool) th.join();
        int cnt = 0;
        for (int r = 0; r < n_reps; ++r) { cnt += ext_flag[static_cast<size_t>(r)]; inv_ok = inv_ok && rep_inv[static_cast<size_t>(r)]; }
        double p = static_cast<double>(cnt) / n_reps, z = 1.95996398454;
        double den = 1.0 + z * z / n_reps, cen = (p + z * z / (2.0 * n_reps)) / den;
        double half = z * std::sqrt((p * (1.0 - p) + z * z / (4.0 * n_reps)) / n_reps) / den;
        p_ext[si] = p; lo[si] = std::max(0.0, cen - half); hi[si] = std::min(1.0, cen + half);
        std::fprintf(fex, "%d,%d,%d,%.6f,%.6f,%.6f\n", W, n_reps, cnt, p, lo[si], hi[si]);
    }
    std::fclose(fex);

    write_svg_timeseries(fig_dir + "/day008b_timeseries.svg", F, S);
    write_svg_phase(fig_dir + "/day008b_phase_portrait.svg", F, S, lvF, lvS);
    write_svg_snapshots(fig_dir + "/day008b_snapshots.svg", snaps, snap_steps, def_p.W, def_p.H);
    write_svg_extinction(fig_dir + "/day008b_extinction_size.svg", sizes, p_ext, lo, hi);
    std::printf("Wa-Tor done: sim_period=%d (acf=%.3f), lv_period=%.2f, inv_ok=%d\n", sim_period, acf_peak, lv_period, inv_ok ? 1 : 0);
    return 0;
}
