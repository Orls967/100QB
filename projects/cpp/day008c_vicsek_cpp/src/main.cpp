#include "vicsek.hpp"
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>

static double now_sec() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void median_mad(std::vector<double> v, double &med, double &mad) {
    std::sort(v.begin(), v.end());
    med = v[v.size() / 2];
    for (double &x : v) x = std::fabs(x - med);
    std::sort(v.begin(), v.end());
    mad = v[v.size() / 2];
}

static void angle_to_rgb(double th, char *buf) {
    double h = std::fmod(th + vicsek::PI, vicsek::TWO_PI) / vicsek::TWO_PI * 6.0;
    int i = static_cast<int>(h) % 6;
    double f = h - std::floor(h);
    int q = static_cast<int>(255 * (1.0 - f)), t = static_cast<int>(255 * f);
    int r = 0, g = 0, b = 0;
    if (i == 0) { r = 255; g = t; } else if (i == 1) { r = q; g = 255; }
    else if (i == 2) { g = 255; b = t; } else if (i == 3) { g = q; b = 255; }
    else if (i == 4) { r = t; b = 255; } else { r = 255; b = q; }
    std::snprintf(buf, 16, "#%02x%02x%02x", r, g, b);
}

static void write_svg_phi_eta(const std::string &path, const std::vector<double> &etas,
                              const double phi[3][25], const double err[3][25]) {
    FILE *fp = std::fopen(path.c_str(), "w");
    const int Ls[3] = {16, 32, 64};
    const char *cols[3] = {"#1f77b4", "#ff7f0e", "#2ca02c"};
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"740\" height=\"420\" viewBox=\"0 0 740 420\">\n");
    std::fprintf(fp, "<rect width=\"740\" height=\"420\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"370\" y=\"24\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Model Vicsek (rho=4): Polarisasi phi vs Derau eta</text>\n");
    std::fprintf(fp, "<rect x=\"65\" y=\"40\" width=\"640\" height=\"325\" fill=\"none\" stroke=\"#cccccc\"/>\n");
    for (int li = 0; li < 3; ++li) {
        std::fprintf(fp, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.8\" points=\"", cols[li]);
        for (int ei = 0; ei < 25; ++ei) {
            double px = 65.0 + 640.0 * etas[ei] / vicsek::TWO_PI;
            double py = 365.0 - 315.0 * phi[li][ei];
            std::fprintf(fp, "%.1f,%.1f ", px, py);
        }
        std::fprintf(fp, "\"/>\n");
        for (int ei = 0; ei < 25; ++ei) {
            double px = 65.0 + 640.0 * etas[ei] / vicsek::TWO_PI;
            double py = 365.0 - 315.0 * phi[li][ei];
            double dy = 315.0 * err[li][ei];
            std::fprintf(fp, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\"/>\n", px, py - dy, px, py + dy, cols[li]);
            std::fprintf(fp, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"3\" fill=\"%s\"/>\n", px, py, cols[li]);
        }
        std::fprintf(fp, "<text x=\"600\" y=\"%d\" font-family=\"sans-serif\" font-size=\"12\" fill=\"%s\">L=%d</text>\n", 65 + li * 18, cols[li], Ls[li]);
    }
    std::fprintf(fp, "<text x=\"385\" y=\"398\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">Derau eta [0, 2*pi]</text>\n</svg>\n");
    std::fclose(fp);
}

static void write_svg_scaling(const std::string &path, const int Ns[8], const double bf_ns[8], const double cl_ns[8]) {
    FILE *fp = std::fopen(path.c_str(), "w");
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"420\" viewBox=\"0 0 720 420\">\n");
    std::fprintf(fp, "<rect width=\"720\" height=\"420\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"360\" y=\"24\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"bold\">Skala Waktu Log-Log: Cell List O(N) vs Brute Force O(N^2)</text>\n");
    std::fprintf(fp, "<rect x=\"75\" y=\"40\" width=\"610\" height=\"325\" fill=\"none\" stroke=\"#cccccc\"/>\n");
    double lx0 = std::log10(500.0), lx1 = std::log10(256000.0), ly0 = 0.5, ly1 = 5.0;
    auto mapx = [&](int N) { return 75.0 + 610.0 * (std::log10(static_cast<double>(N)) - lx0) / (lx1 - lx0); };
    auto mapy = [&](double ns) { return 365.0 - 325.0 * (std::log10(ns) - ly0) / (ly1 - ly0); };
    std::fprintf(fp, "<polyline fill=\"none\" stroke=\"#d62728\" stroke-width=\"2\" points=\"");
    for (int i = 0; i < 6; ++i) std::fprintf(fp, "%.1f,%.1f ", mapx(Ns[i]), mapy(bf_ns[i]));
    std::fprintf(fp, "\"/>\n<polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"2\" points=\"");
    for (int i = 0; i < 8; ++i) std::fprintf(fp, "%.1f,%.1f ", mapx(Ns[i]), mapy(cl_ns[i]));
    std::fprintf(fp, "\"/>\n");
    for (int i = 0; i < 8; ++i) {
        std::fprintf(fp, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"3.5\" fill=\"#1f77b4\"/>\n", mapx(Ns[i]), mapy(cl_ns[i]));
        if (i < 6) std::fprintf(fp, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"3.5\" fill=\"#d62728\"/>\n", mapx(Ns[i]), mapy(bf_ns[i]));
    }
    std::fprintf(fp, "<text x=\"110\" y=\"65\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#d62728\">Brute Force (~N^2)</text>\n");
    std::fprintf(fp, "<text x=\"110\" y=\"85\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1f77b4\">Cell List (~N)</text>\n");
    std::fprintf(fp, "<text x=\"380\" y=\"398\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">log10(N) [500 .. 256000]</text>\n</svg>\n");
    std::fclose(fp);
}

static void write_svg_snapshots(const std::string &path, const std::vector<vicsek::Sim> &snaps, const double etas[3]) {
    FILE *fp = std::fopen(path.c_str(), "w");
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"840\" height=\"295\" viewBox=\"0 0 840 295\">\n");
    std::fprintf(fp, "<rect width=\"840\" height=\"295\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"420\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\" font-weight=\"bold\">Snapshot Kawanan Vicsek L=64 (N=16384, Diwarnai Menurut Arah)</text>\n");
    char col[16];
    for (int k = 0; k < 3; ++k) {
        int ox = 20 + k * 270, oy = 45;
        std::fprintf(fp, "<text x=\"%d\" y=\"38\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">eta = %.2f (phi=%.3f)</text>\n", ox + 120, etas[k], snaps[k].polarization());
        std::fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"240\" height=\"240\" fill=\"#fafafa\" stroke=\"#999999\"/>\n", ox, oy);
        for (int i = 0; i < snaps[k].N; i += 8) {
            double px = ox + 240.0 * snaps[k].x[i] / snaps[k].L;
            double py = oy + 240.0 * snaps[k].y[i] / snaps[k].L;
            double dx = 4.5 * std::cos(snaps[k].theta[i]), dy = 4.5 * std::sin(snaps[k].theta[i]);
            angle_to_rgb(snaps[k].theta[i], col);
            std::fprintf(fp, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"1.1\"/>\n", px, py, px + dx, py + dy, col);
        }
    }
    std::fprintf(fp, "</svg>\n");
    std::fclose(fp);
}

static void write_svg_timeseries(const std::string &path, const std::vector<double> &phi_t, double eta) {
    FILE *fp = std::fopen(path.c_str(), "w");
    std::fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"740\" height=\"320\" viewBox=\"0 0 740 320\">\n");
    std::fprintf(fp, "<rect width=\"740\" height=\"320\" fill=\"#ffffff\"/>\n");
    std::fprintf(fp, "<text x=\"370\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\" font-weight=\"bold\">Deret Waktu Polarisasi phi(t) Dekat Transisi (L=32, eta=%.2f)</text>\n", eta);
    std::fprintf(fp, "<rect x=\"60\" y=\"35\" width=\"650\" height=\"240\" fill=\"none\" stroke=\"#cccccc\"/>\n");
    std::fprintf(fp, "<polyline fill=\"none\" stroke=\"#9467bd\" stroke-width=\"1.2\" points=\"");
    for (size_t i = 0; i < phi_t.size(); i += 4) {
        double px = 60.0 + 650.0 * static_cast<double>(i) / static_cast<double>(phi_t.size() - 1);
        double py = 275.0 - 235.0 * phi_t[i];
        std::fprintf(fp, "%.1f,%.1f ", px, py);
    }
    std::fprintf(fp, "\"/>\n<text x=\"385\" y=\"302\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">Langkah</text>\n</svg>\n");
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

    /* 1. Parity N=1000, 500 steps */
    vicsek::Sim p_bf(1000, std::sqrt(250.0), 2.0, vicsek::MASTER_SEED), p_cl(1000, std::sqrt(250.0), 2.0, vicsek::MASTER_SEED);
    for (int t = 0; t < 500; ++t) { p_bf.step_bruteforce(); p_cl.step_cell_list(true); }
    int parity_diffs = 0;
    for (int i = 0; i < 1000; ++i)
        if (p_bf.x[i] != p_cl.x[i] || p_bf.y[i] != p_cl.y[i] || p_bf.theta[i] != p_cl.theta[i]) ++parity_diffs;
    std::printf("Parity N=1000, 500 steps diffs = %d\n", parity_diffs);

    /* 2. Order parameter phi(eta) for L in {16, 32, 64}, rho=4 */
    const int Ls[3] = {16, 32, 64};
    std::vector<double> etas(25);
    for (int ei = 0; ei < 25; ++ei) etas[ei] = vicsek::TWO_PI * static_cast<double>(ei) / 24.0;
    double phi_mean[3][25]{}, phi_err[3][25]{}, u4_mean[3][25]{};

    struct Task { int li, ei, rep; double p1, p2, p4; };
    std::vector<Task> tasks;
    for (int li = 0; li < 3; ++li)
        for (int ei = 0; ei < 25; ++ei)
            for (int r = 0; r < 4; ++r) tasks.push_back({li, ei, r, 0, 0, 0});

    std::atomic<size_t> next_t{0};
    std::vector<std::thread> pool;
    for (int th = 0; th < n_threads; ++th) {
        pool.emplace_back([&]() {
            while (true) {
                size_t idx = next_t.fetch_add(1);
                if (idx >= tasks.size()) break;
                Task &tk = tasks[idx];
                int L = Ls[tk.li];
                uint64_t sd = vicsek::splitmix64_mix(vicsek::MASTER_SEED ^ static_cast<uint64_t>(L * 1000 + tk.ei * 10 + tk.rep));
                vicsek::Sim sim(4 * L * L, static_cast<double>(L), etas[tk.ei], sd);
                for (int t = 0; t < 2000; ++t) sim.step_cell_list(false);
                double s1 = 0, s2 = 0, s4 = 0;
                for (int t = 0; t < 4000; ++t) {
                    sim.step_cell_list(false);
                    double ph = sim.polarization();
                    s1 += ph; s2 += ph * ph; s4 += ph * ph * ph * ph;
                }
                tk.p1 = s1 / 4000.0; tk.p2 = s2 / 4000.0; tk.p4 = s4 / 4000.0;
            }
        });
    }
    for (auto &th : pool) th.join();

    FILE *fphi = std::fopen((out_dir + "/phi_eta.csv").c_str(), "w");
    std::fprintf(fphi, "L,N,eta,phi_mean,phi_err,binder_u4\n");
    for (int li = 0; li < 3; ++li) {
        for (int ei = 0; ei < 25; ++ei) {
            double vals[4], m2 = 0, m4 = 0;
            for (int r = 0; r < 4; ++r) {
                const auto &tk = tasks[static_cast<size_t>((li * 25 + ei) * 4 + r)];
                vals[r] = tk.p1; m2 += tk.p2 * 0.25; m4 += tk.p4 * 0.25;
            }
            double m = 0.25 * (vals[0] + vals[1] + vals[2] + vals[3]), v = 0;
            for (int r = 0; r < 4; ++r) v += (vals[r] - m) * (vals[r] - m);
            phi_mean[li][ei] = m;
            phi_err[li][ei] = std::sqrt(v / (4.0 * 3.0));
            u4_mean[li][ei] = 1.0 - m4 / (3.0 * m2 * m2);
            std::fprintf(fphi, "%d,%d,%.8f,%.8f,%.8f,%.8f\n", Ls[li], 4 * Ls[li] * Ls[li], etas[ei], m, phi_err[li][ei], u4_mean[li][ei]);
        }
    }
    std::fclose(fphi);

    /* 3. Scaling N in {500..256000} */
    const int Ns[8] = {500, 1000, 2000, 4000, 8000, 16000, 64000, 256000};
    double bf_med[8]{}, bf_mad[8]{}, cl_med[8]{}, cl_mad[8]{};
    FILE *fsc = std::fopen((out_dir + "/scaling.csv").c_str(), "w");
    std::fprintf(fsc, "N,bf_ns_med,bf_ns_mad,cl_ns_med,cl_ns_mad\n");
    for (int ni = 0; ni < 8; ++ni) {
        int N = Ns[ni];
        double L = std::sqrt(static_cast<double>(N) / 4.0);
        int steps = (N <= 16000) ? 12 : 4;
        std::vector<double> vb, vc;
        for (int r = 0; r < 7; ++r) {
            vicsek::Sim scl(N, L, 2.0, vicsek::MASTER_SEED + static_cast<uint64_t>(r));
            double t0 = now_sec();
            for (int s = 0; s < steps; ++s) scl.step_cell_list(false);
            vc.push_back((now_sec() - t0) * 1e9 / (static_cast<double>(steps) * N));
            if (ni < 6) {
                vicsek::Sim sbf(N, L, 2.0, vicsek::MASTER_SEED + static_cast<uint64_t>(r));
                int bsteps = (N <= 4000) ? 10 : 3;
                t0 = now_sec();
                for (int s = 0; s < bsteps; ++s) sbf.step_bruteforce();
                vb.push_back((now_sec() - t0) * 1e9 / (static_cast<double>(bsteps) * N));
            }
        }
        median_mad(vc, cl_med[ni], cl_mad[ni]);
        if (ni < 6) median_mad(vb, bf_med[ni], bf_mad[ni]);
        std::fprintf(fsc, "%d,%.6f,%.6f,%.6f,%.6f\n", N, bf_med[ni], bf_mad[ni], cl_med[ni], cl_mad[ni]);
    }
    std::fclose(fsc);

    /* 4. Snapshots L=64 and transition time series L=32 */
    const double snap_etas[3] = {0.3, 2.5, 5.5};
    std::vector<vicsek::Sim> snaps;
    for (int k = 0; k < 3; ++k) {
        snaps.emplace_back(16384, 64.0, snap_etas[k], vicsek::MASTER_SEED + static_cast<uint64_t>(k));
        for (int t = 0; t < 2000; ++t) snaps.back().step_cell_list(false);
    }
    vicsek::Sim ts_sim(4096, 32.0, 2.5, vicsek::MASTER_SEED + 77);
    std::vector<double> phi_t(4000);
    for (int t = 0; t < 4000; ++t) { ts_sim.step_cell_list(false); phi_t[t] = ts_sim.polarization(); }

    FILE *fp = std::fopen((out_dir + "/summary.csv").c_str(), "w");
    std::fprintf(fp, "parity_diffs\n%d\n", parity_diffs);
    std::fclose(fp);

    write_svg_phi_eta(fig_dir + "/day008c_phi_vs_eta.svg", etas, phi_mean, phi_err);
    write_svg_scaling(fig_dir + "/day008c_scaling.svg", Ns, bf_med, cl_med);
    write_svg_snapshots(fig_dir + "/day008c_snapshots.svg", snaps, snap_etas);
    write_svg_timeseries(fig_dir + "/day008c_phi_timeseries.svg", phi_t, 2.5);
    std::printf("Vicsek simulation and SVGs complete.\n");
    return 0;
}
