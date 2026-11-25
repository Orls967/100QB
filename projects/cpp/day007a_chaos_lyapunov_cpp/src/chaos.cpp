#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace chaos {

struct Vec3 {
    double x, y, z;
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

inline double dot3(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm3(const Vec3& a) { return std::sqrt(dot3(a, a)); }

template <typename Deriv>
auto rk4_step(const auto& state, double t, double dt, Deriv&& f) {
    auto k1 = f(state, t);
    auto k2 = f(state + k1 * (0.5 * dt), t + 0.5 * dt);
    auto k3 = f(state + k2 * (0.5 * dt), t + 0.5 * dt);
    auto k4 = f(state + k3 * dt, t + dt);
    return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

inline double logistic_lyapunov(double r, int n_trans, long long n_iter, double x0 = 0.3141592653589793) {
    double x = x0;
    for (int i = 0; i < n_trans; ++i) x = r * x * (1.0 - x);
    double acc = 0.0;
    for (long long i = 0; i < n_iter; ++i) {
        double df = std::abs(r * (1.0 - 2.0 * x));
        acc += std::log(df > 1e-300 ? df : 1e-300);
        x = r * x * (1.0 - x);
    }
    return acc / static_cast<double>(n_iter);
}

struct FeigenbaumRow { int n; long long period; double R_n; double delta_n; };

inline std::vector<FeigenbaumRow> solve_feigenbaum(int n_max = 8) {
    std::vector<FeigenbaumRow> rows;
    for (int n = 1; n <= n_max; ++n) {
        long long p = 1LL << n;
        double r = (n == 1) ? 3.2 : (n == 2) ? 3.4985 :
            rows[n - 2].R_n + (rows[n - 2].R_n - rows[n - 3].R_n) /
            ((n == 3) ? 4.669 : rows[n - 2].delta_n);
        for (int it = 0; it < 25; ++it) {
            double x = 0.5, dx = 0.0;
            for (long long k = 0; k < p; ++k) {
                dx = x * (1.0 - x) + r * (1.0 - 2.0 * x) * dx;
                x = r * x * (1.0 - x);
            }
            double dr = (x - 0.5) / dx;
            r -= dr;
            if (std::abs(dr) < 1e-15) break;
        }
        double d_n = (n >= 3) ? (rows[n - 2].R_n - rows[n - 3].R_n) / (r - rows[n - 2].R_n) : 0.0;
        rows.push_back({n, p, r, d_n});
    }
    return rows;
}

struct LorenzState12 {
    std::array<double, 12> v;
    LorenzState12 operator+(const LorenzState12& o) const {
        LorenzState12 r{};
        for (int i = 0; i < 12; ++i) r.v[i] = v[i] + o.v[i];
        return r;
    }
    LorenzState12 operator*(double s) const {
        LorenzState12 r{};
        for (int i = 0; i < 12; ++i) r.v[i] = v[i] * s;
        return r;
    }
};

inline Vec3 lorenz_rhs(const Vec3& u, double sig = 10.0, double rho = 28.0, double beta = 8.0 / 3.0) {
    return {sig * (u.y - u.x), u.x * (rho - u.z) - u.y, u.x * u.y - beta * u.z};
}

inline LorenzState12 lorenz_var_rhs(const LorenzState12& s, double sig = 10.0, double rho = 28.0, double beta = 8.0 / 3.0) {
    LorenzState12 d{};
    double x = s.v[0], y = s.v[1], z = s.v[2];
    d.v[0] = sig * (y - x); d.v[1] = x * (rho - z) - y; d.v[2] = x * y - beta * z;
    for (int c = 0; c < 3; ++c) {
        double wx = s.v[3 + 3 * c], wy = s.v[4 + 3 * c], wz = s.v[5 + 3 * c];
        d.v[3 + 3 * c] = sig * (wy - wx);
        d.v[4 + 3 * c] = (rho - z) * wx - wy - x * wz;
        d.v[5 + 3 * c] = y * wx + x * wy - beta * wz;
    }
    return d;
}

struct LyapunovResult {
    std::array<double, 3> lam;
    double sum_lam, d_ky;
    std::vector<std::array<double, 4>> history;
    Vec3 post_transient;
};

inline LyapunovResult compute_lorenz_lyapunov(double dt, double t_trans, double t_run, int gs_every = 10) {
    Vec3 u{1.0, 1.0, 1.0};
    int n_tr = static_cast<int>( std::llround(t_trans / dt) );
    auto f3 = [](const Vec3& st, double) { return lorenz_rhs(st); };
    for (int i = 0; i < n_tr; ++i) u = rk4_step(u, 0.0, dt, f3);
    LorenzState12 s{};
    s.v = {u.x, u.y, u.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::array<double, 3> lsum{0.0, 0.0, 0.0};
    long long n_steps = std::llround(t_run / dt);
    auto f12 = [](const LorenzState12& st, double) { return lorenz_var_rhs(st); };
    std::vector<std::array<double, 4>> hist;
    for (long long step = 1; step <= n_steps; ++step) {
        s = rk4_step(s, 0.0, dt, f12);
        if (step % gs_every == 0) {
            Vec3 w[3] = {{s.v[3], s.v[4], s.v[5]}, {s.v[6], s.v[7], s.v[8]}, {s.v[9], s.v[10], s.v[11]}};
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < j; ++k) w[j] = w[j] - w[k] * dot3(w[j], w[k]);
                double nrm = norm3(w[j]);
                lsum[j] += std::log(nrm);
                w[j] = w[j] * (1.0 / nrm);
                s.v[3 + 3 * j] = w[j].x; s.v[4 + 3 * j] = w[j].y; s.v[5 + 3 * j] = w[j].z;
            }
            if (step % 400 == 0) {
                double t_cur = step * dt;
                hist.push_back({t_cur, lsum[0] / t_cur, lsum[1] / t_cur, lsum[2] / t_cur});
            }
        }
    }
    double T_tot = n_steps * dt;
    std::array<double, 3> lam{lsum[0] / T_tot, lsum[1] / T_tot, lsum[2] / T_tot};
    double d_ky = 2.0 + (lam[0] + lam[1]) / std::abs(lam[2]);
    return {lam, lam[0] + lam[1] + lam[2], d_ky, hist, u};
}

inline std::pair<double, double> verify_lorenz_rk4_order(double dt_base = 0.01, double t_end = 1.0) {
    auto integrate = [&](double h) {
        Vec3 u{1.0, 1.0, 1.0};
        int n = static_cast<int>(std::llround(t_end / h));
        auto f = [](const Vec3& st, double) { return lorenz_rhs(st); };
        for (int i = 0; i < n; ++i) u = rk4_step(u, 0.0, h, f);
        return u;
    };
    Vec3 u_ref = integrate(dt_base / 16.0);
    double e1 = norm3(integrate(dt_base) - u_ref);
    double e2 = norm3(integrate(dt_base / 2.0) - u_ref);
    double e4 = norm3(integrate(dt_base / 4.0) - u_ref);
    return {e1 / e2, e2 / e4};
}

} // namespace chaos

#ifndef CHAOS_NO_MAIN
int main(int argc, char** argv) {
    using namespace chaos;
    std::string out_dir = (argc > 1) ? argv[1] : "../../../outputs/day007a";
    std::string fig_dir = (argc > 2) ? argv[2] : "../../../outputs/figures";

    const int n_r = 1600;
    std::vector<double> r_grid(n_r), lam_grid(n_r);
    std::vector<std::pair<double, double>> bif_pts;
    bif_pts.reserve(n_r * 200);
    for (int i = 0; i < n_r; ++i) {
        double r = 2.5 + 1.5 * i / (n_r - 1.0);
        r_grid[i] = r;
        lam_grid[i] = logistic_lyapunov(r, 1000, 200000);
        double x = 0.3141592653589793;
        for (int k = 0; k < 1000; ++k) x = r * x * (1.0 - x);
        for (int k = 0; k < 200; ++k) { x = r * x * (1.0 - x); bif_pts.push_back({r, x}); }
    }
    double lam_r4 = logistic_lyapunov(4.0, 1000, 10000000LL);
    double lam_r4_dyadic = logistic_lyapunov(4.0, 10, 1000, 0.5);
    auto feig = solve_feigenbaum(8);
    double delta_8 = feig.back().delta_n;
    double r_inf = feig[7].R_n + (feig[7].R_n - feig[6].R_n) / (delta_8 - 1.0);

    auto lres = compute_lorenz_lyapunov(0.005, 50.0, 2000.0, 10);
    auto [ratio12, ratio24] = verify_lorenz_rk4_order(0.01, 1.0);

    // Lorenz trajectory for attractor SVG + separation growth
    auto f3 = [](const Vec3& st, double) { return lorenz_rhs(st); };
    Vec3 ua = lres.post_transient, ub = ua + Vec3{1e-9, 0.0, 0.0};
    std::vector<std::array<double, 4>> attr_pts;
    std::vector<std::pair<double, double>> sep_pts;
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int sn = 0;
    for (int i = 0; i <= 7000; ++i) {
        double t = i * 0.005;
        if (i <= 6000) attr_pts.push_back({t, ua.x, ua.y, ua.z});
        double lnd = std::log(std::max(1e-30, norm3(ub - ua)));
        if (i % 4 == 0) sep_pts.push_back({t, lnd});
        if (t >= 2.0 && t <= 15.0) { sx += t; sy += lnd; sxx += t * t; sxy += t * lnd; ++sn; }
        ua = rk4_step(ua, 0.0, 0.005, f3);
        ub = rk4_step(ub, 0.0, 0.005, f3);
    }
    double sep_slope = (sn * sxy - sx * sy) / (sn * sxx - sx * sx);

    { // Write summary & CSVs
        std::ofstream fs(out_dir + "/summary.csv");
        fs << std::setprecision(16) << "metric,value\n"
           << "lam_r4_1e7," << lam_r4 << "\nlam_r4_dyadic_x05," << lam_r4_dyadic << "\n"
           << "feigenbaum_delta8," << delta_8 << "\nfeigenbaum_r_inf," << r_inf << "\n"
           << "lorenz_lam1," << lres.lam[0] << "\nlorenz_lam2," << lres.lam[1] << "\n"
           << "lorenz_lam3," << lres.lam[2] << "\nlorenz_sum_lam," << lres.sum_lam << "\n"
           << "lorenz_d_ky," << lres.d_ky << "\nrk4_ratio_12," << ratio12 << "\n"
           << "rk4_ratio_24," << ratio24 << "\nsep_slope," << sep_slope << "\n";
        std::ofstream ff(out_dir + "/feigenbaum.csv");
        ff << std::setprecision(16) << "n,period,R_n,delta_n\n";
        for (auto& r : feig) ff << r.n << "," << r.period << "," << r.R_n << "," << r.delta_n << "\n";
        std::ofstream fl(out_dir + "/logistic_lyapunov.csv");
        fl << std::setprecision(10) << "r,lambda\n";
        for (int i = 0; i < n_r; ++i) fl << r_grid[i] << "," << lam_grid[i] << "\n";
    }

    auto svg_head = [](std::ofstream& o, int w, int h, const std::string& title) {
        o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
          << w << "\" height=\"" << h << "\" viewBox=\"0 0 " << w << " " << h << "\">\n"
          << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
          << "<text x=\"" << w / 2 << "\" y=\"24\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\" font-weight=\"bold\">" << title << "</text>\n";
    };

    { // SVG 1: Bifurcation + Lyapunov
        std::ofstream o(fig_dir + "/day007a_bifurcation.svg");
        svg_head(o, 900, 620, "Logistic Map Bifurcation Diagram and Lyapunov Exponent lambda(r)");
        o << std::fixed << std::setprecision(2)
          << "<rect x=\"65\" y=\"40\" width=\"810\" height=\"310\" fill=\"none\" stroke=\"#334155\"/>\n"
          << "<rect x=\"65\" y=\"375\" width=\"810\" height=\"205\" fill=\"none\" stroke=\"#334155\"/>\n";
        double y0_line = 375.0 + (1.0 - 0.0) / 3.5 * 205.0;
        o << "<line x1=\"65\" y1=\"" << y0_line << "\" x2=\"875\" y2=\"" << y0_line << "\" stroke=\"#dc2626\" stroke-dasharray=\"4,4\"/>\n"
          << "<g fill=\"#0f172a\" fill-opacity=\"0.18\">\n";
        for (size_t i = 0; i < bif_pts.size(); i += 4) {
            double px = 65.0 + (bif_pts[i].first - 2.5) / 1.5 * 810.0;
            double py = 350.0 - bif_pts[i].second * 310.0;
            o << "<rect x=\"" << px << "\" y=\"" << py << "\" width=\"0.9\" height=\"0.9\"/>\n";
        }
        o << "</g>\n<polyline fill=\"none\" stroke=\"#1d4ed8\" stroke-width=\"1.2\" points=\"";
        for (int i = 0; i < n_r; ++i) {
            double px = 65.0 + (r_grid[i] - 2.5) / 1.5 * 810.0;
            double py = 375.0 + (1.0 - std::clamp(lam_grid[i], -2.5, 1.0)) / 3.5 * 205.0;
            o << px << "," << py << " ";
        }
        o << "\"/>\n<text x=\"470\" y=\"605\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">r in [2.5, 4.0]</text>\n</svg>\n";
    }
    { // SVG 2: Lorenz Attractor (x,z) and (x,y)
        std::ofstream o(fig_dir + "/day007a_lorenz_attractor.svg");
        svg_head(o, 920, 440, "Lorenz Attractor Projections: (x, z) and (x, y) Colored by Time");
        o << std::fixed << std::setprecision(2)
          << "<rect x=\"45\" y=\"40\" width=\"400\" height=\"360\" fill=\"none\" stroke=\"#334155\"/>\n"
          << "<rect x=\"485\" y=\"40\" width=\"400\" height=\"360\" fill=\"none\" stroke=\"#334155\"/>\n";
        for (size_t i = 1; i < attr_pts.size(); ++i) {
            int r_c = static_cast<int>(30 + 190.0 * i / attr_pts.size());
            int b_c = static_cast<int>(220 - 170.0 * i / attr_pts.size());
            double x1 = 45 + (attr_pts[i - 1][1] + 25) / 50.0 * 400, z1 = 400 - attr_pts[i - 1][3] / 55.0 * 360;
            double x2 = 45 + (attr_pts[i][1] + 25) / 50.0 * 400, z2 = 400 - attr_pts[i][3] / 55.0 * 360;
            double u1 = 485 + (attr_pts[i - 1][1] + 25) / 50.0 * 400, y1 = 400 - (attr_pts[i - 1][2] + 30) / 60.0 * 360;
            double u2 = 485 + (attr_pts[i][1] + 25) / 50.0 * 400, y2 = 400 - (attr_pts[i][2] + 30) / 60.0 * 360;
            o << "<line x1=\"" << x1 << "\" y1=\"" << z1 << "\" x2=\"" << x2 << "\" y2=\"" << z2
              << "\" stroke=\"rgb(" << r_c << ",70," << b_c << ")\" stroke-width=\"0.8\"/>\n"
              << "<line x1=\"" << u1 << "\" y1=\"" << y1 << "\" x2=\"" << u2 << "\" y2=\"" << y2
              << "\" stroke=\"rgb(" << r_c << ",70," << b_c << ")\" stroke-width=\"0.8\"/>\n";
        }
        o << "<text x=\"245\" y=\"422\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">(x, z)</text>\n"
          << "<text x=\"685\" y=\"422\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">(x, y)</text>\n</svg>\n";
    }
    { // SVG 3: Lorenz Lyapunov Convergence
        std::ofstream o(fig_dir + "/day007a_lorenz_lyapunov_convergence.svg");
        svg_head(o, 860, 440, "Lorenz Running Lyapunov Exponents and Liouville Sum (-13.6667)");
        o << std::fixed << std::setprecision(2) << "<rect x=\"65\" y=\"40\" width=\"760\" height=\"355\" fill=\"none\" stroke=\"#334155\"/>\n";
        auto map_y = [](double v) { return 40.0 + (2.0 - v) / 18.0 * 355.0; };
        o << "<line x1=\"65\" y1=\"" << map_y(-13.6667) << "\" x2=\"825\" y2=\"" << map_y(-13.6667) << "\" stroke=\"#dc2626\" stroke-dasharray=\"5,4\"/>\n"
          << "<line x1=\"65\" y1=\"" << map_y(0.0) << "\" x2=\"825\" y2=\"" << map_y(0.0) << "\" stroke=\"#94a3b8\" stroke-dasharray=\"3,3\"/>\n";
        const char* cols[4] = {"#1d4ed8", "#059669", "#7c3aed", "#b91c1c"};
        for (int c = 0; c < 4; ++c) {
            o << "<polyline fill=\"none\" stroke=\"" << cols[c] << "\" stroke-width=\"1.5\" points=\"";
            for (auto& h : lres.history) {
                double val = (c < 3) ? h[c + 1] : (h[1] + h[2] + h[3]);
                o << (65.0 + h[0] / 2000.0 * 760.0) << "," << map_y(val) << " ";
            }
            o << "\"/>\n";
        }
        o << "<text x=\"445\" y=\"420\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">Time t in [0, 2000]</text>\n</svg>\n";
    }
    { // SVG 4: Separation Growth
        std::ofstream o(fig_dir + "/day007a_separation_growth.svg");
        svg_head(o, 860, 420, "Trajectory Divergence ln||delta x(t)|| vs Reference Slope lambda_1");
        o << std::fixed << std::setprecision(2) << "<rect x=\"65\" y=\"40\" width=\"760\" height=\"335\" fill=\"none\" stroke=\"#334155\"/>\n";
        auto my = [](double v) { return 40.0 + (5.0 - std::clamp(v, -23.0, 5.0)) / 28.0 * 335.0; };
        o << "<polyline fill=\"none\" stroke=\"#1d4ed8\" stroke-width=\"1.4\" points=\"";
        for (auto& p : sep_pts) o << (65.0 + p.first / 35.0 * 760.0) << "," << my(p.second) << " ";
        o << "\"/>\n<line x1=\"65\" y1=\"" << my(-20.723) << "\" x2=\"" << (65.0 + 24.0 / 35.0 * 760.0)
          << "\" y2=\"" << my(-20.723 + lres.lam[0] * 24.0) << "\" stroke=\"#dc2626\" stroke-width=\"1.5\" stroke-dasharray=\"5,4\"/>\n"
          << "<text x=\"445\" y=\"402\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">Time t in [0, 35]</text>\n</svg>\n";
    }
    return 0;
}
#endif
