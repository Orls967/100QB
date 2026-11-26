#include "lbm/analysis.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace lbm {

OlsResult ols_fit(const std::vector<double>& xs, const std::vector<double>& ys) {
    OlsResult res{};
    const std::size_t n = xs.size();
    if (n < 2 || ys.size() != n) return res;

    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mx += xs[i];
        my += ys[i];
    }
    mx /= static_cast<double>(n);
    my /= static_cast<double>(n);

    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = xs[i] - mx;
        const double dy = ys[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (sxx <= 1e-30) return res;

    res.slope = sxy / sxx;
    res.intercept = my - res.slope * mx;
    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double r = ys[i] - (res.intercept + res.slope * xs[i]);
        sse += r * r;
    }
    res.r2 = (syy > 1e-30) ? (1.0 - sse / syy) : 1.0;
    if (n > 2) {
        const double s2 = sse / static_cast<double>(n - 2);
        res.slope_se = std::sqrt(std::max(0.0, s2 / sxx));
    }
    return res;
}

MedianMad median_mad(std::vector<double> xs) {
    MedianMad res{};
    if (xs.empty()) return res;
    std::sort(xs.begin(), xs.end());
    const std::size_t n = xs.size();
    res.median = (n & 1ULL) ? xs[n / 2] : 0.5 * (xs[n / 2 - 1] + xs[n / 2]);
    std::vector<double> devs(n);
    for (std::size_t i = 0; i < n; ++i) devs[i] = std::abs(xs[i] - res.median);
    std::sort(devs.begin(), devs.end());
    res.mad = (n & 1ULL) ? devs[n / 2] : 0.5 * (devs[n / 2 - 1] + devs[n / 2]);
    return res;
}

VortexCenter analyze_cavity_vortices(const Grid2D& g) {
    const int nx = g.nx;
    const int ny = g.ny;
    std::vector<double> ux(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0.0);
    std::vector<double> uy(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0.0);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const auto mc = g.cell_macro(x, y);
            const std::size_t c = static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                                  static_cast<std::size_t>(x);
            ux[c] = mc.ux;
            uy[c] = mc.uy;
        }
    }

    // Centerline summary (vertical centerline x=nx/2 for u_x, horizontal centerline y=ny/2 for u_y)
    VortexCenter vc{};
    vc.u_center_min = 1e30;
    vc.v_center_min = 1e30;
    vc.v_center_max = -1e30;
    const int xc = nx / 2;
    const int yc = ny / 2;
    for (int y = 0; y < ny; ++y) {
        vc.u_center_min = std::min(vc.u_center_min, ux[static_cast<std::size_t>(y) * nx + xc]);
    }
    for (int x = 0; x < nx; ++x) {
        const double v = uy[static_cast<std::size_t>(yc) * nx + x];
        vc.v_center_min = std::min(vc.v_center_min, v);
        vc.v_center_max = std::max(vc.v_center_max, v);
    }

    // Streamfunction integration: d(psi)/dy = ux (with psi = 0 on bottom boundary y=-0.5)
    std::vector<double> psi(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0.0);
    const double dx = 1.0 / static_cast<double>(nx);
    for (int x = 0; x < nx; ++x) {
        double p = 0.0;
        for (int y = 0; y < ny; ++y) {
            const double u_here = ux[static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                                     static_cast<std::size_t>(x)];
            p += (y == 0 ? 0.5 * u_here : 0.5 * (ux[static_cast<std::size_t>(y - 1) * nx + x] + u_here)) * dx;
            psi[static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                static_cast<std::size_t>(x)] = p;
        }
    }

    // Locate global interior extremum (negative peak of psi for clockwise top-lid flow)
    int best_x = nx / 2, best_y = ny / 2;
    double best_val = 1e30;
    for (int y = ny / 5; y < ny - 2; ++y) {
        for (int x = 2; x < nx - 2; ++x) {
            const double v = psi[static_cast<std::size_t>(y) * nx + x];
            if (v < best_val) {
                best_val = v;
                best_x = x;
                best_y = y;
            }
        }
    }

    // Parabolic sub-grid quadratic interpolation around (best_x, best_y)
    double sub_x = static_cast<double>(best_x) + 0.5;
    double sub_y = static_cast<double>(best_y) + 0.5;
    if (best_x > 1 && best_x < nx - 2 && best_y > 1 && best_y < ny - 2) {
        const double fm = psi[static_cast<std::size_t>(best_y) * nx + (best_x - 1)];
        const double f0 = best_val;
        const double fp = psi[static_cast<std::size_t>(best_y) * nx + (best_x + 1)];
        const double denom_x = fm - 2.0 * f0 + fp;
        if (std::abs(denom_x) > 1e-18) sub_x += 0.5 * (fm - fp) / denom_x;

        const double gm = psi[static_cast<std::size_t>(best_y - 1) * nx + best_x];
        const double gp = psi[static_cast<std::size_t>(best_y + 1) * nx + best_x];
        const double denom_y = gm - 2.0 * f0 + gp;
        if (std::abs(denom_y) > 1e-18) sub_y += 0.5 * (gm - gp) / denom_y;
    }

    vc.x = sub_x / static_cast<double>(nx);
    vc.y = sub_y / static_cast<double>(ny);
    vc.psi_ext = best_val;

    // Count distinct positive-psi secondary corner vortices in BR, BL, TL quadrants
    auto has_pos_vortex = [&](int x0, int x1, int y0, int y1) -> bool {
        double pmax = 0.0;
        for (int y = y0; y < y1; ++y) {
            for (int x = x0; x < x1; ++x) {
                pmax = std::max(pmax, psi[static_cast<std::size_t>(y) * nx + x]);
            }
        }
        return pmax > 1e-6 * std::abs(best_val);
    };
    int eddies = 0;
    if (has_pos_vortex(nx / 2, nx - 1, 1, ny / 2)) ++eddies;
    if (has_pos_vortex(1, nx / 2, 1, ny / 2)) ++eddies;
    if (has_pos_vortex(1, nx / 3, (2 * ny) / 3, ny - 1)) ++eddies;
    vc.corner_eddies = eddies;
    return vc;
}

double measure_recirculation_length(const Grid2D& g, double cx, double cy, double diameter) {
    const int y_mid = static_cast<int>(std::llround(cy));
    const int x_back = static_cast<int>(std::floor(cx + 0.5 * diameter)) + 1;
    double prev_x = cx + 0.5 * diameter;
    double prev_u = -1e-6;
    for (int x = x_back; x < g.nx - 2; ++x) {
        if (g.solid[static_cast<std::size_t>(y_mid) * static_cast<std::size_t>(g.nx) +
                    static_cast<std::size_t>(x)]) {
            continue;
        }
        const double u = g.cell_macro(x, y_mid).ux;
        if (x > x_back && prev_u <= 0.0 && u > 0.0) {
            const double frac = -prev_u / (u - prev_u);
            const double x_zero = prev_x + frac;
            return (x_zero - (cx + 0.5 * diameter)) / diameter;
        }
        if (u < 0.0 || x == x_back) {
            prev_u = u;
            prev_x = static_cast<double>(x);
        }
    }
    return 0.0;
}

StrouhalStats analyze_strouhal(const std::vector<double>& ts, double diameter, double u_inf) {
    StrouhalStats st{};
    const std::size_t n = ts.size();
    if (n < 16) return st;

    double mean = 0.0;
    for (double v : ts) mean += v;
    mean /= static_cast<double>(n);

    double cmin = ts[0], cmax = ts[0];
    std::vector<double> up_crossings;
    for (std::size_t i = 1; i < n; ++i) {
        cmin = std::min(cmin, ts[i]);
        cmax = std::max(cmax, ts[i]);
        if (ts[i - 1] - mean <= 0.0 && ts[i] - mean > 0.0) {
            const double denom = ts[i] - ts[i - 1];
            const double frac = (denom > 1e-18) ? (mean - ts[i - 1]) / denom : 0.5;
            up_crossings.push_back(static_cast<double>(i - 1) + frac);
        }
    }
    st.cl_amp = 0.5 * (cmax - cmin);
    st.zero_crossings = static_cast<int>(up_crossings.size());
    if (up_crossings.size() >= 2) {
        const double avg_period =
            (up_crossings.back() - up_crossings.front()) / static_cast<double>(up_crossings.size() - 1);
        if (avg_period > 0.0) st.freq_zc = 1.0 / avg_period;
    }

    // Naive DFT power spectrum peak + parabolic sub-bin interpolation for Strouhal cross-validation
    double best_p = -1.0;
    std::size_t best_k = 1;
    const std::size_t kmax = n / 2;
    constexpr double PI = 3.14159265358979323846;
    std::vector<double> pws(kmax, 0.0);
    for (std::size_t k = 1; k < kmax; ++k) {
        double re = 0.0, im = 0.0;
        const double ang0 = 2.0 * PI * static_cast<double>(k) / static_cast<double>(n);
        for (std::size_t t = 0; t < n; ++t) {
            const double val = ts[t] - mean;
            re += val * std::cos(ang0 * static_cast<double>(t));
            im -= val * std::sin(ang0 * static_cast<double>(t));
        }
        const double pw = re * re + im * im;
        pws[k] = pw;
        if (pw > best_p) {
            best_p = pw;
            best_k = k;
        }
    }
    double k_interp = static_cast<double>(best_k);
    if (best_k > 1 && best_k + 1 < kmax) {
        const double pm = pws[best_k - 1], p0 = pws[best_k], pp = pws[best_k + 1];
        const double den = pm - 2.0 * p0 + pp;
        if (std::abs(den) > 1e-18) k_interp += 0.5 * (pm - pp) / den;
    }
    st.freq_dft = k_interp / static_cast<double>(n);
    const double f_use = st.freq_dft;
    st.strouhal = f_use * diameter / u_inf;
    return st;
}

double exact_square_duct_ux(double y_pos, double z_pos, double height, double fx, double rho0, double nu) {
    // Domain y_pos, z_pos in [0, H]. With half-width a = H/2, 16*a^2 = 4*H^2.
    constexpr double PI = 3.14159265358979323846;
    const double yc = y_pos - 0.5 * height;
    const double zc = z_pos - 0.5 * height;
    const double coeff = 4.0 * fx * height * height / (PI * PI * PI * rho0 * nu);
    double sum = 0.0;
    for (int n = 1; n <= 49; n += 2) {
        const double nf = static_cast<double>(n);
        const double sign = (((n - 1) / 2) % 2 == 0) ? 1.0 : -1.0;
        const double arg_z = nf * PI * zc / height;
        const double arg_h = 0.5 * nf * PI;
        // Stable evaluation of cosh(arg_z) / cosh(arg_h)
        const double ratio = std::cosh(arg_z) / std::cosh(arg_h);
        const double term = (sign / (nf * nf * nf)) * (1.0 - ratio) * std::cos(nf * PI * yc / height);
        sum += term;
    }
    return coeff * sum;
}

double measure_copy_bandwidth_gbps(std::size_t num_doubles, int repeats) {
    std::vector<double> src(num_doubles, 1.25);
    std::vector<double> dst(num_doubles, 0.0);
    // Warmup
    std::copy(src.begin(), src.end(), dst.begin());
    std::vector<double> bws;
    bws.reserve(static_cast<std::size_t>(repeats));
    for (int r = 0; r < repeats; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < num_doubles; ++i) {
            dst[i] = src[i];
        }
        const auto t1 = std::chrono::steady_clock::now();
        volatile double sink = dst[num_doubles / 2];
        (void)sink;
        const double sec = std::chrono::duration<double>(t1 - t0).count();
        // 1 read + 1 write = 16 bytes per double copied
        const double gb = static_cast<double>(num_doubles) * 16.0 / 1.0e9;
        bws.push_back(gb / std::max(sec, 1e-9));
    }
    return median_mad(bws).median;
}

RleMetrics compress_delta_rle(const std::vector<double>& field, double quant_step) {
    RleMetrics met{};
    if (field.empty()) return met;
    met.raw_bytes = field.size() * sizeof(double);

    std::vector<std::int64_t> q(field.size());
    for (std::size_t i = 0; i < field.size(); ++i) {
        q[i] = static_cast<std::int64_t>(std::llround(field[i] / quant_step));
    }

    std::vector<std::int64_t> deltas(q.size());
    deltas[0] = q[0];
    for (std::size_t i = 1; i < q.size(); ++i) {
        deltas[i] = q[i] - q[i - 1];
    }

    // RLE encode deltas as (value, run_length) pairs
    struct Run {
        std::int64_t val;
        std::uint32_t len;
    };
    std::vector<Run> runs;
    std::int64_t cur = deltas[0];
    std::uint32_t cnt = 1;
    for (std::size_t i = 1; i < deltas.size(); ++i) {
        if (deltas[i] == cur) {
            ++cnt;
        } else {
            runs.push_back({cur, cnt});
            cur = deltas[i];
            cnt = 1;
        }
    }
    runs.push_back({cur, cnt});

    // Variable-byte storage cost model (1 byte tag + small varint 1/2/4/8B)
    std::size_t cbytes = 8;
    for (const auto& run : runs) {
        const std::int64_t av = std::abs(run.val);
        const std::size_t vbytes = (av < 64) ? 1 : ((av < 16384) ? 2 : 4);
        const std::size_t lbytes = (run.len < 128) ? 1 : 2;
        cbytes += vbytes + lbytes;
    }

    // Exact roundtrip verification on integer quantized sequence
    std::vector<std::int64_t> decoded;
    decoded.reserve(q.size());
    std::int64_t acc = 0;
    bool first = true;
    for (const auto& run : runs) {
        for (std::uint32_t k = 0; k < run.len; ++k) {
            if (first) {
                acc = run.val;
                first = false;
            } else {
                acc += run.val;
            }
            decoded.push_back(acc);
        }
    }
    std::int64_t max_err = 0;
    for (std::size_t i = 0; i < q.size(); ++i) {
        max_err = std::max(max_err, std::abs(decoded[i] - q[i]));
    }

    met.compressed_bytes = cbytes;
    met.ratio = static_cast<double>(met.raw_bytes) / static_cast<double>(std::max<std::size_t>(1, cbytes));
    met.max_roundtrip_q_err = max_err;
    return met;
}

}  // namespace lbm
