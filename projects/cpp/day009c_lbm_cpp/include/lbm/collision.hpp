#pragma once

#include "lbm/equilibrium.hpp"
#include "lbm/lattice.hpp"
#include "lbm/moments.hpp"
#include <array>
#include <cstddef>

namespace lbm {

enum class CollisionOp { BGK, TRT, MRT, REGULARIZED };

inline constexpr const char* collision_name(CollisionOp op) noexcept {
    switch (op) {
        case CollisionOp::BGK: return "BGK";
        case CollisionOp::TRT: return "TRT";
        case CollisionOp::MRT: return "MRT";
        case CollisionOp::REGULARIZED: return "Regularized";
    }
    return "BGK";
}

struct CollisionParams {
    CollisionOp op = CollisionOp::BGK;
    double tau = 0.8;
    double trt_lambda = 3.0 / 16.0;  // Magic parameter (tau_plus - 0.5)*(tau_minus - 0.5)
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    constexpr double nu() const noexcept { return (tau - 0.5) / 3.0; }
    constexpr double omega() const noexcept { return 1.0 / tau; }
    constexpr double tau_minus() const noexcept { return 0.5 + trt_lambda / (tau - 0.5); }
    constexpr double omega_minus() const noexcept { return 1.0 / tau_minus(); }
};

// Unified 2D D2Q9 collision (with optional Guo body forcing) operating in-place on 9 populations
inline void collide_d2q9_cell(double* f, const CollisionParams& cp) noexcept {
    double rho = 0.0, mx = 0.0, my = 0.0;
    for (int i = 0; i < 9; ++i) {
        rho += f[i];
        mx += static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]) * f[i];
        my += static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]) * f[i];
    }
    const bool has_force = (cp.fx != 0.0 || cp.fy != 0.0);
    const double inv_r = (rho > 1e-30) ? (1.0 / rho) : 0.0;
    const double ux = (mx + 0.5 * cp.fx) * inv_r;
    const double uy = (my + 0.5 * cp.fy) * inv_r;
    const double omega_p = cp.omega();

    // Compute Guo forcing per direction if active (exact symmetric+antisymmetric splitting for TRT)
    std::array<double, 9> gforce{};
    if (has_force) {
        std::array<double, 9> raw_g{};
        for (int i = 0; i < 9; ++i) {
            const double cx = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]);
            const double cy = static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]);
            const double cu = cx * ux + cy * uy;
            const double term_x = 3.0 * (cx - ux) + 9.0 * cu * cx;
            const double term_y = 3.0 * (cy - uy) + 9.0 * cu * cy;
            raw_g[static_cast<std::size_t>(i)] =
                D2Q9::w[static_cast<std::size_t>(i)] * (term_x * cp.fx + term_y * cp.fy);
        }
        if (cp.op == CollisionOp::TRT) {
            const double fac_p = 1.0 - 0.5 * omega_p;
            const double fac_m = 1.0 - 0.5 * cp.omega_minus();
            for (int i = 0; i < 9; ++i) {
                const int j = D2Q9::opp[static_cast<std::size_t>(i)];
                const double gp = 0.5 * (raw_g[static_cast<std::size_t>(i)] + raw_g[static_cast<std::size_t>(j)]);
                const double gm = 0.5 * (raw_g[static_cast<std::size_t>(i)] - raw_g[static_cast<std::size_t>(j)]);
                gforce[static_cast<std::size_t>(i)] = fac_p * gp + fac_m * gm;
            }
        } else {
            const double fac = 1.0 - 0.5 * omega_p;
            for (int i = 0; i < 9; ++i) gforce[static_cast<std::size_t>(i)] = fac * raw_g[static_cast<std::size_t>(i)];
        }
    }

    if (cp.op == CollisionOp::BGK) {
        const auto feq = compute_feq<D2Q9>(rho, ux, uy);
        for (int i = 0; i < 9; ++i) {
            f[i] += omega_p * (feq[static_cast<std::size_t>(i)] - f[i]) + gforce[static_cast<std::size_t>(i)];
        }
        return;
    }

    if (cp.op == CollisionOp::TRT) {
        const auto feq = compute_feq<D2Q9>(rho, ux, uy);
        const double omega_m = cp.omega_minus();
        std::array<double, 9> fpost{};
        for (int i = 0; i < 9; ++i) {
            const int j = D2Q9::opp[static_cast<std::size_t>(i)];
            const double f_plus = 0.5 * (f[i] + f[j]);
            const double f_minus = 0.5 * (f[i] - f[j]);
            const double feq_plus = 0.5 * (feq[static_cast<std::size_t>(i)] + feq[static_cast<std::size_t>(j)]);
            const double feq_minus = 0.5 * (feq[static_cast<std::size_t>(i)] - feq[static_cast<std::size_t>(j)]);
            fpost[static_cast<std::size_t>(i)] =
                f[i] - omega_p * (f_plus - feq_plus) - omega_m * (f_minus - feq_minus) +
                gforce[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < 9; ++i) f[i] = fpost[static_cast<std::size_t>(i)];
        return;
    }

    if (cp.op == CollisionOp::REGULARIZED) {
        const auto feq = compute_feq<D2Q9>(rho, ux, uy);
        double pxx = 0.0, pxy = 0.0, pyy = 0.0;
        for (int i = 0; i < 9; ++i) {
            const double fneq = f[i] - feq[static_cast<std::size_t>(i)];
            const double cx = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]);
            const double cy = static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]);
            pxx += cx * cx * fneq;
            pxy += cx * cy * fneq;
            pyy += cy * cy * fneq;
        }
        for (int i = 0; i < 9; ++i) {
            const double cx = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]);
            const double cy = static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]);
            const double qxx = cx * cx - D2Q9::CS2;
            const double qyy = cy * cy - D2Q9::CS2;
            const double qxy = cx * cy;
            const double f1 = 4.5 * D2Q9::w[static_cast<std::size_t>(i)] *
                              (qxx * pxx + 2.0 * qxy * pxy + qyy * pyy);
            f[i] = feq[static_cast<std::size_t>(i)] + (1.0 - omega_p) * f1 + gforce[static_cast<std::size_t>(i)];
        }
        return;
    }

    // MRT (Lallemand-Luo 2000)
    const double u2 = ux * ux + uy * uy;
    std::array<double, 9> m{};
    for (int r = 0; r < 9; ++r) {
        double s = 0.0;
        for (int c = 0; c < 9; ++c) {
            s += D2Q9_MRT_M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] * f[c];
        }
        m[static_cast<std::size_t>(r)] = s;
    }
    const std::array<double, 9> meq = {
        rho,
        rho * (-2.0 + 3.0 * u2),
        rho * (1.0 - 3.0 * u2),
        rho * ux,
        -rho * ux,
        rho * uy,
        -rho * uy,
        rho * (ux * ux - uy * uy),
        rho * ux * uy
    };
    // Lallemand & Luo (2000) tunable stability relaxation rates (keeps nu = (tau-0.5)/3 exact)
    const std::array<double, 9> s_rates = {0.0, 1.64, 1.54, 0.0, 1.70, 0.0, 1.70, omega_p, omega_p};
    for (int r = 0; r < 9; ++r) {
        m[static_cast<std::size_t>(r)] -= s_rates[static_cast<std::size_t>(r)] *
                                         (m[static_cast<std::size_t>(r)] - meq[static_cast<std::size_t>(r)]);
    }
    for (int c = 0; c < 9; ++c) {
        double val = 0.0;
        for (int r = 0; r < 9; ++r) {
            val += D2Q9_MRT_M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] *
                   (m[static_cast<std::size_t>(r)] / D2Q9_MRT_DIAG_MMT[static_cast<std::size_t>(r)]);
        }
        f[c] = val + gforce[static_cast<std::size_t>(c)];
    }
}

void collide_d3q19_cell(double* f, const CollisionParams& cp) noexcept;

}  // namespace lbm
