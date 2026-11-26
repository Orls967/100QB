#pragma once

#include "lbm/equilibrium.hpp"
#include "lbm/lattice.hpp"
#include <array>
#include <cmath>
#include <span>

namespace lbm {

struct Macro2D {
    double rho{};
    double ux{};
    double uy{};
    double pxx_neq{};
    double pxy_neq{};
    double pyy_neq{};
};

struct Macro3D {
    double rho{};
    double ux{};
    double uy{};
    double uz{};
};

inline Macro2D compute_macro_2d(std::span<const double, 9> f, double fx = 0.0, double fy = 0.0) noexcept {
    double rho = 0.0, mx = 0.0, my = 0.0;
    for (int i = 0; i < 9; ++i) {
        const double fi = f[static_cast<std::size_t>(i)];
        rho += fi;
        mx += static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]) * fi;
        my += static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]) * fi;
    }
    const double inv_r = (rho > 1e-30) ? (1.0 / rho) : 0.0;
    const double ux = (mx + 0.5 * fx) * inv_r;
    const double uy = (my + 0.5 * fy) * inv_r;
    const auto feq = compute_feq<D2Q9>(rho, ux, uy);
    double pxx = 0.0, pxy = 0.0, pyy = 0.0;
    for (int i = 0; i < 9; ++i) {
        const double fneq = f[static_cast<std::size_t>(i)] - feq[static_cast<std::size_t>(i)];
        const double cx = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]);
        const double cy = static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]);
        pxx += cx * cx * fneq;
        pxy += cx * cy * fneq;
        pyy += cy * cy * fneq;
    }
    return {rho, ux, uy, pxx, pxy, pyy};
}

template <LatticeModel L>
inline Macro3D compute_macro_3d(const double* f, double fx = 0.0, double fy = 0.0, double fz = 0.0) noexcept {
    double rho = 0.0, mx = 0.0, my = 0.0, mz = 0.0;
    for (int i = 0; i < L::Q; ++i) {
        const double fi = f[i];
        rho += fi;
        mx += static_cast<double>(L::cx[static_cast<std::size_t>(i)]) * fi;
        my += static_cast<double>(L::cy[static_cast<std::size_t>(i)]) * fi;
        mz += static_cast<double>(L::cz[static_cast<std::size_t>(i)]) * fi;
    }
    const double inv_r = (rho > 1e-30) ? (1.0 / rho) : 0.0;
    return {rho, (mx + 0.5 * fx) * inv_r, (my + 0.5 * fy) * inv_r, (mz + 0.5 * fz) * inv_r};
}

// Lallemand & Luo (2000) orthogonal MRT matrix M (9x9) for D2Q9
inline constexpr std::array<std::array<double, 9>, 9> D2Q9_MRT_M = {{
    {1, 1, 1, 1, 1, 1, 1, 1, 1},
    {-4, -1, -1, -1, -1, 2, 2, 2, 2},
    {4, -2, -2, -2, -2, 1, 1, 1, 1},
    {0, 1, 0, -1, 0, 1, -1, -1, 1},
    {0, -2, 0, 2, 0, 1, -1, -1, 1},
    {0, 0, 1, 0, -1, 1, 1, -1, -1},
    {0, 0, -2, 0, 2, 1, 1, -1, -1},
    {0, 1, -1, 1, -1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, -1, 1, -1}
}};

inline constexpr std::array<double, 9> D2Q9_MRT_DIAG_MMT = {9.0, 36.0, 36.0, 6.0, 12.0, 6.0, 12.0, 4.0, 4.0};

// d'Humieres et al. (2002) orthogonal MRT matrix M (19x19) evaluated from D3Q19 velocity vectors
consteval std::array<std::array<double, 19>, 19> build_d3q19_mrt_m() {
    std::array<std::array<double, 19>, 19> m{};
    for (int i = 0; i < 19; ++i) {
        const auto ii = static_cast<std::size_t>(i);
        const double x = static_cast<double>(D3Q19::cx[ii]);
        const double y = static_cast<double>(D3Q19::cy[ii]);
        const double z = static_cast<double>(D3Q19::cz[ii]);
        const double c2 = x * x + y * y + z * z;
        m[0][ii] = 1.0;
        m[1][ii] = 19.0 * c2 - 30.0;
        m[2][ii] = (21.0 * c2 * c2 - 53.0 * c2 + 24.0) / 2.0;
        m[3][ii] = x;
        m[4][ii] = (5.0 * c2 - 9.0) * x;
        m[5][ii] = y;
        m[6][ii] = (5.0 * c2 - 9.0) * y;
        m[7][ii] = z;
        m[8][ii] = (5.0 * c2 - 9.0) * z;
        m[9][ii] = 3.0 * x * x - c2;
        m[10][ii] = (3.0 * c2 - 5.0) * (3.0 * x * x - c2);
        m[11][ii] = y * y - z * z;
        m[12][ii] = (3.0 * c2 - 5.0) * (y * y - z * z);
        m[13][ii] = x * y;
        m[14][ii] = y * z;
        m[15][ii] = x * z;
        m[16][ii] = (y * y - z * z) * x;
        m[17][ii] = (z * z - x * x) * y;
        m[18][ii] = (x * x - y * y) * z;
    }
    return m;
}

inline constexpr std::array<std::array<double, 19>, 19> D3Q19_MRT_M = build_d3q19_mrt_m();

consteval std::array<double, 19> build_d3q19_diag_mmt() {
    std::array<double, 19> d{};
    for (int r = 0; r < 19; ++r) {
        double s = 0.0;
        for (int c = 0; c < 19; ++c) {
            const double v = D3Q19_MRT_M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
            s += v * v;
        }
        d[static_cast<std::size_t>(r)] = s;
    }
    return d;
}

inline constexpr std::array<double, 19> D3Q19_MRT_DIAG_MMT = build_d3q19_diag_mmt();

struct OrthoReport {
    double max_off_diag_err = 0.0;
    double max_inv_err = 0.0;
};

OrthoReport verify_d2q9_mrt_orthogonality();
OrthoReport verify_d3q19_mrt_orthogonality();

}  // namespace lbm
