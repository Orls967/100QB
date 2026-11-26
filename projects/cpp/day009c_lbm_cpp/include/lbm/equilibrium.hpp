#pragma once

#include "lbm/lattice.hpp"
#include <array>

namespace lbm {

template <LatticeModel L>
inline double compute_feq_i(int i, double rho, double ux, double uy, double uz = 0.0) noexcept {
    const double cu = static_cast<double>(L::cx[static_cast<std::size_t>(i)]) * ux +
                      static_cast<double>(L::cy[static_cast<std::size_t>(i)]) * uy +
                      static_cast<double>(L::cz[static_cast<std::size_t>(i)]) * uz;
    const double u2 = ux * ux + uy * uy + uz * uz;
    return L::w[static_cast<std::size_t>(i)] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
}

template <LatticeModel L>
inline std::array<double, L::Q> compute_feq(double rho, double ux, double uy, double uz = 0.0) noexcept {
    std::array<double, L::Q> feq{};
    const double u2 = ux * ux + uy * uy + uz * uz;
    for (int i = 0; i < L::Q; ++i) {
        const double cu = static_cast<double>(L::cx[static_cast<std::size_t>(i)]) * ux +
                          static_cast<double>(L::cy[static_cast<std::size_t>(i)]) * uy +
                          static_cast<double>(L::cz[static_cast<std::size_t>(i)]) * uz;
        feq[static_cast<std::size_t>(i)] =
            L::w[static_cast<std::size_t>(i)] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
    }
    return feq;
}

// Chapman-Enskog 1st-order initialization from strain-rate tensor S_ab = 0.5*(du_a/dx_b + du_b/dx_a):
// f_i^(1) = -tau * rho * c_s^2 * w_i / c_s^4 * sum_{a,b} (c_{ia} c_{ib} - c_s^2 delta_{ab}) * du_a/dx_b
// Note: Pi_ab^(1) = -2 * rho * c_s^2 * tau * S_ab, so for incompressible velocity gradient grad_u:
inline std::array<double, 9> compute_chapman_enskog_init_2d(
    double rho, double ux, double uy,
    double dux_dx, double dux_dy, double duy_dx, double duy_dy, double tau) noexcept {
    const double pxx = -2.0 * rho * D2Q9::CS2 * tau * dux_dx;
    const double pyy = -2.0 * rho * D2Q9::CS2 * tau * duy_dy;
    const double pxy = -rho * D2Q9::CS2 * tau * (dux_dy + duy_dx);
    auto feq = compute_feq<D2Q9>(rho, ux, uy);
    for (int i = 0; i < 9; ++i) {
        const double cx = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(i)]);
        const double cy = static_cast<double>(D2Q9::cy[static_cast<std::size_t>(i)]);
        const double qxx = cx * cx - D2Q9::CS2;
        const double qyy = cy * cy - D2Q9::CS2;
        const double qxy = cx * cy;
        const double f1 = 4.5 * D2Q9::w[static_cast<std::size_t>(i)] *
                          (qxx * pxx + 2.0 * qxy * pxy + qyy * pyy);
        feq[static_cast<std::size_t>(i)] += f1;
    }
    return feq;
}

}  // namespace lbm
