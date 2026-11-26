#include "lbm/lattice.hpp"
#include "lbm/collision.hpp"
#include "lbm/moments.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace lbm {

namespace {

template <LatticeModel L>
IsotropyReport check_lattice_isotropy_impl() {
    IsotropyReport rep{};
    double wsum = 0.0;
    for (int i = 0; i < L::Q; ++i) wsum += L::w[static_cast<std::size_t>(i)];
    rep.weight_sum_err = std::abs(wsum - 1.0);

    // Check opposite table involution & sign flip
    for (int i = 0; i < L::Q; ++i) {
        const auto ii = static_cast<std::size_t>(i);
        const int j = L::opp[ii];
        const auto jj = static_cast<std::size_t>(j);
        const int sx = L::cx[ii] + L::cx[jj];
        const int sy = L::cy[ii] + L::cy[jj];
        const int sz = L::cz[ii] + L::cz[jj];
        rep.opp_check_err = std::max(rep.opp_check_err,
                                     static_cast<double>(std::abs(sx) + std::abs(sy) + std::abs(sz)));
        if (L::opp[jj] != i) rep.opp_check_err = 1.0;
    }

    auto vel = [](int i, int axis) -> double {
        const auto ii = static_cast<std::size_t>(i);
        return static_cast<double>(axis == 0 ? L::cx[ii] : (axis == 1 ? L::cy[ii] : L::cz[ii]));
    };

    const double cs2 = L::CS2;
    const double cs4 = cs2 * cs2;
    const int ndim = L::D;

    // Order 1
    for (int a = 0; a < ndim; ++a) {
        double s = 0.0;
        for (int i = 0; i < L::Q; ++i) s += L::w[static_cast<std::size_t>(i)] * vel(i, a);
        rep.order1_err = std::max(rep.order1_err, std::abs(s));
    }

    // Order 2: sum w_i c_ia c_ib == cs2 * delta_ab
    for (int a = 0; a < ndim; ++a) {
        for (int b = 0; b < ndim; ++b) {
            double s = 0.0;
            for (int i = 0; i < L::Q; ++i)
                s += L::w[static_cast<std::size_t>(i)] * vel(i, a) * vel(i, b);
            const double target = (a == b) ? cs2 : 0.0;
            rep.order2_err = std::max(rep.order2_err, std::abs(s - target));
        }
    }

    // Order 3: sum w_i c_ia c_ib c_ic == 0
    for (int a = 0; a < ndim; ++a) {
        for (int b = 0; b < ndim; ++b) {
            for (int c = 0; c < ndim; ++c) {
                double s = 0.0;
                for (int i = 0; i < L::Q; ++i)
                    s += L::w[static_cast<std::size_t>(i)] * vel(i, a) * vel(i, b) * vel(i, c);
                rep.order3_err = std::max(rep.order3_err, std::abs(s));
            }
        }
    }

    // Order 4: sum w_i c_ia c_ib c_ic c_id == cs4 * (d_ab d_cd + d_ac d_bd + d_ad d_bc)
    for (int a = 0; a < ndim; ++a) {
        for (int b = 0; b < ndim; ++b) {
            for (int c = 0; c < ndim; ++c) {
                for (int d = 0; d < ndim; ++d) {
                    double s = 0.0;
                    for (int i = 0; i < L::Q; ++i)
                        s += L::w[static_cast<std::size_t>(i)] * vel(i, a) * vel(i, b) *
                             vel(i, c) * vel(i, d);
                    const double delta_part =
                        ((a == b && c == d) ? 1.0 : 0.0) +
                        ((a == c && b == d) ? 1.0 : 0.0) +
                        ((a == d && b == c) ? 1.0 : 0.0);
                    rep.order4_err = std::max(rep.order4_err, std::abs(s - cs4 * delta_part));
                }
            }
        }
    }

    // Order 5: odd moments == 0
    for (int a = 0; a < ndim; ++a) {
        double s = 0.0;
        for (int i = 0; i < L::Q; ++i) {
            const double v = vel(i, a);
            s += L::w[static_cast<std::size_t>(i)] * v * v * v * v * v;
        }
        rep.order5_err = std::max(rep.order5_err, std::abs(s));
    }

    return rep;
}

}  // namespace

IsotropyReport check_d2q9_isotropy() { return check_lattice_isotropy_impl<D2Q9>(); }
IsotropyReport check_d3q19_isotropy() { return check_lattice_isotropy_impl<D3Q19>(); }
IsotropyReport check_d3q27_isotropy() { return check_lattice_isotropy_impl<D3Q27>(); }

OrthoReport verify_d2q9_mrt_orthogonality() {
    OrthoReport rep{};
    for (int r1 = 0; r1 < 9; ++r1) {
        for (int r2 = 0; r2 < 9; ++r2) {
            double dot = 0.0;
            for (int c = 0; c < 9; ++c) {
                dot += D2Q9_MRT_M[static_cast<std::size_t>(r1)][static_cast<std::size_t>(c)] *
                       D2Q9_MRT_M[static_cast<std::size_t>(r2)][static_cast<std::size_t>(c)];
            }
            if (r1 != r2) {
                rep.max_off_diag_err = std::max(rep.max_off_diag_err, std::abs(dot));
            } else {
                rep.max_inv_err =
                    std::max(rep.max_inv_err,
                             std::abs(dot - D2Q9_MRT_DIAG_MMT[static_cast<std::size_t>(r1)]));
            }
        }
    }
    return rep;
}

OrthoReport verify_d3q19_mrt_orthogonality() {
    OrthoReport rep{};
    for (int r1 = 0; r1 < 19; ++r1) {
        for (int r2 = 0; r2 < 19; ++r2) {
            double dot = 0.0;
            for (int c = 0; c < 19; ++c) {
                dot += D3Q19_MRT_M[static_cast<std::size_t>(r1)][static_cast<std::size_t>(c)] *
                       D3Q19_MRT_M[static_cast<std::size_t>(r2)][static_cast<std::size_t>(c)];
            }
            if (r1 != r2) {
                rep.max_off_diag_err = std::max(rep.max_off_diag_err, std::abs(dot));
            } else {
                rep.max_inv_err =
                    std::max(rep.max_inv_err,
                             std::abs(dot - D3Q19_MRT_DIAG_MMT[static_cast<std::size_t>(r1)]));
            }
        }
    }
    return rep;
}

void collide_d3q19_cell(double* f, const CollisionParams& cp) noexcept {
    double rho = 0.0, mx = 0.0, my = 0.0, mz = 0.0;
    for (int i = 0; i < 19; ++i) {
        const auto ii = static_cast<std::size_t>(i);
        rho += f[i];
        mx += static_cast<double>(D3Q19::cx[ii]) * f[i];
        my += static_cast<double>(D3Q19::cy[ii]) * f[i];
        mz += static_cast<double>(D3Q19::cz[ii]) * f[i];
    }
    const bool has_force = (cp.fx != 0.0 || cp.fy != 0.0 || cp.fz != 0.0);
    const double inv_r = (rho > 1e-30) ? (1.0 / rho) : 0.0;
    const double ux = (mx + 0.5 * cp.fx) * inv_r;
    const double uy = (my + 0.5 * cp.fy) * inv_r;
    const double uz = (mz + 0.5 * cp.fz) * inv_r;
    const double omega_p = cp.omega();

    std::array<double, 19> gforce{};
    if (has_force) {
        std::array<double, 19> raw_g{};
        for (int i = 0; i < 19; ++i) {
            const auto ii = static_cast<std::size_t>(i);
            const double cx = static_cast<double>(D3Q19::cx[ii]);
            const double cy = static_cast<double>(D3Q19::cy[ii]);
            const double cz = static_cast<double>(D3Q19::cz[ii]);
            const double cu = cx * ux + cy * uy + cz * uz;
            raw_g[ii] = D3Q19::w[ii] * ((3.0 * (cx - ux) + 9.0 * cu * cx) * cp.fx +
                                        (3.0 * (cy - uy) + 9.0 * cu * cy) * cp.fy +
                                        (3.0 * (cz - uz) + 9.0 * cu * cz) * cp.fz);
        }
        if (cp.op == CollisionOp::TRT) {
            const double fac_p = 1.0 - 0.5 * omega_p;
            const double fac_m = 1.0 - 0.5 * cp.omega_minus();
            for (int i = 0; i < 19; ++i) {
                const auto ii = static_cast<std::size_t>(i);
                const auto jj = static_cast<std::size_t>(D3Q19::opp[ii]);
                const double gp = 0.5 * (raw_g[ii] + raw_g[jj]);
                const double gm = 0.5 * (raw_g[ii] - raw_g[jj]);
                gforce[ii] = fac_p * gp + fac_m * gm;
            }
        } else {
            const double fac = 1.0 - 0.5 * omega_p;
            for (int i = 0; i < 19; ++i) {
                gforce[static_cast<std::size_t>(i)] = fac * raw_g[static_cast<std::size_t>(i)];
            }
        }
    }

    const auto feq = compute_feq<D3Q19>(rho, ux, uy, uz);

    if (cp.op == CollisionOp::BGK) {
        for (int i = 0; i < 19; ++i) {
            f[i] += omega_p * (feq[static_cast<std::size_t>(i)] - f[i]) +
                    gforce[static_cast<std::size_t>(i)];
        }
        return;
    }

    if (cp.op == CollisionOp::TRT) {
        const double omega_m = cp.omega_minus();
        std::array<double, 19> fpost{};
        for (int i = 0; i < 19; ++i) {
            const auto ii = static_cast<std::size_t>(i);
            const auto jj = static_cast<std::size_t>(D3Q19::opp[ii]);
            const double fp = 0.5 * (f[i] + f[jj]);
            const double fm = 0.5 * (f[i] - f[jj]);
            const double ep = 0.5 * (feq[ii] + feq[jj]);
            const double em = 0.5 * (feq[ii] - feq[jj]);
            fpost[ii] = f[i] - omega_p * (fp - ep) - omega_m * (fm - em) + gforce[ii];
        }
        for (int i = 0; i < 19; ++i) f[i] = fpost[static_cast<std::size_t>(i)];
        return;
    }

    if (cp.op == CollisionOp::REGULARIZED) {
        double pxx = 0.0, pyy = 0.0, pzz = 0.0, pxy = 0.0, pxz = 0.0, pyz = 0.0;
        for (int i = 0; i < 19; ++i) {
            const auto ii = static_cast<std::size_t>(i);
            const double fneq = f[i] - feq[ii];
            const double cx = static_cast<double>(D3Q19::cx[ii]);
            const double cy = static_cast<double>(D3Q19::cy[ii]);
            const double cz = static_cast<double>(D3Q19::cz[ii]);
            pxx += cx * cx * fneq;
            pyy += cy * cy * fneq;
            pzz += cz * cz * fneq;
            pxy += cx * cy * fneq;
            pxz += cx * cz * fneq;
            pyz += cy * cz * fneq;
        }
        for (int i = 0; i < 19; ++i) {
            const auto ii = static_cast<std::size_t>(i);
            const double cx = static_cast<double>(D3Q19::cx[ii]);
            const double cy = static_cast<double>(D3Q19::cy[ii]);
            const double cz = static_cast<double>(D3Q19::cz[ii]);
            const double f1 = 4.5 * D3Q19::w[ii] *
                              ((cx * cx - D3Q19::CS2) * pxx +
                               (cy * cy - D3Q19::CS2) * pyy +
                               (cz * cz - D3Q19::CS2) * pzz +
                               2.0 * (cx * cy * pxy + cx * cz * pxz + cy * cz * pyz));
            f[i] = feq[ii] + (1.0 - omega_p) * f1 + gforce[ii];
        }
        return;
    }

    // D3Q19 MRT (d'Humieres 2002)
    std::array<double, 19> m{};
    std::array<double, 19> meq{};
    for (int r = 0; r < 19; ++r) {
        double sm = 0.0, seq = 0.0;
        const auto rr = static_cast<std::size_t>(r);
        for (int c = 0; c < 19; ++c) {
            const auto cc = static_cast<std::size_t>(c);
            const double mij = D3Q19_MRT_M[rr][cc];
            sm += mij * f[c];
            seq += mij * feq[cc];
        }
        m[rr] = sm;
        meq[rr] = seq;
    }
    const std::array<double, 19> s_rates = {
        0.0, 1.19, 1.4, 0.0, 1.2, 0.0, 1.2, 0.0, 1.2,
        omega_p, 1.4, omega_p, 1.4, omega_p, omega_p, omega_p,
        1.2, 1.2, 1.2
    };
    for (int r = 0; r < 19; ++r) {
        const auto rr = static_cast<std::size_t>(r);
        m[rr] -= s_rates[rr] * (m[rr] - meq[rr]);
    }
    for (int c = 0; c < 19; ++c) {
        const auto cc = static_cast<std::size_t>(c);
        double val = 0.0;
        for (int r = 0; r < 19; ++r) {
            const auto rr = static_cast<std::size_t>(r);
            val += D3Q19_MRT_M[rr][cc] * (m[rr] / D3Q19_MRT_DIAG_MMT[rr]);
        }
        f[c] = val + gforce[cc];
    }
}

}  // namespace lbm
