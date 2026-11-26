#include "lbm/stages.hpp"
#include "lbm/analysis.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/grid.hpp"
#include "lbm/lattice.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

namespace lbm {

bool run_selftest() {
    std::printf("[selftest] Starting LBM C++20 self-test suite...\n");
    int failures = 0;

    // 1. D2Q9 & D3Q19 weights & 4th-order isotropy
    const auto r2 = check_d2q9_isotropy();
    const auto r3 = check_d3q19_isotropy();
    const double max_iso = std::max({r2.weight_sum_err, r2.order2_err, r2.order4_err,
                                     r3.weight_sum_err, r3.order2_err, r3.order4_err});
    if (max_iso > 1e-13) {
        std::printf("  [FAIL] 1. Lattice weights & isotropy max_err=%.3e\n", max_iso);
        ++failures;
    } else {
        std::printf("  [OK] 1. D2Q9 & D3Q19 weights & 4th-order isotropy: max_err=%.2e\n", max_iso);
    }

    // 2. Hand-computed 1-cell equilibrium & moment recovery
    const double rho0 = 1.2, ux0 = 0.1, uy0 = -0.05;
    const auto feq = compute_feq<D2Q9>(rho0, ux0, uy0);
    // Exact u2 = 0.0125 => f_0 = (4/9)*1.2*(1 - 1.5*0.0125) = 0.5233333333333333...
    const double hand_f0 = (4.0 / 9.0) * 1.2 * (1.0 - 1.5 * 0.0125);
    const double hand_f1 = (1.0 / 9.0) * 1.2 * (1.0 + 0.3 + 0.045 - 0.01875);
    const auto mc = compute_macro_2d(std::span<const double, 9>(feq));
    const double cell_err = std::max({std::abs(feq[0] - hand_f0), std::abs(feq[1] - hand_f1),
                                      std::abs(mc.rho - rho0), std::abs(mc.ux - ux0),
                                      std::abs(mc.uy - uy0)});
    if (cell_err > 1e-14) {
        std::printf("  [FAIL] 2. Single-cell equilibrium & moments err=%.3e\n", cell_err);
        ++failures;
    } else {
        std::printf("  [OK] 2. Hand-computed 1-cell equilibrium & moments: max_err=%.2e\n", cell_err);
    }

    // 3. Periodic 100-step mass & momentum conservation
    Grid2D g(16, 16, Scheme::PULL_2GRID, MemLayout::SOA);
    for (int y = 0; y < 16; ++y) {
        for (int x = 0; x < 16; ++x) {
            const double ang = 2.0 * 3.14159265358979323846 * x / 16.0;
            g.set_cell_array(x, y, compute_feq<D2Q9>(1.0 + 0.01 * std::cos(ang), 0.02 * std::sin(ang), 0.01));
        }
    }
    double m0, px0, py0, ek0, rmin0, rmax0, umax0;
    bool nan0 = false;
    g.totals(m0, px0, py0, ek0, rmin0, rmax0, umax0, nan0);
    CollisionParams cp{CollisionOp::BGK, 0.8, 3.0 / 16.0, 0.0, 0.0, 0.0};
    for (int t = 0; t < 100; ++t) step_periodic_2d(g, cp, 1);
    double m1, px1, py1, ek1, rmin1, rmax1, umax1;
    bool nan1 = false;
    g.totals(m1, px1, py1, ek1, rmin1, rmax1, umax1, nan1);
    const double dm = std::abs(m1 - m0) / m0;
    const double dp = std::hypot(px1 - px0, py1 - py0) / m0;
    if (nan1 || dm > 1e-13 || dp > 1e-13) {
        std::printf("  [FAIL] 3. 100-step periodic conservation dm=%.3e dp=%.3e\n", dm, dp);
        ++failures;
    } else {
        std::printf("  [OK] 3. 100-step periodic mass/momentum drift: dm=%.2e, dp=%.2e\n", dm, dp);
    }

    // 4. Shear wave N=16 vs nu = (tau - 0.5) / 3 (tau = 0.8 -> nu = 0.1)
    Grid2D sw(16, 16, Scheme::PULL_2GRID, MemLayout::SOA);
    const double tau = 0.8;
    const double k = 2.0 * 3.14159265358979323846 / 16.0;
    const double u0 = 0.01;
    for (int y = 0; y < 16; ++y) {
        const double ux = u0 * std::sin(k * (y + 0.5));
        const double dux_dy = u0 * k * std::cos(k * (y + 0.5));
        const auto f_init = compute_chapman_enskog_init_2d(1.0, ux, 0.0, 0.0, dux_dy, 0.0, 0.0, tau);
        for (int x = 0; x < 16; ++x) sw.set_cell_array(x, y, f_init);
    }
    CollisionParams cp_sw{CollisionOp::BGK, tau, 3.0 / 16.0, 0.0, 0.0, 0.0};
    std::vector<double> ts, ln_amp;
    for (int t = 1; t <= 80; ++t) {
        step_periodic_2d(sw, cp_sw, 1);
        if (t >= 10 && t % 5 == 0) {
            double num = 0.0;
            for (int y = 0; y < 16; ++y) {
                num += sw.cell_macro(0, y).ux * std::sin(k * (y + 0.5));
            }
            const double amp = std::abs(2.0 * num / 16.0);
            ts.push_back(static_cast<double>(t));
            ln_amp.push_back(std::log(amp));
        }
    }
    const auto fit = ols_fit(ts, ln_amp);
    const double nu_meas = -fit.slope / (k * k);
    const double nu_th = (tau - 0.5) / 3.0;
    const double rel_nu = std::abs(nu_meas - nu_th) / nu_th;
    // At coarse N=16 (k = pi/8), O(k^2) spatial truncation is ~0.83%; tolerance 1.5%
    if (rel_nu > 0.015) {
        std::printf("  [FAIL] 4. N=16 shear wave nu=%.6f vs %.6f (rel=%.3e)\n", nu_meas, nu_th, rel_nu);
        ++failures;
    } else {
        std::printf("  [OK] 4. N=16 shear wave nu=%.6f vs %.6f (rel=%.3e <= 0.015)\n",
                    nu_meas, nu_th, rel_nu);
    }

    std::printf("[selftest] Finished with %d failure(s).\n", failures);
    return failures == 0;
}

}  // namespace lbm
