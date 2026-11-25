#define CHAOS_NO_MAIN
#include "../src/chaos.cpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using namespace chaos;

    // Tes 1: RK4 pada osilator harmonik x'' = -x (T = 2.0, rasio galat orde 4 ~ 16)
    auto ho_err = [](double dt) {
        Vec3 u{1.0, 0.0, 0.0};
        int n = static_cast<int>(std::llround(2.0 / dt));
        auto f = [](const Vec3& s, double) -> Vec3 { return {s.y, -s.x, 0.0}; };
        for (int i = 0; i < n; ++i) u = rk4_step(u, 0.0, dt, f);
        return std::hypot(u.x - std::cos(2.0), u.y + std::sin(2.0));
    };
    double ratio_ho = ho_err(0.05) / ho_err(0.025);
    assert(ratio_ho > 15.0 && ratio_ho < 17.0);

    // Tes 2: Kontraksi volume Liouville Lorenz pada T pendek (T = 60)
    auto lshort = compute_lorenz_lyapunov(0.005, 15.0, 60.0, 10);
    double exact_sum = -(10.0 + 1.0 + 8.0 / 3.0);
    assert(std::abs(lshort.sum_lam - exact_sum) < 0.02);
    assert(lshort.lam[0] > 0.7 && lshort.lam[0] < 1.1);

    // Tes 3: Titik tetap stabil peta logistik r = 2.5 -> x* = 1 - 1/r = 0.6, f'(x*) = 2 - r = -0.5 -> lambda = ln(0.5)
    double lam_25 = logistic_lyapunov(2.5, 2000, 50000);
    assert(std::abs(lam_25 - std::log(0.5)) < 1e-6);

    // Tes 4: Siklus stabil periode-2 pada r = 3.2 (lambda < 0) dan titik superstabil R_1 = 1 + sqrt(5)
    double lam_32 = logistic_lyapunov(3.2, 2000, 50000);
    assert(lam_32 < -0.1);
    auto feig = solve_feigenbaum(4);
    assert(std::abs(feig[0].R_n - (1.0 + std::sqrt(5.0))) < 1e-12);

    std::cout << "ALL 4 CHAOS TESTS PASSED (HO_ratio=" << ratio_ho
              << ", Lorenz_sum=" << lshort.sum_lam << ", lam_2.5=" << lam_25 << ")\n";
    return 0;
}
