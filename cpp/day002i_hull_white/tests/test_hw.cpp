#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hw/bermudan.hpp"
#include "hw/calibration.hpp"
#include "hw/curve.hpp"
#include "hw/hull_white.hpp"
#include "hw/montecarlo.hpp"
#include "hw/trinomial_tree.hpp"

#include <cmath>
#include <vector>

namespace {

hw::NelsonSiegelCurve default_curve() {
    return hw::NelsonSiegelCurve{0.045, -0.02, 0.015, 2.0, 0.0};
}

hw::HwParams default_params() {
    return hw::HwParams{0.08, 0.012};
}

}  // namespace

TEST_CASE("Trinomial tree reproduces initial curve and moments to 1e-12") {
    auto curve = default_curve();
    auto params = default_params();
    hw::TrinomialTree tree(curve, params, 12, 5);

    CHECK(tree.max_curve_error() < 1e-12);

    double dt = 1.0 / 12.0;
    double v = params.sigma * params.sigma * dt;
    double dx = std::sqrt(3.0 * v);

    for (int j = -3; j <= 3; ++j) {
        auto b = tree.branching(j);
        double p_sum = b.pu + b.pm + b.pd;
        CHECK(std::abs(p_sum - 1.0) < 1e-12);
        CHECK(b.pu >= 0.0);
        CHECK(b.pm >= 0.0);
        CHECK(b.pd >= 0.0);

        double x = static_cast<double>(j) * dx;
        double mean_target = x - params.a * x * dt;
        double mean_actual = b.pu * (b.k_mid + 1) * dx + b.pm * b.k_mid * dx + b.pd * (b.k_mid - 1) * dx;
        CHECK(std::abs(mean_actual - mean_target) < 1e-12);

        double m2_target = v + mean_target * mean_target;
        double x_u = (b.k_mid + 1) * dx;
        double x_m = b.k_mid * dx;
        double x_d = (b.k_mid - 1) * dx;
        double m2_actual = b.pu * x_u * x_u + b.pm * x_m * x_m + b.pd * x_d * x_d;
        CHECK(std::abs(m2_actual - m2_target) < 1e-12);
    }
}

TEST_CASE("ZBP, ZBC, put-call parity and 5y cap match analytic formulas") {
    auto curve = default_curve();
    auto params = default_params();
    double t_opt = 1.0;
    double s_bond = 5.0;
    double k_fwd = curve.discount(s_bond) / curve.discount(t_opt);

    double zbp_ana = hw::hw_zbp(curve, params, t_opt, s_bond, k_fwd);
    double zbc_ana = hw::hw_zbc(curve, params, t_opt, s_bond, k_fwd);
    double parity = zbc_ana - zbp_ana - (curve.discount(s_bond) - k_fwd * curve.discount(t_opt));
    CHECK(std::abs(parity) < 1e-12);

    hw::TrinomialTree tree(curve, params, 96, 5);
    double zbp_tree = tree.price_zbp(t_opt, s_bond, k_fwd);
    double zbc_tree = tree.price_zbc(t_opt, s_bond, k_fwd);
    CHECK(std::abs(zbp_tree - zbp_ana) < 2e-5);
    CHECK(std::abs(zbc_tree - zbc_ana) < 2e-5);

    double cap_ana = hw::hw_cap_5y(curve, params, 0.045);
    double cap_tree = tree.price_cap_5y(0.045);
    CHECK(std::abs(cap_tree - cap_ana) < 5e-5);
}

TEST_CASE("Jamshidian swaption matches tree and Monte Carlo") {
    auto curve = default_curve();
    auto params = default_params();
    hw::TrinomialTree tree(curve, params, 48, 10);

    for (double t0 : {1.0, 5.0}) {
        double p_jam = hw::hw_payer_swaption_jamshidian(curve, params, t0, 10.0, 0.045);
        double p_tree = tree.price_payer_swaption(t0, 10.0, 0.045);
        CHECK(std::abs(p_tree - p_jam) / p_jam < 5e-3);
    }

    auto mc = hw::run_mc_swaptions(
        curve, params, {1.0, 5.0}, 10.0, 0.045, {1.0, 5.0}, 20'000, 52, 42ULL, 2);
    for (const auto& mt : mc.martingale_tests) {
        CHECK(mt.abs_error < 4.0 * mt.std_error + 1e-4);
    }
    for (const auto& sw : mc.swaptions) {
        double p_jam = hw::hw_payer_swaption_jamshidian(curve, params, sw.expiry_years, 10.0, 0.045);
        CHECK(std::abs(sw.price - p_jam) < 4.0 * sw.std_error + 5e-4);
    }
}

TEST_CASE("Noise-free calibration recovers true parameters and noisy calibration is thread-invariant") {
    auto curve = default_curve();
    auto params = default_params();
    std::vector<double> expiries{1.0, 2.0, 3.0, 5.0, 7.0, 9.0};
    auto quotes = hw::generate_coterminal_quotes(curve, params, expiries, 10.0, 0.045);

    auto cal1 = hw::calibrate_swaptions(curve, quotes, hw::HwParams{0.2, 0.02});
    CHECK(std::abs(cal1.fitted.a - params.a) < 1e-5);
    CHECK(std::abs(cal1.fitted.sigma - params.sigma) < 1e-5);
    CHECK(cal1.price_rmse < 1e-9);

    auto s1 = hw::run_noisy_calibrations(curve, quotes, params, 0.005, 8, 123ULL, 1);
    auto s4 = hw::run_noisy_calibrations(curve, quotes, params, 0.005, 8, 123ULL, 4);
    CHECK(std::abs(s1.mean_a - s4.mean_a) < 1e-12);
    CHECK(std::abs(s1.mean_sigma - s4.mean_sigma) < 1e-12);
}

TEST_CASE("Bermudan swaption satisfies European bounds and DV01 matches Jamshidian") {
    auto curve = default_curve();
    auto params = default_params();

    double max_euro = 0.0;
    double sum_euro = 0.0;
    for (int yr = 1; yr <= 9; ++yr) {
        double pe = hw::hw_payer_swaption_jamshidian(
            curve, params, static_cast<double>(yr), 10.0, 0.045);
        max_euro = std::max(max_euro, pe);
        sum_euro += pe;
    }

    hw::TrinomialTree tree(curve, params, 48, 10);
    double berm_tree = tree.price_bermudan_payer(10, 0.045);
    CHECK(berm_tree >= max_euro);
    CHECK(berm_tree <= sum_euro);

    auto lsm = hw::price_bermudan_lsm(curve, params, 1, 9, 10, 0.045, 20'000, 12, 99ULL, 2);
    CHECK(std::abs(lsm.price - berm_tree) / berm_tree < 0.03);

    auto dv_jam = hw::dv01_euro_jamshidian(curve, params, 5.0, 10.0, 0.045);
    auto dv_tree = hw::dv01_euro_tree(curve, params, 5.0, 10.0, 0.045, 96);
    CHECK(std::abs(dv_tree.dv01_per_1bp - dv_jam.dv01_per_1bp) / std::abs(dv_jam.dv01_per_1bp) < 0.04);
}
