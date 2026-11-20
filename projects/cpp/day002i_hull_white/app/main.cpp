#include "hw/bermudan.hpp"
#include "hw/calibration.hpp"
#include "hw/curve.hpp"
#include "hw/hull_white.hpp"
#include "hw/montecarlo.hpp"
#include "hw/svg.hpp"
#include "hw/trinomial_tree.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

namespace {

double log_log_slope(const std::vector<int>& ms, const std::vector<double>& errs) {
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    double n = static_cast<double>(ms.size());
    for (std::size_t i = 0; i < ms.size(); ++i) {
        double x = std::log(static_cast<double>(ms[i]));
        double y = std::log(std::max(1e-15, errs[i]));
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

}  // namespace

int main() {
    std::filesystem::create_directories("outputs/day002i");
    std::filesystem::create_directories("outputs/figures");

    const hw::NelsonSiegelCurve curve{0.045, -0.02, 0.015, 2.0, 0.0};
    const hw::HwParams params{0.08, 0.012};
    unsigned int n_threads = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));

    // --- Bagian 1: Trinomial tree, Arrow-Debreu, ZBP/ZBC, Cap ---
    const std::vector<int> m_values{12, 24, 48, 96, 192};
    const double t_opt = 1.0;
    const double s_bond = 5.0;
    const double k_fwd = curve.discount(s_bond) / curve.discount(t_opt);
    const std::vector<double> strikes{0.90 * k_fwd, k_fwd, 1.10 * k_fwd};
    const double cap_strike = 0.045;

    std::vector<double> zbp_ana(3), zbc_ana(3);
    for (std::size_t i = 0; i < 3; ++i) {
        zbp_ana[i] = hw::hw_zbp(curve, params, t_opt, s_bond, strikes[i]);
        zbc_ana[i] = hw::hw_zbc(curve, params, t_opt, s_bond, strikes[i]);
    }
    const double cap_ana = hw::hw_cap_5y(curve, params, cap_strike);

    struct TreeRow {
        int m{0};
        double max_curve_err{0.0};
        std::vector<double> zbp{};
        std::vector<double> zbc{};
        double cap{0.0};
        double zbp_fwd_err{0.0};
        double zbc_fwd_err{0.0};
        double cap_err{0.0};
    };
    std::vector<TreeRow> tree_rows;
    std::vector<double> zbp_fwd_errs, zbc_fwd_errs, cap_errs;

    for (int m : m_values) {
        hw::TrinomialTree tree(curve, params, m, 5);
        TreeRow row;
        row.m = m;
        row.max_curve_err = tree.max_curve_error();
        for (std::size_t i = 0; i < 3; ++i) {
            row.zbp.push_back(tree.price_zbp(t_opt, s_bond, strikes[i]));
            row.zbc.push_back(tree.price_zbc(t_opt, s_bond, strikes[i]));
        }
        row.cap = tree.price_cap_5y(cap_strike);
        row.zbp_fwd_err = std::abs(row.zbp[1] - zbp_ana[1]);
        row.zbc_fwd_err = std::abs(row.zbc[1] - zbc_ana[1]);
        row.cap_err = std::abs(row.cap - cap_ana);

        zbp_fwd_errs.push_back(row.zbp_fwd_err);
        zbc_fwd_errs.push_back(row.zbc_fwd_err);
        cap_errs.push_back(row.cap_err);
        tree_rows.push_back(std::move(row));
    }

    double slope_zbp = log_log_slope(m_values, zbp_fwd_errs);
    double slope_zbc = log_log_slope(m_values, zbc_fwd_errs);
    double slope_cap = log_log_slope(m_values, cap_errs);

    // First-order Richardson extrapolation from m=96 and m=192: 2*V(2m) - V(m)
    double zbp_rich = 2.0 * tree_rows[4].zbp[1] - tree_rows[3].zbp[1];
    double zbc_rich = 2.0 * tree_rows[4].zbc[1] - tree_rows[3].zbc[1];
    double cap_rich = 2.0 * tree_rows[4].cap - tree_rows[3].cap;

    // --- Bagian 2: Jamshidian vs Tree vs Monte Carlo ---
    const std::vector<double> expiries{1.0, 2.0, 3.0, 5.0, 7.0, 9.0};
    const double tenor_end = 10.0;
    const double swap_rate = 0.045;

    std::vector<double> jam_prices(expiries.size());
    std::vector<double> tree10_prices(expiries.size());
    hw::TrinomialTree tree10_192(curve, params, 192, 10);
    for (std::size_t i = 0; i < expiries.size(); ++i) {
        jam_prices[i] = hw::hw_payer_swaption_jamshidian(curve, params, expiries[i], tenor_end, swap_rate);
        tree10_prices[i] = tree10_192.price_payer_swaption(expiries[i], tenor_end, swap_rate);
    }

    auto mc_suite = hw::run_mc_swaptions(
        curve, params, expiries, tenor_end, swap_rate, {1.0, 5.0, 10.0}, 200'000, 52, 20261120ULL, n_threads);

    std::vector<double> mc_prices;
    for (const auto& sw : mc_suite.swaptions) {
        mc_prices.push_back(sw.price);
    }

    // --- Bagian 3: Calibration ---
    auto quotes = hw::generate_coterminal_quotes(curve, params, expiries, tenor_end, swap_rate);
    auto cal_start1 = hw::calibrate_swaptions(curve, quotes, hw::HwParams{0.2, 0.02});
    auto cal_start2 = hw::calibrate_swaptions(curve, quotes, hw::HwParams{0.02, 0.005});
    auto noisy_cal = hw::run_noisy_calibrations(
        curve, quotes, params, 0.005, 100, 987654321ULL, n_threads);

    // --- Bagian 4: Bermudan Swaption & LSM ---
    std::vector<double> berm_tree_prices;
    for (int m : m_values) {
        if (m == 192) {
            berm_tree_prices.push_back(tree10_192.price_bermudan_payer(10, swap_rate));
        } else {
            hw::TrinomialTree tm(curve, params, m, 10);
            berm_tree_prices.push_back(tm.price_bermudan_payer(10, swap_rate));
        }
    }
    double berm_rich = 2.0 * berm_tree_prices[4] - berm_tree_prices[3];
    double max_euro = 0.0;
    double sum_euro = 0.0;
    for (int yr = 1; yr <= 9; ++yr) {
        double pe = hw::hw_payer_swaption_jamshidian(curve, params, static_cast<double>(yr), 10.0, swap_rate);
        max_euro = std::max(max_euro, pe);
        sum_euro += pe;
    }
    double early_exercise_premium = berm_rich - max_euro;

    auto berm_lsm = hw::price_bermudan_lsm(
        curve, params, 1, 9, 10, swap_rate, 100'000, 24, 7777777ULL, n_threads);

    // --- Bagian 5: Parallel +-1bp DV01 ---
    auto dv01_5x5_jam = hw::dv01_euro_jamshidian(curve, params, 5.0, 10.0, swap_rate, 1.0);
    auto dv01_5x5_tree = hw::dv01_euro_tree(curve, params, 5.0, 10.0, swap_rate, 192, 1.0);
    auto dv01_berm = hw::dv01_bermudan_tree(curve, params, 1, 9, 10, swap_rate, 192, 1.0);

    // --- SVG Figures ---
    hw::svg::write_tree_convergence_svg(
        "outputs/figures/day002i_tree-convergence.svg",
        m_values, zbp_fwd_errs, zbc_fwd_errs, cap_errs, slope_zbp, slope_zbc);

    std::vector<double> mats, zrs, fwds;
    for (int i = 1; i <= 40; ++i) {
        double t = 0.25 * static_cast<double>(i);
        mats.push_back(t);
        zrs.push_back(curve.zero_rate(t));
        fwds.push_back(curve.inst_forward(t));
    }
    hw::svg::write_curve_and_swaptions_svg(
        "outputs/figures/day002i_curve-and-swaptions.svg",
        mats, zrs, fwds, expiries, jam_prices, tree10_prices, mc_prices);

    std::vector<double> a_samps, sig_samps;
    for (const auto& s : noisy_cal.samples) {
        a_samps.push_back(s.a);
        sig_samps.push_back(s.sigma);
    }
    hw::svg::write_calib_scatter_svg(
        "outputs/figures/day002i_calib-scatter.svg",
        params.a, params.sigma, a_samps, sig_samps,
        noisy_cal.corr_a_sigma, noisy_cal.jacobian_cond_number);

    hw::svg::write_bermudan_bounds_dv01_svg(
        "outputs/figures/day002i_bermudan-bounds-dv01.svg",
        m_values, berm_tree_prices, berm_rich, berm_lsm.price, berm_lsm.std_error,
        max_euro, sum_euro, dv01_5x5_jam.dv01_per_1bp, dv01_berm.dv01_per_1bp);

    // --- Write results.json (17 significant digits) ---
    std::ofstream js("outputs/day002i/results.json");
    js << std::setprecision(17);
    js << "{\n";
    js << "  \"model\": {\"a\": " << params.a << ", \"sigma\": " << params.sigma << "},\n";
    js << "  \"bagian1_tree_and_options\": {\n";
    js << "    \"k_fwd\": " << k_fwd << ",\n";
    js << "    \"analytic\": {\"zbp\": [" << zbp_ana[0] << ", " << zbp_ana[1] << ", " << zbp_ana[2] << "], "
       << "\"zbc\": [" << zbc_ana[0] << ", " << zbc_ana[1] << ", " << zbc_ana[2] << "], "
       << "\"cap_5y\": " << cap_ana << "},\n";
    js << "    \"richardson_192_96\": {\"zbp_fwd\": " << zbp_rich << ", \"zbc_fwd\": " << zbc_rich
       << ", \"cap_5y\": " << cap_rich << "},\n";
    js << "    \"loglog_slopes\": {\"zbp_fwd\": " << slope_zbp << ", \"zbc_fwd\": " << slope_zbc
       << ", \"cap_5y\": " << slope_cap << "},\n";
    js << "    \"convergence\": [\n";
    for (std::size_t i = 0; i < tree_rows.size(); ++i) {
        const auto& r = tree_rows[i];
        js << "      {\"m\": " << r.m << ", \"max_curve_err\": " << r.max_curve_err
           << ", \"zbp_fwd\": " << r.zbp[1] << ", \"zbc_fwd\": " << r.zbc[1]
           << ", \"cap_5y\": " << r.cap << ", \"zbp_fwd_err\": " << r.zbp_fwd_err
           << ", \"zbc_fwd_err\": " << r.zbc_fwd_err << ", \"cap_err\": " << r.cap_err << "}"
           << (i + 1 < tree_rows.size() ? "," : "") << "\n";
    }
    js << "    ]\n  },\n";

    js << "  \"bagian2_swaptions_and_mc\": {\n";
    js << "    \"martingale_tests\": [\n";
    for (std::size_t i = 0; i < mc_suite.martingale_tests.size(); ++i) {
        const auto& mt = mc_suite.martingale_tests[i];
        js << "      {\"T\": " << mt.maturity_years << ", \"market_df\": " << mt.market_df
           << ", \"mc_df\": " << mt.mc_df << ", \"se\": " << mt.std_error
           << ", \"abs_err\": " << mt.abs_error << "}"
           << (i + 1 < mc_suite.martingale_tests.size() ? "," : "") << "\n";
    }
    js << "    ],\n    \"coterminal_10y\": [\n";
    for (std::size_t i = 0; i < expiries.size(); ++i) {
        const auto& sw = mc_suite.swaptions[i];
        js << "      {\"T0\": " << expiries[i] << ", \"jamshidian\": " << jam_prices[i]
           << ", \"tree_m192\": " << tree10_prices[i]
           << ", \"mc_cv\": " << sw.price << ", \"mc_cv_se\": " << sw.std_error
           << ", \"mc_raw\": " << sw.raw_price << ", \"mc_raw_se\": " << sw.raw_std_error << "}"
           << (i + 1 < expiries.size() ? "," : "") << "\n";
    }
    js << "    ]\n  },\n";

    js << "  \"bagian3_calibration\": {\n";
    js << "    \"start_high\": {\"a\": " << cal_start1.fitted.a << ", \"sigma\": " << cal_start1.fitted.sigma
       << ", \"rmse\": " << cal_start1.price_rmse << ", \"iters\": " << cal_start1.iterations << "},\n";
    js << "    \"start_low\": {\"a\": " << cal_start2.fitted.a << ", \"sigma\": " << cal_start2.fitted.sigma
       << ", \"rmse\": " << cal_start2.price_rmse << ", \"iters\": " << cal_start2.iterations << "},\n";
    js << "    \"noisy_100\": {\"mean_a\": " << noisy_cal.mean_a << ", \"sd_a\": " << noisy_cal.sd_a
       << ", \"mean_sigma\": " << noisy_cal.mean_sigma << ", \"sd_sigma\": " << noisy_cal.sd_sigma
       << ", \"corr_a_sigma\": " << noisy_cal.corr_a_sigma
       << ", \"jacobian_cond\": " << noisy_cal.jacobian_cond_number << "}\n";
    js << "  },\n";

    js << "  \"bagian4_bermudan\": {\n";
    js << "    \"max_euro\": " << max_euro << ", \"sum_euro\": " << sum_euro << ",\n";
    js << "    \"tree_by_m\": [";
    for (std::size_t i = 0; i < m_values.size(); ++i) {
        js << "{\"m\": " << m_values[i] << ", \"price\": " << berm_tree_prices[i] << "}"
           << (i + 1 < m_values.size() ? ", " : "");
    }
    js << "],\n";
    js << "    \"tree_richardson\": " << berm_rich << ",\n";
    js << "    \"early_exercise_premium\": " << early_exercise_premium << ",\n";
    js << "    \"lsm\": {\"price\": " << berm_lsm.price << ", \"se\": " << berm_lsm.std_error
       << ", \"paths_per_pass\": " << berm_lsm.train_paths << "}\n";
    js << "  },\n";

    js << "  \"bagian5_dv01\": {\n";
    js << "    \"euro_5x5_jamshidian\": " << dv01_5x5_jam.dv01_per_1bp << ",\n";
    js << "    \"euro_5x5_tree_m192\": " << dv01_5x5_tree.dv01_per_1bp << ",\n";
    js << "    \"euro_5x5_rel_diff\": "
       << std::abs(dv01_5x5_tree.dv01_per_1bp - dv01_5x5_jam.dv01_per_1bp) / std::abs(dv01_5x5_jam.dv01_per_1bp)
       << ",\n";
    js << "    \"bermudan_10y_tree_m192\": " << dv01_berm.dv01_per_1bp << "\n";
    js << "  }\n}\n";

    std::cout << "day002i complete. Results written to outputs/day002i/results.json\n";
    std::cout << "Max tree curve err (m=192): " << tree_rows.back().max_curve_err << "\n";
    std::cout << "ZBP/ZBC slopes: " << slope_zbp << " / " << slope_zbc << "\n";
    std::cout << "Bermudan Richardson: " << berm_rich << " | LSM: " << berm_lsm.price
              << " +- " << berm_lsm.std_error << "\n";
    std::cout << "DV01 5x5 Jam=" << dv01_5x5_jam.dv01_per_1bp
              << " Tree=" << dv01_5x5_tree.dv01_per_1bp
              << " | Berm=" << dv01_berm.dv01_per_1bp << "\n";
    return 0;
}
