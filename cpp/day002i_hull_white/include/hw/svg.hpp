#pragma once

#include <string>
#include <vector>

namespace hw::svg {

void write_tree_convergence_svg(const std::string& path,
                                const std::vector<int>& m_values,
                                const std::vector<double>& zbp_errors,
                                const std::vector<double>& zbc_errors,
                                const std::vector<double>& cap_errors,
                                double slope_zbp,
                                double slope_zbc);

void write_curve_and_swaptions_svg(const std::string& path,
                                   const std::vector<double>& maturities,
                                   const std::vector<double>& zero_rates,
                                   const std::vector<double>& fwd_rates,
                                   const std::vector<double>& expiries,
                                   const std::vector<double>& jamshidian_prices,
                                   const std::vector<double>& tree_prices,
                                   const std::vector<double>& mc_prices);

void write_calib_scatter_svg(const std::string& path,
                             double true_a,
                             double true_sigma,
                             const std::vector<double>& a_samples,
                             const std::vector<double>& sigma_samples,
                             double corr_a_sigma,
                             double cond_number);

void write_bermudan_bounds_dv01_svg(const std::string& path,
                                    const std::vector<int>& m_values,
                                    const std::vector<double>& bermudan_tree,
                                    double bermudan_rich,
                                    double lsm_price,
                                    double lsm_se,
                                    double max_euro,
                                    double sum_euro,
                                    double dv01_euro,
                                    double dv01_berm);

}  // namespace hw::svg
