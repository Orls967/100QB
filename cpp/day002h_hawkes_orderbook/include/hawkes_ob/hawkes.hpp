#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hob {

struct HawkesParams {
  double mu{0.5};
  double alpha{0.6};
  double beta{1.2};
};

struct HawkesMleFit {
  HawkesParams params;
  double log_lik{0.0};
  double se_mu{0.0};
  double se_alpha{0.0};
  double se_beta{0.0};
  double corr_alpha_beta{0.0};
};

struct HawkesMleMcSummary {
  std::size_t r_reps{0};
  double mean_mu{0.0}, sd_mu{0.0}, bias_mu{0.0}, mean_se_mu{0.0}, cov_mu{0.0};
  double mean_alpha{0.0}, sd_alpha{0.0}, bias_alpha{0.0}, mean_se_alpha{0.0}, cov_alpha{0.0};
  double mean_beta{0.0}, sd_beta{0.0}, bias_beta{0.0}, mean_se_beta{0.0}, cov_beta{0.0};
  double corr_alpha_beta{0.0};
  std::vector<HawkesParams> estimates;
};

std::vector<double> simulate_hawkes_ogata(const HawkesParams& p, double t_end, std::uint64_t seed);

std::vector<double> simulate_hawkes_dassios_zhao(const HawkesParams& p, double t_end,
                                                 std::uint64_t seed);

std::vector<double> inter_event_times(const std::vector<double>& times);

std::vector<double> time_rescaling_residuals(const std::vector<double>& times,
                                             const HawkesParams& p);

double hawkes_log_likelihood(const std::vector<double>& times, double t_end,
                             const HawkesParams& p);

HawkesMleFit fit_hawkes_mle(const std::vector<double>& times, double t_end);

HawkesMleMcSummary run_hawkes_mle_mc(const HawkesParams& true_p, std::size_t r_reps, double t_end,
                                     std::uint64_t master_seed, std::size_t n_threads = 0);

} // namespace hob
