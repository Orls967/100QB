#include "asian/options.hpp"
#include "asian/svg.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;

namespace {

fs::path find_repo_root() {
  fs::path dir = fs::current_path();
  while (true) {
    if (fs::exists(dir / "src" / "qb100")) {
      return dir;
    }
    if (dir == dir.parent_path()) {
      throw std::runtime_error("Repository root not found");
    }
    dir = dir.parent_path();
  }
}

double elapsed_ms(std::chrono::steady_clock::time_point start) {
  const auto end = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count();
}

} // namespace

int main() {
  const asian::OptionParams params{
      .s0 = 100.0,
      .k = 100.0,
      .r = 0.05,
      .sigma = 0.20,
      .t = 1.0,
      .n_steps = 252,
  };
  constexpr std::size_t n_paths = 200000;
  constexpr std::uint64_t seed = 42;

  const auto t_total_start = std::chrono::steady_clock::now();

  auto t0 = std::chrono::steady_clock::now();
  const double bs_euro = asian::bs_call(params);
  const auto mc_euro = asian::price_european_mc(params, n_paths, seed);
  const double t_sanity = elapsed_ms(t0);
  const double diff_euro = mc_euro.price - bs_euro;
  const double z_euro = std::abs(diff_euro) / mc_euro.se;

  t0 = std::chrono::steady_clock::now();
  const auto res_naive = asian::price_asian_naive(params, n_paths, seed);
  const double t_naive = elapsed_ms(t0);

  t0 = std::chrono::steady_clock::now();
  const auto res_anti = asian::price_asian_antithetic(params, n_paths, seed);
  const double t_anti = elapsed_ms(t0);

  t0 = std::chrono::steady_clock::now();
  const auto res_cv = asian::price_asian_control_variate(params, n_paths, seed);
  const double t_cv = elapsed_ms(t0);

  t0 = std::chrono::steady_clock::now();
  const auto delta_pw = asian::delta_asian_pathwise(params, n_paths, seed);
  const double t_pw = elapsed_ms(t0);

  t0 = std::chrono::steady_clock::now();
  const auto delta_bmp = asian::delta_asian_bump(params, n_paths, 0.01, seed);
  const double t_bmp = elapsed_ms(t0);

  const double core_ms = elapsed_ms(t_total_start);

  const double delta_diff = delta_pw.delta - delta_bmp.delta;
  const double delta_z = std::abs(delta_diff) / std::hypot(delta_pw.se, delta_bmp.se);

  const double var_naive = res_naive.se * res_naive.se;
  const std::vector<asian::EstimatorSummary> estimators{
      {"Naive MC (200k)", res_naive.price, res_naive.se, res_naive.ci_low, res_naive.ci_high, 1.0},
      {"Antithetic (100k pasang)", res_anti.price, res_anti.se, res_anti.ci_low, res_anti.ci_high,
       var_naive / (res_anti.se * res_anti.se)},
      {"Control Variate (Geom. Asia)", res_cv.price, res_cv.se, res_cv.ci_low, res_cv.ci_high,
       var_naive / (res_cv.se * res_cv.se)},
  };

  t0 = std::chrono::steady_clock::now();
  const std::vector<std::size_t> n_grid{1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000};
  std::vector<asian::ConvergencePoint> conv_points;
  std::vector<double> se_n;
  std::vector<double> se_a;
  std::vector<double> se_c;
  for (const auto n : n_grid) {
    const auto rn = asian::price_asian_naive(params, n, seed);
    const auto ra = asian::price_asian_antithetic(params, n, seed);
    const auto rc = asian::price_asian_control_variate(params, n, seed);
    se_n.push_back(rn.se);
    se_a.push_back(ra.se);
    se_c.push_back(rc.se);
    conv_points.push_back({n, rn.se, ra.se, rc.se});
  }
  const double t_conv = elapsed_ms(t0);

  const asian::Slopes slopes{
      .naive = asian::fit_loglog_slope(n_grid, se_n),
      .antithetic = asian::fit_loglog_slope(n_grid, se_a),
      .control_variate = asian::fit_loglog_slope(n_grid, se_c),
  };

  const fs::path root = find_repo_root();
  const fs::path out_dir = root / "outputs" / "day001g";
  const fs::path fig_dir = root / "outputs" / "figures";
  fs::create_directories(out_dir);
  fs::create_directories(fig_dir);

  std::ostringstream js;
  js << std::setprecision(15);
  js << "{\n"
     << "  \"params\": {\n"
     << "    \"s0\": " << params.s0 << ", \"k\": " << params.k << ", \"r\": " << params.r
     << ", \"sigma\": " << params.sigma << ", \"t\": " << params.t
     << ", \"n_steps\": " << params.n_steps << "\n  },\n"
     << "  \"n_paths\": " << n_paths << ",\n"
     << "  \"seed\": " << seed << ",\n"
     << "  \"sanity_check\": {\n"
     << "    \"bs_call\": " << bs_euro << ",\n"
     << "    \"mc_call\": {\"price\": " << mc_euro.price << ", \"se\": " << mc_euro.se
     << ", \"ci_low\": " << mc_euro.ci_low << ", \"ci_high\": " << mc_euro.ci_high << "},\n"
     << "    \"diff\": " << diff_euro << ",\n"
     << "    \"z_score\": " << z_euro << "\n  },\n"
     << "  \"estimators\": [\n";
  for (std::size_t i = 0; i < estimators.size(); ++i) {
    const auto &e = estimators[i];
    js << "    {\"name\": \"" << e.name << "\", \"price\": " << e.price << ", \"se\": " << e.se
       << ", \"ci_low\": " << e.ci_low << ", \"ci_high\": " << e.ci_high << ", \"vrf\": " << e.vrf;
    if (i == 2) {
      js << ", \"beta\": " << res_cv.beta << ", \"y_exact\": " << res_cv.y_exact;
    }
    js << "}" << (i + 1 < estimators.size() ? "," : "") << "\n";
  }
  js << "  ],\n"
     << "  \"delta\": {\n"
     << "    \"pathwise\": {\"delta\": " << delta_pw.delta << ", \"se\": " << delta_pw.se
     << ", \"ci_low\": " << delta_pw.ci_low << ", \"ci_high\": " << delta_pw.ci_high << "},\n"
     << "    \"bump_1pct\": {\"delta\": " << delta_bmp.delta << ", \"se\": " << delta_bmp.se
     << ", \"ci_low\": " << delta_bmp.ci_low << ", \"ci_high\": " << delta_bmp.ci_high << "},\n"
     << "    \"diff\": " << delta_diff << ",\n"
     << "    \"z_score\": " << delta_z << "\n  },\n"
     << "  \"convergence\": [\n";
  for (std::size_t i = 0; i < conv_points.size(); ++i) {
    const auto &pt = conv_points[i];
    js << "    {\"n\": " << pt.n << ", \"se_naive\": " << pt.se_naive
       << ", \"se_antithetic\": " << pt.se_antithetic
       << ", \"se_control_variate\": " << pt.se_control_variate << "}"
       << (i + 1 < conv_points.size() ? "," : "") << "\n";
  }
  js << "  ],\n"
     << "  \"slopes\": {\"naive\": " << slopes.naive << ", \"antithetic\": " << slopes.antithetic
     << ", \"control_variate\": " << slopes.control_variate << "},\n"
     << "  \"timings_ms\": {\n"
     << "    \"sanity_check\": " << t_sanity << ", \"naive\": " << t_naive
     << ", \"antithetic\": " << t_anti << ", \"control_variate\": " << t_cv
     << ", \"delta_pathwise\": " << t_pw << ", \"delta_bump\": " << t_bmp
     << ", \"convergence\": " << t_conv << ", \"total_200k_core\": " << core_ms << "\n  }\n}\n";

  {
    std::ofstream f(out_dir / "results.json");
    f << js.str();
  }
  {
    std::ofstream f(fig_dir / "day001g_asian-convergence-loglog.svg");
    f << asian::render_convergence_svg(conv_points, slopes);
  }
  {
    std::ofstream f(fig_dir / "day001g_asian-estimators-ci.svg");
    f << asian::render_estimators_svg(estimators);
  }

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "BS Eropa: " << bs_euro << " | MC Eropa: " << mc_euro.price << " ± " << mc_euro.se
            << " (" << std::setprecision(2) << z_euro << " SE)\n";
  std::cout << std::setprecision(6) << "Geom. Asian exact: " << res_cv.y_exact
            << " | beta: " << res_cv.beta << "\n";
  for (const auto &e : estimators) {
    std::cout << std::left << std::setw(28) << e.name << ": " << e.price << " ± " << e.se
              << " (CI [" << e.ci_low << ", " << e.ci_high << "]) | VRF = " << std::setprecision(2)
              << e.vrf << "x\n"
              << std::setprecision(6);
  }
  std::cout << "Delta PW: " << delta_pw.delta << " ± " << delta_pw.se
            << " | Bump: " << delta_bmp.delta << " ± " << delta_bmp.se << " | diff = "
            << (delta_diff >= 0 ? "+" : "") << delta_diff << " (" << std::setprecision(2) << delta_z
            << " SE)\n";
  std::cout << std::setprecision(4) << "Slopes -> Naive: " << slopes.naive
            << ", Antithetic: " << slopes.antithetic << ", CV: " << slopes.control_variate << "\n";
  std::cout << std::setprecision(1) << "Total 200k core time: " << core_ms << " ms\n";
  return 0;
}
