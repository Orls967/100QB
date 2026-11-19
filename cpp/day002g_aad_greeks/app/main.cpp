#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "aad/dual.hpp"
#include "aad/mc_engine.hpp"
#include "aad/models.hpp"
#include "aad/number.hpp"
#include "aad/rng.hpp"
#include "aad/svg.hpp"

using namespace aad;

namespace {

std::filesystem::path find_repo_root() {
  auto dir = std::filesystem::current_path();
  while (!std::filesystem::is_directory(dir / "src" / "qb100")) {
    if (dir == dir.parent_path()) {
      throw std::runtime_error("Cannot find 100QB repo root (src/qb100)");
    }
    dir = dir.parent_path();
  }
  return dir;
}

double now_sec() noexcept {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

void write_est_json(std::ostream& os, const McEstimate& e) {
  os << "{\"mean\": " << e.mean << ", \"se\": " << e.se << ", \"ci_low\": " << e.ci_low
     << ", \"ci_high\": " << e.ci_high << "}";
}

} // namespace

int main() {
  const auto root = find_repo_root();
  const auto out_dir = root / "outputs" / "day002g";
  const auto fig_dir = root / "outputs" / "figures";
  std::filesystem::create_directories(out_dir);
  std::filesystem::create_directories(fig_dir);

  // Bagian 1: AAD vs Dual vs FD
  double max_aad_dual = 0.0, max_aad_fd = 0.0;
  {
    Xoshiro256StarStar rng(7);
    for (int i = 0; i < 10; ++i) {
      const double xv = 0.5 + rng.next_uniform();
      const double yv = 0.5 + rng.next_uniform();
      auto expr = [](const auto& x, const auto& y) {
        return math::exp(x * y) + math::normal_cdf(x - y) + math::sqrt(x + y);
      };
      Tape tape;
      TapeGuard guard(tape);
      Number x(xv), y(yv);
      Number r = expr(x, y);
      tape.backprop(r.index());
      const auto dx = expr(Dual<double>(xv, 1.0), Dual<double>(yv, 0.0));
      max_aad_dual = std::max(max_aad_dual, std::abs(x.adjoint() - dx.dot));
      const double fd = (expr(xv + 1e-6, yv) - expr(xv - 1e-6, yv)) / 2e-6;
      max_aad_fd = std::max(max_aad_fd, std::abs(x.adjoint() - fd) / std::abs(fd));
    }
  }

  // Bagian 2: BS 5 gradients + Forward-over-Reverse Gamma + timing ratio
  const auto an = bs_analytic_greeks(100.0, 100.0, 0.03, 0.20, 1.0);
  double bs_grad_ratio = 0.0;
  double gamma_for = 0.0;
  {
    constexpr int reps = 50'000;
    const double t0 = now_sec();
    volatile double sink = 0.0;
    for (int i = 0; i < reps; ++i) {
      sink += bs_call(100.0, 100.0, 0.03, 0.20, 1.0);
    }
    const double t_val = now_sec() - t0;

    Tape tape;
    TapeGuard guard(tape);
    const double t1 = now_sec();
    for (int i = 0; i < reps; ++i) {
      tape.clear();
      Number s(100.0), k(100.0), r(0.03), sig(0.20), t(1.0);
      Number c = bs_call(s, k, r, sig, t);
      tape.backprop(c.index());
      sink += s.adjoint();
    }
    const double t_grad = now_sec() - t1;
    bs_grad_ratio = t_grad / std::max(t_val, 1e-9);
    (void)sink;

    tape.clear();
    Number s_base(100.0), k(100.0), r(0.03), sig(0.20), t(1.0);
    Dual<Number> dc = bs_call(Dual<Number>(s_base, Number(1.0)), Dual<Number>(k, Number(0.0)),
                              Dual<Number>(r, Number(0.0)), Dual<Number>(sig, Number(0.0)),
                              Dual<Number>(t, Number(0.0)));
    tape.backprop(dc.dot.index());
    gamma_for = s_base.adjoint();
  }

  // Bagian 3: Monte Carlo 200,000 paths
  const auto gbm = run_gbm_1step_aad(200'000, 42);
  const auto asian = run_asian_aad(200'000, 42);
  const auto digital = run_digital_aad(200'000, 42);
  const auto arith_bask = run_basket_10_arith(200'000, 42, 0, true);
  const auto geom_bask = run_basket_10_geom(200'000, 42);

  svg::write_aad_vs_analytic_svg((fig_dir / "day002g_aad-vs-analytic.svg").string(), gbm,
                                 an.delta, an.vega, an.rho, an.theta);
  svg::write_digital_bias_var_svg((fig_dir / "day002g_digital-bias-variance.svg").string(), digital);
  svg::write_basket_speedup_svg((fig_dir / "day002g_basket-speedup.svg").string(), arith_bask);
  svg::write_tape_growth_svg((fig_dir / "day002g_tape-growth.svg").string(), 10,
                             arith_bask.nodes_per_path, asian.nodes_per_path);

  std::ofstream jf(out_dir / "results.json");
  jf << std::setprecision(17);
  jf << "{\n"
     << "  \"max_aad_vs_dual\": " << max_aad_dual << ",\n"
     << "  \"max_aad_vs_fd_rel\": " << max_aad_fd << ",\n"
     << "  \"bs_grad_to_val_ratio\": " << bs_grad_ratio << ",\n"
     << "  \"bs_gamma_for\": " << gamma_for << ",\n"
     << "  \"bs_gamma_analytic\": " << an.gamma << ",\n"
     << "  \"gbm_1step\": {\"delta\": ";
  write_est_json(jf, gbm.delta);
  jf << ", \"vega\": ";
  write_est_json(jf, gbm.vega);
  jf << ", \"rho\": ";
  write_est_json(jf, gbm.rho);
  jf << ", \"theta\": ";
  write_est_json(jf, gbm.theta);
  jf << ", \"time_sec\": " << gbm.time_sec << "},\n"
     << "  \"asian\": {\"price\": ";
  write_est_json(jf, asian.price);
  jf << ", \"delta\": ";
  write_est_json(jf, asian.delta);
  jf << ", \"delta_formula\": ";
  write_est_json(jf, asian.delta_formula);
  jf << ", \"vega\": ";
  write_est_json(jf, asian.vega);
  jf << ", \"vega_bump\": ";
  write_est_json(jf, asian.vega_bump);
  jf << ", \"rho\": ";
  write_est_json(jf, asian.rho);
  jf << ", \"nodes_per_path\": " << asian.nodes_per_path << ", \"paths_per_sec\": "
     << asian.paths_per_sec << ", \"time_sec\": " << asian.time_sec << "},\n"
     << "  \"digital\": {\"analytic_delta\": " << digital.analytic_delta
     << ", \"raw_pathwise\": " << digital.pathwise_raw_delta << ", \"optimal_eps\": "
     << digital.optimal_eps << ", \"smoothed\": [";
  for (std::size_t i = 0; i < digital.smoothed.size(); ++i) {
    const auto& r = digital.smoothed[i];
    if (i > 0) jf << ", ";
    jf << "{\"eps\": " << r.eps << ", \"delta\": ";
    write_est_json(jf, r.delta);
    jf << ", \"bias\": " << r.bias << ", \"mse\": " << r.mse << "}";
  }
  jf << "]},\n"
     << "  \"basket_arith\": {\"max_diff_se\": " << arith_bask.max_diff_se
     << ", \"violations_3se\": " << arith_bask.violations_3se
     << ", \"violations_4se\": " << arith_bask.violations_4se
     << ", \"nodes_per_path\": " << arith_bask.nodes_per_path
     << ", \"peak_tape_bytes\": " << arith_bask.peak_tape_bytes
     << ", \"aad_time_sec\": " << arith_bask.aad_time_sec
     << ", \"fd_time_sec\": " << arith_bask.fd_time_sec
     << ", \"speedup_vs_fd\": " << arith_bask.speedup_vs_fd
     << ", \"grad_to_val_ratio\": " << arith_bask.grad_to_val_ratio
     << ", \"paths_per_sec\": " << arith_bask.paths_per_sec << "},\n"
     << "  \"basket_geom\": {\"max_diff_se\": " << geom_bask.max_diff_se
     << ", \"violations_3se\": " << geom_bask.violations_3se
     << ", \"violations_4se\": " << geom_bask.violations_4se
     << ", \"nodes_per_path\": " << geom_bask.nodes_per_path
     << ", \"aad_time_sec\": " << geom_bask.aad_time_sec << "}\n"
     << "}\n";

  std::cout << "BS Gamma FOR=" << gamma_for << " vs Analytic=" << an.gamma
            << " | grad/val=" << bs_grad_ratio << "x\n"
            << "Asian price=" << asian.price.mean << " delta=" << asian.delta.mean
            << " vega=" << asian.vega.mean << " (" << asian.paths_per_sec << " paths/s)\n"
            << "Basket Arith (67 grads): max_z=" << arith_bask.max_diff_se
            << " viol_3se=" << arith_bask.violations_3se
            << " speedup=" << arith_bask.speedup_vs_fd << "x\n"
            << "Basket Geom (67 grads): max_z=" << geom_bask.max_diff_se
            << " viol_3se=" << geom_bask.violations_3se << "\n";
  return 0;
}
