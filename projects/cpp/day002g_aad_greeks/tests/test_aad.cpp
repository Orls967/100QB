#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <cmath>
#include <vector>

#include "aad/cholesky.hpp"
#include "aad/dual.hpp"
#include "aad/mc_engine.hpp"
#include "aad/models.hpp"
#include "aad/number.hpp"
#include "aad/rng.hpp"

using namespace aad;

namespace {

template <typename Scalar>
Scalar random_expr(const Scalar& x, const Scalar& y, const Scalar& z) {
  const Scalar t1 = math::exp(Scalar(0.3) * x) * math::sqrt(math::abs(y) + Scalar(1.0));
  const Scalar t2 = math::log(Scalar(2.5) + z * z) / (Scalar(1.2) + x * x);
  const Scalar t3 = math::normal_cdf(x - Scalar(0.5) * y) + math::pow(Scalar(1.5) + z * z, Scalar(0.7));
  return math::max(t1 + t2, t3) - Scalar(0.2) * x * y;
}

} // namespace

TEST_CASE("AAD vs Dual vs Central Finite Differences on random expressions") {
  Xoshiro256StarStar rng(12345);
  for (int trial = 0; trial < 5; ++trial) {
    const double xv = rng.next_uniform() * 1.5 - 0.5;
    const double yv = rng.next_uniform() * 1.5 - 0.5;
    const double zv = rng.next_uniform() * 1.5 - 0.5;

    Tape tape;
    TapeGuard guard(tape);
    const Number x(xv), y(yv), z(zv);
    const Number res = random_expr(x, y, z);
    tape.backprop(res.index());

    const auto dx = random_expr(Dual<double>(xv, 1.0), Dual<double>(yv, 0.0), Dual<double>(zv, 0.0));
    const auto dy = random_expr(Dual<double>(xv, 0.0), Dual<double>(yv, 1.0), Dual<double>(zv, 0.0));
    const auto dz = random_expr(Dual<double>(xv, 0.0), Dual<double>(yv, 0.0), Dual<double>(zv, 1.0));

    CHECK(std::abs(res.value() - dx.val) < 1e-13);
    CHECK(std::abs(x.adjoint() - dx.dot) < 1e-12);
    CHECK(std::abs(y.adjoint() - dy.dot) < 1e-12);
    CHECK(std::abs(z.adjoint() - dz.dot) < 1e-12);

    constexpr double h = 1e-6;
    const double fdx = (random_expr(xv + h, yv, zv) - random_expr(xv - h, yv, zv)) / (2.0 * h);
    CHECK(std::abs(x.adjoint() - fdx) / std::max(1.0, std::abs(fdx)) < 1e-6);
  }
}

TEST_CASE("Black-Scholes 5 AAD gradients and Forward-over-Reverse Gamma") {
  const auto an = bs_analytic_greeks(100.0, 100.0, 0.03, 0.20, 1.0);

  Tape tape;
  TapeGuard guard(tape);
  const Number s(100.0), k(100.0), r(0.03), sig(0.20), t(1.0);
  const Number c = bs_call(s, k, r, sig, t);
  tape.backprop(c.index());

  CHECK(std::abs(c.value() - an.price) < 1e-12);
  CHECK(std::abs(s.adjoint() - an.delta) / an.delta < 1e-12);
  CHECK(std::abs(k.adjoint() - an.d_k) / std::abs(an.d_k) < 1e-12);
  CHECK(std::abs(r.adjoint() - an.rho) / an.rho < 1e-12);
  CHECK(std::abs(sig.adjoint() - an.vega) / an.vega < 1e-12);
  CHECK(std::abs(t.adjoint() - an.theta) / an.theta < 1e-12);

  // Forward-over-Reverse Gamma
  tape.clear();
  const Number s_base(100.0), k2(100.0), r2(0.03), sig2(0.20), t2(1.0);
  const Dual<Number> ds(s_base, Number(1.0));
  const Dual<Number> dc = bs_call(ds, Dual<Number>(k2, Number(0.0)), Dual<Number>(r2, Number(0.0)),
                                  Dual<Number>(sig2, Number(0.0)), Dual<Number>(t2, Number(0.0)));
  tape.backprop(dc.dot.index());
  CHECK(std::abs(s_base.adjoint() - an.gamma) / an.gamma < 1e-9);
}

TEST_CASE("Monte Carlo AAD small-n statistical checks, seed reproducibility and thread invariance") {
  const auto an = bs_analytic_greeks(100.0, 100.0, 0.03, 0.20, 1.0);
  const auto gbm1 = run_gbm_1step_aad(12'000, 42, 1);
  const auto gbm4 = run_gbm_1step_aad(12'000, 42, 4);

  CHECK(gbm1.delta.mean == doctest::Approx(gbm4.delta.mean).epsilon(1e-15));
  CHECK(std::abs(gbm4.delta.mean - an.delta) < 3.0 * gbm4.delta.se);
  CHECK(std::abs(gbm4.vega.mean - an.vega) < 3.0 * gbm4.vega.se);

  const auto dig = run_digital_aad(12'000, 42, 4);
  CHECK(dig.pathwise_raw_delta == 0.0);

  const auto geom = run_basket_10_geom(16'000, 42, 4);
  CHECK(geom.violations_4se == 0);
  CHECK(geom.violations_3se <= 2);
}
