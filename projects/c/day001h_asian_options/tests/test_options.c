#include "asian/options.h"
#include "asian/rng.h"

#include <math.h>
#include <stdio.h>

static int g_failures = 0;

#define TEST_ASSERT(cond)                                                                          \
  do {                                                                                             \
    if (!(cond)) {                                                                                 \
      fprintf(stderr, "FAIL %s (%s:%d): %s\n", __func__, __FILE__, __LINE__, #cond);               \
      g_failures += 1;                                                                             \
      return;                                                                                      \
    }                                                                                              \
  } while (0)

static asian_params_t default_params(void) {
  asian_params_t p = {100.0, 100.0, 0.05, 0.20, 1.0, 252U};
  return p;
}

static void test_bs_put_call_parity(void) {
  const double strikes[3] = {90.0, 100.0, 105.0};
  for (size_t i = 0; i < 3U; ++i) {
    asian_params_t p = default_params();
    p.k = strikes[i];
    double c = 0.0;
    double put = 0.0;
    TEST_ASSERT(asian_bs_call(&p, &c) == ASIAN_OK);
    TEST_ASSERT(asian_bs_put(&p, &put) == ASIAN_OK);
    const double rhs = p.s0 - p.k * exp(-p.r * p.t);
    TEST_ASSERT(fabs((c - put) - rhs) < 1e-10);
  }
}

static void test_geometric_asian_matches_large_mc(void) {
  const asian_params_t p = default_params();
  double exact = 0.0;
  TEST_ASSERT(asian_geometric_call(&p, &exact) == ASIAN_OK);

  const size_t n_paths = 40000U;
  const double dt = p.t / (double)p.n_steps;
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * sqrt(dt);
  const double disc = exp(-p.r * p.t);
  const double log_s0 = log(p.s0);
  const double inv_n = 1.0 / (double)p.n_steps;

  asian_rng_t rng;
  asian_rng_seed(&rng, 42ULL);
  double sum = 0.0;
  double sum_sq = 0.0;
  for (size_t i = 0; i < n_paths; ++i) {
    double log_s = log_s0;
    double sum_log_s = 0.0;
    for (size_t step = 0; step < p.n_steps; ++step) {
      log_s += drift + vol * asian_rng_normal(&rng);
      sum_log_s += log_s;
    }
    const double geom_avg = exp(sum_log_s * inv_n);
    const double pay = (geom_avg > p.k) ? (geom_avg - p.k) * disc : 0.0;
    sum += pay;
    sum_sq += pay * pay;
  }
  const double n = (double)n_paths;
  const double mean = sum / n;
  const double var = (sum_sq - n * mean * mean) / (n - 1.0);
  const double se = sqrt(var / n);
  TEST_ASSERT(fabs(mean - exact) < 3.0 * se);
}

static void test_control_variate_beats_naive_se(void) {
  const asian_params_t p = default_params();
  asian_estimator_t naive;
  asian_cv_result_t cv;
  TEST_ASSERT(asian_price_naive(&p, 20000U, 42ULL, 0U, &naive) == ASIAN_OK);
  TEST_ASSERT(asian_price_control_variate(&p, 20000U, 42ULL, 0U, &cv) == ASIAN_OK);
  TEST_ASSERT(cv.se < naive.se / 10.0);
}

static void test_deterministic_seed_and_threads(void) {
  asian_params_t p = default_params();
  p.n_steps = 50U;
  asian_cv_result_t a;
  asian_cv_result_t b;
  TEST_ASSERT(asian_price_control_variate(&p, 10000U, 7ULL, 1U, &a) == ASIAN_OK);
  TEST_ASSERT(asian_price_control_variate(&p, 10000U, 7ULL, 4U, &b) == ASIAN_OK);
  TEST_ASSERT(a.price == b.price);
  TEST_ASSERT(a.se == b.se);
}

static void test_delta_pathwise_matches_bump(void) {
  const asian_params_t p = default_params();
  asian_delta_t pw;
  asian_delta_t bmp;
  TEST_ASSERT(asian_delta_pathwise(&p, 50000U, 42ULL, 0U, &pw) == ASIAN_OK);
  TEST_ASSERT(asian_delta_bump(&p, 50000U, 0.01, 42ULL, 0U, &bmp) == ASIAN_OK);
  const double se_comb = hypot(pw.se, bmp.se);
  TEST_ASSERT(fabs(pw.delta - bmp.delta) < 3.0 * se_comb);
}

static void test_rng_normal_mean_variance(void) {
  asian_rng_t rng;
  asian_rng_seed(&rng, 42ULL);
  const size_t n = 40000U;
  double sum = 0.0;
  double sum_sq = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double z = asian_rng_normal(&rng);
    sum += z;
    sum_sq += z * z;
  }
  const double mean = sum / (double)n;
  const double var = sum_sq / (double)n - mean * mean;
  TEST_ASSERT(fabs(mean) < 0.02);
  TEST_ASSERT(fabs(var - 1.0) < 0.03);
}

static void test_invalid_inputs_return_error(void) {
  asian_params_t bad = default_params();
  bad.n_steps = 0U;
  asian_estimator_t est;
  asian_delta_t dlt;
  double out = 0.0;

  TEST_ASSERT(asian_bs_call(NULL, &out) == ASIAN_ERR_INVALID_ARG);
  TEST_ASSERT(asian_geometric_call(&bad, &out) == ASIAN_ERR_INVALID_ARG);
  TEST_ASSERT(asian_price_naive(&bad, 1000U, 42ULL, 0U, &est) == ASIAN_ERR_INVALID_ARG);

  const asian_params_t ok = default_params();
  TEST_ASSERT(asian_price_antithetic(&ok, 101U, 42ULL, 0U, &est) == ASIAN_ERR_INVALID_ARG);
  TEST_ASSERT(asian_delta_bump(&ok, 100U, 0.0, 42ULL, 0U, &dlt) == ASIAN_ERR_INVALID_ARG);
  TEST_ASSERT(asian_fit_loglog_slope(NULL, NULL, 1U, &out) == ASIAN_ERR_INVALID_ARG);
}

int main(void) {
  test_bs_put_call_parity();
  test_geometric_asian_matches_large_mc();
  test_control_variate_beats_naive_se();
  test_deterministic_seed_and_threads();
  test_delta_pathwise_matches_bump();
  test_rng_normal_mean_variance();
  test_invalid_inputs_return_error();

  if (g_failures > 0) {
    fprintf(stderr, "%d test(s) failed.\n", g_failures);
    return 1;
  }
  printf("All 7 C17 tests passed.\n");
  return 0;
}
