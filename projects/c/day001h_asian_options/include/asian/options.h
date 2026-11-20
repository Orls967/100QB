#ifndef ASIAN_OPTIONS_H
#define ASIAN_OPTIONS_H

#include <stddef.h>
#include <stdint.h>

typedef enum {
  ASIAN_OK = 0,
  ASIAN_ERR_INVALID_ARG = 1,
  ASIAN_ERR_ALLOC = 2,
  ASIAN_ERR_ZERO_VAR = 3,
  ASIAN_ERR_THREAD = 4,
  ASIAN_ERR_IO = 5
} asian_status_t;

typedef struct {
  double s0;
  double k;
  double r;
  double sigma;
  double t;
  size_t n_steps;
} asian_params_t;

typedef struct {
  double price;
  double se;
  double ci_low;
  double ci_high;
} asian_estimator_t;

typedef struct {
  double price;
  double se;
  double ci_low;
  double ci_high;
  double beta;
  double y_exact;
} asian_cv_result_t;

typedef struct {
  double delta;
  double se;
  double ci_low;
  double ci_high;
} asian_delta_t;

double asian_norm_cdf(double x);
asian_status_t asian_bs_call(const asian_params_t *p, double *out);
asian_status_t asian_bs_put(const asian_params_t *p, double *out);
asian_status_t asian_geometric_call(const asian_params_t *p, double *out);

asian_status_t asian_price_european_mc(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                       size_t workers, asian_estimator_t *out);
asian_status_t asian_price_naive(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                 size_t workers, asian_estimator_t *out);
asian_status_t asian_price_antithetic(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                      size_t workers, asian_estimator_t *out);
asian_status_t asian_price_control_variate(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                           size_t workers, asian_cv_result_t *out);

asian_status_t asian_delta_pathwise(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                    size_t workers, asian_delta_t *out);
asian_status_t asian_delta_bump(const asian_params_t *p, size_t n_paths, double h_rel,
                                uint64_t seed, size_t workers, asian_delta_t *out);

asian_status_t asian_fit_loglog_slope(const size_t *ns, const double *ses, size_t len,
                                      double *out_slope);

#endif /* ASIAN_OPTIONS_H */
