#ifndef ASIAN_SVG_H
#define ASIAN_SVG_H

#include "asian/options.h"
#include <stddef.h>

typedef struct {
  size_t n;
  double se_naive;
  double se_antithetic;
  double se_control_variate;
} asian_conv_point_t;

typedef struct {
  const char *name;
  double price;
  double se;
  double ci_low;
  double ci_high;
  double vrf;
} asian_est_summary_t;

typedef struct {
  double naive;
  double antithetic;
  double control_variate;
} asian_slopes_t;

asian_status_t asian_write_convergence_svg(const char *path, const asian_conv_point_t *points,
                                           size_t len, const asian_slopes_t *slopes);
asian_status_t asian_write_estimators_svg(const char *path, const asian_est_summary_t *estimators,
                                          size_t len);

#endif /* ASIAN_SVG_H */
