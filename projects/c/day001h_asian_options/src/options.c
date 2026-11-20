#include "asian/options.h"
#include "asian/rng.h"

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

#define ASIAN_CHUNK_SIZE 4000U
#define ASIAN_CI_Z95 1.96
#define ASIAN_INV_SQRT2 0.70710678118654752440

typedef enum {
  TASK_EUROPEAN = 0,
  TASK_NAIVE = 1,
  TASK_ANTITHETIC = 2,
  TASK_CV = 3,
  TASK_DELTA_PW = 4,
  TASK_DELTA_BUMP = 5
} asian_task_t;

typedef struct {
  double sum;
  double sum_sq;
  double sx;
  double sy;
  double sxx;
  double syy;
  double sxy;
} chunk_acc_t;

typedef struct {
  size_t wid;
  size_t n_workers;
  size_t n_chunks;
  size_t n_items;
  asian_task_t task;
  const asian_params_t *p;
  uint64_t seed;
  double h_rel;
  chunk_acc_t *results;
} worker_arg_t;

static size_t resolve_workers(size_t workers) {
  if (workers > 0) {
    return workers;
  }
  const long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
  return (ncpu > 0) ? (size_t)ncpu : 1U;
}

static asian_status_t validate_mc(const asian_params_t *p, size_t n_paths) {
  if (p == NULL || p->n_steps < 1U || p->s0 <= 0.0 || p->k <= 0.0 || n_paths < 2U) {
    return ASIAN_ERR_INVALID_ARG;
  }
  return ASIAN_OK;
}

static asian_estimator_t summarize_moments(size_t count, double sum, double sum_sq) {
  const double n = (double)count;
  const double mean = sum / n;
  double var = (sum_sq - n * mean * mean) / (n - 1.0);
  if (var < 0.0) {
    var = 0.0;
  }
  const double se = sqrt(var / n);
  asian_estimator_t res;
  res.price = mean;
  res.se = se;
  res.ci_low = mean - ASIAN_CI_Z95 * se;
  res.ci_high = mean + ASIAN_CI_Z95 * se;
  return res;
}

static chunk_acc_t run_single_chunk(asian_task_t task, const asian_params_t *p, size_t chunk_idx,
                                    size_t len, uint64_t seed, double h_rel) {
  asian_rng_t rng;
  asian_rng_seed(&rng, asian_derive_chunk_seed(seed, chunk_idx));

  const double dt = p->t / (double)p->n_steps;
  const double drift = (p->r - 0.5 * p->sigma * p->sigma) * dt;
  const double vol = p->sigma * sqrt(dt);
  const double disc = exp(-p->r * p->t);
  const double log_s0 = log(p->s0);
  const double inv_n = 1.0 / (double)p->n_steps;

  chunk_acc_t acc = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if (task == TASK_CV) {
    for (size_t i = 0; i < len; ++i) {
      double log_s = log_s0;
      double sum_s = 0.0;
      double sum_log_s = 0.0;
      for (size_t step = 0; step < p->n_steps; ++step) {
        log_s += drift + vol * asian_rng_normal(&rng);
        sum_s += exp(log_s);
        sum_log_s += log_s;
      }
      const double arith_avg = sum_s * inv_n;
      const double geom_avg = exp(sum_log_s * inv_n);
      const double x = (arith_avg > p->k) ? (arith_avg - p->k) * disc : 0.0;
      const double y = (geom_avg > p->k) ? (geom_avg - p->k) * disc : 0.0;
      acc.sx += x;
      acc.sy += y;
      acc.sxx += x * x;
      acc.syy += y * y;
      acc.sxy += x * y;
    }
    return acc;
  }

  if (task == TASK_ANTITHETIC) {
    for (size_t i = 0; i < len; ++i) {
      double log_pos = log_s0;
      double log_neg = log_s0;
      double sum_pos = 0.0;
      double sum_neg = 0.0;
      for (size_t step = 0; step < p->n_steps; ++step) {
        const double inc = vol * asian_rng_normal(&rng);
        log_pos += drift + inc;
        log_neg += drift - inc;
        sum_pos += exp(log_pos);
        sum_neg += exp(log_neg);
      }
      const double avg_pos = sum_pos * inv_n;
      const double avg_neg = sum_neg * inv_n;
      const double p_pos = (avg_pos > p->k) ? (avg_pos - p->k) : 0.0;
      const double p_neg = (avg_neg > p->k) ? (avg_neg - p->k) : 0.0;
      const double pair_avg = 0.5 * (p_pos + p_neg) * disc;
      acc.sum += pair_avg;
      acc.sum_sq += pair_avg * pair_avg;
    }
    return acc;
  }

  const double h = h_rel * p->s0;
  const double scale_up = (p->s0 + h) / p->s0;
  const double scale_dn = (p->s0 - h) / p->s0;

  for (size_t i = 0; i < len; ++i) {
    double log_s = log_s0;
    double sum_s = 0.0;
    for (size_t step = 0; step < p->n_steps; ++step) {
      log_s += drift + vol * asian_rng_normal(&rng);
      if (task != TASK_EUROPEAN) {
        sum_s += exp(log_s);
      }
    }
    double val = 0.0;
    if (task == TASK_EUROPEAN) {
      const double s_t = exp(log_s);
      val = (s_t > p->k) ? (s_t - p->k) * disc : 0.0;
    } else if (task == TASK_NAIVE) {
      const double avg = sum_s * inv_n;
      val = (avg > p->k) ? (avg - p->k) * disc : 0.0;
    } else if (task == TASK_DELTA_PW) {
      const double avg = sum_s * inv_n;
      val = (avg > p->k) ? disc * (avg / p->s0) : 0.0;
    } else {
      const double avg = sum_s * inv_n;
      const double pay_up = (avg * scale_up > p->k) ? (avg * scale_up - p->k) : 0.0;
      const double pay_dn = (avg * scale_dn > p->k) ? (avg * scale_dn - p->k) : 0.0;
      val = disc * (pay_up - pay_dn) / (2.0 * h);
    }
    acc.sum += val;
    acc.sum_sq += val * val;
  }
  return acc;
}

static void *worker_entry(void *raw_arg) {
  const worker_arg_t *arg = (const worker_arg_t *)raw_arg;
  for (size_t c = arg->wid; c < arg->n_chunks; c += arg->n_workers) {
    const size_t start = c * ASIAN_CHUNK_SIZE;
    const size_t rem = arg->n_items - start;
    const size_t len = (rem < ASIAN_CHUNK_SIZE) ? rem : ASIAN_CHUNK_SIZE;
    arg->results[c] = run_single_chunk(arg->task, arg->p, c, len, arg->seed, arg->h_rel);
  }
  return NULL;
}

static asian_status_t execute_parallel(asian_task_t task, const asian_params_t *p, size_t n_items,
                                       uint64_t seed, double h_rel, size_t workers,
                                       chunk_acc_t *out_total) {
  const size_t n_chunks = (n_items + ASIAN_CHUNK_SIZE - 1U) / ASIAN_CHUNK_SIZE;
  chunk_acc_t *results = calloc(n_chunks, sizeof(chunk_acc_t));
  if (results == NULL) {
    return ASIAN_ERR_ALLOC;
  }

  size_t w = resolve_workers(workers);
  if (w > n_chunks) {
    w = n_chunks;
  }

  asian_status_t status = ASIAN_OK;
  if (w <= 1U) {
    for (size_t c = 0; c < n_chunks; ++c) {
      const size_t start = c * ASIAN_CHUNK_SIZE;
      const size_t rem = n_items - start;
      const size_t len = (rem < ASIAN_CHUNK_SIZE) ? rem : ASIAN_CHUNK_SIZE;
      results[c] = run_single_chunk(task, p, c, len, seed, h_rel);
    }
  } else {
    pthread_t *threads = calloc(w, sizeof(pthread_t));
    worker_arg_t *args = calloc(w, sizeof(worker_arg_t));
    if (threads == NULL || args == NULL) {
      free(threads);
      free(args);
      free(results);
      return ASIAN_ERR_ALLOC;
    }

    size_t spawned = 0U;
    for (size_t wid = 0; wid < w; ++wid) {
      args[wid].wid = wid;
      args[wid].n_workers = w;
      args[wid].n_chunks = n_chunks;
      args[wid].n_items = n_items;
      args[wid].task = task;
      args[wid].p = p;
      args[wid].seed = seed;
      args[wid].h_rel = h_rel;
      args[wid].results = results;
      if (pthread_create(&threads[wid], NULL, worker_entry, &args[wid]) != 0) {
        status = ASIAN_ERR_THREAD;
        break;
      }
      spawned += 1U;
    }
    for (size_t wid = 0; wid < spawned; ++wid) {
      (void)pthread_join(threads[wid], NULL);
    }
    free(threads);
    free(args);
  }

  if (status == ASIAN_OK) {
    chunk_acc_t total = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (size_t c = 0; c < n_chunks; ++c) {
      total.sum += results[c].sum;
      total.sum_sq += results[c].sum_sq;
      total.sx += results[c].sx;
      total.sy += results[c].sy;
      total.sxx += results[c].sxx;
      total.syy += results[c].syy;
      total.sxy += results[c].sxy;
    }
    *out_total = total;
  }

  free(results);
  return status;
}

double asian_norm_cdf(double x) {
  return 0.5 * erfc(-x * ASIAN_INV_SQRT2);
}

asian_status_t asian_bs_call(const asian_params_t *p, double *out) {
  if (p == NULL || out == NULL || p->s0 <= 0.0 || p->k <= 0.0) {
    return ASIAN_ERR_INVALID_ARG;
  }
  if (p->t <= 0.0 || p->sigma <= 0.0) {
    const double val = p->s0 - p->k * exp(-p->r * p->t);
    *out = (val > 0.0) ? val : 0.0;
    return ASIAN_OK;
  }
  const double vol_sqrt_t = p->sigma * sqrt(p->t);
  const double d1 = (log(p->s0 / p->k) + (p->r + 0.5 * p->sigma * p->sigma) * p->t) / vol_sqrt_t;
  const double d2 = d1 - vol_sqrt_t;
  *out = p->s0 * asian_norm_cdf(d1) - p->k * exp(-p->r * p->t) * asian_norm_cdf(d2);
  return ASIAN_OK;
}

asian_status_t asian_bs_put(const asian_params_t *p, double *out) {
  if (p == NULL || out == NULL || p->s0 <= 0.0 || p->k <= 0.0) {
    return ASIAN_ERR_INVALID_ARG;
  }
  if (p->t <= 0.0 || p->sigma <= 0.0) {
    const double val = p->k * exp(-p->r * p->t) - p->s0;
    *out = (val > 0.0) ? val : 0.0;
    return ASIAN_OK;
  }
  const double vol_sqrt_t = p->sigma * sqrt(p->t);
  const double d1 = (log(p->s0 / p->k) + (p->r + 0.5 * p->sigma * p->sigma) * p->t) / vol_sqrt_t;
  const double d2 = d1 - vol_sqrt_t;
  *out = p->k * exp(-p->r * p->t) * asian_norm_cdf(-d2) - p->s0 * asian_norm_cdf(-d1);
  return ASIAN_OK;
}

asian_status_t asian_geometric_call(const asian_params_t *p, double *out) {
  if (p == NULL || out == NULL || p->n_steps < 1U || p->s0 <= 0.0 || p->k <= 0.0) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const double n = (double)p->n_steps;
  const double mu_ln_g = log(p->s0) + (p->r - 0.5 * p->sigma * p->sigma) * p->t * (n + 1.0) / (2.0 * n);
  const double var_ln_g = (p->sigma * p->sigma * p->t) * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
  const double sigma_ln_g = sqrt(var_ln_g);
  const double expected_g = exp(mu_ln_g + 0.5 * var_ln_g);
  const double d1 = (mu_ln_g - log(p->k) + var_ln_g) / sigma_ln_g;
  const double d2 = d1 - sigma_ln_g;
  *out = exp(-p->r * p->t) * (expected_g * asian_norm_cdf(d1) - p->k * asian_norm_cdf(d2));
  return ASIAN_OK;
}

asian_status_t asian_price_european_mc(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                       size_t workers, asian_estimator_t *out) {
  if (out == NULL) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const asian_status_t v = validate_mc(p, n_paths);
  if (v != ASIAN_OK) {
    return v;
  }
  chunk_acc_t total;
  const asian_status_t st = execute_parallel(TASK_EUROPEAN, p, n_paths, seed, 0.01, workers, &total);
  if (st != ASIAN_OK) {
    return st;
  }
  *out = summarize_moments(n_paths, total.sum, total.sum_sq);
  return ASIAN_OK;
}

asian_status_t asian_price_naive(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                 size_t workers, asian_estimator_t *out) {
  if (out == NULL) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const asian_status_t v = validate_mc(p, n_paths);
  if (v != ASIAN_OK) {
    return v;
  }
  chunk_acc_t total;
  const asian_status_t st = execute_parallel(TASK_NAIVE, p, n_paths, seed, 0.01, workers, &total);
  if (st != ASIAN_OK) {
    return st;
  }
  *out = summarize_moments(n_paths, total.sum, total.sum_sq);
  return ASIAN_OK;
}

asian_status_t asian_price_antithetic(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                      size_t workers, asian_estimator_t *out) {
  if (out == NULL || p == NULL || p->n_steps < 1U || p->s0 <= 0.0 || p->k <= 0.0 ||
      n_paths < 4U || (n_paths % 2U) != 0U) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const size_t n_pairs = n_paths / 2U;
  chunk_acc_t total;
  const asian_status_t st = execute_parallel(TASK_ANTITHETIC, p, n_pairs, seed, 0.01, workers, &total);
  if (st != ASIAN_OK) {
    return st;
  }
  *out = summarize_moments(n_pairs, total.sum, total.sum_sq);
  return ASIAN_OK;
}

asian_status_t asian_price_control_variate(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                           size_t workers, asian_cv_result_t *out) {
  if (out == NULL) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const asian_status_t v = validate_mc(p, n_paths);
  if (v != ASIAN_OK) {
    return v;
  }
  double y_exact = 0.0;
  const asian_status_t gst = asian_geometric_call(p, &y_exact);
  if (gst != ASIAN_OK) {
    return gst;
  }
  chunk_acc_t total;
  const asian_status_t st = execute_parallel(TASK_CV, p, n_paths, seed, 0.01, workers, &total);
  if (st != ASIAN_OK) {
    return st;
  }

  const double n = (double)n_paths;
  const double mx = total.sx / n;
  const double my = total.sy / n;
  const double cov_xy = (total.sxy - n * mx * my) / (n - 1.0);
  const double var_y = (total.syy - n * my * my) / (n - 1.0);
  if (var_y <= 0.0) {
    return ASIAN_ERR_ZERO_VAR;
  }
  const double beta = cov_xy / var_y;
  const double z_mean = mx - beta * (my - y_exact);
  const double var_x = (total.sxx - n * mx * mx) / (n - 1.0);
  double var_z = var_x + beta * beta * var_y - 2.0 * beta * cov_xy;
  if (var_z < 0.0) {
    var_z = 0.0;
  }
  const double se = sqrt(var_z / n);

  out->price = z_mean;
  out->se = se;
  out->ci_low = z_mean - ASIAN_CI_Z95 * se;
  out->ci_high = z_mean + ASIAN_CI_Z95 * se;
  out->beta = beta;
  out->y_exact = y_exact;
  return ASIAN_OK;
}

asian_status_t asian_delta_pathwise(const asian_params_t *p, size_t n_paths, uint64_t seed,
                                    size_t workers, asian_delta_t *out) {
  if (out == NULL) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const asian_status_t v = validate_mc(p, n_paths);
  if (v != ASIAN_OK) {
    return v;
  }
  chunk_acc_t total;
  const asian_status_t st = execute_parallel(TASK_DELTA_PW, p, n_paths, seed, 0.01, workers, &total);
  if (st != ASIAN_OK) {
    return st;
  }
  const asian_estimator_t est = summarize_moments(n_paths, total.sum, total.sum_sq);
  out->delta = est.price;
  out->se = est.se;
  out->ci_low = est.ci_low;
  out->ci_high = est.ci_high;
  return ASIAN_OK;
}

asian_status_t asian_delta_bump(const asian_params_t *p, size_t n_paths, double h_rel,
                                uint64_t seed, size_t workers, asian_delta_t *out) {
  if (out == NULL || h_rel <= 0.0) {
    return ASIAN_ERR_INVALID_ARG;
  }
  const asian_status_t v = validate_mc(p, n_paths);
  if (v != ASIAN_OK) {
    return v;
  }
  chunk_acc_t total;
  const asian_status_t st = execute_parallel(TASK_DELTA_BUMP, p, n_paths, seed, h_rel, workers, &total);
  if (st != ASIAN_OK) {
    return st;
  }
  const asian_estimator_t est = summarize_moments(n_paths, total.sum, total.sum_sq);
  out->delta = est.price;
  out->se = est.se;
  out->ci_low = est.ci_low;
  out->ci_high = est.ci_high;
  return ASIAN_OK;
}

asian_status_t asian_fit_loglog_slope(const size_t *ns, const double *ses, size_t len,
                                      double *out_slope) {
  if (ns == NULL || ses == NULL || out_slope == NULL || len < 2U) {
    return ASIAN_ERR_INVALID_ARG;
  }
  double *xs = calloc(len, sizeof(double));
  double *ys = calloc(len, sizeof(double));
  if (xs == NULL || ys == NULL) {
    free(xs);
    free(ys);
    return ASIAN_ERR_ALLOC;
  }

  double sum_x = 0.0;
  double sum_y = 0.0;
  for (size_t i = 0; i < len; ++i) {
    xs[i] = log10((double)ns[i]);
    ys[i] = log10(ses[i]);
    sum_x += xs[i];
    sum_y += ys[i];
  }
  const double n = (double)len;
  const double mx = sum_x / n;
  const double my = sum_y / n;
  double num = 0.0;
  double den = 0.0;
  for (size_t i = 0; i < len; ++i) {
    const double dx = xs[i] - mx;
    num += dx * (ys[i] - my);
    den += dx * dx;
  }
  free(xs);
  free(ys);
  if (den <= 0.0) {
    return ASIAN_ERR_ZERO_VAR;
  }
  *out_slope = num / den;
  return ASIAN_OK;
}
