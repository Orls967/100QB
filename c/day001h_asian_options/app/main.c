#include "asian/options.h"
#include "asian/svg.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

static double now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

static int ensure_dir(const char *path) {
  if (mkdir(path, 0755) == 0 || errno == EEXIST) {
    return 0;
  }
  return -1;
}

static int find_repo_root(char *buf, size_t size) {
  if (getcwd(buf, size) == NULL) {
    return -1;
  }
  char probe[4096];
  while (1) {
    snprintf(probe, sizeof(probe), "%s/src/qb100", buf);
    struct stat st;
    if (stat(probe, &st) == 0 && S_ISDIR(st.st_mode)) {
      return 0;
    }
    char *slash = strrchr(buf, '/');
    if (slash == NULL || slash == buf) {
      return -1;
    }
    *slash = '\0';
  }
}

static void write_json_escaped(FILE *f, const char *s) {
  fputc('"', f);
  for (const unsigned char *p = (const unsigned char *)s; *p != '\0'; ++p) {
    if (*p == '"' || *p == '\\') {
      fputc('\\', f);
      fputc(*p, f);
    } else if (*p == '\n') {
      fputs("\\n", f);
    } else {
      fputc(*p, f);
    }
  }
  fputc('"', f);
}

int main(void) {
  const asian_params_t params = {100.0, 100.0, 0.05, 0.20, 1.0, 252U};
  const size_t n_paths = 200000U;
  const uint64_t seed = 42ULL;

  const double t_total_start = now_ms();

  double t0 = now_ms();
  double bs_euro = 0.0;
  asian_estimator_t mc_euro;
  if (asian_bs_call(&params, &bs_euro) != ASIAN_OK ||
      asian_price_european_mc(&params, n_paths, seed, 0U, &mc_euro) != ASIAN_OK) {
    fprintf(stderr, "Error in sanity check\n");
    return 1;
  }
  const double t_sanity = now_ms() - t0;
  const double diff_euro = mc_euro.price - bs_euro;
  const double z_euro = fabs(diff_euro) / mc_euro.se;

  t0 = now_ms();
  asian_estimator_t res_naive;
  if (asian_price_naive(&params, n_paths, seed, 0U, &res_naive) != ASIAN_OK) {
    fprintf(stderr, "Error in naive MC\n");
    return 1;
  }
  const double t_naive = now_ms() - t0;

  t0 = now_ms();
  asian_estimator_t res_anti;
  if (asian_price_antithetic(&params, n_paths, seed, 0U, &res_anti) != ASIAN_OK) {
    fprintf(stderr, "Error in antithetic MC\n");
    return 1;
  }
  const double t_anti = now_ms() - t0;

  t0 = now_ms();
  asian_cv_result_t res_cv;
  if (asian_price_control_variate(&params, n_paths, seed, 0U, &res_cv) != ASIAN_OK) {
    fprintf(stderr, "Error in control variate MC\n");
    return 1;
  }
  const double t_cv = now_ms() - t0;

  t0 = now_ms();
  asian_delta_t delta_pw;
  if (asian_delta_pathwise(&params, n_paths, seed, 0U, &delta_pw) != ASIAN_OK) {
    fprintf(stderr, "Error in pathwise delta\n");
    return 1;
  }
  const double t_pw = now_ms() - t0;

  t0 = now_ms();
  asian_delta_t delta_bmp;
  if (asian_delta_bump(&params, n_paths, 0.01, seed, 0U, &delta_bmp) != ASIAN_OK) {
    fprintf(stderr, "Error in bump delta\n");
    return 1;
  }
  const double t_bmp = now_ms() - t0;

  const double core_ms = now_ms() - t_total_start;

  const double delta_diff = delta_pw.delta - delta_bmp.delta;
  const double delta_z = fabs(delta_diff) / hypot(delta_pw.se, delta_bmp.se);

  const double var_naive = res_naive.se * res_naive.se;
  const asian_est_summary_t estimators[3] = {
      {"Naive MC (200k)", res_naive.price, res_naive.se, res_naive.ci_low, res_naive.ci_high, 1.0},
      {"Antithetic (100k pasang)", res_anti.price, res_anti.se, res_anti.ci_low, res_anti.ci_high,
       var_naive / (res_anti.se * res_anti.se)},
      {"Control Variate (Geom. Asia)", res_cv.price, res_cv.se, res_cv.ci_low, res_cv.ci_high,
       var_naive / (res_cv.se * res_cv.se)},
  };

  t0 = now_ms();
  const size_t n_grid[8] = {1000U, 2000U, 5000U, 10000U, 20000U, 50000U, 100000U, 200000U};
  asian_conv_point_t conv_pts[8];
  double se_n[8];
  double se_a[8];
  double se_c[8];
  for (size_t i = 0; i < 8U; ++i) {
    asian_estimator_t rn;
    asian_estimator_t ra;
    asian_cv_result_t rc;
    if (asian_price_naive(&params, n_grid[i], seed, 0U, &rn) != ASIAN_OK ||
        asian_price_antithetic(&params, n_grid[i], seed, 0U, &ra) != ASIAN_OK ||
        asian_price_control_variate(&params, n_grid[i], seed, 0U, &rc) != ASIAN_OK) {
      fprintf(stderr, "Error in convergence loop\n");
      return 1;
    }
    se_n[i] = rn.se;
    se_a[i] = ra.se;
    se_c[i] = rc.se;
    conv_pts[i].n = n_grid[i];
    conv_pts[i].se_naive = rn.se;
    conv_pts[i].se_antithetic = ra.se;
    conv_pts[i].se_control_variate = rc.se;
  }
  const double t_conv = now_ms() - t0;

  asian_slopes_t slopes;
  if (asian_fit_loglog_slope(n_grid, se_n, 8U, &slopes.naive) != ASIAN_OK ||
      asian_fit_loglog_slope(n_grid, se_a, 8U, &slopes.antithetic) != ASIAN_OK ||
      asian_fit_loglog_slope(n_grid, se_c, 8U, &slopes.control_variate) != ASIAN_OK) {
    fprintf(stderr, "Error fitting slopes\n");
    return 1;
  }

  char root[4096];
  if (find_repo_root(root, sizeof(root)) != 0) {
    fprintf(stderr, "Could not find repo root\n");
    return 1;
  }
  char dir_buf[4096];
  snprintf(dir_buf, sizeof(dir_buf), "%s/outputs", root);
  (void)ensure_dir(dir_buf);
  snprintf(dir_buf, sizeof(dir_buf), "%s/outputs/day001h", root);
  (void)ensure_dir(dir_buf);
  snprintf(dir_buf, sizeof(dir_buf), "%s/outputs/figures", root);
  (void)ensure_dir(dir_buf);

  char json_path[4096];
  snprintf(json_path, sizeof(json_path), "%s/outputs/day001h/results.json", root);
  FILE *jf = fopen(json_path, "w");
  if (jf == NULL) {
    fprintf(stderr, "Cannot open %s\n", json_path);
    return 1;
  }

  fprintf(jf, "{\n");
  fprintf(jf,
          "  \"params\": {\"s0\": %.17g, \"k\": %.17g, \"r\": %.17g, \"sigma\": %.17g, \"t\": "
          "%.17g, \"n_steps\": %zu},\n",
          params.s0, params.k, params.r, params.sigma, params.t, params.n_steps);
  fprintf(jf, "  \"n_paths\": %zu,\n", n_paths);
  fprintf(jf, "  \"seed\": %llu,\n", (unsigned long long)seed);
  fprintf(jf, "  \"sanity_check\": {\n");
  fprintf(jf, "    \"bs_call\": %.17g,\n", bs_euro);
  fprintf(jf,
          "    \"mc_call\": {\"price\": %.17g, \"se\": %.17g, \"ci_low\": %.17g, \"ci_high\": "
          "%.17g},\n",
          mc_euro.price, mc_euro.se, mc_euro.ci_low, mc_euro.ci_high);
  fprintf(jf, "    \"diff\": %.17g,\n", diff_euro);
  fprintf(jf, "    \"z_score\": %.17g\n  },\n", z_euro);

  fprintf(jf, "  \"estimators\": [\n");
  for (size_t i = 0; i < 3U; ++i) {
    fprintf(jf, "    {\"name\": ");
    write_json_escaped(jf, estimators[i].name);
    fprintf(jf,
            ", \"price\": %.17g, \"se\": %.17g, \"ci_low\": %.17g, \"ci_high\": %.17g, \"vrf\": "
            "%.17g",
            estimators[i].price, estimators[i].se, estimators[i].ci_low, estimators[i].ci_high,
            estimators[i].vrf);
    if (i == 2U) {
      fprintf(jf, ", \"beta\": %.17g, \"y_exact\": %.17g", res_cv.beta, res_cv.y_exact);
    }
    fprintf(jf, "}%s\n", (i + 1U < 3U) ? "," : "");
  }
  fprintf(jf, "  ],\n");

  fprintf(jf, "  \"delta\": {\n");
  fprintf(jf,
          "    \"pathwise\": {\"delta\": %.17g, \"se\": %.17g, \"ci_low\": %.17g, \"ci_high\": "
          "%.17g},\n",
          delta_pw.delta, delta_pw.se, delta_pw.ci_low, delta_pw.ci_high);
  fprintf(jf,
          "    \"bump_1pct\": {\"delta\": %.17g, \"se\": %.17g, \"ci_low\": %.17g, \"ci_high\": "
          "%.17g},\n",
          delta_bmp.delta, delta_bmp.se, delta_bmp.ci_low, delta_bmp.ci_high);
  fprintf(jf, "    \"diff\": %.17g,\n", delta_diff);
  fprintf(jf, "    \"z_score\": %.17g\n  },\n", delta_z);

  fprintf(jf, "  \"convergence\": [\n");
  for (size_t i = 0; i < 8U; ++i) {
    fprintf(jf,
            "    {\"n\": %zu, \"se_naive\": %.17g, \"se_antithetic\": %.17g, "
            "\"se_control_variate\": %.17g}%s\n",
            conv_pts[i].n, conv_pts[i].se_naive, conv_pts[i].se_antithetic,
            conv_pts[i].se_control_variate, (i + 1U < 8U) ? "," : "");
  }
  fprintf(jf, "  ],\n");
  fprintf(jf, "  \"slopes\": {\"naive\": %.17g, \"antithetic\": %.17g, \"control_variate\": %.17g},\n",
          slopes.naive, slopes.antithetic, slopes.control_variate);
  fprintf(jf,
          "  \"timings_ms\": {\"sanity_check\": %.3f, \"naive\": %.3f, \"antithetic\": %.3f, "
          "\"control_variate\": %.3f, \"delta_pathwise\": %.3f, \"delta_bump\": %.3f, "
          "\"convergence\": %.3f, \"total_200k_core\": %.3f}\n",
          t_sanity, t_naive, t_anti, t_cv, t_pw, t_bmp, t_conv, core_ms);
  fprintf(jf, "}\n");
  fclose(jf);

  char svg1[4096];
  char svg2[4096];
  snprintf(svg1, sizeof(svg1), "%s/outputs/figures/day001h_asian-convergence-loglog.svg", root);
  snprintf(svg2, sizeof(svg2), "%s/outputs/figures/day001h_asian-estimators-ci.svg", root);
  if (asian_write_convergence_svg(svg1, conv_pts, 8U, &slopes) != ASIAN_OK ||
      asian_write_estimators_svg(svg2, estimators, 3U) != ASIAN_OK) {
    fprintf(stderr, "Failed to write SVG files\n");
    return 1;
  }

  printf("BS Eropa: %.17g | MC Eropa: %.6f ± %.6f (%.2f SE)\n", bs_euro, mc_euro.price, mc_euro.se,
         z_euro);
  printf("Geom. Asian exact: %.17g | beta: %.6f\n", res_cv.y_exact, res_cv.beta);
  for (size_t i = 0; i < 3U; ++i) {
    printf("%-28s: %.6f ± %.6f (CI [%.6f, %.6f]) | VRF = %.2fx\n", estimators[i].name,
           estimators[i].price, estimators[i].se, estimators[i].ci_low, estimators[i].ci_high,
           estimators[i].vrf);
  }
  printf("Delta PW: %.6f ± %.6f | Bump: %.6f ± %.6f | diff = %+.6f (%.2f SE)\n", delta_pw.delta,
         delta_pw.se, delta_bmp.delta, delta_bmp.se, delta_diff, delta_z);
  printf("Slopes -> Naive: %.4f, Antithetic: %.4f, CV: %.4f\n", slopes.naive, slopes.antithetic,
         slopes.control_variate);
  printf("Total 200k core time: %.1f ms\n", core_ms);
  return 0;
}
