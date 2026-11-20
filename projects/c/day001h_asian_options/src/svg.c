#include "asian/svg.h"

#include <math.h>
#include <stdio.h>

static double map_x_log(size_t n, double min_log_x, double max_log_x, int pad_l, int plot_w) {
  return (double)pad_l + ((log10((double)n) - min_log_x) / (max_log_x - min_log_x)) * (double)plot_w;
}

static double map_y_log(double se, double min_log_y, double max_log_y, int pad_t, int plot_h) {
  return (double)pad_t +
         (1.0 - (log10(se) - min_log_y) / (max_log_y - min_log_y)) * (double)plot_h;
}

asian_status_t asian_write_convergence_svg(const char *path, const asian_conv_point_t *points,
                                           size_t len, const asian_slopes_t *slopes) {
  if (path == NULL || points == NULL || slopes == NULL || len == 0U) {
    return ASIAN_ERR_INVALID_ARG;
  }
  FILE *f = fopen(path, "w");
  if (f == NULL) {
    return ASIAN_ERR_IO;
  }

  const int w = 900;
  const int h = 520;
  const int pad_l = 85;
  const int pad_r = 30;
  const int pad_t = 55;
  const int pad_b = 65;
  const int plot_w = w - pad_l - pad_r;
  const int plot_h = h - pad_t - pad_b;

  const double min_log_x = log10(900.0);
  const double max_log_x = log10(230000.0);
  double min_se = points[0].se_control_variate;
  double max_se = points[0].se_naive;
  for (size_t i = 1; i < len; ++i) {
    if (points[i].se_control_variate < min_se) {
      min_se = points[i].se_control_variate;
    }
    if (points[i].se_naive > max_se) {
      max_se = points[i].se_naive;
    }
  }
  const double min_log_y = log10(min_se * 0.7);
  const double max_log_y = log10(max_se * 1.35);

  fprintf(f,
          "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d "
          "%d\" font-family=\"system-ui, -apple-system, sans-serif\">\n",
          w, h, w, h);
  fprintf(f, "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", w, h);
  fprintf(f,
          "<text x=\"%d\" y=\"30\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"bold\" "
          "fill=\"#111827\">Konvergensi Standard Error vs Jumlah Path (Log-Log, C17)</text>\n",
          w / 2);

  const size_t x_ticks[3] = {1000U, 10000U, 100000U};
  for (size_t i = 0; i < 3U; ++i) {
    const double x = map_x_log(x_ticks[i], min_log_x, max_log_x, pad_l, plot_w);
    fprintf(f,
            "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e5e7eb\" "
            "stroke-dasharray=\"4,4\"/>\n",
            x, pad_t, x, pad_t + plot_h);
    fprintf(f,
            "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-size=\"12\" "
            "fill=\"#374151\">%zu</text>\n",
            x, pad_t + plot_h + 22, x_ticks[i]);
  }

  for (int exp_i = (int)ceil(min_log_y); exp_i <= (int)floor(max_log_y); ++exp_i) {
    const double val = pow(10.0, (double)exp_i);
    const double y = map_y_log(val, min_log_y, max_log_y, pad_t, plot_h);
    fprintf(f,
            "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e5e7eb\" "
            "stroke-dasharray=\"4,4\"/>\n",
            pad_l, y, pad_l + plot_w, y);
    fprintf(f,
            "<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-size=\"12\" "
            "fill=\"#374151\">%.0e</text>\n",
            pad_l - 10, y + 4.0, val);
  }

  fprintf(f,
          "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"#9ca3af\"/>\n",
          pad_l, pad_t, plot_w, plot_h);

  const char *colors[3] = {"#1f6feb", "#d1495b", "#2a9d8f"};
  const double slope_vals[3] = {slopes->naive, slopes->antithetic, slopes->control_variate};
  const char *names[3] = {"Naive MC", "Antithetic", "Control Variate"};

  for (int s = 0; s < 3; ++s) {
    fprintf(f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2.5\" points=\"", colors[s]);
    for (size_t i = 0; i < len; ++i) {
      const double se = (s == 0) ? points[i].se_naive
                                 : ((s == 1) ? points[i].se_antithetic : points[i].se_control_variate);
      fprintf(f, "%s%.1f,%.1f", (i > 0U ? " " : ""),
              map_x_log(points[i].n, min_log_x, max_log_x, pad_l, plot_w),
              map_y_log(se, min_log_y, max_log_y, pad_t, plot_h));
    }
    fprintf(f, "\"/>\n");
    for (size_t i = 0; i < len; ++i) {
      const double se = (s == 0) ? points[i].se_naive
                                 : ((s == 1) ? points[i].se_antithetic : points[i].se_control_variate);
      fprintf(f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4\" fill=\"%s\"/>\n",
              map_x_log(points[i].n, min_log_x, max_log_x, pad_l, plot_w),
              map_y_log(se, min_log_y, max_log_y, pad_t, plot_h), colors[s]);
    }
  }

  fprintf(f,
          "<rect x=\"%d\" y=\"%d\" width=\"244\" height=\"72\" rx=\"6\" fill=\"#ffffff\" "
          "fill-opacity=\"0.9\" stroke=\"#d1d5db\"/>\n",
          pad_l + plot_w - 252, pad_t + 8);
  for (int s = 0; s < 3; ++s) {
    const int ly = pad_t + 22 + s * 22;
    const int lx = pad_l + plot_w - 240;
    fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-width=\"2.5\"/>\n",
            lx, ly, lx + 24, ly, colors[s]);
    fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#111827\">%s (slope=%.3f)</text>\n",
            lx + 32, ly + 4, names[s], slope_vals[s]);
  }

  fprintf(f,
          "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\" font-size=\"13\" fill=\"#111827\">"
          "Jumlah Path N (skala log)</text>\n",
          pad_l + plot_w / 2, h - 16);
  fprintf(f,
          "<text x=\"22\" y=\"%d\" text-anchor=\"middle\" font-size=\"13\" fill=\"#111827\" "
          "transform=\"rotate(-90 22 %d)\">Standard Error (SE)</text>\n",
          pad_t + plot_h / 2, pad_t + plot_h / 2);
  fprintf(f, "</svg>\n");

  if (fclose(f) != 0) {
    return ASIAN_ERR_IO;
  }
  return ASIAN_OK;
}

asian_status_t asian_write_estimators_svg(const char *path, const asian_est_summary_t *estimators,
                                          size_t len) {
  if (path == NULL || estimators == NULL || len == 0U) {
    return ASIAN_ERR_INVALID_ARG;
  }
  FILE *f = fopen(path, "w");
  if (f == NULL) {
    return ASIAN_ERR_IO;
  }

  const int w = 900;
  const int h = 440;
  const int pad_l = 210;
  const int pad_r = 60;
  const int pad_t = 60;
  const int pad_b = 60;
  const int plot_w = w - pad_l - pad_r;
  const int plot_h = h - pad_t - pad_b;

  double min_x = estimators[0].ci_low;
  double max_x = estimators[0].ci_high;
  for (size_t i = 1; i < len; ++i) {
    if (estimators[i].ci_low < min_x) {
      min_x = estimators[i].ci_low;
    }
    if (estimators[i].ci_high > max_x) {
      max_x = estimators[i].ci_high;
    }
  }
  const double margin = (max_x - min_x) * 0.15;
  min_x -= margin;
  max_x += margin;

  const char *colors[3] = {"#1f6feb", "#d1495b", "#2a9d8f"};

  fprintf(f,
          "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d "
          "%d\" font-family=\"system-ui, -apple-system, sans-serif\">\n",
          w, h, w, h);
  fprintf(f, "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", w, h);
  fprintf(f,
          "<text x=\"%d\" y=\"32\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"bold\" "
          "fill=\"#111827\">Perbandingan Estimasi Harga + CI 95%% &amp; VRF (200k Path, "
          "C17)</text>\n",
          w / 2);
  fprintf(f,
          "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"#9ca3af\"/>\n",
          pad_l, pad_t, plot_w, plot_h);

  for (size_t i = 0; i < len; ++i) {
    const double cy = (double)pad_t + (((double)i + 0.5) / (double)len) * (double)plot_h;
    const double x_lo = (double)pad_l + ((estimators[i].ci_low - min_x) / (max_x - min_x)) * (double)plot_w;
    const double x_hi = (double)pad_l + ((estimators[i].ci_high - min_x) / (max_x - min_x)) * (double)plot_w;
    const double x_mid = (double)pad_l + ((estimators[i].price - min_x) / (max_x - min_x)) * (double)plot_w;
    const char *col = colors[i % 3U];

    fprintf(f,
            "<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-size=\"12\" "
            "fill=\"#111827\">%s</text>\n",
            pad_l - 12, cy + 4.0, estimators[i].name);
    fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"3\"/>\n",
            x_lo, cy, x_hi, cy, col);
    fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"2\"/>\n",
            x_lo, cy - 6.0, x_lo, cy + 6.0, col);
    fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"2\"/>\n",
            x_hi, cy - 6.0, x_hi, cy + 6.0, col);
    fprintf(f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"5\" fill=\"%s\"/>\n", x_mid, cy, col);
    fprintf(f,
            "<text x=\"%.1f\" y=\"%.1f\" font-size=\"11\" fill=\"#374151\">%.4f (VRF %.1fx)</text>\n",
            x_hi + 8.0, cy + 4.0, estimators[i].price, estimators[i].vrf);
  }

  fprintf(f,
          "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\" font-size=\"13\" fill=\"#111827\">"
          "Harga Opsi Asia</text>\n",
          pad_l + plot_w / 2, h - 18);
  fprintf(f, "</svg>\n");

  if (fclose(f) != 0) {
    return ASIAN_ERR_IO;
  }
  return ASIAN_OK;
}
