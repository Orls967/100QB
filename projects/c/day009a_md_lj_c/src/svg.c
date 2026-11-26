#include "md_svg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void svg_escape(FILE *fp, const char *s) {
    if (!s) return;
    for (const char *p = s; *p; ++p) {
        if (*p == '&') fputs("&amp;", fp);
        else if (*p == '<') fputs("&lt;", fp);
        else if (*p == '>') fputs("&gt;", fp);
        else if (*p == '"') fputs("&quot;", fp);
        else fputc(*p, fp);
    }
}

static void svg_header(FILE *fp, int w, int h, const char *title, const char *subtitle) {
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
            w, h, w, h);
    fprintf(fp, "  <rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#fafbfc\"/>\n", w, h);
    fprintf(fp, "  <rect x=\"1\" y=\"1\" width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"#d8dee9\" stroke-width=\"1\"/>\n",
            w - 2, h - 2);
    fprintf(fp, "  <text x=\"%d\" y=\"28\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"15\" font-weight=\"bold\" fill=\"#1f2937\">",
            w / 2);
    svg_escape(fp, title);
    fprintf(fp, "</text>\n");
    if (subtitle && *subtitle) {
        fprintf(fp, "  <text x=\"%d\" y=\"46\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"11\" fill=\"#4b5563\">",
                w / 2);
        svg_escape(fp, subtitle);
        fprintf(fp, "</text>\n");
    }
}

static double safe_tx(double v, bool is_log) {
    if (!is_log) return v;
    return log10((v > 1e-30) ? v : 1e-30);
}

bool md_svg_write_xy(const char *path, const md_svg_opts *opts,
                     const md_svg_series *series, int n_series) {
    FILE *fp = fopen(path, "w");
    if (!fp) return false;

    int w = (opts->width > 0) ? opts->width : 760;
    int h = (opts->height > 0) ? opts->height : 460;
    int ml = 76, mr = 32, mt = 62, mb = 60;
    int pw = w - ml - mr;
    int ph = h - mt - mb;

    double xmin = 1e300, xmax = -1e300;
    double ymin = 1e300, ymax = -1e300;

    for (int s = 0; s < n_series; ++s) {
        for (int i = 0; i < series[s].n; ++i) {
            double xv = series[s].x[i];
            double yv = series[s].y[i];
            if (!isfinite(xv) || !isfinite(yv)) continue;
            if (opts->log_x && xv <= 0.0) continue;
            if (opts->log_y && yv <= 0.0) continue;
            double tx = safe_tx(xv, opts->log_x);
            double ty = safe_tx(yv, opts->log_y);
            if (tx < xmin) xmin = tx;
            if (tx > xmax) xmax = tx;
            if (ty < ymin) ymin = ty;
            if (ty > ymax) ymax = ty;
            if (series[s].yerr && !opts->log_y) {
                double e = fabs(series[s].yerr[i]);
                if (isfinite(e)) {
                    if (ty - e < ymin) ymin = ty - e;
                    if (ty + e > ymax) ymax = ty + e;
                }
            }
        }
    }
    if (opts->has_ref_y && (!opts->log_y || opts->ref_y > 0.0)) {
        double ry = safe_tx(opts->ref_y, opts->log_y);
        if (ry < ymin) ymin = ry;
        if (ry > ymax) ymax = ry;
    }

    if (!(xmin < xmax)) { xmin -= 1.0; xmax += 1.0; }
    if (!(ymin < ymax)) { ymin -= 1.0; ymax += 1.0; }
    double xpad = 0.04 * (xmax - xmin);
    double ypad = 0.08 * (ymax - ymin);
    xmin -= xpad; xmax += xpad;
    ymin -= ypad; ymax += ypad;

    svg_header(fp, w, h, opts->title, opts->subtitle);

    /* Plot frame */
    fprintf(fp, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"#ffffff\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n",
            ml, mt, pw, ph);

    /* Grid and ticks (5 divisions each axis) */
    for (int t = 0; t <= 5; ++t) {
        double frac = (double)t / 5.0;
        double gx = ml + frac * pw;
        double val_x = xmin + frac * (xmax - xmin);
        double disp_x = opts->log_x ? pow(10.0, val_x) : val_x;
        fprintf(fp, "  <line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n",
                gx, mt, gx, mt + ph);
        fprintf(fp, "  <text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">%.3g</text>\n",
                gx, mt + ph + 16, disp_x);

        double gy = mt + ph - frac * ph;
        double val_y = ymin + frac * (ymax - ymin);
        double disp_y = opts->log_y ? pow(10.0, val_y) : val_y;
        fprintf(fp, "  <line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n",
                ml, gy, ml + pw, gy);
        fprintf(fp, "  <text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">%.3g</text>\n",
                ml - 8, gy + 3.5, disp_y);
    }

    /* Axis labels */
    fprintf(fp, "  <text x=\"%d\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\">",
            ml + pw / 2, h - 14);
    svg_escape(fp, opts->xlabel);
    fprintf(fp, "</text>\n");
    fprintf(fp, "  <text x=\"18\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\" transform=\"rotate(-90 18 %d)\">",
            mt + ph / 2, mt + ph / 2);
    svg_escape(fp, opts->ylabel);
    fprintf(fp, "</text>\n");

    /* Optional horizontal reference line */
    if (opts->has_ref_y && (!opts->log_y || opts->ref_y > 0.0)) {
        double ry = safe_tx(opts->ref_y, opts->log_y);
        double py = mt + ph - ((ry - ymin) / (ymax - ymin)) * ph;
        fprintf(fp, "  <line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#6b7280\" stroke-width=\"1.2\" stroke-dasharray=\"5,4\"/>\n",
                ml, py, ml + pw, py);
        if (opts->ref_label && *opts->ref_label) {
            fprintf(fp, "  <text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#4b5563\">",
                    ml + pw - 6, py - 4.0);
            svg_escape(fp, opts->ref_label);
            fprintf(fp, "</text>\n");
        }
    }

    /* Draw series */
    for (int s = 0; s < n_series; ++s) {
        const md_svg_series *sr = &series[s];
        const char *col = sr->color ? sr->color : "#2563eb";
        if (sr->draw_line && sr->n >= 2) {
            fprintf(fp, "  <polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2\"%s points=\"",
                    col, sr->dashed ? " stroke-dasharray=\"6,4\"" : "");
            for (int i = 0; i < sr->n; ++i) {
                double xv = sr->x[i], yv = sr->y[i];
                if (!isfinite(xv) || !isfinite(yv)) continue;
                if (opts->log_x && xv <= 0.0) continue;
                if (opts->log_y && yv <= 0.0) continue;
                double px = ml + ((safe_tx(xv, opts->log_x) - xmin) / (xmax - xmin)) * pw;
                double py = mt + ph - ((safe_tx(yv, opts->log_y) - ymin) / (ymax - ymin)) * ph;
                fprintf(fp, "%.1f,%.1f ", px, py);
            }
            fprintf(fp, "\"/>\n");
        }
        if (sr->yerr && !opts->log_y) {
            for (int i = 0; i < sr->n; ++i) {
                double xv = sr->x[i], yv = sr->y[i], ev = fabs(sr->yerr[i]);
                if (!isfinite(xv) || !isfinite(yv) || !isfinite(ev) || ev <= 0.0) continue;
                if (opts->log_x && xv <= 0.0) continue;
                double px = ml + ((safe_tx(xv, opts->log_x) - xmin) / (xmax - xmin)) * pw;
                double py1 = mt + ph - ((yv - ev - ymin) / (ymax - ymin)) * ph;
                double py2 = mt + ph - ((yv + ev - ymin) / (ymax - ymin)) * ph;
                fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"1.2\"/>\n",
                        px, py1, px, py2, col);
                fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"1.2\"/>\n",
                        px - 3.0, py1, px + 3.0, py1, col);
                fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"1.2\"/>\n",
                        px - 3.0, py2, px + 3.0, py2, col);
            }
        }
        if (sr->draw_points) {
            for (int i = 0; i < sr->n; ++i) {
                double xv = sr->x[i], yv = sr->y[i];
                if (!isfinite(xv) || !isfinite(yv)) continue;
                if (opts->log_x && xv <= 0.0) continue;
                if (opts->log_y && yv <= 0.0) continue;
                double px = ml + ((safe_tx(xv, opts->log_x) - xmin) / (xmax - xmin)) * pw;
                double py = mt + ph - ((safe_tx(yv, opts->log_y) - ymin) / (ymax - ymin)) * ph;
                fprintf(fp, "  <circle cx=\"%.1f\" cy=\"%.1f\" r=\"3.2\" fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"0.8\"/>\n",
                        px, py, col);
            }
        }
    }

    /* Compact legend inside top-right */
    if (n_series > 0) {
        int lx = ml + pw - 195;
        int ly = mt + 10;
        int lh = 16 * n_series + 10;
        fprintf(fp, "  <rect x=\"%d\" y=\"%d\" width=\"188\" height=\"%d\" rx=\"4\" fill=\"#ffffff\" fill-opacity=\"0.9\" stroke=\"#d1d5db\"/>\n",
                lx, ly, lh);
        for (int s = 0; s < n_series; ++s) {
            const char *col = series[s].color ? series[s].color : "#2563eb";
            int cy = ly + 15 + s * 16;
            fprintf(fp, "  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-width=\"2.2\"%s/>\n",
                    lx + 8, cy - 3, lx + 26, cy - 3, col, series[s].dashed ? " stroke-dasharray=\"4,3\"" : "");
            fprintf(fp, "  <text x=\"%d\" y=\"%d\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#1f2937\">",
                    lx + 32, cy);
            svg_escape(fp, series[s].label ? series[s].label : "series");
            fprintf(fp, "</text>\n");
        }
    }

    fprintf(fp, "</svg>\n");
    fclose(fp);
    return true;
}

bool md_svg_write_bars(const char *path, const char *title, const char *subtitle,
                       const char *ylabel, const char **labels, const double *vals,
                       const double *errs, const char **colors, int n_bars, bool log_y) {
    FILE *fp = fopen(path, "w");
    if (!fp) return false;

    int w = 760, h = 460;
    int ml = 80, mr = 30, mt = 62, mb = 74;
    int pw = w - ml - mr;
    int ph = h - mt - mb;

    double ymin = log_y ? 1e300 : 0.0;
    double ymax = -1e300;
    for (int i = 0; i < n_bars; ++i) {
        double v = vals[i];
        double e = errs ? fabs(errs[i]) : 0.0;
        if (log_y) {
            double lv = safe_tx(v > 1e-12 ? v : 1e-12, true);
            if (lv < ymin) ymin = lv;
            if (lv > ymax) ymax = lv;
        } else {
            if (v - e < ymin) ymin = v - e;
            if (v + e > ymax) ymax = v + e;
        }
    }
    if (log_y) {
        ymin -= 0.3;
        ymax += 0.3;
    } else {
        if (!(ymin < ymax)) ymax = ymin + 1.0;
        ymax += 0.12 * (ymax - ymin);
    }

    svg_header(fp, w, h, title, subtitle);
    fprintf(fp, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"#ffffff\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n",
            ml, mt, pw, ph);

    for (int t = 0; t <= 5; ++t) {
        double frac = (double)t / 5.0;
        double gy = mt + ph - frac * ph;
        double val_y = ymin + frac * (ymax - ymin);
        double disp_y = log_y ? pow(10.0, val_y) : val_y;
        fprintf(fp, "  <line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n",
                ml, gy, ml + pw, gy);
        fprintf(fp, "  <text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">%.3g</text>\n",
                ml - 8, gy + 3.5, disp_y);
    }

    double slot = (double)pw / (double)(n_bars > 0 ? n_bars : 1);
    double bw = slot * 0.58;
    double y0_px = mt + ph - ((0.0 - ymin) / (ymax - ymin)) * ph;
    if (log_y || y0_px > mt + ph) y0_px = mt + ph;

    for (int i = 0; i < n_bars; ++i) {
        double cx = ml + ((double)i + 0.5) * slot;
        double tv = log_y ? safe_tx(vals[i] > 1e-12 ? vals[i] : 1e-12, true) : vals[i];
        double py = mt + ph - ((tv - ymin) / (ymax - ymin)) * ph;
        double top = (py < y0_px) ? py : y0_px;
        double bh = fabs(y0_px - py);
        const char *col = (colors && colors[i]) ? colors[i] : "#2563eb";
        fprintf(fp, "  <rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\" fill-opacity=\"0.85\" stroke=\"#1f2937\" stroke-width=\"0.8\"/>\n",
                cx - 0.5 * bw, top, bw, bh, col);

        if (errs && errs[i] > 0.0 && !log_y) {
            double py1 = mt + ph - ((vals[i] - errs[i] - ymin) / (ymax - ymin)) * ph;
            double py2 = mt + ph - ((vals[i] + errs[i] - ymin) / (ymax - ymin)) * ph;
            fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#111827\" stroke-width=\"1.4\"/>\n",
                    cx, py1, cx, py2);
            fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#111827\" stroke-width=\"1.4\"/>\n",
                    cx - 5.0, py1, cx + 5.0, py1);
            fprintf(fp, "  <line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#111827\" stroke-width=\"1.4\"/>\n",
                    cx - 5.0, py2, cx + 5.0, py2);
        }
        fprintf(fp, "  <text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"9.5\" fill=\"#111827\">%.3g</text>\n",
                cx, top - 5.0, vals[i]);
        fprintf(fp, "  <text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">",
                cx, mt + ph + 18);
        svg_escape(fp, labels[i]);
        fprintf(fp, "</text>\n");
    }

    fprintf(fp, "  <text x=\"18\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\" transform=\"rotate(-90 18 %d)\">",
            mt + ph / 2, mt + ph / 2);
    svg_escape(fp, ylabel);
    fprintf(fp, "</text>\n");

    fprintf(fp, "</svg>\n");
    fclose(fp);
    return true;
}

bool md_svg_write_hist(const char *path, const char *title, const char *subtitle,
                       const char *xlabel, const char *ylabel,
                       const double *bin_centers, const double *counts, int n_bins,
                       const double *theory_y, const char *theory_label) {
    FILE *fp = fopen(path, "w");
    if (!fp) return false;

    int w = 760, h = 460;
    int ml = 76, mr = 32, mt = 62, mb = 60;
    int pw = w - ml - mr;
    int ph = h - mt - mb;

    double xmin = bin_centers[0], xmax = bin_centers[n_bins - 1];
    double dx = (n_bins > 1) ? (bin_centers[1] - bin_centers[0]) : 1.0;
    xmin -= 0.5 * dx;
    xmax += 0.5 * dx;
    double ymax = 0.0;
    for (int i = 0; i < n_bins; ++i) {
        if (counts[i] > ymax) ymax = counts[i];
        if (theory_y && theory_y[i] > ymax) ymax = theory_y[i];
    }
    if (ymax <= 0.0) ymax = 1.0;
    ymax *= 1.12;

    svg_header(fp, w, h, title, subtitle);
    fprintf(fp, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"#ffffff\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n",
            ml, mt, pw, ph);

    for (int t = 0; t <= 5; ++t) {
        double frac = (double)t / 5.0;
        double gx = ml + frac * pw;
        double gy = mt + ph - frac * ph;
        fprintf(fp, "  <line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e5e7eb\"/>\n", gx, mt, gx, mt + ph);
        fprintf(fp, "  <line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e5e7eb\"/>\n", ml, gy, ml + pw, gy);
        fprintf(fp, "  <text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">%.3g</text>\n",
                gx, mt + ph + 16, xmin + frac * (xmax - xmin));
        fprintf(fp, "  <text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">%.3g</text>\n",
                ml - 8, gy + 3.5, frac * ymax);
    }

    double bw = (dx / (xmax - xmin)) * pw * 0.92;
    for (int i = 0; i < n_bins; ++i) {
        double cx = ml + ((bin_centers[i] - xmin) / (xmax - xmin)) * pw;
        double bh = (counts[i] / ymax) * ph;
        fprintf(fp, "  <rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#60a5fa\" fill-opacity=\"0.72\" stroke=\"#2563eb\" stroke-width=\"0.6\"/>\n",
                cx - 0.5 * bw, mt + ph - bh, bw, bh);
    }

    if (theory_y) {
        fprintf(fp, "  <polyline fill=\"none\" stroke=\"#dc2626\" stroke-width=\"2.2\" points=\"");
        for (int i = 0; i < n_bins; ++i) {
            double cx = ml + ((bin_centers[i] - xmin) / (xmax - xmin)) * pw;
            double py = mt + ph - (theory_y[i] / ymax) * ph;
            fprintf(fp, "%.1f,%.1f ", cx, py);
        }
        fprintf(fp, "\"/>\n");
        if (theory_label && *theory_label) {
            fprintf(fp, "  <rect x=\"%d\" y=\"%d\" width=\"170\" height=\"26\" rx=\"4\" fill=\"#ffffff\" fill-opacity=\"0.9\" stroke=\"#d1d5db\"/>\n",
                    ml + pw - 178, mt + 10);
            fprintf(fp, "  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#dc2626\" stroke-width=\"2.2\"/>\n",
                    ml + pw - 168, mt + 23, ml + pw - 148, mt + 23);
            fprintf(fp, "  <text x=\"%d\" y=\"%d\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#1f2937\">",
                    ml + pw - 142, mt + 26);
            svg_escape(fp, theory_label);
            fprintf(fp, "</text>\n");
        }
    }

    fprintf(fp, "  <text x=\"%d\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\">",
            ml + pw / 2, h - 14);
    svg_escape(fp, xlabel);
    fprintf(fp, "</text>\n");
    fprintf(fp, "  <text x=\"18\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\" transform=\"rotate(-90 18 %d)\">",
            mt + ph / 2, mt + ph / 2);
    svg_escape(fp, ylabel);
    fprintf(fp, "</text>\n");

    fprintf(fp, "</svg>\n");
    fclose(fp);
    return true;
}

bool md_svg_write_particles_2d(const char *path, const char *title, const char *subtitle,
                               const md_system *sys, const double *color_scalar) {
    FILE *fp = fopen(path, "w");
    if (!fp) return false;

    int w = 620, h = 620;
    int ml = 60, mr = 40, mt = 66, mb = 44;
    int pw = w - ml - mr;
    int ph = h - mt - mb;

    svg_header(fp, w, h, title, subtitle);
    fprintf(fp, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"#ffffff\" stroke=\"#374151\" stroke-width=\"1.4\"/>\n",
            ml, mt, pw, ph);

    double cmin = 0.0, cmax = 1.0;
    if (color_scalar && sys->n > 0) {
        cmin = cmax = color_scalar[0];
        for (int i = 1; i < sys->n; ++i) {
            if (color_scalar[i] < cmin) cmin = color_scalar[i];
            if (color_scalar[i] > cmax) cmax = color_scalar[i];
        }
        if (cmax <= cmin) cmax = cmin + 1.0;
    }

    double rad = (sys->n > 400) ? 3.2 : 5.0;
    for (int i = 0; i < sys->n; ++i) {
        double px = ml + (sys->pos[i].x * sys->box.inv_lx) * pw;
        double py = mt + ph - (sys->pos[i].y * sys->box.inv_ly) * ph;
        int r = 37, g = 99, b = 235;
        if (color_scalar) {
            double t = (color_scalar[i] - cmin) / (cmax - cmin);
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            r = (int)(37 + t * (220 - 37));
            g = (int)(99 + (1.0 - fabs(2.0 * t - 1.0)) * 80);
            b = (int)(235 - t * (235 - 38));
        }
        fprintf(fp, "  <circle cx=\"%.1f\" cy=\"%.1f\" r=\"%.1f\" fill=\"rgb(%d,%d,%d)\" fill-opacity=\"0.82\"/>\n",
                px, py, rad, r, g, b);
    }

    fprintf(fp, "  <text x=\"%d\" y=\"%d\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"11\" fill=\"#374151\">x in [0, %.2f], y in [0, %.2f] (N=%d)</text>\n",
            w / 2, h - 14, sys->box.lx, sys->box.ly, sys->n);
    fprintf(fp, "</svg>\n");
    fclose(fp);
    return true;
}
