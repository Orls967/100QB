#ifndef MD_SVG_H
#define MD_SVG_H

#include "md_types.h"

#define MD_SVG_MAX_SERIES 8

typedef struct {
    const char *label;
    const char *color;
    const double *x;
    const double *y;
    const double *yerr; /* Optional error bars (may be NULL) */
    int n;
    bool draw_line;
    bool draw_points;
    bool dashed;
} md_svg_series;

typedef struct {
    const char *title;
    const char *subtitle;
    const char *xlabel;
    const char *ylabel;
    int width;
    int height;
    bool log_x;
    bool log_y;
    bool has_ref_y;
    double ref_y;
    const char *ref_label;
} md_svg_opts;

/* Hand-written multi-series line/scatter/error-bar SVG plot (linear or log-log) */
bool md_svg_write_xy(const char *path, const md_svg_opts *opts,
                     const md_svg_series *series, int n_series);

/* Hand-written bar chart with optional MAD/SE whiskers */
bool md_svg_write_bars(const char *path, const char *title, const char *subtitle,
                       const char *ylabel, const char **labels, const double *vals,
                       const double *errs, const char **colors, int n_bars, bool log_y);

/* Hand-written histogram with optional theoretical curve overlay */
bool md_svg_write_hist(const char *path, const char *title, const char *subtitle,
                       const char *xlabel, const char *ylabel,
                       const double *bin_centers, const double *counts, int n_bins,
                       const double *theory_y, const char *theory_label);

/* Hand-written 2D particle snapshot scatter inside periodic box */
bool md_svg_write_particles_2d(const char *path, const char *title, const char *subtitle,
                               const md_system *sys, const double *color_scalar);

#endif /* MD_SVG_H */
