#ifndef MD_STATS_H
#define MD_STATS_H

#include "md_types.h"

/* High-resolution monotonic wall-clock timer in seconds */
double md_wall_time_sec(void);

/* Non-overlapping block averaging (n_blocks fixed a priori, e.g. 10) */
md_block_stats md_block_average(const double *data, int n, int n_blocks, double *block_means_out);

/* Sample mean and variance (unbiased n-1) */
double md_mean(const double *x, int n);
double md_variance(const double *x, int n);

/* Normalized autocorrelation C(lag) and integrated autocorrelation time tau_int */
double md_autocorr_time(const double *x, int n, int max_lag, double *acf_out);

/* Ordinary least squares linear regression y = a + b*x with SE of slope */
md_ols_result md_ols_fit(const double *x, const double *y, int n);

/* Robust location/scale: median and MAD (1.4826 * median(|x - med|) or raw MAD) */
void md_median_mad(const double *x, int n, double *median_out, double *mad_out);

/* Maxwell-Boltzmann speed CDF F(v; T) in 3D */
double md_maxwell_speed_cdf(double v, double temp);

/* One-sample Kolmogorov-Smirnov test against Maxwell-Boltzmann speed CDF */
void md_ks_maxwell(const double *speeds, int n, double temp, double *d_stat, double *p_val);

/* FNV-1a 64-bit hash helpers */
uint64_t md_fnv1a_init(void);
uint64_t md_fnv1a_update(uint64_t h, const void *data, size_t len);
uint64_t md_fnv1a_double(uint64_t h, double val);
uint64_t md_fnv1a_system_state(const md_system *sys);

/* Quantized delta + RLE lossless trajectory compression roundtrip */
typedef struct {
    size_t raw_bytes;
    size_t comp_bytes;
    double ratio;
    bool exact_roundtrip;
    int64_t max_err_q;
} md_compress_report;

md_compress_report md_compress_trajectory_roundtrip(const double *coords, int n_frames,
                                                    int n_vals, double quant_step);

#endif /* MD_STATS_H */
