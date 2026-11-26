#define _POSIX_C_SOURCE 200809L
#include "md_stats.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

double md_wall_time_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

double md_mean(const double *x, int n) {
    if (n <= 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i];
    return s / (double)n;
}

double md_variance(const double *x, int n) {
    if (n <= 1) return 0.0;
    double m = md_mean(x, n);
    double s2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - m;
        s2 += d * d;
    }
    return s2 / (double)(n - 1);
}

md_block_stats md_block_average(const double *data, int n, int n_blocks, double *block_means_out) {
    md_block_stats bs;
    memset(&bs, 0, sizeof(bs));
    if (n <= 0) return bs;
    if (n_blocks < 2) n_blocks = 2;
    if (n_blocks > n) n_blocks = n;

    int bsz = n / n_blocks;
    if (bsz < 1) {
        bsz = 1;
        n_blocks = n;
    }
    double *bm = (double *)malloc((size_t)n_blocks * sizeof(double));
    bs.min_val = data[0];
    bs.max_val = data[0];
    for (int i = 0; i < n; ++i) {
        if (data[i] < bs.min_val) bs.min_val = data[i];
        if (data[i] > bs.max_val) bs.max_val = data[i];
    }

    for (int b = 0; b < n_blocks; ++b) {
        double s = 0.0;
        int off = b * bsz;
        for (int k = 0; k < bsz; ++k) {
            s += data[off + k];
        }
        bm[b] = s / (double)bsz;
        if (block_means_out) block_means_out[b] = bm[b];
    }

    bs.mean = md_mean(bm, n_blocks);
    double var_b = md_variance(bm, n_blocks);
    bs.sd = sqrt(md_variance(data, n_blocks * bsz));
    bs.se = sqrt(var_b / (double)n_blocks);
    bs.n_blocks = n_blocks;
    bs.block_size = bsz;
    free(bm);
    return bs;
}

double md_autocorr_time(const double *x, int n, int max_lag, double *acf_out) {
    if (n <= 2) return 0.5;
    if (max_lag >= n / 2) max_lag = n / 2;
    double m = md_mean(x, n);
    double c0 = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - m;
        c0 += d * d;
    }
    c0 /= (double)n;
    if (c0 <= 1e-30) {
        if (acf_out) {
            for (int k = 0; k < max_lag; ++k) acf_out[k] = (k == 0) ? 1.0 : 0.0;
        }
        return 0.5;
    }
    double tau = 0.5;
    bool stopped = false;
    for (int lag = 0; lag < max_lag; ++lag) {
        double sum = 0.0;
        int cnt = n - lag;
        for (int i = 0; i < cnt; ++i) {
            sum += (x[i] - m) * (x[i + lag] - m);
        }
        double r = (sum / (double)cnt) / c0;
        if (acf_out) acf_out[lag] = r;
        if (lag > 0 && !stopped) {
            if (r <= 0.0) {
                stopped = true;
            } else {
                tau += r;
            }
        }
    }
    return tau;
}

md_ols_result md_ols_fit(const double *x, const double *y, int n) {
    md_ols_result res;
    memset(&res, 0, sizeof(res));
    res.n_pts = n;
    if (n < 2) return res;

    double mx = md_mean(x, n);
    double my = md_mean(y, n);
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - mx;
        double dy = y[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (sxx <= 1e-30) return res;

    res.slope = sxy / sxx;
    res.intercept = my - res.slope * mx;
    double sse = 0.0;
    for (int i = 0; i < n; ++i) {
        double r = y[i] - (res.intercept + res.slope * x[i]);
        sse += r * r;
    }
    res.r2 = (syy > 1e-30) ? (1.0 - sse / syy) : 1.0;
    if (n > 2) {
        double s2 = sse / (double)(n - 2);
        res.slope_se = sqrt(s2 / sxx);
    } else {
        res.slope_se = 0.0;
    }
    return res;
}

static int cmp_double_asc(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

void md_median_mad(const double *x, int n, double *median_out, double *mad_out) {
    if (n <= 0) {
        *median_out = 0.0;
        *mad_out = 0.0;
        return;
    }
    double *tmp = (double *)malloc((size_t)n * sizeof(double));
    memcpy(tmp, x, (size_t)n * sizeof(double));
    qsort(tmp, (size_t)n, sizeof(double), cmp_double_asc);
    double med = (n % 2 == 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    for (int i = 0; i < n; ++i) {
        tmp[i] = fabs(x[i] - med);
    }
    qsort(tmp, (size_t)n, sizeof(double), cmp_double_asc);
    double mad = (n % 2 == 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    free(tmp);
    *median_out = med;
    *mad_out = mad;
}

double md_maxwell_speed_cdf(double v, double temp) {
    if (v <= 0.0) return 0.0;
    double a = sqrt(temp);
    double u = v / (sqrt(2.0) * a);
    return erf(u) - sqrt(2.0 / MD_PI) * (v / a) * exp(-0.5 * (v * v) / temp);
}

static double ks_kolmogorov_p(double d, int n) {
    /* Asymptotic Kolmogorov distribution with Stephens small-sample correction */
    double en = sqrt((double)n);
    double lam = (en + 0.12 + 0.11 / en) * d;
    if (lam < 1e-8) return 1.0;
    double sum = 0.0;
    for (int k = 1; k <= 100; ++k) {
        double term = 2.0 * exp(-2.0 * (double)(k * k) * lam * lam);
        if (k % 2 == 1) sum += term;
        else sum -= term;
        if (term < 1e-12) break;
    }
    if (sum < 0.0) sum = 0.0;
    if (sum > 1.0) sum = 1.0;
    return sum;
}

void md_ks_maxwell(const double *speeds, int n, double temp, double *d_stat, double *p_val) {
    if (n <= 0) {
        *d_stat = 0.0;
        *p_val = 1.0;
        return;
    }
    double *s = (double *)malloc((size_t)n * sizeof(double));
    memcpy(s, speeds, (size_t)n * sizeof(double));
    qsort(s, (size_t)n, sizeof(double), cmp_double_asc);

    double dmax = 0.0;
    for (int i = 0; i < n; ++i) {
        double f = md_maxwell_speed_cdf(s[i], temp);
        double f_lo = (double)i / (double)n;
        double f_hi = (double)(i + 1) / (double)n;
        double d1 = fabs(f - f_lo);
        double d2 = fabs(f_hi - f);
        if (d1 > dmax) dmax = d1;
        if (d2 > dmax) dmax = d2;
    }
    free(s);
    *d_stat = dmax;
    *p_val = ks_kolmogorov_p(dmax, n);
}

uint64_t md_fnv1a_init(void) {
    return 14695981039346656037ULL;
}

uint64_t md_fnv1a_update(uint64_t h, const void *data, size_t len) {
    const uint8_t *p = (const uint8_t *)data;
    for (size_t i = 0; i < len; ++i) {
        h ^= (uint64_t)p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

uint64_t md_fnv1a_double(uint64_t h, double val) {
    /* Normalize -0.0 to +0.0 for canonical bit representation */
    if (val == 0.0) val = 0.0;
    uint64_t bits = 0;
    memcpy(&bits, &val, sizeof(bits));
    return md_fnv1a_update(h, &bits, sizeof(bits));
}

uint64_t md_fnv1a_system_state(const md_system *sys) {
    uint64_t h = md_fnv1a_init();
    for (int i = 0; i < sys->n; ++i) {
        h = md_fnv1a_double(h, sys->pos[i].x);
        h = md_fnv1a_double(h, sys->pos[i].y);
        h = md_fnv1a_double(h, sys->pos[i].z);
        h = md_fnv1a_double(h, sys->vel[i].x);
        h = md_fnv1a_double(h, sys->vel[i].y);
        h = md_fnv1a_double(h, sys->vel[i].z);
        h = md_fnv1a_double(h, sys->force[i].x);
        h = md_fnv1a_double(h, sys->force[i].y);
        h = md_fnv1a_double(h, sys->force[i].z);
    }
    h = md_fnv1a_double(h, sys->pe);
    h = md_fnv1a_double(h, sys->virial);
    return h;
}

static size_t encode_varint_zigzag(int64_t v, uint8_t *out) {
    uint64_t zz = ((uint64_t)v << 1) ^ (uint64_t)(v >> 63);
    size_t n = 0;
    while (zz >= 0x80ULL) {
        out[n++] = (uint8_t)((zz & 0x7FULL) | 0x80ULL);
        zz >>= 7;
    }
    out[n++] = (uint8_t)zz;
    return n;
}

static size_t decode_varint_zigzag(const uint8_t *in, int64_t *v_out) {
    uint64_t zz = 0;
    int shift = 0;
    size_t n = 0;
    for (;;) {
        uint8_t b = in[n++];
        zz |= (uint64_t)(b & 0x7FU) << shift;
        if ((b & 0x80U) == 0) break;
        shift += 7;
    }
    *v_out = (int64_t)((zz >> 1) ^ (uint64_t)(-(int64_t)(zz & 1ULL)));
    return n;
}

md_compress_report md_compress_trajectory_roundtrip(const double *coords, int n_frames,
                                                    int n_vals, double quant_step) {
    md_compress_report rep;
    memset(&rep, 0, sizeof(rep));
    size_t total = (size_t)n_frames * (size_t)n_vals;
    rep.raw_bytes = total * sizeof(double);

    /* Quantize coordinates to int64 grid */
    int64_t *q = (int64_t *)malloc(total * sizeof(int64_t));
    for (size_t i = 0; i < total; ++i) {
        q[i] = llround(coords[i] / quant_step);
    }

    /* Compute frame-to-frame deltas */
    int64_t *deltas = (int64_t *)malloc(total * sizeof(int64_t));
    for (int j = 0; j < n_vals; ++j) {
        deltas[j] = q[j];
    }
    for (int f = 1; f < n_frames; ++f) {
        for (int j = 0; j < n_vals; ++j) {
            size_t idx = (size_t)f * (size_t)n_vals + (size_t)j;
            size_t prev = (size_t)(f - 1) * (size_t)n_vals + (size_t)j;
            deltas[idx] = q[idx] - q[prev];
        }
    }

    /* RLE + ZigZag Varint encode: runs of zeros encoded with sentinel (0, run_len) */
    uint8_t *buf = (uint8_t *)malloc(total * 12 + 64);
    size_t pos = 0;
    for (size_t i = 0; i < total; ) {
        if (deltas[i] == 0) {
            size_t run = 0;
            while (i + run < total && deltas[i + run] == 0 && run < 65535) {
                run++;
            }
            pos += encode_varint_zigzag(0, buf + pos);
            pos += encode_varint_zigzag((int64_t)run, buf + pos);
            i += run;
        } else {
            pos += encode_varint_zigzag(deltas[i], buf + pos);
            i++;
        }
    }
    rep.comp_bytes = pos;
    rep.ratio = (pos > 0) ? ((double)rep.raw_bytes / (double)pos) : 0.0;

    /* Decode and verify exact lossless roundtrip on quantized stream */
    int64_t *dec_d = (int64_t *)malloc(total * sizeof(int64_t));
    size_t rpos = 0;
    size_t out_i = 0;
    while (rpos < pos && out_i < total) {
        int64_t val = 0;
        rpos += decode_varint_zigzag(buf + rpos, &val);
        if (val == 0) {
            int64_t run = 0;
            rpos += decode_varint_zigzag(buf + rpos, &run);
            for (int64_t k = 0; k < run && out_i < total; ++k) {
                dec_d[out_i++] = 0;
            }
        } else {
            dec_d[out_i++] = val;
        }
    }

    /* Reconstruct quantized coordinates */
    int64_t *rec = (int64_t *)malloc(total * sizeof(int64_t));
    for (int j = 0; j < n_vals; ++j) {
        rec[j] = dec_d[j];
    }
    for (int f = 1; f < n_frames; ++f) {
        for (int j = 0; j < n_vals; ++j) {
            size_t idx = (size_t)f * (size_t)n_vals + (size_t)j;
            size_t prev = (size_t)(f - 1) * (size_t)n_vals + (size_t)j;
            rec[idx] = rec[prev] + dec_d[idx];
        }
    }

    int64_t max_err = 0;
    for (size_t i = 0; i < total; ++i) {
        int64_t e = llabs(rec[i] - q[i]);
        if (e > max_err) max_err = e;
    }
    rep.max_err_q = max_err;
    rep.exact_roundtrip = (out_i == total && max_err == 0);

    free(q);
    free(deltas);
    free(buf);
    free(dec_d);
    free(rec);
    return rep;
}
