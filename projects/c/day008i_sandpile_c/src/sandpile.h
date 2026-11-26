#ifndef SANDPILE_H
#define SANDPILE_H

#include <stdint.h>
#include <stdio.h>

#define MASTER_SEED 20261126ULL
__extension__ typedef __int128 i128;

typedef struct { uint64_t s[4]; } Rng;

typedef struct {
    int n;
    int32_t *z, *q_a, *q_b;
    int64_t *odom, sink;
    uint64_t *stamp_a, *stamp_b, cur_stamp;
} Grid;

typedef struct {
    int64_t size;
    int32_t area, dur;
} Aval;

typedef struct {
    int abelian_configs, abelian_orders, abelian_mismatches;
    uint64_t abelian_digest;
    i128 det_val[7];
    int det_ovf[7], safe_det_n;
    int64_t bf_n2, bf_n3;
    int id_recurrent, id_d4_ok, id_ee_ok, id_er_tested, id_er_ok;
    int ord_tested, ord_div_ok, selftest_ok, svgs_written;
    int64_t mass_max_diff, soc_aval_64, soc_aval_128, warmup_64, warmup_128, bench_topples;
    double rad_slope, rad_se, tau_64, tau_64_se, tau_128, tau_128_se;
    double tau_win_lo, tau_win_hi, tail_64, tail_128;
    double h_mean_ctr, h_mean_se, h_mean_full, h_p[4], h_p_se[4], bench_med, bench_mad;
} Metrics;

void rng_seed(Rng *r, uint64_t seed);
uint64_t rng_u64(Rng *r);
int32_t rng_int(Rng *r, int32_t lo, int32_t hi);
uint64_t fnv1a_init(void);
uint64_t fnv1a_bytes(uint64_t h, const void *buf, size_t len);
uint64_t grid_state_hash(const Grid *g);

void grid_init(Grid *g, int n);
void grid_clear(Grid *g);
void grid_free(Grid *g);
void stabilize_det(Grid *g);
Aval stabilize_waves(Grid *g);
void stabilize_rand(Grid *g, Rng *rng);
int is_recurrent(const Grid *g, int32_t *stk, uint8_t *deg);
void compute_identity(int n, int32_t *e_out);
void fmt_i128(i128 v, char *buf);

int cmd_selftest(Metrics *m);
int cmd_abelian(Metrics *m, const char *out_dir);
int cmd_group(Metrics *m, const char *out_dir);
int cmd_radius(Metrics *m, const char *out_dir);
int cmd_soc_heights(Metrics *m, const char *out_dir, int32_t *stat_128_out);
int cmd_pics(Metrics *m, const char *out_dir, const char *fig_dir, const int32_t *stat_128);

#endif
