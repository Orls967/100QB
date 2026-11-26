#ifndef EPI_COMMON_H
#define EPI_COMMON_H

#include <stdint.h>
#include <stdio.h>

#define BP 10000LL
#define MAX_DAYS 200
#define ST_S   0
#define ST_E   1
#define ST_IA  2
#define ST_IS  3
#define ST_H   4
#define ST_ICU 5
#define ST_R   6
#define ST_D   7

enum {
    BR_INF_HH = 0, BR_INF_SW, BR_INF_COMM, BR_FAIL_VAX, BR_FAIL_MASK,
    BR_FAIL_QUAR, BR_PROG_IA, BR_PROG_IS, BR_HOSP_ADM, BR_HOSP_REJ,
    BR_ICU_ADM, BR_ICU_REJ, BR_DEATH, BR_RECOVER, BR_WANING,
    BR_TEST_POS, BR_TEST_NEG, BR_QUAR_TRACE, BR_VAX_DOSE, BR_VAX_OUT,
    BR_VAR_IMPORT, NUM_BRANCHES
};

typedef struct {
    int64_t master_seed, parity_N, parity_days, main_days;
    int64_t beta_bp_init, asymp_prob_bp, asymp_rel_bp;
    int64_t latent_days, inf_days_min, inf_days_max;
    int64_t hosp_beds_per_10k, icu_beds_per_10k;
    int64_t waning_min_days, waning_bp;
    int64_t variant_day, variant_mult_bp;
    int64_t school_close_day, lockdown_day, mask_day, mask_eff_bp;
    int64_t test_day, test_prob_bp, trace_prob_bp;
    int64_t vax_day, vax_daily_bp, vax_eff1_bp, vax_eff2_bp, vax_hosp_bp;
    int64_t init_infected;
} Config;

typedef struct { uint64_t s[4]; } Xoshiro256;

typedef struct {
    int32_t age_grp;      /* 0: 0-17, 1: 18-59, 2: 60+ */
    int32_t hh_id;
    int32_t sw_id;        /* school or workplace id, -1 for 60+ */
    int32_t status;       /* 0..7 */
    int32_t variant;      /* 0 or 1 */
    int32_t state_day;    /* day entered current status */
    int32_t dur_days;     /* planned duration in current state */
    int32_t inf_day;      /* day infected */
    int32_t vax_doses;    /* 0, 1, 2 */
    int32_t vax_day;      /* day of last dose */
    int32_t compliance;   /* basis points 0..10000 */
    int32_t quar_until;   /* day quarantine/isolation ends (-1 if none) */
    int32_t super_flag;   /* 1 if 10% community superspreader */
    int32_t offspring;    /* secondary infections caused */
    int32_t ever_inf;     /* 1 if ever infected */
    int32_t hosp_rej;     /* 1 if rejected from H or ICU */
} Person;

typedef struct {
    int32_t day, src_id, dst_id, layer, variant, gen_int;
} TransRec;

typedef struct {
    int32_t day, pid, from_st, to_st;
} ProgRec;

typedef struct {
    int64_t N, days;
    int64_t beta_bp, compliance_ovr, test_ovr, vax_day_ovr, vax_cov_ovr;
    int64_t school_ovr, lock_ovr, var_day_ovr, mask_ovr, contact_scale_bp, vax_perfect;
    int32_t interv_on;
} RunParams;

typedef struct {
    Config cfg;
    RunParams rp;
    int64_t N, hosp_cap, icu_cap, hosp_occ, icu_occ;
    Person *pop;
    int32_t n_hh, n_sw;
    int32_t *hh_off, *hh_mem;
    int32_t *sw_off, *sw_mem;
    Xoshiro256 rng[5]; /* 0:pop, 1:contact, 2:disease, 3:interv, 4:hosp */
    int32_t *cand_tgt, *cand_src, *cand_lay, *cand_var;
    int32_t *seen_day;
    int32_t n_cand, cap_cand;
    int64_t daily_cnt[MAX_DAYS][8];
    int64_t daily_inc_age[MAX_DAYS][3];
    int64_t daily_deaths_age[MAX_DAYS][3];
    int64_t daily_hosp[MAX_DAYS][4]; /* h_occ, icu_occ, h_rej, icu_rej */
    int64_t daily_tests[MAX_DAYS][2]; /* tests, pos */
    int64_t daily_quar[MAX_DAYS];
    int64_t daily_vax[MAX_DAYS][2]; /* d1, d2 */
    int64_t daily_var[MAX_DAYS][2]; /* v0, v1 new inf */
    int64_t layer_inf[3], age_inf[3], age_tot[3];
    int64_t branches[NUM_BRANCHES];
    int64_t inv_violations, total_trans, total_reinf;
    uint64_t event_hash;
    int32_t log_enabled;
    TransRec *trans; int32_t n_trans, cap_trans;
    ProgRec *progs; int32_t n_progs, cap_progs;
} Sim;

uint64_t splitmix64(uint64_t *x);
uint64_t fnv1a_bytes(const void *buf, size_t n);
uint64_t fnv1a_u32pair(uint32_t a, uint32_t b);
uint64_t fnv1a_i64(uint64_t h, int64_t v);
void x256_seed(Xoshiro256 *st, uint64_t seed);
uint64_t x256_next(Xoshiro256 *st);
int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi);

int load_config(const char *path, Config *c);
void sim_init(Sim *s, const Config *c, const RunParams *rp, int32_t log_en);
void sim_free(Sim *s);
void sim_set_state(Sim *s, int32_t pid, int32_t new_st, int32_t day, int32_t dur);
void build_population(Sim *s);
void sim_contact_step(Sim *s, int32_t day);
void sim_disease_step(Sim *s, int32_t day);
void sim_interv_step(Sim *s, int32_t day);
void sim_run(Sim *s);
int check_invariants(Sim *s, int32_t day);
double analytical_r0(const Config *c, int64_t beta_bp, double c_hh, double c_sw, double c_comm);
double solve_final_size(double r0);

void write_s01_to_s09(const Config *c, const char *out_dir);
void write_s10_to_s37(Sim *s_no, Sim *s_main, const Config *c, const char *out_dir, int64_t beta_cal);
void write_s38_svgs(Sim *s_no, Sim *s_main, const char *out_dir);

#endif
