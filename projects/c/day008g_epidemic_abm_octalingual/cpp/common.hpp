#pragma once
#include <cstdint>
#include <vector>

constexpr int64_t BP = 10000;
constexpr int MAX_DAYS = 200;
constexpr int ST_S = 0, ST_E = 1, ST_IA = 2, ST_IS = 3, ST_H = 4, ST_ICU = 5, ST_R = 6, ST_D = 7;

enum {
    BR_INF_HH = 0, BR_INF_SW, BR_INF_COMM, BR_FAIL_VAX, BR_FAIL_MASK,
    BR_FAIL_QUAR, BR_PROG_IA, BR_PROG_IS, BR_HOSP_ADM, BR_HOSP_REJ,
    BR_ICU_ADM, BR_ICU_REJ, BR_DEATH, BR_RECOVER, BR_WANING,
    BR_TEST_POS, BR_TEST_NEG, BR_QUAR_TRACE, BR_VAX_DOSE, BR_VAX_OUT,
    BR_VAR_IMPORT, NUM_BRANCHES
};

struct Config {
    int64_t master_seed = 0, parity_N = 0, parity_days = 0, main_days = 0;
    int64_t beta_bp_init = 0, asymp_prob_bp = 0, asymp_rel_bp = 0;
    int64_t latent_days = 0, inf_days_min = 0, inf_days_max = 0;
    int64_t hosp_beds_per_10k = 0, icu_beds_per_10k = 0;
    int64_t waning_min_days = 0, waning_bp = 0;
    int64_t variant_day = 0, variant_mult_bp = 0;
    int64_t school_close_day = 0, lockdown_day = 0, mask_day = 0, mask_eff_bp = 0;
    int64_t test_day = 0, test_prob_bp = 0, trace_prob_bp = 0;
    int64_t vax_day = 0, vax_daily_bp = 0, vax_eff1_bp = 0, vax_eff2_bp = 0, vax_hosp_bp = 0;
    int64_t init_infected = 0;
};

struct Xoshiro256 { uint64_t s[4]{}; };

struct Person {
    int32_t age_grp = 0, hh_id = 0, sw_id = -1, status = 0, variant = 0;
    int32_t state_day = 0, dur_days = 0, inf_day = -1, vax_doses = 0, vax_day = -1;
    int32_t compliance = 0, quar_until = -1, super_flag = 0, offspring = 0, ever_inf = 0, hosp_rej = 0;
};

struct TransRec { int32_t day, src_id, dst_id, layer, variant, gen_int; };
struct ProgRec { int32_t day, pid, from_st, to_st; };

struct RunParams {
    int64_t N = 0, days = 0, beta_bp = 0;
    int64_t compliance_ovr = -1, test_ovr = -1, vax_day_ovr = -1, vax_cov_ovr = -1;
    int64_t school_ovr = -1, lock_ovr = -1, var_day_ovr = -1, mask_ovr = -1;
    int64_t contact_scale_bp = 10000, vax_perfect = 0;
    int32_t interv_on = 0;
};

struct Sim {
    Config cfg{};
    RunParams rp{};
    int64_t N = 0, hosp_cap = 0, icu_cap = 0, hosp_occ = 0, icu_occ = 0;
    std::vector<Person> pop;
    int32_t n_hh = 0, n_sw = 0;
    std::vector<int32_t> hh_off, hh_mem, sw_off, sw_mem, seen_day;
    Xoshiro256 rng[5]{};
    std::vector<int32_t> cand_tgt, cand_src, cand_lay, cand_var;
    int64_t daily_cnt[MAX_DAYS][8]{};
    int64_t daily_inc_age[MAX_DAYS][3]{};
    int64_t daily_deaths_age[MAX_DAYS][3]{};
    int64_t daily_hosp[MAX_DAYS][4]{};
    int64_t daily_tests[MAX_DAYS][2]{};
    int64_t daily_quar[MAX_DAYS]{};
    int64_t daily_vax[MAX_DAYS][2]{};
    int64_t daily_var[MAX_DAYS][2]{};
    int64_t layer_inf[3]{}, age_inf[3]{}, age_tot[3]{};
    int64_t branches[NUM_BRANCHES]{};
    int64_t inv_violations = 0, total_trans = 0, total_reinf = 0;
    uint64_t event_hash = 0xcbf29ce484222325ULL;
    int32_t log_enabled = 0;
    std::vector<TransRec> trans;
    std::vector<ProgRec> progs;
};

uint64_t splitmix64(uint64_t *x);
uint64_t fnv1a_bytes(const void *buf, size_t n);
uint64_t fnv1a_u32pair(uint32_t a, uint32_t b);
uint64_t fnv1a_i64(uint64_t h, int64_t v);
void x256_seed(Xoshiro256 *st, uint64_t seed);
uint64_t x256_next(Xoshiro256 *st);
int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi);

bool load_config(const char *path, Config &c);
void sim_init(Sim &s, const Config &c, const RunParams &rp, int32_t log_en);
void sim_set_state(Sim &s, int32_t pid, int32_t new_st, int32_t day, int32_t dur);
void sim_contact_step(Sim &s, int32_t day);
void sim_disease_step(Sim &s, int32_t day);
void sim_interv_step(Sim &s, int32_t day);
void sim_run(Sim &s);
int check_invariants(Sim &s, int32_t day);
double analytical_r0(const Config &c, int64_t beta_bp, double c_hh, double c_sw, double c_comm);
double solve_final_size(double r0);

void write_s01_to_s09(const Config &c, const char *out_dir);
void write_s10_to_s37(Sim &s_no, Sim &s_main, const Config &c, const char *out_dir, int64_t beta_cal);
void write_s38_svgs(Sim &s_no, Sim &s_main, const char *out_dir);
