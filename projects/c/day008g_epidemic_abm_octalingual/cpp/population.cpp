#include "common.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>

bool load_config(const char *path, Config &c) {
    FILE *f = fopen(path, "r");
    if (!f) return false;
    c = Config{};
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        int64_t v = atoll(eq + 1);
        if (!strcmp(line, "master_seed")) c.master_seed = v;
        else if (!strcmp(line, "parity_N")) c.parity_N = v;
        else if (!strcmp(line, "parity_days")) c.parity_days = v;
        else if (!strcmp(line, "main_days")) c.main_days = v;
        else if (!strcmp(line, "beta_bp_init")) c.beta_bp_init = v;
        else if (!strcmp(line, "asymp_prob_bp")) c.asymp_prob_bp = v;
        else if (!strcmp(line, "asymp_rel_bp")) c.asymp_rel_bp = v;
        else if (!strcmp(line, "latent_days")) c.latent_days = v;
        else if (!strcmp(line, "inf_days_min")) c.inf_days_min = v;
        else if (!strcmp(line, "inf_days_max")) c.inf_days_max = v;
        else if (!strcmp(line, "hosp_beds_per_10k")) c.hosp_beds_per_10k = v;
        else if (!strcmp(line, "icu_beds_per_10k")) c.icu_beds_per_10k = v;
        else if (!strcmp(line, "waning_min_days")) c.waning_min_days = v;
        else if (!strcmp(line, "waning_bp")) c.waning_bp = v;
        else if (!strcmp(line, "variant_day")) c.variant_day = v;
        else if (!strcmp(line, "variant_mult_bp")) c.variant_mult_bp = v;
        else if (!strcmp(line, "school_close_day")) c.school_close_day = v;
        else if (!strcmp(line, "lockdown_day")) c.lockdown_day = v;
        else if (!strcmp(line, "mask_day")) c.mask_day = v;
        else if (!strcmp(line, "mask_eff_bp")) c.mask_eff_bp = v;
        else if (!strcmp(line, "test_day")) c.test_day = v;
        else if (!strcmp(line, "test_prob_bp")) c.test_prob_bp = v;
        else if (!strcmp(line, "trace_prob_bp")) c.trace_prob_bp = v;
        else if (!strcmp(line, "vax_day")) c.vax_day = v;
        else if (!strcmp(line, "vax_daily_bp")) c.vax_daily_bp = v;
        else if (!strcmp(line, "vax_eff1_bp")) c.vax_eff1_bp = v;
        else if (!strcmp(line, "vax_eff2_bp")) c.vax_eff2_bp = v;
        else if (!strcmp(line, "vax_hosp_bp")) c.vax_hosp_bp = v;
        else if (!strcmp(line, "init_infected")) c.init_infected = v;
    }
    fclose(f);
    return true;
}

void sim_init(Sim &s, const Config &c, const RunParams &rp, int32_t log_en) {
    s = Sim{}; s.cfg = c; s.rp = rp; s.N = rp.N; s.log_enabled = log_en;
    s.hosp_cap = std::max<int64_t>(1, (rp.N * c.hosp_beds_per_10k) / BP);
    s.icu_cap  = std::max<int64_t>(1, (rp.N * c.icu_beds_per_10k) / BP);
    for (int i = 0; i < 5; ++i) {
        uint64_t x = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(100, static_cast<uint32_t>(i));
        x256_seed(&s.rng[i], splitmix64(&x));
    }
    int32_t N = static_cast<int32_t>(rp.N);
    s.pop.assign(N, Person{});
    s.seen_day.assign(N, -1);
    Xoshiro256 *r = &s.rng[0];
    std::vector<int32_t> hh_cnt(N + 1, 0);
    int32_t pid = 0, hh = 0;
    while (pid < N) {
        int32_t sz = static_cast<int32_t>(rand_int(r, 1, 6));
        if (pid + sz > N) sz = N - pid;
        hh_cnt[hh] = sz;
        for (int32_t k = 0; k < sz; ++k) {
            Person &p = s.pop[pid];
            p.hh_id = hh;
            int64_t u = rand_int(r, 1, 100);
            p.age_grp = (u <= 22) ? 0 : ((u <= 82) ? 1 : 2);
            s.age_tot[p.age_grp]++;
            p.compliance = static_cast<int32_t>(rand_int(r, 4000, 9500));
            p.super_flag = (rand_int(r, 1, 10) == 1) ? 1 : 0;
            pid++;
        }
        hh++;
    }
    s.n_hh = hh;
    s.hh_off.assign(hh + 1, 0);
    for (int32_t i = 0; i < hh; ++i) s.hh_off[i + 1] = s.hh_off[i] + hh_cnt[i];
    s.hh_mem.resize(N);
    std::vector<int32_t> cur(hh, 0);
    for (int32_t i = 0; i < N; ++i) { int32_t h = s.pop[i].hh_id; s.hh_mem[s.hh_off[h] + cur[h]++] = i; }

    int32_t n_schools = static_cast<int32_t>(s.age_tot[0] / 25) + 1;
    int32_t n_works = static_cast<int32_t>(s.age_tot[1] / 15) + 1;
    s.n_sw = n_schools + n_works;
    std::vector<int32_t> sw_cnt(s.n_sw + 1, 0);
    for (int32_t i = 0; i < N; ++i) {
        if (s.pop[i].age_grp == 0) { int32_t sw = static_cast<int32_t>(rand_int(r, 0, n_schools - 1)); s.pop[i].sw_id = sw; sw_cnt[sw]++; }
        else if (s.pop[i].age_grp == 1) { int32_t sw = n_schools + static_cast<int32_t>(rand_int(r, 0, n_works - 1)); s.pop[i].sw_id = sw; sw_cnt[sw]++; }
        else s.pop[i].sw_id = -1;
    }
    s.sw_off.assign(s.n_sw + 1, 0);
    for (int32_t i = 0; i < s.n_sw; ++i) s.sw_off[i + 1] = s.sw_off[i] + sw_cnt[i];
    s.sw_mem.resize(N);
    std::vector<int32_t> sw_cur(s.n_sw, 0);
    for (int32_t i = 0; i < N; ++i) { int32_t sw = s.pop[i].sw_id; if (sw >= 0) s.sw_mem[s.sw_off[sw] + sw_cur[sw]++] = i; }
}
