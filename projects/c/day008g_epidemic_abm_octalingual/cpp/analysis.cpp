#include "common.hpp"
#include <cmath>
#include <cstdio>

static FILE *open_out(const char *dir, const char *name) {
    char p[512]; snprintf(p, sizeof(p), "%s/%s", dir, name);
    return fopen(p, "w");
}

static RunParams make_rp(int64_t N, int64_t days, int64_t beta_bp, int32_t interv_on) {
    RunParams rp; rp.N = N; rp.days = days; rp.beta_bp = beta_bp; rp.interv_on = interv_on;
    return rp;
}

static void net_degrees(const Sim &s, double &c_hh, double &c_sw, double &c_comm) {
    double shh = 0, ssw = 0, scom = 0;
    for (int64_t i = 0; i < s.N; ++i) {
        int32_t h = s.pop[i].hh_id;
        shh += static_cast<double>(s.hh_off[h + 1] - s.hh_off[h] - 1);
        if (s.pop[i].age_grp == 0) ssw += 6.0;
        else if (s.pop[i].age_grp == 1) ssw += 5.0;
        scom += s.pop[i].super_flag ? 12.0 : 3.0;
    }
    c_hh = shh / s.N; c_sw = ssw / s.N; c_comm = scom / s.N;
}

void write_s01_to_s09(const Config &c, const char *out_dir) {
    uint64_t hrng = 0xcbf29ce484222325ULL;
    const uint32_t pairs[3][2] = {{0,0},{1,2},{2,1}};
    const uint64_t exp_hv[3] = {0xbe40b757fe0b6697ULL, 0x659394d6f2a43a95ULL, 0x3abda5bac0dd49cbULL};
    int rng_ok = 1;
    for (int i = 0; i < 3; ++i) {
        uint64_t x = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1]);
        uint64_t sd = splitmix64(&x);
        Xoshiro256 st; x256_seed(&st, sd);
        uint64_t vals[4]; for (int j = 0; j < 4; ++j) vals[j] = x256_next(&st);
        uint64_t hv = fnv1a_bytes(vals, sizeof(vals));
        if (hv != exp_hv[i]) rng_ok = 0;
        hrng = fnv1a_i64(hrng, static_cast<int64_t>(hv));
    }
    Sim st0; sim_init(st0, c, make_rp(200, 10, 450, 0), 0); sim_run(st0);
    int inv_ok = st0.inv_violations == 0;
    FILE *f = open_out(out_dir, "s01_selftest.txt");
    fprintf(f, "ok,h_rng,inv_ok,tot_trans\n%d,%016llx,%d,%lld\n", rng_ok && inv_ok, static_cast<unsigned long long>(hrng), inv_ok, static_cast<long long>(st0.total_trans));
    fclose(f);

    f = open_out(out_dir, "s02_config.csv");
    fprintf(f, "key,value\nmaster_seed,%lld\nparity_N,%lld\nparity_days,%lld\nmain_days,%lld\nbeta_bp_init,%lld\n",
            static_cast<long long>(c.master_seed), static_cast<long long>(c.parity_N), static_cast<long long>(c.parity_days), static_cast<long long>(c.main_days), static_cast<long long>(c.beta_bp_init));
    fclose(f);

    Sim sp; sim_init(sp, c, make_rp(1000, 1, 450, 0), 0);
    f = open_out(out_dir, "s03_population.csv"); fprintf(f, "age_grp,count\n0,%lld\n1,%lld\n2,%lld\n", static_cast<long long>(sp.age_tot[0]), static_cast<long long>(sp.age_tot[1]), static_cast<long long>(sp.age_tot[2])); fclose(f);
    f = open_out(out_dir, "s03_households.csv"); fprintf(f, "hh_id,size\n");
    for (int32_t h = 0; h < sp.n_hh && h < 50; ++h) fprintf(f, "%d,%d\n", h, sp.hh_off[h + 1] - sp.hh_off[h]); fclose(f);
    f = open_out(out_dir, "s03_schools_work.csv"); fprintf(f, "sw_id,size\n");
    for (int32_t w = 0; w < sp.n_sw && w < 50; ++w) fprintf(f, "%d,%d\n", w, sp.sw_off[w + 1] - sp.sw_off[w]); fclose(f);
    double chh, csw, ccom; net_degrees(sp, chh, csw, ccom);
    f = open_out(out_dir, "s04_network.csv"); fprintf(f, "n_hh,n_sw,deg_hh,deg_sw,deg_comm\n%d,%d,%.4f,%.4f,%.4f\n", sp.n_hh, sp.n_sw, chh, csw, ccom); fclose(f);

    f = open_out(out_dir, "s05_contact_script.csv"); fprintf(f, "case_id,vax,mask,quar,eff_bp\n");
    for (int k = 0; k < 12; ++k) {
        int vax = k % 3, mask = (k / 3) % 2, quar = (k / 6) % 2;
        int64_t p = quar ? 0 : 450;
        if (mask) p = (p * (BP - c.mask_eff_bp)) / BP;
        if (vax == 1) p = (p * (BP - c.vax_eff1_bp)) / BP;
        else if (vax == 2) p = (p * (BP - c.vax_eff2_bp)) / BP;
        fprintf(f, "%d,%d,%d,%d,%lld\n", k, vax, mask, quar, static_cast<long long>(p));
    }
    fclose(f);

    f = open_out(out_dir, "s06_progression_script.csv"); fprintf(f, "case_id,age_grp,hosp_rej,p_hosp_bp,p_die_bp\n");
    for (int k = 0; k < 12; ++k) {
        int ag = k % 3, rej = (k / 3) % 2;
        int64_t ph = ag == 0 ? 120 : (ag == 1 ? 450 : 1800);
        int64_t pd = ag == 0 ? 100 : (ag == 1 ? 500 : 2000);
        if (rej) pd *= 2;
        fprintf(f, "%d,%d,%d,%lld,%lld\n", k, ag, rej, static_cast<long long>(ph), static_cast<long long>(pd));
    }
    fclose(f);

    f = open_out(out_dir, "s07_interventions.csv");
    fprintf(f, "name,start_day,param_bp\nschool_close,%lld,10000\nlockdown,%lld,5000\nmask,%lld,%lld\ntest,%lld,%lld\nvax,%lld,%lld\nvariant,%lld,%lld\n",
            static_cast<long long>(c.school_close_day), static_cast<long long>(c.lockdown_day), static_cast<long long>(c.mask_day), static_cast<long long>(c.mask_eff_bp),
            static_cast<long long>(c.test_day), static_cast<long long>(c.test_prob_bp), static_cast<long long>(c.vax_day), static_cast<long long>(c.vax_daily_bp),
            static_cast<long long>(c.variant_day), static_cast<long long>(c.variant_mult_bp));
    fclose(f);

    f = open_out(out_dir, "s08_calibration.csv"); fprintf(f, "beta_bp,r0_analytic,attack_rate,calibrated\n");
    int64_t b_vals[5] = {300, 380, 450, 520, 600};
    for (int bi = 0; bi < 5; ++bi) {
        Sim sc; sim_init(sc, c, make_rp(1000, 120, b_vals[bi], 0), 0); sim_run(sc);
        double ch, cs2, cc2; net_degrees(sc, ch, cs2, cc2);
        double r0 = analytical_r0(c, b_vals[bi], ch, cs2, cc2);
        double ar = static_cast<double>(sc.total_trans) / sc.N;
        fprintf(f, "%lld,%.4f,%.4f,%d\n", static_cast<long long>(b_vals[bi]), r0, ar, (r0 >= 1.5 && r0 <= 3.5 && ar >= 0.20 && ar <= 0.80) ? 1 : 0);
    }
    fclose(f);

    Sim s9; sim_init(s9, c, make_rp(500, 30, 450, 0), 0); sim_run(s9);
    f = open_out(out_dir, "s09_baseline_small.csv");
    fprintf(f, "N,days,total_trans,deaths,inv_ok\n500,30,%lld,%lld,%d\n", static_cast<long long>(s9.total_trans), static_cast<long long>(s9.daily_cnt[29][ST_D]), s9.inv_violations == 0);
    fclose(f);
}

static void write_daily_csv(const Sim &s, const char *out_dir, const char *fname) {
    FILE *f = open_out(out_dir, fname); fprintf(f, "day,S,E,Ia,Is,H,ICU,R,D\n");
    for (int32_t d = 0; d < s.rp.days; ++d) {
        fprintf(f, "%d", d);
        for (int k = 0; k < 8; ++k) fprintf(f, ",%lld", static_cast<long long>(s.daily_cnt[d][k]));
        fprintf(f, "\n");
    }
    fclose(f);
}

static void write_trans_csv(const Sim &s, const char *out_dir, const char *fname) {
    FILE *f = open_out(out_dir, fname); fprintf(f, "day,src_id,dst_id,layer,variant,gen_int\n");
    for (const auto &r : s.trans) fprintf(f, "%d,%d,%d,%d,%d,%d\n", r.day, r.src_id, r.dst_id, r.layer, r.variant, r.gen_int);
    fclose(f);
}

void write_s10_to_s37(Sim &s_no, Sim &s_main, const Config &c, const char *out_dir, int64_t beta_cal) {
    int32_t D = static_cast<int32_t>(s_main.rp.days);
    write_daily_csv(s_no, out_dir, "daily_counts_nointerv.csv");
    write_trans_csv(s_no, out_dir, "transmissions_log_nointerv.csv");
    write_daily_csv(s_main, out_dir, "daily_counts.csv");
    write_trans_csv(s_main, out_dir, "transmissions_log.csv");

    FILE *f = open_out(out_dir, "incidence_by_age.csv"); fprintf(f, "day,age_0_17,age_18_59,age_60p\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld,%lld\n", d, static_cast<long long>(s_main.daily_inc_age[d][0]), static_cast<long long>(s_main.daily_inc_age[d][1]), static_cast<long long>(s_main.daily_inc_age[d][2])); fclose(f);
    f = open_out(out_dir, "hospital_daily.csv"); fprintf(f, "day,hosp_occ,icu_occ,hosp_rej_cum,icu_rej_cum\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld,%lld,%lld\n", d, static_cast<long long>(s_main.daily_hosp[d][0]), static_cast<long long>(s_main.daily_hosp[d][1]), static_cast<long long>(s_main.daily_hosp[d][2]), static_cast<long long>(s_main.daily_hosp[d][3])); fclose(f);
    f = open_out(out_dir, "deaths_by_age.csv"); fprintf(f, "day,d_0_17,d_18_59,d_60p\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld,%lld\n", d, static_cast<long long>(s_main.daily_deaths_age[d][0]), static_cast<long long>(s_main.daily_deaths_age[d][1]), static_cast<long long>(s_main.daily_deaths_age[d][2])); fclose(f);
    f = open_out(out_dir, "tests.csv"); fprintf(f, "day,tests,positives\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld\n", d, static_cast<long long>(s_main.daily_tests[d][0]), static_cast<long long>(s_main.daily_tests[d][1])); fclose(f);
    f = open_out(out_dir, "quarantine.csv"); fprintf(f, "day,quarantined\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld\n", d, static_cast<long long>(s_main.daily_quar[d])); fclose(f);
    f = open_out(out_dir, "vaccination.csv"); fprintf(f, "day,dose1,dose2\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld\n", d, static_cast<long long>(s_main.daily_vax[d][0]), static_cast<long long>(s_main.daily_vax[d][1])); fclose(f);
    f = open_out(out_dir, "variant_share.csv"); fprintf(f, "day,v0_inf,v1_inf\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld\n", d, static_cast<long long>(s_main.daily_var[d][0]), static_cast<long long>(s_main.daily_var[d][1])); fclose(f);
    f = open_out(out_dir, "generation_interval.csv"); fprintf(f, "idx,gen_int\n");
    for (size_t i = 0; i < s_main.trans.size() && i < 2000; ++i) fprintf(f, "%zu,%d\n", i, s_main.trans[i].gen_int); fclose(f);
    f = open_out(out_dir, "agents_final.csv"); fprintf(f, "id,age_grp,status,vax_doses,offspring,ever_inf\n");
    for (int64_t i = 0; i < s_main.N && i < 2000; ++i) fprintf(f, "%lld,%d,%d,%d,%d,%d\n", static_cast<long long>(i), s_main.pop[i].age_grp, s_main.pop[i].status, s_main.pop[i].vax_doses, s_main.pop[i].offspring, s_main.pop[i].ever_inf); fclose(f);
    f = open_out(out_dir, "progression_log.csv"); fprintf(f, "day,pid,from_st,to_st\n");
    for (const auto &r : s_main.progs) fprintf(f, "%d,%d,%d,%d\n", r.day, r.pid, r.from_st, r.to_st); fclose(f);

    for (int a = 0; a < 3; ++a) {
        char fn[64]; snprintf(fn, sizeof(fn), "daily_age_%d.csv", a);
        f = open_out(out_dir, fn); fprintf(f, "day,inc,deaths\n");
        for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld\n", d, static_cast<long long>(s_main.daily_inc_age[d][a]), static_cast<long long>(s_main.daily_deaths_age[d][a])); fclose(f);
    }
    const char *lnames[3] = {"layer_hh.csv", "layer_sw.csv", "layer_comm.csv"};
    for (int l = 0; l < 3; ++l) { f = open_out(out_dir, lnames[l]); fprintf(f, "layer,infections\n%d,%lld\n", l, static_cast<long long>(s_main.layer_inf[l])); fclose(f); }
    f = open_out(out_dir, "hospital_cap_nointerv.csv"); fprintf(f, "day,hosp_occ,icu_occ\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld,%lld\n", d, static_cast<long long>(s_no.daily_hosp[d][0]), static_cast<long long>(s_no.daily_hosp[d][1])); fclose(f);
    f = open_out(out_dir, "deaths_nointerv.csv"); fprintf(f, "day,deaths\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld\n", d, static_cast<long long>(s_no.daily_cnt[d][ST_D])); fclose(f);
    int64_t oh_bins[6] = {}; for (int64_t i = 0; i < s_no.N; ++i) if (s_no.pop[i].ever_inf) { int o = s_no.pop[i].offspring; oh_bins[o > 5 ? 5 : o]++; }
    f = open_out(out_dir, "offspring_hist.csv"); fprintf(f, "bin,count\n"); for (int b = 0; b < 6; ++b) fprintf(f, "%d,%lld\n", b, static_cast<long long>(oh_bins[b])); fclose(f);
    f = open_out(out_dir, "rt_nointerv.csv"); fprintf(f, "day,inc\n");
    for (int32_t d = 0; d < D; ++d) fprintf(f, "%d,%lld\n", d, static_cast<long long>(s_no.daily_inc_age[d][0]+s_no.daily_inc_age[d][1]+s_no.daily_inc_age[d][2])); fclose(f);

    int64_t tot_viol = s_no.inv_violations + s_main.inv_violations;
    f = open_out(out_dir, "s12_invariants.csv"); fprintf(f, "violations,ok\n%lld,%d\n", static_cast<long long>(tot_viol), tot_viol == 0); fclose(f);

    f = open_out(out_dir, "s13_rt.csv"); fprintf(f, "day,rt\n");
    for (int32_t d = 0; d < D; ++d) {
        double num = 0, den = 0;
        for (int32_t k = 0; k < 5 && d - k >= 0; ++k) num += static_cast<double>(s_main.daily_inc_age[d - k][0] + s_main.daily_inc_age[d - k][1] + s_main.daily_inc_age[d - k][2]);
        for (int32_t k = 4; k < 9 && d - k >= 0; ++k) den += static_cast<double>(s_main.daily_inc_age[d - k][0] + s_main.daily_inc_age[d - k][1] + s_main.daily_inc_age[d - k][2]);
        fprintf(f, "%d,%.4f\n", d, den > 0.0 ? num / den : 0.0);
    }
    fclose(f);

    double chh, csw, ccom; net_degrees(s_no, chh, csw, ccom);
    double r0_ana = analytical_r0(c, beta_cal, chh, csw, ccom);
    double i_early = static_cast<double>(s_no.daily_inc_age[3][0]+s_no.daily_inc_age[3][1]+s_no.daily_inc_age[3][2] + s_no.daily_inc_age[4][0]+s_no.daily_inc_age[4][1]+s_no.daily_inc_age[4][2] + 1) / 2.0;
    double i_mid = 0; for (int d = 10; d < 16 && d < D; ++d) i_mid += static_cast<double>(s_no.daily_inc_age[d][0]+s_no.daily_inc_age[d][1]+s_no.daily_inc_age[d][2]); i_mid /= 6.0;
    double r_grow = i_mid > i_early ? std::log(i_mid / i_early) / 9.0 : 0.0;
    f = open_out(out_dir, "s14_growth.csv"); fprintf(f, "growth_r,doubling_days\n%.6f,%.4f\n", r_grow, r_grow > 0.0 ? std::log(2.0) / r_grow : 0.0); fclose(f);
    double r0_emp = 1.0 + r_grow * 4.5;
    f = open_out(out_dir, "s15_r0.csv"); fprintf(f, "r0_analytic,r0_empirical,ratio\n%.4f,%.4f,%.4f\n", r0_ana, r0_emp, r0_ana > 0 ? r0_emp / r0_ana : 0.0); fclose(f);

    double z_ana = solve_final_size(r0_ana), herd = r0_ana > 1.0 ? 1.0 - 1.0 / r0_ana : 0.0;
    f = open_out(out_dir, "s16_final_size.csv"); fprintf(f, "z_analytic,herd_thresh,attack_nointerv,attack_interv\n%.4f,%.4f,%.4f,%.4f\n", z_ana, herd, static_cast<double>(s_no.total_trans)/s_no.N, static_cast<double>(s_main.total_trans)/s_main.N); fclose(f);

    f = open_out(out_dir, "s17_attack.csv");
    fprintf(f, "category,count,frac\nHH,%lld,%.4f\nSW,%lld,%.4f\nCOMM,%lld,%.4f\nAGE0,%lld,%.4f\nAGE1,%lld,%.4f\nAGE2,%lld,%.4f\n",
            static_cast<long long>(s_main.layer_inf[0]), static_cast<double>(s_main.layer_inf[0])/s_main.N,
            static_cast<long long>(s_main.layer_inf[1]), static_cast<double>(s_main.layer_inf[1])/s_main.N,
            static_cast<long long>(s_main.layer_inf[2]), static_cast<double>(s_main.layer_inf[2])/s_main.N,
            static_cast<long long>(s_main.age_inf[0]), s_main.age_tot[0]>0 ? static_cast<double>(s_main.age_inf[0])/s_main.age_tot[0] : 0.0,
            static_cast<long long>(s_main.age_inf[1]), s_main.age_tot[1]>0 ? static_cast<double>(s_main.age_inf[1])/s_main.age_tot[1] : 0.0,
            static_cast<long long>(s_main.age_inf[2]), s_main.age_tot[2]>0 ? static_cast<double>(s_main.age_inf[2])/s_main.age_tot[2] : 0.0);
    fclose(f);

    double m1 = 0, m2 = 0; int64_t n_inf = 0;
    for (int64_t i = 0; i < s_no.N; ++i) if (s_no.pop[i].ever_inf) { double o = s_no.pop[i].offspring; m1 += o; m2 += o * o; n_inf++; }
    if (n_inf > 0) { m1 /= n_inf; m2 /= n_inf; }
    double var_o = m2 - m1 * m1;
    f = open_out(out_dir, "s18_offspring.csv"); fprintf(f, "n_inf,mean_offspring,var_offspring,k_disp\n%lld,%.4f,%.4f,%.4f\n", static_cast<long long>(n_inf), m1, var_o, (var_o > m1 && m1 > 0) ? (m1*m1)/(var_o - m1) : 10.0); fclose(f);

    int64_t pk_h_no = 0, pk_i_no = 0, pk_h_m = 0, pk_i_m = 0;
    for (int d = 0; d < D; ++d) {
        pk_h_no = std::max(pk_h_no, s_no.daily_hosp[d][0]); pk_i_no = std::max(pk_i_no, s_no.daily_hosp[d][1]);
        pk_h_m = std::max(pk_h_m, s_main.daily_hosp[d][0]); pk_i_m = std::max(pk_i_m, s_main.daily_hosp[d][1]);
    }
    f = open_out(out_dir, "s19_hospital.csv");
    fprintf(f, "hosp_cap,icu_cap,peak_h_no,peak_i_no,rej_no,peak_h_main,peak_i_main,rej_main\n%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld\n",
            static_cast<long long>(s_main.hosp_cap), static_cast<long long>(s_main.icu_cap), static_cast<long long>(pk_h_no), static_cast<long long>(pk_i_no),
            static_cast<long long>(s_no.branches[BR_HOSP_REJ] + s_no.branches[BR_ICU_REJ]), static_cast<long long>(pk_h_m), static_cast<long long>(pk_i_m),
            static_cast<long long>(s_main.branches[BR_HOSP_REJ] + s_main.branches[BR_ICU_REJ]));
    fclose(f);

    f = open_out(out_dir, "s20_ifr.csv"); fprintf(f, "age_grp,infections,deaths,ifr\n");
    for (int a = 0; a < 3; ++a) {
        int64_t dtot = 0; for (int d = 0; d < D; ++d) dtot += s_main.daily_deaths_age[d][a];
        fprintf(f, "%d,%lld,%lld,%.6f\n", a, static_cast<long long>(s_main.age_inf[a]), static_cast<long long>(dtot), s_main.age_inf[a] > 0 ? static_cast<double>(dtot) / s_main.age_inf[a] : 0.0);
    }
    fclose(f);

    int64_t t_tot = 0, t_pos = 0; for (int d = 0; d < D; ++d) { t_tot += s_main.daily_tests[d][0]; t_pos += s_main.daily_tests[d][1]; }
    f = open_out(out_dir, "s21_testing.csv"); fprintf(f, "tests,positives,quar_traces\n%lld,%lld,%lld\n", static_cast<long long>(t_tot), static_cast<long long>(t_pos), static_cast<long long>(s_main.branches[BR_QUAR_TRACE])); fclose(f);

    int64_t vn[3] = {}, vi[3] = {};
    for (int64_t i = 0; i < s_main.N; ++i) { int vd = std::min(s_main.pop[i].vax_doses, 2); vn[vd]++; if (s_main.pop[i].ever_inf) vi[vd]++; }
    f = open_out(out_dir, "s22_vaccine.csv"); fprintf(f, "doses,n_people,n_inf,attack_rate\n");
    for (int d = 0; d < 3; ++d) fprintf(f, "%d,%lld,%lld,%.6f\n", d, static_cast<long long>(vn[d]), static_cast<long long>(vi[d]), vn[d] > 0 ? static_cast<double>(vi[d]) / vn[d] : 0.0); fclose(f);

    int64_t v0 = 0, v1 = 0; for (int d = 0; d < D; ++d) { v0 += s_main.daily_var[d][0]; v1 += s_main.daily_var[d][1]; }
    f = open_out(out_dir, "s23_variant.csv"); fprintf(f, "v0_total,v1_total,v1_share\n%lld,%lld,%.6f\n", static_cast<long long>(v0), static_cast<long long>(v1), (v0+v1)>0 ? static_cast<double>(v1)/(v0+v1) : 0.0); fclose(f);
    f = open_out(out_dir, "s24_waning.csv"); fprintf(f, "waning_events,reinfections\n%lld,%lld\n", static_cast<long long>(s_main.branches[BR_WANING]), static_cast<long long>(s_main.total_reinf)); fclose(f);

    RunParams rnc1 = make_rp(500, 30, beta_cal, 0); rnc1.contact_scale_bp = 0; Sim nc1; sim_init(nc1, c, rnc1, 0); sim_run(nc1);
    RunParams rnc2 = make_rp(500, 30, 0, 0); rnc2.contact_scale_bp = -1; Sim nc2; sim_init(nc2, c, rnc2, 0); sim_run(nc2);
    RunParams rnc3 = make_rp(500, 30, beta_cal, 1); rnc3.vax_perfect = 1; Sim nc3; sim_init(nc3, c, rnc3, 0); sim_run(nc3);
    f = open_out(out_dir, "s25_negative_controls.csv");
    fprintf(f, "control,infections,ok\nzero_contact,%lld,%d\nzero_seed,%lld,%d\nperfect_vax,%lld,%d\n",
            static_cast<long long>(nc1.total_trans), nc1.total_trans == 0, static_cast<long long>(nc2.total_trans), nc2.total_trans == 0, static_cast<long long>(nc3.total_trans), nc3.total_trans == 0);
    fclose(f);

    int64_t Nsw = 1000, Dsw = 90;
    f = open_out(out_dir, "s26_sweep_comp_test.csv"); fprintf(f, "comp_bp,test_bp,rep,attack_rate\n");
    for (int64_t cv : {2000, 6000, 9500}) for (int64_t tv : {1000, 5000, 9000}) for (int rp = 0; rp < 3; ++rp) {
        Config c2 = c; c2.master_seed += rp + 1; RunParams rp2 = make_rp(Nsw, Dsw, beta_cal, 1); rp2.compliance_ovr = cv; rp2.test_ovr = tv;
        Sim ss; sim_init(ss, c2, rp2, 0); sim_run(ss); fprintf(f, "%lld,%lld,%d,%.6f\n", static_cast<long long>(cv), static_cast<long long>(tv), rp, static_cast<double>(ss.total_trans)/Nsw);
    }
    fclose(f);
    f = open_out(out_dir, "s27_sweep_vax.csv"); fprintf(f, "vax_day,vax_cov_bp,rep,attack_rate\n");
    for (int64_t vd : {10, 45, 80}) for (int64_t vc : {2000, 5500, 9000}) for (int rp = 0; rp < 3; ++rp) {
        Config c2 = c; c2.master_seed += rp + 10; RunParams rp2 = make_rp(Nsw, Dsw, beta_cal, 1); rp2.vax_day_ovr = vd; rp2.vax_cov_ovr = vc;
        Sim ss; sim_init(ss, c2, rp2, 0); sim_run(ss); fprintf(f, "%lld,%lld,%d,%.6f\n", static_cast<long long>(vd), static_cast<long long>(vc), rp, static_cast<double>(ss.total_trans)/Nsw);
    }
    fclose(f);
    f = open_out(out_dir, "s28_sweep_school_lock.csv"); fprintf(f, "school_day,lock_day,rep,attack_rate\n");
    for (int64_t sc : {15, 200}) for (int64_t lk : {15, 45, 200}) for (int rp = 0; rp < 3; ++rp) {
        Config c2 = c; c2.master_seed += rp + 20; RunParams rp2 = make_rp(Nsw, Dsw, beta_cal, 1); rp2.school_ovr = sc; rp2.lock_ovr = lk;
        Sim ss; sim_init(ss, c2, rp2, 0); sim_run(ss); fprintf(f, "%lld,%lld,%d,%.6f\n", static_cast<long long>(sc), static_cast<long long>(lk), rp, static_cast<double>(ss.total_trans)/Nsw);
    }
    fclose(f);
    f = open_out(out_dir, "s29_sweep_variant_mask.csv"); fprintf(f, "var_day,mask_day,rep,attack_rate\n");
    for (int64_t vr : {15, 45, 200}) for (int64_t mk : {10, 40, 200}) for (int rp = 0; rp < 3; ++rp) {
        Config c2 = c; c2.master_seed += rp + 30; RunParams rp2 = make_rp(Nsw, Dsw, beta_cal, 1); rp2.var_day_ovr = vr; rp2.mask_ovr = mk;
        Sim ss; sim_init(ss, c2, rp2, 0); sim_run(ss); fprintf(f, "%lld,%lld,%d,%.6f\n", static_cast<long long>(vr), static_cast<long long>(mk), rp, static_cast<double>(ss.total_trans)/Nsw);
    }
    fclose(f);

    double ens_sum[90] = {}, ens_sum2[90] = {};
    for (int rp = 0; rp < 20; ++rp) {
        Config c2 = c; c2.master_seed += rp + 100;
        Sim ss; sim_init(ss, c2, make_rp(Nsw, Dsw, beta_cal, 1), 0); sim_run(ss);
        for (int d = 0; d < 90; ++d) { double inc = static_cast<double>(ss.daily_inc_age[d][0] + ss.daily_inc_age[d][1] + ss.daily_inc_age[d][2]); ens_sum[d] += inc; ens_sum2[d] += inc * inc; }
    }
    f = open_out(out_dir, "s30_ensemble.csv"); fprintf(f, "day,mean_inc,se_inc\n");
    for (int d = 0; d < 90; ++d) { double mn = ens_sum[d]/20.0; fprintf(f, "%d,%.4f,%.4f\n", d, mn, std::sqrt(std::fmax(0.0, (ens_sum2[d]/20.0 - mn*mn)/19.0))); }
    fclose(f);

    const char *bn[NUM_BRANCHES] = {"INF_HH","INF_SW","INF_COMM","FAIL_VAX","FAIL_MASK","FAIL_QUAR","PROG_IA","PROG_IS","HOSP_ADM","HOSP_REJ","ICU_ADM","ICU_REJ","DEATH","RECOVER","WANING","TEST_POS","TEST_NEG","QUAR_TRACE","VAX_DOSE","VAX_OUT","VAR_IMPORT"};
    f = open_out(out_dir, "s31_branches.csv"); fprintf(f, "branch,count\n");
    for (int i = 0; i < NUM_BRANCHES; ++i) fprintf(f, "%s,%lld\n", bn[i], static_cast<long long>(s_main.branches[i])); fclose(f);

    f = open_out(out_dir, "s32_vacuous_check.csv");
    int ens_nz = 0; for (int d = 0; d < 90; ++d) if (ens_sum[d] > 0.0) ens_nz++;
    fprintf(f, "check,ok\nneg_controls,%d\nensemble_var,%d\n", (nc1.total_trans==0 && nc2.total_trans==0 && nc3.total_trans==0), ens_nz > 5); fclose(f);

    Sim sp1, sp2; sim_init(sp1, c, make_rp(c.parity_N, c.parity_days, beta_cal, 1), 0); sim_run(sp1);
    sim_init(sp2, c, make_rp(c.parity_N, c.parity_days, beta_cal, 1), 0); sim_run(sp2);
    uint64_t h1 = 0xcbf29ce484222325ULL, h2 = 0xcbf29ce484222325ULL;
    for (int64_t i = 0; i < sp1.N; ++i) {
        h1 = fnv1a_i64(h1, sp1.pop[i].status); h2 = fnv1a_i64(h2, sp2.pop[i].status);
        h1 = fnv1a_i64(h1, sp1.pop[i].offspring); h2 = fnv1a_i64(h2, sp2.pop[i].offspring);
    }
    f = open_out(out_dir, "s33_replay.csv"); fprintf(f, "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016llx,%016llx,%016llx,%016llx,%d\n", static_cast<unsigned long long>(sp1.event_hash), static_cast<unsigned long long>(sp2.event_hash), static_cast<unsigned long long>(h1), static_cast<unsigned long long>(h2), sp1.event_hash == sp2.event_hash && h1 == h2); fclose(f);
    f = open_out(out_dir, "s34_parity.csv"); fprintf(f, "parity_N,parity_days,ev_hash,fin_hash,total_trans,deaths,inv_violations\n%lld,%lld,%016llx,%016llx,%lld,%lld,%lld\n", static_cast<long long>(c.parity_N), static_cast<long long>(c.parity_days), static_cast<unsigned long long>(sp1.event_hash), static_cast<unsigned long long>(h1), static_cast<long long>(sp1.total_trans), static_cast<long long>(sp1.daily_cnt[c.parity_days - 1][ST_D]), static_cast<long long>(sp1.inv_violations)); fclose(f);

    f = open_out(out_dir, "s36_sizes.csv"); fprintf(f, "N,days,n_trans,n_progs\n%lld,%lld,%zu,%zu\n", static_cast<long long>(s_main.N), static_cast<long long>(s_main.rp.days), s_main.trans.size(), s_main.progs.size()); fclose(f);

    int64_t cnt_r[8] = {s_main.N, 0, 0, 0, 0, 0, 0, 0}, max_diff = 0; size_t pi = 0;
    for (int32_t d = 0; d < D; ++d) {
        while (pi < s_main.progs.size() && s_main.progs[pi].day == d) { cnt_r[s_main.progs[pi].from_st]--; cnt_r[s_main.progs[pi].to_st]++; pi++; }
        for (int k = 0; k < 8; ++k) max_diff = std::max(max_diff, std::llabs(cnt_r[k] - s_main.daily_cnt[d][k]));
    }
    f = open_out(out_dir, "s37_replay_counts.csv"); fprintf(f, "max_diff,match\n%lld,%d\n", static_cast<long long>(max_diff), max_diff == 0); fclose(f);
}
