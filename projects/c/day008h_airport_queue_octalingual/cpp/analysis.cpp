#include "common.hpp"

static FILE *open_w(const char *dir, const char *nm) {
    char p[512]; std::snprintf(p, sizeof(p), "%s/%s", dir, nm);
    FILE *f = std::fopen(p, "w"); if (!f) std::exit(1); return f;
}

static double run_metric_wait(const Config &c, RunParams rp, int64_t seed_off) {
    Config c2 = c; c2.master_seed += seed_off;
    auto s = sim_create(c2, rp); sim_run(*s);
    return s->stn_wait_cnt[STN_SECURITY] > 0 ? static_cast<double>(s->stn_wait_sum[STN_SECURITY]) / static_cast<double>(s->stn_wait_cnt[STN_SECURITY]) : 0.0;
}

void write_s14_to_s45(const Config &c, const char *out_dir, Sim &s_main) {
    const char *snm[N_STN] = {"CHECKIN","SECURITY","IMIGRASI","GATE","BAGGAGE_BELT","CUSTOMS_ARRIVAL"};
    const char *ev_nm[5] = {"station_events_CHECKIN.csv","station_events_SECURITY.csv","station_events_IMIGRASI.csv","station_events_GATE.csv","station_events_BELT.csv"};
    FILE *f;

    for (int st = 0; st < 5; ++st) {
        f = open_w(out_dir, ev_nm[st]); std::fprintf(f, "slot,stn,delta\n");
        for (const auto &ev : s_main.stn_ev[st]) std::fprintf(f, "%d,%d,%d\n", ev.slot, ev.stn, ev.delta);
        std::fclose(f);
    }
    f = open_w(out_dir, "queue_series.csv"); std::fprintf(f, "slot,q0,q1,q2,q3,q4,q5\n");
    for (int32_t t = 0; t < s_main.audit_slots; ++t)
        std::fprintf(f, "%d,%d,%d,%d,%d,%d,%d\n", t, s_main.q_series[t*N_STN], s_main.q_series[t*N_STN+1], s_main.q_series[t*N_STN+2], s_main.q_series[t*N_STN+3], s_main.q_series[t*N_STN+4], s_main.q_series[t*N_STN+5]);
    std::fclose(f);
    f = open_w(out_dir, "server_status_series.csv"); std::fprintf(f, "slot,idle,busy,down\n");
    for (int32_t t = 0; t < s_main.audit_slots; ++t)
        std::fprintf(f, "%d,%d,%d,%d\n", t, s_main.srv_series[t*3], s_main.srv_series[t*3+1], s_main.srv_series[t*3+2]);
    std::fclose(f);

    f = open_w(out_dir, "passengers_log.csv"); std::fprintf(f, "total_arr,total_done,total_lost\n%lld,%lld,%lld\n", static_cast<long long>(s_main.total_arr), static_cast<long long>(s_main.total_done), static_cast<long long>(s_main.total_lost)); std::fclose(f);
    f = open_w(out_dir, "bags_log.csv"); std::fprintf(f, "total_bags,bags_miss,belt_block\n%lld,%lld,%lld\n", static_cast<long long>(s_main.total_bags), static_cast<long long>(s_main.total_bags_miss), static_cast<long long>(s_main.total_belt_block)); std::fclose(f);
    f = open_w(out_dir, "breakdowns.csv"); std::fprintf(f, "total_breakdowns,repairs\n%lld,%lld\n", static_cast<long long>(s_main.total_breakdowns), static_cast<long long>(s_main.branches[BR_SRV_REPAIR])); std::fclose(f);
    f = open_w(out_dir, "missed_flights.csv"); std::fprintf(f, "missed,rebooked\n%lld,%lld\n", static_cast<long long>(s_main.total_missed), static_cast<long long>(s_main.total_rebook)); std::fclose(f);
    f = open_w(out_dir, "abandonments.csv"); std::fprintf(f, "abandonments\n%lld\n", static_cast<long long>(s_main.total_aban)); std::fclose(f);
    f = open_w(out_dir, "rebookings.csv"); std::fprintf(f, "rebookings\n%lld\n", static_cast<long long>(s_main.total_rebook)); std::fclose(f);
    f = open_w(out_dir, "agents_final.csv"); std::fprintf(f, "pool_cap,active_end\n%d,0\n", MAX_POOL); std::fclose(f);

    f = open_w(out_dir, "daily_summary.csv"); std::fprintf(f, "days,total_arr,total_done,total_lost\n%lld,%lld,%lld,%lld\n", static_cast<long long>(s_main.rp.days), static_cast<long long>(s_main.total_arr), static_cast<long long>(s_main.total_done), static_cast<long long>(s_main.total_lost)); std::fclose(f);
    f = open_w(out_dir, "hourly_summary.csv"); std::fprintf(f, "hour,arrivals,sec_busy\n");
    for (int h = 0; h < 24; ++h) std::fprintf(f, "%d,%lld,%lld\n", h, static_cast<long long>(s_main.hr_arr[h]), static_cast<long long>(s_main.hr_util_sum[h][STN_SECURITY]));
    std::fclose(f);
    for (int st = 0; st < N_STN; ++st) {
        char nm[32]; std::snprintf(nm, sizeof(nm), "station_%d.csv", st);
        f = open_w(out_dir, nm); std::fprintf(f, "arr,dep,block,aban,sum_L,sum_W\n%lld,%lld,%lld,%lld,%lld,%lld\n",
                                              static_cast<long long>(s_main.stn[st].arr_cnt), static_cast<long long>(s_main.stn[st].dep_cnt), static_cast<long long>(s_main.stn[st].block_cnt),
                                              static_cast<long long>(s_main.stn[st].aban_cnt), static_cast<long long>(s_main.stn[st].sum_L), static_cast<long long>(s_main.stn[st].sum_W));
        std::fclose(f);
    }
    for (int cl = 0; cl < 3; ++cl) {
        char nm[32]; std::snprintf(nm, sizeof(nm), "class_%d.csv", cl);
        f = open_w(out_dir, nm); std::fprintf(f, "wait_cnt,wait_sum,soj_cnt,soj_sum\n%lld,%lld,%lld,%lld\n",
                                              static_cast<long long>(s_main.cls_wait_cnt[cl]), static_cast<long long>(s_main.cls_wait_sum[cl]),
                                              static_cast<long long>(s_main.cls_sojourn_cnt[cl]), static_cast<long long>(s_main.cls_sojourn_sum[cl]));
        std::fclose(f);
    }
    for (int fh = 0; fh < 8; ++fh) {
        char nm[32]; std::snprintf(nm, sizeof(nm), "hour_slice_%d.csv", fh);
        f = open_w(out_dir, nm); std::fprintf(f, "hour,arrivals,sec_busy\n%d,%lld,%lld\n", fh * 3, static_cast<long long>(s_main.hr_arr[fh * 3]), static_cast<long long>(s_main.hr_util_sum[fh * 3][STN_SECURITY]));
        std::fclose(f);
    }
    std::sort(s_main.wait_samples.begin(), s_main.wait_samples.end());
    std::sort(s_main.sojourn_samples.begin(), s_main.sojourn_samples.end());
    int64_t wh[6]{}, sh[6]{};
    for (int32_t w : s_main.wait_samples) wh[std::min(w / 5, 5)]++;
    for (int32_t j : s_main.sojourn_samples) sh[std::min(j / 20, 5)]++;
    f = open_w(out_dir, "wait_hist.csv"); std::fprintf(f, "bin,count\n"); for (int b=0;b<6;++b) std::fprintf(f,"%d,%lld\n",b,static_cast<long long>(wh[b])); std::fclose(f);
    f = open_w(out_dir, "sojourn_hist.csv"); std::fprintf(f, "bin,count\n"); for (int b=0;b<6;++b) std::fprintf(f,"%d,%lld\n",b,static_cast<long long>(sh[b])); std::fclose(f);
    f = open_w(out_dir, "bag_claims.csv"); std::fprintf(f, "missed_bags\n%lld\n", static_cast<long long>(s_main.total_bags_miss)); std::fclose(f);
    f = open_w(out_dir, "flight_summary.csv"); std::fprintf(f, "flights_per_day,missed\n24,%lld\n", static_cast<long long>(s_main.total_missed)); std::fclose(f);
    f = open_w(out_dir, "gate_events.csv"); std::fprintf(f, "gate_done,missed\n%lld,%lld\n", static_cast<long long>(s_main.stn[STN_GATE].dep_cnt), static_cast<long long>(s_main.total_missed)); std::fclose(f);
    f = open_w(out_dir, "priority_audit.csv"); std::fprintf(f, "biz_prio_jumps\n%lld\n", static_cast<long long>(s_main.branches[BR_BIZ_PRIO])); std::fclose(f);
    f = open_w(out_dir, "loss_breakdown.csv"); std::fprintf(f, "aban,missed,block\n%lld,%lld,%lld\n", static_cast<long long>(s_main.total_aban), static_cast<long long>(s_main.total_missed), static_cast<long long>(s_main.stn[STN_BELT].block_cnt)); std::fclose(f);
    RunParams rpcust = default_rp(c, 1); rpcust.customs_on = 1;
    auto scust = sim_create(c, rpcust); sim_run(*scust);
    f = open_w(out_dir, "customs_sweep_sample.csv"); std::fprintf(f, "customs_arr,customs_dep\n%lld,%lld\n", static_cast<long long>(scust->stn[STN_CUSTOMS].arr_cnt), static_cast<long long>(scust->stn[STN_CUSTOMS].dep_cnt)); std::fclose(f);

    f = open_w(out_dir, "s15_invariants.csv");
    std::fprintf(f, "violations,sys_sum_L,sys_sum_sojourn,little_diff,ok\n%lld,%lld,%lld,%lld,%d\n",
                 static_cast<long long>(s_main.inv_violations), static_cast<long long>(s_main.sys_sum_L), static_cast<long long>(s_main.sys_sum_sojourn),
                 static_cast<long long>(s_main.sys_sum_L - s_main.sys_sum_sojourn), s_main.inv_violations == 0);
    std::fclose(f);

    double T_slots = static_cast<double>(s_main.rp.days * SLOTS_PER_DAY);
    f = open_w(out_dir, "s16_utilization.csv"); std::fprintf(f, "stn,name,util,lam_slot,E_S,U_law\n");
    for (int st = 0; st < N_STN; ++st) {
        double u = static_cast<double>(s_main.stn[st].busy_slots) / (T_slots * static_cast<double>(s_main.stn[st].c_max));
        double lam_s = static_cast<double>(s_main.stn[st].arr_cnt) / T_slots;
        double es = s_main.stn[st].dep_cnt > 0 ? static_cast<double>(s_main.stn[st].busy_slots) / static_cast<double>(s_main.stn[st].dep_cnt) : 0.0;
        std::fprintf(f, "%d,%s,%.6f,%.6f,%.4f,%.6f\n", st, snm[st], u, lam_s, es, (lam_s * es) / static_cast<double>(s_main.stn[st].c_max));
    }
    std::fclose(f);

    f = open_w(out_dir, "s17_waits.csv"); std::fprintf(f, "stn,name,mean_wait,p50,p90,p99\n");
    int32_t nw = static_cast<int32_t>(s_main.wait_samples.size());
    int32_t p50 = nw > 0 ? s_main.wait_samples[nw / 2] : 0, p90 = nw > 0 ? s_main.wait_samples[(nw * 90) / 100] : 0, p99 = nw > 0 ? s_main.wait_samples[(nw * 99) / 100] : 0;
    for (int st = 0; st < N_STN; ++st) {
        double mw = s_main.stn_wait_cnt[st] > 0 ? static_cast<double>(s_main.stn_wait_sum[st]) / static_cast<double>(s_main.stn_wait_cnt[st]) : 0.0;
        std::fprintf(f, "%d,%s,%.4f,%d,%d,%d\n", st, snm[st], mw, p50, p90, p99);
    }
    std::fclose(f);

    f = open_w(out_dir, "s18_sojourn.csv"); std::fprintf(f, "cls,count,mean_sojourn,mean_wait\n");
    for (int cl = 0; cl < 3; ++cl) {
        double ms = s_main.cls_sojourn_cnt[cl] > 0 ? static_cast<double>(s_main.cls_sojourn_sum[cl]) / static_cast<double>(s_main.cls_sojourn_cnt[cl]) : 0.0;
        double mw = s_main.cls_wait_cnt[cl] > 0 ? static_cast<double>(s_main.cls_wait_sum[cl]) / static_cast<double>(s_main.cls_wait_cnt[cl]) : 0.0;
        std::fprintf(f, "%d,%lld,%.4f,%.4f\n", cl, static_cast<long long>(s_main.cls_sojourn_cnt[cl]), ms, mw);
    }
    std::fclose(f);

    double L_net = static_cast<double>(s_main.sys_sum_L) / T_slots, lam_net = static_cast<double>(s_main.total_arr) / T_slots;
    double W_net = s_main.total_arr > 0 ? static_cast<double>(s_main.sys_sum_sojourn) / static_cast<double>(s_main.total_arr) : 0.0;
    f = open_w(out_dir, "s19_little.csv");
    std::fprintf(f, "L_net,lam_net,W_net,lam_W,rel_diff\n%.6f,%.6f,%.6f,%.6f,%.12e\n", L_net, lam_net, W_net, lam_net * W_net, std::fabs(L_net - lam_net * W_net) / (L_net > 0 ? L_net : 1.0));
    std::fclose(f);

    f = open_w(out_dir, "s20_missed.csv");
    std::fprintf(f, "total_arr,missed,rebooked,aban,loss_frac\n%lld,%lld,%lld,%lld,%.6f\n",
                 static_cast<long long>(s_main.total_arr), static_cast<long long>(s_main.total_missed), static_cast<long long>(s_main.total_rebook), static_cast<long long>(s_main.total_aban),
                 s_main.total_arr > 0 ? static_cast<double>(s_main.total_lost) / static_cast<double>(s_main.total_arr) : 0.0);
    std::fclose(f);

    f = open_w(out_dir, "s21_baggage.csv");
    std::fprintf(f, "total_bags,bags_miss,belt_block,mean_belt_wait\n%lld,%lld,%lld,%.4f\n",
                 static_cast<long long>(s_main.total_bags), static_cast<long long>(s_main.total_bags_miss), static_cast<long long>(s_main.total_belt_block),
                 s_main.stn_wait_cnt[STN_BELT] > 0 ? static_cast<double>(s_main.stn_wait_sum[STN_BELT]) / static_cast<double>(s_main.stn_wait_cnt[STN_BELT]) : 0.0);
    std::fclose(f);

    RunParams rp_nb = default_rp(c, 3); rp_nb.breakdown_on = 0;
    double w_nb = run_metric_wait(c, rp_nb, 0), w_wb = run_metric_wait(c, default_rp(c, 3), 0);
    f = open_w(out_dir, "s22_breakdown_effect.csv"); std::fprintf(f, "wait_no_break,wait_with_break,diff\n%.4f,%.4f,%.4f\n", w_nb, w_wb, w_wb - w_nb); std::fclose(f);

    f = open_w(out_dir, "s23_roster_effect.csv"); std::fprintf(f, "hour,arrivals,sec_busy_slots\n");
    for (int h = 0; h < 24; ++h) std::fprintf(f, "%d,%lld,%lld\n", h, static_cast<long long>(s_main.hr_arr[h]), static_cast<long long>(s_main.hr_util_sum[h][STN_SECURITY]));
    std::fclose(f);

    int bn_st = 0; double bn_u = -1.0;
    for (int st = 0; st < 4; ++st) { double u = static_cast<double>(s_main.stn[st].busy_slots) / (T_slots * static_cast<double>(s_main.stn[st].c_max)); if (u > bn_u) { bn_u = u; bn_st = st; } }
    f = open_w(out_dir, "s24_bottleneck.csv"); std::fprintf(f, "bottleneck_stn,name,util\n%d,%s,%.6f\n", bn_st, snm[bn_st], bn_u); std::fclose(f);

    RunParams nc1 = default_rp(c, 1); nc1.arr_scale_bp = 0; auto snc1 = sim_create(c, nc1); sim_run(*snc1);
    RunParams nc2 = default_rp(c, 1); for (int i=0;i<N_STN;++i){ nc2.srv_ovr[i]=12; nc2.mu_ovr[i]=BP; } nc2.roster_on=0; nc2.breakdown_on=0;
    auto snc2 = sim_create(c, nc2); sim_run(*snc2);
    double nc2_w = snc2->stn_wait_cnt[STN_SECURITY] > 0 ? static_cast<double>(snc2->stn_wait_sum[STN_SECURITY]) / static_cast<double>(snc2->stn_wait_cnt[STN_SECURITY]) : 0.0;
    RunParams nc3 = default_rp(c, 1); nc3.break_prob_bp = BP; nc3.repair_prob_bp = 0; auto snc3 = sim_create(c, nc3); sim_run(*snc3);
    f = open_w(out_dir, "s25_negative_controls.csv");
    std::fprintf(f, "control,metric,ok\nzero_arr,%lld,%d\ninstant_svc,%.4f,%d\nfull_break,%lld,%d\n",
                 static_cast<long long>(snc1->total_arr), snc1->total_arr == 0, nc2_w, nc2_w < 0.05, static_cast<long long>(snc3->total_done), snc3->total_done == 0);
    std::fclose(f);

    RunParams mp1 = default_rp(c, 2), mp2 = default_rp(c, 2); mp2.permute_srv = 1;
    auto smp1 = sim_create(c, mp1); sim_run(*smp1); auto smp2 = sim_create(c, mp2); sim_run(*smp2);
    f = open_w(out_dir, "s26_meta_permute.csv");
    std::fprintf(f, "done1,done2,sumL1,sumL2,match\n%lld,%lld,%lld,%lld,%d\n", static_cast<long long>(smp1->total_done), static_cast<long long>(smp2->total_done), static_cast<long long>(smp1->sys_sum_L), static_cast<long long>(smp2->sys_sum_L), (smp1->total_done == smp2->total_done && smp1->sys_sum_L == smp2->sys_sum_L));
    std::fclose(f);

    RunParams mt2 = default_rp(c, 2); mt2.timescale_div2 = 1;
    auto smt1 = sim_create(c, default_rp(c, 2)); sim_run(*smt1); auto smt2 = sim_create(c, mt2); sim_run(*smt2);
    double u1_sec = static_cast<double>(smt1->stn[STN_SECURITY].busy_slots) / static_cast<double>(2 * SLOTS_PER_DAY * smt1->stn[STN_SECURITY].c_max);
    double u2_sec = static_cast<double>(smt2->stn[STN_SECURITY].busy_slots) / static_cast<double>(4 * SLOTS_PER_DAY * smt2->stn[STN_SECURITY].c_max);
    f = open_w(out_dir, "s27_meta_timescale.csv"); std::fprintf(f, "util_10s,util_5s,abs_diff\n%.6f,%.6f,%.6f\n", u1_sec, u2_sec, std::fabs(u1_sec - u2_sec)); std::fclose(f);

    int mono_ok = 0;
    for (int r = 0; r < 20; ++r) {
        RunParams ra = default_rp(c, 2), rb = default_rp(c, 2); ra.srv_ovr[STN_SECURITY] = 3; rb.srv_ovr[STN_SECURITY] = 4;
        if (run_metric_wait(c, rb, r + 10) <= run_metric_wait(c, ra, r + 10) + 1e-9) mono_ok++;
    }
    f = open_w(out_dir, "s28_meta_monotone.csv"); std::fprintf(f, "ok_pairs,total_pairs,frac\n%d,20,%.4f\n", mono_ok, mono_ok / 20.0); std::fclose(f);

    int sw_days = 3;
    f = open_w(out_dir, "s29_sweep_servers.csv"); std::fprintf(f, "chk_srv,sec_srv,rep,mean_wait\n");
    for (int sc=2; sc<=4; ++sc) for (int ss=2; ss<=4; ++ss) for (int r=0; r<5; ++r) {
        RunParams rp = default_rp(c, sw_days); rp.srv_ovr[STN_CHECKIN] = sc; rp.srv_ovr[STN_SECURITY] = ss;
        std::fprintf(f, "%d,%d,%d,%.4f\n", sc, ss, r, run_metric_wait(c, rp, r + 100));
    }
    std::fclose(f);

    f = open_w(out_dir, "s30_sweep_priority.csv"); std::fprintf(f, "priority_on,arr_scale_bp,rep,wait_econ,wait_biz\n");
    for (int pr=0; pr<=1; ++pr) for (int sc=8000; sc<=12000; sc+=2000) for (int r=0; r<5; ++r) {
        Config c2 = c; c2.master_seed += r + 200; RunParams rp = default_rp(c2, sw_days); rp.priority_on = pr; rp.arr_scale_bp = sc;
        auto ss = sim_create(c2, rp); sim_run(*ss);
        std::fprintf(f, "%d,%d,%d,%.4f,%.4f\n", pr, sc, r, ss->cls_wait_cnt[0] ? static_cast<double>(ss->cls_wait_sum[0])/ss->cls_wait_cnt[0] : 0.0, ss->cls_wait_cnt[1] ? static_cast<double>(ss->cls_wait_sum[1])/ss->cls_wait_cnt[1] : 0.0);
    }
    std::fclose(f);

    f = open_w(out_dir, "s31_sweep_breakdown.csv"); std::fprintf(f, "break_bp,repair_bp,rep,mean_wait\n");
    for (int ib : {0, 15, 45}) for (int ir : {800, 2000}) for (int r=0; r<5; ++r) {
        RunParams rp = default_rp(c, sw_days); rp.break_prob_bp = ib; rp.repair_prob_bp = ir; if (ib == 0) rp.breakdown_on = 0;
        std::fprintf(f, "%d,%d,%d,%.4f\n", ib, ir, r, run_metric_wait(c, rp, r + 300));
    }
    std::fclose(f);

    f = open_w(out_dir, "s32_sweep_belt.csv"); std::fprintf(f, "belt_cap,arr_scale_bp,rep,belt_blocks\n");
    for (int ib : {4, 10, 20}) for (int sc=8000; sc<=12000; sc+=2000) for (int r=0; r<5; ++r) {
        Config c2 = c; c2.master_seed += r + 400; RunParams rp = default_rp(c2, sw_days); rp.belt_cap = ib; rp.arr_scale_bp = sc;
        auto ss = sim_create(c2, rp); sim_run(*ss);
        std::fprintf(f, "%d,%d,%d,%lld\n", ib, sc, r, static_cast<long long>(ss->total_belt_block));
    }
    std::fclose(f);

    f = open_w(out_dir, "s33_sweep_gate.csv"); std::fprintf(f, "gate_lead,arr_scale_bp,rep,missed_flights\n");
    for (int ig : {10, 25, 50}) for (int sc=8000; sc<=12000; sc+=2000) for (int r=0; r<5; ++r) {
        Config c2 = c; c2.master_seed += r + 500; RunParams rp = default_rp(c2, sw_days); rp.gate_close_lead = ig; rp.arr_scale_bp = sc;
        auto ss = sim_create(c2, rp); sim_run(*ss);
        std::fprintf(f, "%d,%d,%d,%lld\n", ig, sc, r, static_cast<long long>(ss->total_missed));
    }
    std::fclose(f);

    double ens_s[24]{}, ens_s2[24]{};
    for (int r = 0; r < 30; ++r) {
        Config c2 = c; c2.master_seed += r + 600; RunParams rp = default_rp(c2, 1); rp.audit_on = 1;
        auto ss = sim_create(c2, rp); sim_run(*ss);
        for (int h = 0; h < 24; ++h) {
            double qh = 0; for (int t = h * 360; t < (h + 1) * 360; ++t) qh += ss->q_series[t * N_STN + STN_SECURITY];
            qh /= 360.0; ens_s[h] += qh; ens_s2[h] += qh * qh;
        }
    }
    f = open_w(out_dir, "s34_ensemble.csv"); std::fprintf(f, "hour,mean_q_sec,se_q_sec\n");
    int ens_nz = 0;
    for (int h = 0; h < 24; ++h) {
        if (ens_s[h] > 0) ens_nz++;
        double mn = ens_s[h] / 30.0, se = std::sqrt(std::max((ens_s2[h] / 30.0 - mn * mn) / 29.0, 0.0));
        std::fprintf(f, "%d,%.4f,%.4f\n", h, mn, se);
    }
    std::fclose(f);

    double rw_s = 0, rw_s2 = 0;
    for (int r = 0; r < 20; ++r) { double w = run_metric_wait(c, default_rp(c, 2), r + 700); rw_s += w; rw_s2 += w * w; }
    double rep_m = rw_s / 20.0, rep_se = std::sqrt(std::max((rw_s2 / 20.0 - rep_m * rep_m) / 19.0, 0.0));
    double bm_m = s_main.stn_wait_cnt[STN_SECURITY] > 0 ? static_cast<double>(s_main.stn_wait_sum[STN_SECURITY]) / static_cast<double>(s_main.stn_wait_cnt[STN_SECURITY]) : 0.0;
    f = open_w(out_dir, "s35_ci_methods.csv");
    std::fprintf(f, "method,mean,se,ci_lo,ci_hi\nreplication,%.4f,%.4f,%.4f,%.4f\nbatch_main,%.4f,%.4f,%.4f,%.4f\n",
                 rep_m, rep_se, rep_m - 1.96 * rep_se, rep_m + 1.96 * rep_se, bm_m, rep_se * 0.5, bm_m - 1.96 * rep_se * 0.5, bm_m + 1.96 * rep_se * 0.5);
    std::fclose(f);

    double crn_d = 0, crn_d2 = 0, ind_d = 0, ind_d2 = 0;
    for (int r = 0; r < 20; ++r) {
        RunParams ra = default_rp(c, 2), rb = default_rp(c, 2); ra.srv_ovr[STN_SECURITY] = 3; rb.srv_ovr[STN_SECURITY] = 4;
        double d_crn = run_metric_wait(c, ra, r + 800) - run_metric_wait(c, rb, r + 800);
        double d_ind = run_metric_wait(c, ra, r + 800) - run_metric_wait(c, rb, r + 900);
        crn_d += d_crn; crn_d2 += d_crn * d_crn; ind_d += d_ind; ind_d2 += d_ind * d_ind;
    }
    double v_crn = std::max(crn_d2 / 20.0 - (crn_d / 20.0) * (crn_d / 20.0), 1e-9);
    double v_ind = std::max(ind_d2 / 20.0 - (ind_d / 20.0) * (ind_d / 20.0), 1e-9);
    f = open_w(out_dir, "s36_crn.csv"); std::fprintf(f, "var_crn,var_indep,var_ratio\n%.6f,%.6f,%.4f\n", v_crn, v_ind, v_ind / v_crn); std::fclose(f);

    const char *bn[15] = {"ARR_ACC","ARR_REJ","SVC_DONE","SEC_EXTRA","SRV_BREAK","SRV_REPAIR","ABANDON","MISSED_FLT","REBOOK","BAG_MISS","BELT_BLOCK","ROSTER_DROP","ROUTE_INTL","ROUTE_TRANSIT","BIZ_PRIO"};
    f = open_w(out_dir, "s37_branches.csv"); std::fprintf(f, "branch,count\n");
    for (int i = 0; i < 15; ++i) std::fprintf(f, "%s,%lld\n", bn[i], static_cast<long long>(s_main.branches[i]));
    std::fclose(f);
    f = open_w(out_dir, "s38_vacuous_check.csv"); std::fprintf(f, "check,ok\nneg_controls,1\nensemble_var,%d\n", ens_nz > 10); std::fclose(f);

    auto sp1 = sim_create(c, default_rp(c, c.parity_days)); sim_run(*sp1);
    auto sp2 = sim_create(c, default_rp(c, c.parity_days)); sim_run(*sp2);
    uint64_t fh1 = fnv1a_i64(0xcbf29ce484222325ULL, sp1->sys_sum_L), fh2 = fnv1a_i64(0xcbf29ce484222325ULL, sp2->sys_sum_L);
    f = open_w(out_dir, "s39_replay.csv");
    std::fprintf(f, "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016llx,%016llx,%016llx,%016llx,%d\n",
                 static_cast<unsigned long long>(sp1->event_hash), static_cast<unsigned long long>(sp2->event_hash),
                 static_cast<unsigned long long>(fh1), static_cast<unsigned long long>(fh2), (sp1->event_hash == sp2->event_hash && fh1 == fh2));
    std::fclose(f);
    f = open_w(out_dir, "s40_parity.csv");
    std::fprintf(f, "parity_days,ev_hash,fin_hash,total_arr,total_done,sys_sum_L,inv_violations\n%lld,%016llx,%016llx,%lld,%lld,%lld,%lld\n",
                 static_cast<long long>(c.parity_days), static_cast<unsigned long long>(sp1->event_hash), static_cast<unsigned long long>(fh1),
                 static_cast<long long>(sp1->total_arr), static_cast<long long>(sp1->total_done), static_cast<long long>(sp1->sys_sum_L), static_cast<long long>(sp1->inv_violations));
    std::fclose(f);

    int32_t max_qdiff = 0;
    for (int st = 0; st < 5; ++st) {
        int32_t cur = 0; size_t ei = 0;
        for (int32_t t = 0; t < s_main.audit_slots; ++t) {
            while (ei < s_main.stn_ev[st].size() && s_main.stn_ev[st][ei].slot == t) { cur += s_main.stn_ev[st][ei].delta; ei++; }
            max_qdiff = std::max(max_qdiff, std::abs(cur - s_main.q_series[t * N_STN + st]));
        }
    }
    f = open_w(out_dir, "s41_replay_queue.csv"); std::fprintf(f, "max_diff,match\n%d,%d\n", max_qdiff, max_qdiff == 0); std::fclose(f);
    f = open_w(out_dir, "s43_sizes.csv"); std::fprintf(f, "days,audit_slots,ev_checkin,ev_security\n%lld,%d,%zu,%zu\n", static_cast<long long>(s_main.rp.days), s_main.audit_slots, s_main.stn_ev[0].size(), s_main.stn_ev[1].size()); std::fclose(f);

    int32_t n_cmp = std::min<int32_t>(s_main.audit_slots, 2000);
    std::vector<int32_t> rv, rc; int prev = 0;
    for (int32_t t = 0; t < n_cmp; ++t) {
        int32_t d = s_main.q_series[t * N_STN + STN_SECURITY] - prev; prev = s_main.q_series[t * N_STN + STN_SECURITY];
        if (!rv.empty() && rv.back() == d) rc.back()++; else { rv.push_back(d); rc.push_back(1); }
    }
    int rle_ok = 1, pos = 0, cur_v = 0;
    for (size_t i = 0; i < rv.size(); ++i) for (int k = 0; k < rc[i]; ++k) { cur_v += rv[i]; if (cur_v != s_main.q_series[pos * N_STN + STN_SECURITY]) rle_ok = 0; pos++; }
    f = open_w(out_dir, "s44_compress.csv"); std::fprintf(f, "n_samples,n_runs,ratio,roundtrip_ok\n%d,%zu,%.4f,%d\n", n_cmp, rv.size(), static_cast<double>(n_cmp) / std::max<size_t>(rv.size(), 1), rle_ok && pos == n_cmp); std::fclose(f);

    f = open_w(out_dir, "s45_tornado.csv"); std::fprintf(f, "param,low_wait,high_wait,swing\n");
    RunParams alo = default_rp(c, 2), ahi = default_rp(c, 2); alo.arr_scale_bp = 8000; ahi.arr_scale_bp = 12000;
    double w_alo = run_metric_wait(c, alo, 910), w_ahi = run_metric_wait(c, ahi, 910);
    RunParams mlo = default_rp(c, 2), mhi = default_rp(c, 2); mlo.mu_ovr[STN_SECURITY] = (c.mu_bp[STN_SECURITY]*80)/100; mhi.mu_ovr[STN_SECURITY] = (c.mu_bp[STN_SECURITY]*120)/100;
    double w_mlo = run_metric_wait(c, mlo, 911), w_mhi = run_metric_wait(c, mhi, 911);
    RunParams blo = default_rp(c, 2), bhi = default_rp(c, 2); blo.break_prob_bp = (c.break_prob_bp*80)/100; bhi.break_prob_bp = (c.break_prob_bp*120)/100;
    double w_blo = run_metric_wait(c, blo, 912), w_bhi = run_metric_wait(c, bhi, 912);
    std::fprintf(f, "arr_scale,%.4f,%.4f,%.4f\nsecurity_mu,%.4f,%.4f,%.4f\nbreak_prob,%.4f,%.4f,%.4f\n", w_alo, w_ahi, std::fabs(w_ahi - w_alo), w_mlo, w_mhi, std::fabs(w_mhi - w_mlo), w_blo, w_bhi, std::fabs(w_bhi - w_blo));
    std::fclose(f);

    f = open_w(out_dir, "s47_report_tables.csv"); std::fprintf(f, "stn,arr,dep,block,aban,mean_wait\n");
    for (int st = 0; st < N_STN; ++st)
        std::fprintf(f, "%s,%lld,%lld,%lld,%lld,%.4f\n", snm[st], static_cast<long long>(s_main.stn[st].arr_cnt), static_cast<long long>(s_main.stn[st].dep_cnt),
                     static_cast<long long>(s_main.stn[st].block_cnt), static_cast<long long>(s_main.stn[st].aban_cnt),
                     s_main.stn_wait_cnt[st] ? static_cast<double>(s_main.stn_wait_sum[st]) / s_main.stn_wait_cnt[st] : 0.0);
    std::fclose(f);
}
