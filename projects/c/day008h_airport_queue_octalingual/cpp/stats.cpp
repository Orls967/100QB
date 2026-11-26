#include "common.hpp"

struct SimValRes { double p_block, L, L_q, U, W, W_q; };

static SimValRes run_single_sim(Xoshiro256 &rng, int c, int K, int64_t lam_bp, int64_t mu_bp, int64_t slots) {
    int q = 0; int64_t arr = 0, blk = 0, sum_L = 0, sum_Lq = 0, sum_busy = 0;
    for (int64_t t = 0; t < slots; ++t) {
        int a = bernoulli(&rng, lam_bp), m = std::min(q, c), d = 0;
        for (int s = 0; s < m; ++s) if (bernoulli(&rng, mu_bp)) d++;
        if (q < K) { q += a; if (a) arr++; } else if (a) blk++;
        q -= d;
        if (t >= slots / 5) { int mb = std::min(q, c); sum_L += q; sum_Lq += (q - mb); sum_busy += mb; }
    }
    double n_s = static_cast<double>(slots - slots / 5), lam_eff = static_cast<double>(arr) / static_cast<double>(slots);
    SimValRes r{};
    r.p_block = static_cast<double>(blk) / ((arr + blk) > 0 ? static_cast<double>(arr + blk) : 1.0);
    r.L = static_cast<double>(sum_L) / n_s; r.L_q = static_cast<double>(sum_Lq) / n_s;
    r.U = static_cast<double>(sum_busy) / (n_s * static_cast<double>(c));
    r.W = lam_eff > 0.0 ? r.L / lam_eff : 0.0; r.W_q = lam_eff > 0.0 ? r.L_q / lam_eff : 0.0;
    return r;
}

static void run_tandem_sim(Xoshiro256 &rng, int K1, int K2, int64_t lam_bp, int64_t mu1_bp, int64_t mu2_bp, int64_t slots,
                           double &L1, double &L2, double &U1, double &U2) {
    int q1 = 0, q2 = 0; int64_t sL1 = 0, sL2 = 0, sU1 = 0, sU2 = 0;
    for (int64_t t = 0; t < slots; ++t) {
        int a = bernoulli(&rng, lam_bp), d1 = (q1 > 0) ? bernoulli(&rng, mu1_bp) : 0, d2 = (q2 > 0) ? bernoulli(&rng, mu2_bp) : 0;
        if (q1 < K1) q1 += a;
        q1 -= d1; q2 -= d2;
        if (d1 && q2 < K2) q2++;
        if (t >= slots / 5) { sL1 += q1; sL2 += q2; if (q1 > 0) sU1++; if (q2 > 0) sU2++; }
    }
    double ns = static_cast<double>(slots - slots / 5);
    L1 = static_cast<double>(sL1) / ns; L2 = static_cast<double>(sL2) / ns;
    U1 = static_cast<double>(sU1) / ns; U2 = static_cast<double>(sU2) / ns;
}

static FILE *open_w(const char *dir, const char *nm) {
    char p[512]; std::snprintf(p, sizeof(p), "%s/%s", dir, nm);
    FILE *f = std::fopen(p, "w"); if (!f) std::exit(1); return f;
}

void write_s01_to_s13(const Config &c, const char *out_dir) {
    uint64_t hrng = 0xcbf29ce484222325ULL;
    const uint32_t pairs[3][2] = {{0,0},{1,2},{2,1}};
    const uint64_t exp_hv[3] = {0xbe40b757fe0b6697ULL, 0x659394d6f2a43a95ULL, 0x3abda5bac0dd49cbULL};
    int rng_ok = 1;
    for (int i = 0; i < 3; ++i) {
        uint64_t x = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1]);
        Xoshiro256 st = xoshiro_seed(splitmix64(&x)); uint64_t buf[4];
        for (int j = 0; j < 4; ++j) buf[j] = xoshiro_next(&st);
        uint64_t hv = fnv1a_bytes(buf, 32); if (hv != exp_hv[i]) rng_ok = 0;
        hrng = fnv1a_i64(hrng, static_cast<int64_t>(hv));
    }
    auto s0 = sim_create(c, default_rp(c, 1)); sim_run(*s0);
    int inv_ok = (s0->inv_violations == 0);
    FILE *f = open_w(out_dir, "s01_selftest.txt");
    std::fprintf(f, "ok,h_rng,inv_ok,tot_arr,tot_done\n%d,%016llx,%d,%lld,%lld\n",
                 rng_ok && inv_ok, static_cast<unsigned long long>(hrng), inv_ok, static_cast<long long>(s0->total_arr), static_cast<long long>(s0->total_done));
    std::fclose(f);

    f = open_w(out_dir, "s02_config.csv");
    std::fprintf(f, "key,value\nmaster_seed,%lld\nslots_per_day,%lld\nparity_days,%lld\narr_scale_bp,%lld\nbelt_cap,%lld\n",
                 static_cast<long long>(c.master_seed), static_cast<long long>(c.slots_per_day), static_cast<long long>(c.parity_days),
                 static_cast<long long>(c.arr_scale_bp), static_cast<long long>(c.belt_cap));
    std::fclose(f);

    f = open_w(out_dir, "s03_flights.csv"); std::fprintf(f, "flight_id,dep_slot,gate_close_slot,seats\n");
    for (int h = 0; h < 24; ++h) std::fprintf(f, "%d,%d,%d,180\n", h, h * 360 + 300, h * 360 + 300 - static_cast<int>(c.gate_close_lead));
    std::fclose(f);
    f = open_w(out_dir, "s03_arrival_profile.csv"); std::fprintf(f, "hour,arr_bp\n");
    for (int h = 0; h < 24; ++h) std::fprintf(f, "%d,%d\n", h, s0->arr_prof[h]);
    std::fclose(f);

    f = open_w(out_dir, "s04_stations.csv"); std::fprintf(f, "stn_id,name,c_max,mu_bp,cap,has_prio\n");
    const char *snm[N_STN] = {"CHECKIN","SECURITY","IMIGRASI","GATE","BAGGAGE_BELT","CUSTOMS_ARRIVAL"};
    for (int i = 0; i < N_STN; ++i) std::fprintf(f, "%d,%s,%d,%lld,%d,%d\n", i, snm[i], s0->stn[i].c_max, static_cast<long long>(c.mu_bp[i]), s0->stn[i].cap, s0->stn[i].has_prio);
    std::fclose(f);
    f = open_w(out_dir, "s04_roster.csv"); std::fprintf(f, "hour,checkin,security,imigrasi,gate,belt,customs\n");
    for (int h = 0; h < 24; ++h) std::fprintf(f, "%d,%d,%d,%d,%d,%d,%d\n", h, s0->roster[h][0], s0->roster[h][1], s0->roster[h][2], s0->roster[h][3], s0->roster[h][4], s0->roster[h][5]);
    std::fclose(f);

    f = open_w(out_dir, "s05_queue_script.csv"); std::fprintf(f, "case_id,mode,q_init,cap,prio,breakdown,gate_closed,accepted,served,lost\n");
    for (int k = 0; k < 12; ++k) {
        int q0=(k%4)+1, cap=(k>=6)?3:10, brk=(k%3==2)?1:0, gc=(k==11)?1:0, acc=(q0<=cap)?q0:cap, srv=gc?0:(brk?acc-1:acc), lost=gc?acc:(q0-acc);
        std::fprintf(f, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", k, k%3, q0, cap, k%2, brk, gc, acc, srv, lost);
    }
    std::fclose(f);
    f = open_w(out_dir, "s06_bag_script.csv"); std::fprintf(f, "case_id,n_bags,belt_len,belt_cap,miss_bp,accepted,blocked,missed\n");
    for (int k = 0; k < 12; ++k) {
        int nb=k%3, bl=(k/3)*4, cap=10, acc=(bl+nb<=cap)?nb:(cap-bl>0?cap-bl:0), blk=nb-acc, ms=(k%4==3&&acc>0)?1:0;
        std::fprintf(f, "%d,%d,%d,%d,%lld,%d,%d,%d\n", k, nb, bl, cap, static_cast<long long>(c.bag_miss_bp), acc, blk, ms);
    }
    std::fclose(f);

    f = open_w(out_dir, "s07_dtmc_single.csv"); std::fprintf(f, "c,K,lam,mu,p_block,L,L_q,U,W,W_q\n");
    double lam = 0.08, mu = 0.12; const int cs[3] = {1,2,3}, Ks[3] = {5,10,60};
    for (int ic=0; ic<3; ++ic) for (int ik=0; ik<3; ++ik) {
        DtmcRes dr = solve_dtmc_single(cs[ic], Ks[ik], lam, mu);
        std::fprintf(f, "%d,%d,%.4f,%.4f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", cs[ic], Ks[ik], lam, mu, dr.p_block, dr.L, dr.L_q, dr.U, dr.W, dr.W_q);
    }
    std::fclose(f);

    f = open_w(out_dir, "s08_dtmc_tandem.csv"); std::fprintf(f, "K1,K2,lam,mu1,mu2,p_blk1,L1,U1,p_blk2,L2,U2,p_blk_tnd\n");
    DtmcRes t1, t2; double pbt = 0.0; solve_dtmc_tandem(6, 6, 0.08, 0.14, 0.12, t1, t2, pbt);
    std::fprintf(f, "6,6,0.0800,0.1400,0.1200,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", t1.p_block, t1.L, t1.U, t2.p_block, t2.L, t2.U, pbt);
    std::fclose(f);

    f = open_w(out_dir, "s09_geo_closed_form.csv");
    DtmcRes cf = solve_geo_closed_form(lam, mu), d60 = solve_dtmc_single(1, 60, lam, mu);
    std::fprintf(f, "lam,mu,L_closed,W_closed,U_closed,L_dtmc60,rel_err_L\n%.4f,%.4f,%.12f,%.12f,%.12f,%.12f,%.12e\n",
                 lam, mu, cf.L, cf.W, cf.U, d60.L, std::fabs(cf.L - d60.L) / cf.L);
    std::fclose(f);

    f = open_w(out_dir, "s10_single_validation.csv"); std::fprintf(f, "c,K,metric,dtmc,sim_mean,sim_se,z_score,within_3se\n");
    int64_t val_slots = 20000;
    for (int ic = 0; ic < 3; ++ic) {
        int cc = cs[ic], KK = (cc == 1) ? 10 : 5;
        DtmcRes dr = solve_dtmc_single(cc, KK, lam, mu);
        double sL=0,sL2=0,sU=0,sU2=0,sP=0,sP2=0;
        for (int r = 0; r < 50; ++r) {
            uint64_t x = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(500 + static_cast<uint32_t>(ic), static_cast<uint32_t>(r));
            Xoshiro256 rng = xoshiro_seed(splitmix64(&x));
            SimValRes vr = run_single_sim(rng, cc, KK, 800, 1200, val_slots);
            sL+=vr.L; sL2+=vr.L*vr.L; sU+=vr.U; sU2+=vr.U*vr.U; sP+=vr.p_block; sP2+=vr.p_block*vr.p_block;
        }
        double mL=sL/50, seL=std::sqrt(std::max((sL2/50-mL*mL)/49,1e-18));
        double mU=sU/50, seU=std::sqrt(std::max((sU2/50-mU*mU)/49,1e-18));
        double mP=sP/50, seP=std::sqrt(std::max((sP2/50-mP*mP)/49,1e-18));
        std::fprintf(f, "%d,%d,L,%.6f,%.6f,%.6f,%.3f,%d\n", cc, KK, dr.L, mL, seL, (mL-dr.L)/seL, std::fabs((mL-dr.L)/seL)<=3.0);
        std::fprintf(f, "%d,%d,U,%.6f,%.6f,%.6f,%.3f,%d\n", cc, KK, dr.U, mU, seU, (mU-dr.U)/seU, std::fabs((mU-dr.U)/seU)<=3.0);
        std::fprintf(f, "%d,%d,p_block,%.6f,%.6f,%.6f,%.3f,%d\n", cc, KK, dr.p_block, mP, seP, (mP-dr.p_block)/seP, std::fabs((mP-dr.p_block)/seP)<=3.0);
    }
    std::fclose(f);

    f = open_w(out_dir, "s11_tandem_validation.csv"); std::fprintf(f, "metric,dtmc,sim_mean,sim_se,z_score,within_3se\n");
    double sL1=0,sL12=0,sL2=0,sL22=0,sU1=0,sU12=0,sU2=0,sU22=0;
    for (int r = 0; r < 50; ++r) {
        uint64_t x = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(600, static_cast<uint32_t>(r));
        Xoshiro256 rng = xoshiro_seed(splitmix64(&x));
        double l1, l2, u1, u2; run_tandem_sim(rng, 6, 6, 800, 1400, 1200, val_slots, l1, l2, u1, u2);
        sL1+=l1; sL12+=l1*l1; sL2+=l2; sL22+=l2*l2; sU1+=u1; sU12+=u1*u1; sU2+=u2; sU22+=u2*u2;
    }
    double mL1=sL1/50, seL1=std::sqrt(std::max((sL12/50-mL1*mL1)/49,1e-18));
    double mL2=sL2/50, seL2=std::sqrt(std::max((sL22/50-mL2*mL2)/49,1e-18));
    double mU1=sU1/50, seU1=std::sqrt(std::max((sU12/50-mU1*mU1)/49,1e-18));
    double mU2=sU2/50, seU2=std::sqrt(std::max((sU22/50-mU2*mU2)/49,1e-18));
    std::fprintf(f, "L1,%.6f,%.6f,%.6f,%.3f,%d\nL2,%.6f,%.6f,%.6f,%.3f,%d\nU1,%.6f,%.6f,%.6f,%.3f,%d\nU2,%.6f,%.6f,%.6f,%.3f,%d\n",
                 t1.L, mL1, seL1, (mL1-t1.L)/seL1, std::fabs((mL1-t1.L)/seL1)<=3.0,
                 t2.L, mL2, seL2, (mL2-t2.L)/seL2, std::fabs((mL2-t2.L)/seL2)<=3.0,
                 t1.U, mU1, seU1, (mU1-t1.U)/seU1, std::fabs((mU1-t1.U)/seU1)<=3.0,
                 t2.U, mU2, seU2, (mU2-t2.U)/seU2, std::fabs((mU2-t2.U)/seU2)<=3.0);
    std::fclose(f);

    RunParams rpw = default_rp(c, 2); rpw.audit_on = 1;
    auto sw = sim_create(c, rpw); sim_run(*sw);
    int n_blk = sw->audit_slots / 5, best_d = 0; double best_mser = 1e30;
    std::vector<double> blk_y(n_blk);
    for (int b = 0; b < n_blk; ++b) {
        double sm = 0; for (int k = 0; k < 5; ++k) sm += sw->q_series[(b * 5 + k) * N_STN + STN_SECURITY];
        blk_y[b] = sm / 5.0;
    }
    for (int d = 0; d < n_blk / 2; d += 10) {
        int m = n_blk - d; double mn = 0, s2 = 0;
        for (int i = d; i < n_blk; ++i) mn += blk_y[i];
        mn /= static_cast<double>(m);
        for (int i = d; i < n_blk; ++i) { double df = blk_y[i] - mn; s2 += df * df; }
        double mser = s2 / (static_cast<double>(m) * static_cast<double>(m));
        if (mser < best_mser) { best_mser = mser; best_d = d * 5; }
    }
    f = open_w(out_dir, "s12_warmup.csv"); std::fprintf(f, "warmup_slots,mser_min\n%d,%.6f\n", best_d, best_mser); std::fclose(f);

    f = open_w(out_dir, "s13_calibration.csv"); std::fprintf(f, "arr_scale_bp,max_util,loss_frac,calibrated\n");
    for (int scale : {15000, 18000, 22000, 26000, 30000}) {
        RunParams rpc = default_rp(c, 2); rpc.arr_scale_bp = scale;
        auto sc = sim_create(c, rpc); sim_run(*sc);
        double max_u = 0.0;
        for (int st = 0; st < 4; ++st) max_u = std::max(max_u, static_cast<double>(sc->stn[st].busy_slots) / static_cast<double>(2 * SLOTS_PER_DAY * sc->stn[st].c_max));
        double lf = sc->total_arr > 0 ? static_cast<double>(sc->total_lost) / static_cast<double>(sc->total_arr) : 0.0;
        std::fprintf(f, "%d,%.4f,%.4f,%d\n", scale, max_u, lf, (max_u >= 0.55 && max_u <= 0.92 && lf < 0.15) ? 1 : 0);
    }
    std::fclose(f);
}
