import math
import struct
from engine import analytical_r0, solve_final_size, sim_run
from population import BP, ST_D, RunParams, Sim
from rng import FNV_OFFSET, MASK64, Xoshiro256, fnv1a_bytes, fnv1a_i64, fnv1a_u32pair, splitmix64


def net_degrees(s: Sim):
    shh = sum(s.hh_off[p.hh_id + 1] - s.hh_off[p.hh_id] - 1 for p in s.pop)
    ssw = sum(6.0 if p.age_grp == 0 else (5.0 if p.age_grp == 1 else 0.0) for p in s.pop)
    scom = sum(12.0 if p.super_flag != 0 else 3.0 for p in s.pop)
    nf = float(s.n)
    return shh / nf, ssw / nf, scom / nf


def write_s01_to_s09(c: dict, out_dir: str):
    hrng = FNV_OFFSET
    pairs = [(0, 0), (1, 2), (2, 1)]
    exp_hv = [0xBE40B757FE0B6697, 0x659394D6F2A43A95, 0x3ABDA5BAC0DD49CB]
    rng_ok = True
    for i in range(3):
        x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1])
        _, sd = splitmix64(x)
        st = Xoshiro256(sd)
        buf = b"".join(struct.pack("<Q", st.next_u64()) for _ in range(4))
        hv = fnv1a_bytes(buf)
        if hv != exp_hv[i]:
            rng_ok = False
        hrng = fnv1a_i64(hrng, hv)
    st0 = Sim(c, RunParams(200, 10, 450, interv_on=0), False)
    sim_run(st0)
    inv_ok = st0.inv_violations == 0

    with open(f"{out_dir}/s01_selftest.txt", "w") as f:
        f.write(f"ok,h_rng,inv_ok,tot_trans\n{int(rng_ok and inv_ok)},{hrng:016x},{int(inv_ok)},{st0.total_trans}\n")
    with open(f"{out_dir}/s02_config.csv", "w") as f:
        f.write(f"key,value\nmaster_seed,{c['master_seed']}\nparity_N,{c['parity_N']}\nparity_days,{c['parity_days']}\nmain_days,{c['main_days']}\nbeta_bp_init,{c['beta_bp_init']}\n")

    sp = Sim(c, RunParams(1000, 1, 450, interv_on=0), False)
    with open(f"{out_dir}/s03_population.csv", "w") as f:
        f.write(f"age_grp,count\n0,{sp.age_tot[0]}\n1,{sp.age_tot[1]}\n2,{sp.age_tot[2]}\n")
    with open(f"{out_dir}/s03_households.csv", "w") as f:
        f.write("hh_id,size\n")
        for h in range(min(sp.n_hh, 50)):
            f.write(f"{h},{sp.hh_off[h + 1] - sp.hh_off[h]}\n")
    with open(f"{out_dir}/s03_schools_work.csv", "w") as f:
        f.write("sw_id,size\n")
        for w in range(min(sp.n_sw, 50)):
            f.write(f"{w},{sp.sw_off[w + 1] - sp.sw_off[w]}\n")
    chh, csw, ccom = net_degrees(sp)
    with open(f"{out_dir}/s04_network.csv", "w") as f:
        f.write(f"n_hh,n_sw,deg_hh,deg_sw,deg_comm\n{sp.n_hh},{sp.n_sw},{chh:.4f},{csw:.4f},{ccom:.4f}\n")

    with open(f"{out_dir}/s05_contact_script.csv", "w") as f:
        f.write("case_id,vax,mask,quar,eff_bp\n")
        for k in range(12):
            vax, mask, quar = k % 3, (k // 3) % 2, (k // 6) % 2
            p = 0 if quar != 0 else 450
            if mask != 0:
                p = (p * (BP - c["mask_eff_bp"])) // BP
            if vax == 1:
                p = (p * (BP - c["vax_eff1_bp"])) // BP
            elif vax == 2:
                p = (p * (BP - c["vax_eff2_bp"])) // BP
            f.write(f"{k},{vax},{mask},{quar},{p}\n")
    with open(f"{out_dir}/s06_progression_script.csv", "w") as f:
        f.write("case_id,age_grp,hosp_rej,p_hosp_bp,p_die_bp\n")
        for k in range(12):
            ag, rej = k % 3, (k // 3) % 2
            ph = 120 if ag == 0 else (450 if ag == 1 else 1800)
            pd = 100 if ag == 0 else (500 if ag == 1 else 2000)
            if rej != 0:
                pd *= 2
            f.write(f"{k},{ag},{rej},{ph},{pd}\n")
    with open(f"{out_dir}/s07_interventions.csv", "w") as f:
        f.write(f"name,start_day,param_bp\nschool_close,{c['school_close_day']},10000\nlockdown,{c['lockdown_day']},5000\nmask,{c['mask_day']},{c['mask_eff_bp']}\ntest,{c['test_day']},{c['test_prob_bp']}\nvax,{c['vax_day']},{c['vax_daily_bp']}\nvariant,{c['variant_day']},{c['variant_mult_bp']}\n")

    with open(f"{out_dir}/s08_calibration.csv", "w") as f:
        f.write("beta_bp,r0_analytic,attack_rate,calibrated\n")
        for bv in (300, 380, 450, 520, 600):
            sc = Sim(c, RunParams(1000, 120, bv, interv_on=0), False)
            sim_run(sc)
            ch, cs2, cc2 = net_degrees(sc)
            r0 = analytical_r0(c, bv, ch, cs2, cc2)
            ar = sc.total_trans / sc.n
            f.write(f"{bv},{r0:.4f},{ar:.4f},{int(1.5 <= r0 <= 3.5 and 0.20 <= ar <= 0.80)}\n")
    s9 = Sim(c, RunParams(500, 30, 450, interv_on=0), False)
    sim_run(s9)
    with open(f"{out_dir}/s09_baseline_small.csv", "w") as f:
        f.write(f"N,days,total_trans,deaths,inv_ok\n500,30,{s9.total_trans},{s9.daily_cnt[29][ST_D]},{int(s9.inv_violations == 0)}\n")


def write_daily_csv(s: Sim, out_dir: str, fname: str):
    with open(f"{out_dir}/{fname}", "w") as f:
        f.write("day,S,E,Ia,Is,H,ICU,R,D\n")
        for d in range(s.rp.days):
            f.write(f"{d}," + ",".join(str(x) for x in s.daily_cnt[d]) + "\n")


def write_trans_csv(s: Sim, out_dir: str, fname: str):
    with open(f"{out_dir}/{fname}", "w") as f:
        f.write("day,src_id,dst_id,layer,variant,gen_int\n")
        for r in s.trans:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]}\n")


def write_s10_to_s37(s_no: Sim, s_main: Sim, c: dict, out_dir: str, beta_cal: int):
    d_cnt = s_main.rp.days
    write_daily_csv(s_no, out_dir, "daily_counts_nointerv.csv")
    write_trans_csv(s_no, out_dir, "transmissions_log_nointerv.csv")
    write_daily_csv(s_main, out_dir, "daily_counts.csv")
    write_trans_csv(s_main, out_dir, "transmissions_log.csv")

    for fname, hdr, rows in [
        ("incidence_by_age.csv", "day,age_0_17,age_18_59,age_60p", (f"{d},{s_main.daily_inc_age[d][0]},{s_main.daily_inc_age[d][1]},{s_main.daily_inc_age[d][2]}" for d in range(d_cnt))),
        ("hospital_daily.csv", "day,hosp_occ,icu_occ,hosp_rej_cum,icu_rej_cum", (f"{d},{s_main.daily_hosp[d][0]},{s_main.daily_hosp[d][1]},{s_main.daily_hosp[d][2]},{s_main.daily_hosp[d][3]}" for d in range(d_cnt))),
        ("deaths_by_age.csv", "day,d_0_17,d_18_59,d_60p", (f"{d},{s_main.daily_deaths_age[d][0]},{s_main.daily_deaths_age[d][1]},{s_main.daily_deaths_age[d][2]}" for d in range(d_cnt))),
        ("tests.csv", "day,tests,positives", (f"{d},{s_main.daily_tests[d][0]},{s_main.daily_tests[d][1]}" for d in range(d_cnt))),
        ("quarantine.csv", "day,quarantined", (f"{d},{s_main.daily_quar[d]}" for d in range(d_cnt))),
        ("vaccination.csv", "day,dose1,dose2", (f"{d},{s_main.daily_vax[d][0]},{s_main.daily_vax[d][1]}" for d in range(d_cnt))),
        ("variant_share.csv", "day,v0_inf,v1_inf", (f"{d},{s_main.daily_var[d][0]},{s_main.daily_var[d][1]}" for d in range(d_cnt))),
        ("generation_interval.csv", "idx,gen_int", (f"{i},{r[5]}" for i, r in enumerate(s_main.trans[:2000]))),
        ("agents_final.csv", "id,age_grp,status,vax_doses,offspring,ever_inf", (f"{i},{p.age_grp},{p.status},{p.vax_doses},{p.offspring},{p.ever_inf}" for i, p in enumerate(s_main.pop[:2000]))),
        ("progression_log.csv", "day,pid,from_st,to_st", (f"{r[0]},{r[1]},{r[2]},{r[3]}" for r in s_main.progs)),
        ("hospital_cap_nointerv.csv", "day,hosp_occ,icu_occ", (f"{d},{s_no.daily_hosp[d][0]},{s_no.daily_hosp[d][1]}" for d in range(d_cnt))),
        ("deaths_nointerv.csv", "day,deaths", (f"{d},{s_no.daily_cnt[d][ST_D]}" for d in range(d_cnt))),
        ("rt_nointerv.csv", "day,inc", (f"{d},{sum(s_no.daily_inc_age[d])}" for d in range(d_cnt))),
    ]:
        with open(f"{out_dir}/{fname}", "w") as f:
            f.write(hdr + "\n")
            for ln in rows:
                f.write(ln + "\n")

    for a in range(3):
        with open(f"{out_dir}/daily_age_{a}.csv", "w") as f:
            f.write("day,inc,deaths\n")
            for d in range(d_cnt):
                f.write(f"{d},{s_main.daily_inc_age[d][a]},{s_main.daily_deaths_age[d][a]}\n")
    for l, nm in enumerate(("layer_hh.csv", "layer_sw.csv", "layer_comm.csv")):
        with open(f"{out_dir}/{nm}", "w") as f:
            f.write(f"layer,infections\n{l},{s_main.layer_inf[l]}\n")
    oh_bins = [0] * 6
    for p in s_no.pop:
        if p.ever_inf != 0:
            oh_bins[min(p.offspring, 5)] += 1
    with open(f"{out_dir}/offspring_hist.csv", "w") as f:
        f.write("bin,count\n")
        for b in range(6):
            f.write(f"{b},{oh_bins[b]}\n")

    tot_viol = s_no.inv_violations + s_main.inv_violations
    with open(f"{out_dir}/s12_invariants.csv", "w") as f:
        f.write(f"violations,ok\n{tot_viol},{int(tot_viol == 0)}\n")
    with open(f"{out_dir}/s13_rt.csv", "w") as f:
        f.write("day,rt\n")
        for d in range(d_cnt):
            num = sum(sum(s_main.daily_inc_age[d - k]) for k in range(5) if d >= k)
            den = sum(sum(s_main.daily_inc_age[d - k]) for k in range(4, 9) if d >= k)
            f.write(f"{d},{(num / den if den > 0 else 0.0):.4f}\n")

    chh, csw, ccom = net_degrees(s_no)
    r0_ana = analytical_r0(c, beta_cal, chh, csw, ccom)
    i_early = (sum(s_no.daily_inc_age[3]) + sum(s_no.daily_inc_age[4]) + 1) / 2.0
    i_mid = sum(sum(s_no.daily_inc_age[d]) for d in range(10, min(16, d_cnt))) / 6.0
    r_grow = math.log(i_mid / i_early) / 9.0 if i_mid > i_early else 0.0
    with open(f"{out_dir}/s14_growth.csv", "w") as f:
        f.write(f"growth_r,doubling_days\n{r_grow:.6f},{(math.log(2.0) / r_grow if r_grow > 0 else 0.0):.4f}\n")
    r0_emp = 1.0 + r_grow * 4.5
    with open(f"{out_dir}/s15_r0.csv", "w") as f:
        f.write(f"r0_analytic,r0_empirical,ratio\n{r0_ana:.4f},{r0_emp:.4f},{(r0_emp / r0_ana if r0_ana > 0 else 0.0):.4f}\n")
    z_ana = solve_final_size(r0_ana)
    herd = 1.0 - 1.0 / r0_ana if r0_ana > 1.0 else 0.0
    with open(f"{out_dir}/s16_final_size.csv", "w") as f:
        f.write(f"z_analytic,herd_thresh,attack_nointerv,attack_interv\n{z_ana:.4f},{herd:.4f},{s_no.total_trans / s_no.n:.4f},{s_main.total_trans / s_main.n:.4f}\n")
    with open(f"{out_dir}/s17_attack.csv", "w") as f:
        f.write(f"category,count,frac\nHH,{s_main.layer_inf[0]},{s_main.layer_inf[0] / s_main.n:.4f}\nSW,{s_main.layer_inf[1]},{s_main.layer_inf[1] / s_main.n:.4f}\nCOMM,{s_main.layer_inf[2]},{s_main.layer_inf[2] / s_main.n:.4f}\nAGE0,{s_main.age_inf[0]},{s_main.age_inf[0] / max(s_main.age_tot[0], 1):.4f}\nAGE1,{s_main.age_inf[1]},{s_main.age_inf[1] / max(s_main.age_tot[1], 1):.4f}\nAGE2,{s_main.age_inf[2]},{s_main.age_inf[2] / max(s_main.age_tot[2], 1):.4f}\n")

    m1, m2, n_inf = 0.0, 0.0, 0
    for p in s_no.pop:
        if p.ever_inf != 0:
            o = float(p.offspring)
            m1 += o
            m2 += o * o
            n_inf += 1
    if n_inf > 0:
        m1 /= n_inf
        m2 /= n_inf
    var_o = m2 - m1 * m1
    k_disp = (m1 * m1) / (var_o - m1) if var_o > m1 and m1 > 0 else 10.0
    with open(f"{out_dir}/s18_offspring.csv", "w") as f:
        f.write(f"n_inf,mean_offspring,var_offspring,k_disp\n{n_inf},{m1:.4f},{var_o:.4f},{k_disp:.4f}\n")

    pk_h_no = max(s_no.daily_hosp[d][0] for d in range(d_cnt))
    pk_i_no = max(s_no.daily_hosp[d][1] for d in range(d_cnt))
    pk_h_m = max(s_main.daily_hosp[d][0] for d in range(d_cnt))
    pk_i_m = max(s_main.daily_hosp[d][1] for d in range(d_cnt))
    with open(f"{out_dir}/s19_hospital.csv", "w") as f:
        f.write(f"hosp_cap,icu_cap,peak_h_no,peak_i_no,rej_no,peak_h_main,peak_i_main,rej_main\n{s_main.hosp_cap},{s_main.icu_cap},{pk_h_no},{pk_i_no},{s_no.branches[9] + s_no.branches[11]},{pk_h_m},{pk_i_m},{s_main.branches[9] + s_main.branches[11]}\n")

    with open(f"{out_dir}/s20_ifr.csv", "w") as f:
        f.write("age_grp,infections,deaths,ifr\n")
        for a in range(3):
            dt = sum(s_main.daily_deaths_age[d][a] for d in range(d_cnt))
            f.write(f"{a},{s_main.age_inf[a]},{dt},{(dt / s_main.age_inf[a] if s_main.age_inf[a] > 0 else 0.0):.6f}\n")
    t_tot = sum(s_main.daily_tests[d][0] for d in range(d_cnt))
    t_pos = sum(s_main.daily_tests[d][1] for d in range(d_cnt))
    with open(f"{out_dir}/s21_testing.csv", "w") as f:
        f.write(f"tests,positives,quar_traces\n{t_tot},{t_pos},{s_main.branches[17]}\n")

    vn, vi = [0, 0, 0], [0, 0, 0]
    for p in s_main.pop:
        vd = min(p.vax_doses, 2)
        vn[vd] += 1
        if p.ever_inf != 0:
            vi[vd] += 1
    with open(f"{out_dir}/s22_vaccine.csv", "w") as f:
        f.write("doses,n_people,n_inf,attack_rate\n")
        for d in range(3):
            f.write(f"{d},{vn[d]},{vi[d]},{(vi[d] / vn[d] if vn[d] > 0 else 0.0):.6f}\n")

    v0 = sum(s_main.daily_var[d][0] for d in range(d_cnt))
    v1 = sum(s_main.daily_var[d][1] for d in range(d_cnt))
    with open(f"{out_dir}/s23_variant.csv", "w") as f:
        f.write(f"v0_total,v1_total,v1_share\n{v0},{v1},{(v1 / (v0 + v1) if v0 + v1 > 0 else 0.0):.6f}\n")
    with open(f"{out_dir}/s24_waning.csv", "w") as f:
        f.write(f"waning_events,reinfections\n{s_main.branches[14]},{s_main.total_reinf}\n")

    nc1 = Sim(c, RunParams(500, 30, beta_cal, contact_scale_bp=0, interv_on=0), False); sim_run(nc1)
    nc2 = Sim(c, RunParams(500, 30, 0, contact_scale_bp=-1, interv_on=0), False); sim_run(nc2)
    nc3 = Sim(c, RunParams(500, 30, beta_cal, vax_perfect=1, interv_on=1), False); sim_run(nc3)
    with open(f"{out_dir}/s25_negative_controls.csv", "w") as f:
        f.write(f"control,infections,ok\nzero_contact,{nc1.total_trans},{int(nc1.total_trans == 0)}\nzero_seed,{nc2.total_trans},{int(nc2.total_trans == 0)}\nperfect_vax,{nc3.total_trans},{int(nc3.total_trans == 0)}\n")

    nsw, dsw = 500, 90
    with open(f"{out_dir}/s26_sweep_comp_test.csv", "w") as f:
        f.write("comp_bp,test_bp,rep,attack_rate\n")
        for cv in (2000, 6000, 9500):
            for tv in (1000, 5000, 9000):
                for rp in range(3):
                    c2 = dict(c); c2["master_seed"] += rp + 1
                    ss = Sim(c2, RunParams(nsw, dsw, beta_cal, compliance_ovr=cv, test_ovr=tv, interv_on=1), False); sim_run(ss)
                    f.write(f"{cv},{tv},{rp},{ss.total_trans / nsw:.6f}\n")
    with open(f"{out_dir}/s27_sweep_vax.csv", "w") as f:
        f.write("vax_day,vax_cov_bp,rep,attack_rate\n")
        for vd in (10, 45, 80):
            for vc in (2000, 5500, 9000):
                for rp in range(3):
                    c2 = dict(c); c2["master_seed"] += rp + 10
                    ss = Sim(c2, RunParams(nsw, dsw, beta_cal, vax_day_ovr=vd, vax_cov_ovr=vc, interv_on=1), False); sim_run(ss)
                    f.write(f"{vd},{vc},{rp},{ss.total_trans / nsw:.6f}\n")
    with open(f"{out_dir}/s28_sweep_school_lock.csv", "w") as f:
        f.write("school_day,lock_day,rep,attack_rate\n")
        for sc in (15, 200):
            for lk in (15, 45, 200):
                for rp in range(3):
                    c2 = dict(c); c2["master_seed"] += rp + 20
                    ss = Sim(c2, RunParams(nsw, dsw, beta_cal, school_ovr=sc, lock_ovr=lk, interv_on=1), False); sim_run(ss)
                    f.write(f"{sc},{lk},{rp},{ss.total_trans / nsw:.6f}\n")
    with open(f"{out_dir}/s29_sweep_variant_mask.csv", "w") as f:
        f.write("var_day,mask_day,rep,attack_rate\n")
        for vr in (15, 45, 200):
            for mk in (10, 40, 200):
                for rp in range(3):
                    c2 = dict(c); c2["master_seed"] += rp + 30
                    ss = Sim(c2, RunParams(nsw, dsw, beta_cal, var_day_ovr=vr, mask_ovr=mk, interv_on=1), False); sim_run(ss)
                    f.write(f"{vr},{mk},{rp},{ss.total_trans / nsw:.6f}\n")

    ens_s, ens_s2 = [0.0] * 90, [0.0] * 90
    for rp in range(20):
        c2 = dict(c); c2["master_seed"] += rp + 100
        ss = Sim(c2, RunParams(nsw, dsw, beta_cal, interv_on=1), False); sim_run(ss)
        for d in range(90):
            inc = float(sum(ss.daily_inc_age[d]))
            ens_s[d] += inc; ens_s2[d] += inc * inc
    with open(f"{out_dir}/s30_ensemble.csv", "w") as f:
        f.write("day,mean_inc,se_inc\n")
        for d in range(90):
            mn = ens_s[d] / 20.0
            f.write(f"{d},{mn:.4f},{math.sqrt(max((ens_s2[d] / 20.0 - mn * mn) / 19.0, 0.0)):.4f}\n")

    bn = ["INF_HH", "INF_SW", "INF_COMM", "FAIL_VAX", "FAIL_MASK", "FAIL_QUAR", "PROG_IA", "PROG_IS", "HOSP_ADM", "HOSP_REJ", "ICU_ADM", "ICU_REJ", "DEATH", "RECOVER", "WANING", "TEST_POS", "TEST_NEG", "QUAR_TRACE", "VAX_DOSE", "VAX_OUT", "VAR_IMPORT"]
    with open(f"{out_dir}/s31_branches.csv", "w") as f:
        f.write("branch,count\n")
        for i, nm in enumerate(bn):
            f.write(f"{nm},{s_main.branches[i]}\n")
    with open(f"{out_dir}/s32_vacuous_check.csv", "w") as f:
        f.write(f"check,ok\nneg_controls,{int(nc1.total_trans == 0 and nc2.total_trans == 0 and nc3.total_trans == 0)}\nensemble_var,{int(sum(1 for x in ens_s if x > 0) > 5)}\n")

    sp1 = Sim(c, RunParams(c["parity_N"], c["parity_days"], beta_cal, interv_on=1), False); sim_run(sp1)
    sp2 = Sim(c, RunParams(c["parity_N"], c["parity_days"], beta_cal, interv_on=1), False); sim_run(sp2)
    h1, h2 = FNV_OFFSET, FNV_OFFSET
    for i in range(sp1.n):
        h1 = fnv1a_i64(h1, sp1.pop[i].status); h2 = fnv1a_i64(h2, sp2.pop[i].status)
        h1 = fnv1a_i64(h1, sp1.pop[i].offspring); h2 = fnv1a_i64(h2, sp2.pop[i].offspring)
    with open(f"{out_dir}/s33_replay.csv", "w") as f:
        f.write(f"ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n{sp1.event_hash:016x},{sp2.event_hash:016x},{h1:016x},{h2:016x},{int(sp1.event_hash == sp2.event_hash and h1 == h2)}\n")
    with open(f"{out_dir}/s34_parity.csv", "w") as f:
        f.write(f"parity_N,parity_days,ev_hash,fin_hash,total_trans,deaths,inv_violations\n{c['parity_N']},{c['parity_days']},{sp1.event_hash:016x},{h1:016x},{sp1.total_trans},{sp1.daily_cnt[c['parity_days'] - 1][ST_D]},{sp1.inv_violations}\n")
    with open(f"{out_dir}/s36_sizes.csv", "w") as f:
        f.write(f"N,days,n_trans,n_progs\n{s_main.n},{s_main.rp.days},{len(s_main.trans)},{len(s_main.progs)}\n")

    cnt_r = [s_main.n, 0, 0, 0, 0, 0, 0, 0]
    pi, max_diff = 0, 0
    for d in range(d_cnt):
        while pi < len(s_main.progs) and s_main.progs[pi][0] == d:
            cnt_r[s_main.progs[pi][2]] -= 1
            cnt_r[s_main.progs[pi][3]] += 1
            pi += 1
        for k in range(8):
            max_diff = max(max_diff, abs(cnt_r[k] - s_main.daily_cnt[d][k]))
    with open(f"{out_dir}/s37_replay_counts.csv", "w") as f:
        f.write(f"max_diff,match\n{max_diff},{int(max_diff == 0)}\n")
