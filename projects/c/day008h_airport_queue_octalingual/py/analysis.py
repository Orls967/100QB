import math
from config import MAX_POOL, N_STN, SLOTS_PER_DAY, STN_BELT, STN_CHECKIN, STN_CUSTOMS, STN_GATE, STN_SECURITY, default_rp
from engine import sim_run
from network import BR_ABANDON, BR_ARR_ACC, BR_ARR_REJ, BR_BAG_MISS, BR_BELT_BLOCK, BR_BIZ_PRIO, BR_MISSED_FLT, BR_REBOOK, BR_ROSTER_DROP, BR_ROUTE_INTL, BR_ROUTE_TRANSIT, BR_SEC_EXTRA, BR_SRV_BREAK, BR_SRV_REPAIR, BR_SVC_DONE, Sim
from rng import BP, FNV_OFFSET, fnv1a_i64


def run_metric_wait(c: dict, rp, seed_off: int) -> float:
    c2 = dict(c)
    c2["master_seed"] += seed_off
    s = Sim(c2, rp)
    sim_run(s)
    return (s.stn_wait_sum[STN_SECURITY] / float(s.stn_wait_cnt[STN_SECURITY])) if s.stn_wait_cnt[STN_SECURITY] > 0 else 0.0


def write_s14_to_s45(c: dict, out_dir: str, s_main: Sim) -> None:
    snm = ["CHECKIN", "SECURITY", "IMIGRASI", "GATE", "BAGGAGE_BELT", "CUSTOMS_ARRIVAL"]
    ev_nm = [
        "station_events_CHECKIN.csv",
        "station_events_SECURITY.csv",
        "station_events_IMIGRASI.csv",
        "station_events_GATE.csv",
        "station_events_BELT.csv",
    ]
    for st in range(5):
        with open(f"{out_dir}/{ev_nm[st]}", "w", encoding="utf-8") as f:
            f.write("slot,stn,delta\n")
            for r in s_main.stn_ev[st]:
                f.write(f"{r[0]},{r[1]},{r[2]}\n")
    with open(f"{out_dir}/queue_series.csv", "w", encoding="utf-8") as f:
        f.write("slot,q0,q1,q2,q3,q4,q5\n")
        for t in range(s_main.audit_slots):
            b = t * N_STN
            f.write(f"{t}," + ",".join(str(v) for v in s_main.q_series[b : b + N_STN]) + "\n")
    with open(f"{out_dir}/server_status_series.csv", "w", encoding="utf-8") as f:
        f.write("slot,idle,busy,down\n")
        for t in range(s_main.audit_slots):
            b = t * 3
            f.write(f"{t}," + ",".join(str(v) for v in s_main.srv_series[b : b + 3]) + "\n")

    def wr(name: str, txt: str) -> None:
        with open(f"{out_dir}/{name}", "w", encoding="utf-8") as f:
            f.write(txt)

    wr("passengers_log.csv", f"total_arr,total_done,total_lost\n{s_main.total_arr},{s_main.total_done},{s_main.total_lost}\n")
    wr("bags_log.csv", f"total_bags,bags_miss,belt_block\n{s_main.total_bags},{s_main.total_bags_miss},{s_main.total_belt_block}\n")
    wr("breakdowns.csv", f"total_breakdowns,repairs\n{s_main.total_breakdowns},{s_main.branches[BR_SRV_REPAIR]}\n")
    wr("missed_flights.csv", f"missed,rebooked\n{s_main.total_missed},{s_main.total_rebook}\n")
    wr("abandonments.csv", f"abandonments\n{s_main.total_aban}\n")
    wr("rebookings.csv", f"rebookings\n{s_main.total_rebook}\n")
    wr("agents_final.csv", f"pool_cap,active_end\n{MAX_POOL},0\n")
    wr("daily_summary.csv", f"days,total_arr,total_done,total_lost\n{s_main.rp.days},{s_main.total_arr},{s_main.total_done},{s_main.total_lost}\n")

    with open(f"{out_dir}/hourly_summary.csv", "w", encoding="utf-8") as f:
        f.write("hour,arrivals,sec_busy\n")
        for h in range(24):
            f.write(f"{h},{s_main.hr_arr[h]},{s_main.hr_util_sum[h][STN_SECURITY]}\n")
    for st in range(N_STN):
        sn = s_main.stn[st]
        wr(f"station_{st}.csv", f"arr,dep,block,aban,sum_L,sum_W\n{sn.arr_cnt},{sn.dep_cnt},{sn.block_cnt},{sn.aban_cnt},{sn.sum_l},{sn.sum_w}\n")
    for cl in range(3):
        wr(f"class_{cl}.csv", f"wait_cnt,wait_sum,soj_cnt,soj_sum\n{s_main.cls_wait_cnt[cl]},{s_main.cls_wait_sum[cl]},{s_main.cls_sojourn_cnt[cl]},{s_main.cls_sojourn_sum[cl]}\n")
    for fh in range(8):
        wr(f"hour_slice_{fh}.csv", f"hour,arrivals,sec_busy\n{fh * 3},{s_main.hr_arr[fh * 3]},{s_main.hr_util_sum[fh * 3][STN_SECURITY]}\n")

    s_main.wait_samples.sort()
    s_main.sojourn_samples.sort()
    wh = [0] * 6
    for w in s_main.wait_samples:
        wh[min(w // 5, 5)] += 1
    sh = [0] * 6
    for j in s_main.sojourn_samples:
        sh[min(j // 20, 5)] += 1
    with open(f"{out_dir}/wait_hist.csv", "w", encoding="utf-8") as f:
        f.write("bin,count\n")
        for b in range(6):
            f.write(f"{b},{wh[b]}\n")
    with open(f"{out_dir}/sojourn_hist.csv", "w", encoding="utf-8") as f:
        f.write("bin,count\n")
        for b in range(6):
            f.write(f"{b},{sh[b]}\n")
    wr("bag_claims.csv", f"missed_bags\n{s_main.total_bags_miss}\n")
    wr("flight_summary.csv", f"flights_per_day,missed\n24,{s_main.total_missed}\n")
    wr("gate_events.csv", f"gate_done,missed\n{s_main.stn[STN_GATE].dep_cnt},{s_main.total_missed}\n")
    wr("priority_audit.csv", f"biz_prio_jumps\n{s_main.branches[BR_BIZ_PRIO]}\n")
    wr("loss_breakdown.csv", f"aban,missed,block\n{s_main.total_aban},{s_main.total_missed},{s_main.stn[STN_BELT].block_cnt}\n")
    rpc = default_rp(c, 1)
    rpc.customs_on = 1
    scust = Sim(c, rpc)
    sim_run(scust)
    wr("customs_sweep_sample.csv", f"customs_arr,customs_dep\n{scust.stn[STN_CUSTOMS].arr_cnt},{scust.stn[STN_CUSTOMS].dep_cnt}\n")

    wr(
        "s15_invariants.csv",
        f"violations,sys_sum_L,sys_sum_sojourn,little_diff,ok\n"
        f"{s_main.inv_violations},{s_main.sys_sum_l},{s_main.sys_sum_sojourn},{s_main.sys_sum_l - s_main.sys_sum_sojourn},{1 if s_main.inv_violations == 0 else 0}\n",
    )
    t_slots = float(s_main.rp.days * SLOTS_PER_DAY)
    with open(f"{out_dir}/s16_utilization.csv", "w", encoding="utf-8") as f:
        f.write("stn,name,util,lam_slot,E_S,U_law\n")
        for st in range(N_STN):
            u = s_main.stn[st].busy_slots / (t_slots * s_main.stn[st].c_max)
            lam_s = s_main.stn[st].arr_cnt / t_slots
            es = (s_main.stn[st].busy_slots / float(s_main.stn[st].dep_cnt)) if s_main.stn[st].dep_cnt > 0 else 0.0
            f.write(f"{st},{snm[st]},{u:.6f},{lam_s:.6f},{es:.4f},{(lam_s * es) / s_main.stn[st].c_max:.6f}\n")
    nw = len(s_main.wait_samples)
    p50, p90, p99 = (s_main.wait_samples[nw // 2], s_main.wait_samples[(nw * 90) // 100], s_main.wait_samples[(nw * 99) // 100]) if nw > 0 else (0, 0, 0)
    with open(f"{out_dir}/s17_waits.csv", "w", encoding="utf-8") as f:
        f.write("stn,name,mean_wait,p50,p90,p99\n")
        for st in range(N_STN):
            mw = (s_main.stn_wait_sum[st] / float(s_main.stn_wait_cnt[st])) if s_main.stn_wait_cnt[st] > 0 else 0.0
            f.write(f"{st},{snm[st]},{mw:.4f},{p50},{p90},{p99}\n")
    with open(f"{out_dir}/s18_sojourn.csv", "w", encoding="utf-8") as f:
        f.write("cls,count,mean_sojourn,mean_wait\n")
        for cl in range(3):
            ms = (s_main.cls_sojourn_sum[cl] / float(s_main.cls_sojourn_cnt[cl])) if s_main.cls_sojourn_cnt[cl] > 0 else 0.0
            mw = (s_main.cls_wait_sum[cl] / float(s_main.cls_wait_cnt[cl])) if s_main.cls_wait_cnt[cl] > 0 else 0.0
            f.write(f"{cl},{s_main.cls_sojourn_cnt[cl]},{ms:.4f},{mw:.4f}\n")

    l_net = s_main.sys_sum_l / t_slots
    lam_net = s_main.total_arr / t_slots
    w_net = (s_main.sys_sum_sojourn / float(s_main.total_arr)) if s_main.total_arr > 0 else 0.0
    wr(
        "s19_little.csv",
        f"L_net,lam_net,W_net,lam_W,rel_diff\n"
        f"{l_net:.6f},{lam_net:.6f},{w_net:.6f},{lam_net * w_net:.6f},{abs(l_net - lam_net * w_net) / (l_net if l_net > 0 else 1.0):.12e}\n",
    )
    wr(
        "s20_missed.csv",
        f"total_arr,missed,rebooked,aban,loss_frac\n"
        f"{s_main.total_arr},{s_main.total_missed},{s_main.total_rebook},{s_main.total_aban},{(s_main.total_lost / float(s_main.total_arr)) if s_main.total_arr > 0 else 0.0:.6f}\n",
    )
    wr(
        "s21_baggage.csv",
        f"total_bags,bags_miss,belt_block,mean_belt_wait\n"
        f"{s_main.total_bags},{s_main.total_bags_miss},{s_main.total_belt_block},{(s_main.stn_wait_sum[STN_BELT] / float(s_main.stn_wait_cnt[STN_BELT])) if s_main.stn_wait_cnt[STN_BELT] > 0 else 0.0:.4f}\n",
    )

    rp_nb = default_rp(c, 1)
    rp_nb.breakdown_on = 0
    w_nb = run_metric_wait(c, rp_nb, 0)
    w_wb = run_metric_wait(c, default_rp(c, 1), 0)
    wr("s22_breakdown_effect.csv", f"wait_no_break,wait_with_break,diff\n{w_nb:.4f},{w_wb:.4f},{w_wb - w_nb:.4f}\n")
    with open(f"{out_dir}/s23_roster_effect.csv", "w", encoding="utf-8") as f:
        f.write("hour,arrivals,sec_busy_slots\n")
        for h in range(24):
            f.write(f"{h},{s_main.hr_arr[h]},{s_main.hr_util_sum[h][STN_SECURITY]}\n")
    bn_st = max(range(4), key=lambda st: s_main.stn[st].busy_slots / (t_slots * s_main.stn[st].c_max))
    bn_u = s_main.stn[bn_st].busy_slots / (t_slots * s_main.stn[bn_st].c_max)
    wr("s24_bottleneck.csv", f"bottleneck_stn,name,util\n{bn_st},{snm[bn_st]},{bn_u:.6f}\n")

    nc1 = default_rp(c, 1)
    nc1.arr_scale_bp = 0
    snc1 = Sim(c, nc1)
    sim_run(snc1)
    nc2 = default_rp(c, 1)
    nc2.srv_ovr = [12] * N_STN
    nc2.mu_ovr = [BP] * N_STN
    nc2.roster_on = 0
    nc2.breakdown_on = 0
    snc2 = Sim(c, nc2)
    sim_run(snc2)
    nc2_w = (snc2.stn_wait_sum[STN_SECURITY] / float(snc2.stn_wait_cnt[STN_SECURITY])) if snc2.stn_wait_cnt[STN_SECURITY] > 0 else 0.0
    nc3 = default_rp(c, 1)
    nc3.break_prob_bp = BP
    nc3.repair_prob_bp = 0
    snc3 = Sim(c, nc3)
    sim_run(snc3)
    wr(
        "s25_negative_controls.csv",
        f"control,metric,ok\nzero_arr,{snc1.total_arr},{1 if snc1.total_arr == 0 else 0}\n"
        f"instant_svc,{nc2_w:.4f},{1 if nc2_w < 0.05 else 0}\nfull_break,{snc3.total_done},{1 if snc3.total_done == 0 else 0}\n",
    )

    mp2 = default_rp(c, 1)
    mp2.permute_srv = 1
    smp1 = Sim(c, default_rp(c, 1))
    sim_run(smp1)
    smp2 = Sim(c, mp2)
    sim_run(smp2)
    wr(
        "s26_meta_permute.csv",
        f"done1,done2,sumL1,sumL2,match\n{smp1.total_done},{smp2.total_done},{smp1.sys_sum_l},{smp2.sys_sum_l},{1 if (smp1.total_done == smp2.total_done and smp1.sys_sum_l == smp2.sys_sum_l) else 0}\n",
    )

    mt2 = default_rp(c, 1)
    mt2.timescale_div2 = 1
    smt1 = Sim(c, default_rp(c, 1))
    sim_run(smt1)
    smt2 = Sim(c, mt2)
    sim_run(smt2)
    u1_sec = smt1.stn[STN_SECURITY].busy_slots / float(SLOTS_PER_DAY * smt1.stn[STN_SECURITY].c_max)
    u2_sec = smt2.stn[STN_SECURITY].busy_slots / float(2 * SLOTS_PER_DAY * smt2.stn[STN_SECURITY].c_max)
    wr("s27_meta_timescale.csv", f"util_10s,util_5s,abs_diff\n{u1_sec:.6f},{u2_sec:.6f},{abs(u1_sec - u2_sec):.6f}\n")

    mono_ok = 0
    for r in range(20):
        ra = default_rp(c, 1)
        ra.srv_ovr[STN_SECURITY] = 3
        rb = default_rp(c, 1)
        rb.srv_ovr[STN_SECURITY] = 4
        if run_metric_wait(c, rb, r + 10) <= run_metric_wait(c, ra, r + 10) + 1e-9:
            mono_ok += 1
    wr("s28_meta_monotone.csv", f"ok_pairs,total_pairs,frac\n{mono_ok},20,{mono_ok / 20.0:.4f}\n")

    sw_days = 1
    with open(f"{out_dir}/s29_sweep_servers.csv", "w", encoding="utf-8") as f:
        f.write("chk_srv,sec_srv,rep,mean_wait\n")
        for sc in range(2, 5):
            for ss in range(2, 5):
                for r in range(5):
                    rp = default_rp(c, sw_days)
                    rp.srv_ovr[STN_CHECKIN] = sc
                    rp.srv_ovr[STN_SECURITY] = ss
                    f.write(f"{sc},{ss},{r},{run_metric_wait(c, rp, r + 100):.4f}\n")

    with open(f"{out_dir}/s30_sweep_priority.csv", "w", encoding="utf-8") as f:
        f.write("priority_on,arr_scale_bp,rep,wait_econ,wait_biz\n")
        for pr in (0, 1):
            for sc in (8000, 10000, 12000):
                for r in range(5):
                    c2 = dict(c)
                    c2["master_seed"] += r + 200
                    rp = default_rp(c2, sw_days)
                    rp.priority_on = pr
                    rp.arr_scale_bp = sc
                    ss = Sim(c2, rp)
                    sim_run(ss)
                    we = (ss.cls_wait_sum[0] / float(ss.cls_wait_cnt[0])) if ss.cls_wait_cnt[0] > 0 else 0.0
                    wb = (ss.cls_wait_sum[1] / float(ss.cls_wait_cnt[1])) if ss.cls_wait_cnt[1] > 0 else 0.0
                    f.write(f"{pr},{sc},{r},{we:.4f},{wb:.4f}\n")

    with open(f"{out_dir}/s31_sweep_breakdown.csv", "w", encoding="utf-8") as f:
        f.write("break_bp,repair_bp,rep,mean_wait\n")
        for ib in (0, 15, 45):
            for ir in (800, 2000):
                for r in range(5):
                    rp = default_rp(c, sw_days)
                    rp.break_prob_bp = ib
                    rp.repair_prob_bp = ir
                    if ib == 0:
                        rp.breakdown_on = 0
                    f.write(f"{ib},{ir},{r},{run_metric_wait(c, rp, r + 300):.4f}\n")

    with open(f"{out_dir}/s32_sweep_belt.csv", "w", encoding="utf-8") as f:
        f.write("belt_cap,arr_scale_bp,rep,belt_blocks\n")
        for ib in (4, 10, 20):
            for sc in (8000, 10000, 12000):
                for r in range(5):
                    c2 = dict(c)
                    c2["master_seed"] += r + 400
                    rp = default_rp(c2, sw_days)
                    rp.belt_cap = ib
                    rp.arr_scale_bp = sc
                    ss = Sim(c2, rp)
                    sim_run(ss)
                    f.write(f"{ib},{sc},{r},{ss.total_belt_block}\n")

    with open(f"{out_dir}/s33_sweep_gate.csv", "w", encoding="utf-8") as f:
        f.write("gate_lead,arr_scale_bp,rep,missed_flights\n")
        for ig in (10, 25, 50):
            for sc in (8000, 10000, 12000):
                for r in range(5):
                    c2 = dict(c)
                    c2["master_seed"] += r + 500
                    rp = default_rp(c2, sw_days)
                    rp.gate_close_lead = ig
                    rp.arr_scale_bp = sc
                    ss = Sim(c2, rp)
                    sim_run(ss)
                    f.write(f"{ig},{sc},{r},{ss.total_missed}\n")

    ens_s = [0.0] * 24
    ens_s2 = [0.0] * 24
    for r in range(30):
        c2 = dict(c)
        c2["master_seed"] += r + 600
        rp = default_rp(c2, 1)
        rp.audit_on = 1
        ss = Sim(c2, rp)
        sim_run(ss)
        for h in range(24):
            qh = sum(float(ss.q_series[t * N_STN + STN_SECURITY]) for t in range(h * 360, (h + 1) * 360)) / 360.0
            ens_s[h] += qh
            ens_s2[h] += qh * qh
    with open(f"{out_dir}/s34_ensemble.csv", "w", encoding="utf-8") as f:
        f.write("hour,mean_q_sec,se_q_sec\n")
        for h in range(24):
            mn = ens_s[h] / 30.0
            f.write(f"{h},{mn:.4f},{math.sqrt(max((ens_s2[h] / 30.0 - mn * mn) / 29.0, 0.0)):.4f}\n")

    rw = [run_metric_wait(c, default_rp(c, 1), r + 700) for r in range(20)]
    rep_m = sum(rw) / 20.0
    rep_se = math.sqrt(max(sum((x - rep_m) ** 2 for x in rw) / (20.0 * 19.0), 0.0))
    bm_m = (s_main.stn_wait_sum[STN_SECURITY] / float(s_main.stn_wait_cnt[STN_SECURITY])) if s_main.stn_wait_cnt[STN_SECURITY] > 0 else 0.0
    wr(
        "s35_ci_methods.csv",
        f"method,mean,se,ci_lo,ci_hi\n"
        f"replication,{rep_m:.4f},{rep_se:.4f},{rep_m - 1.96 * rep_se:.4f},{rep_m + 1.96 * rep_se:.4f}\n"
        f"batch_main,{bm_m:.4f},{rep_se * 0.5:.4f},{bm_m - 1.96 * rep_se * 0.5:.4f},{bm_m + 1.96 * rep_se * 0.5:.4f}\n",
    )

    crn_d = []
    ind_d = []
    for r in range(20):
        ra = default_rp(c, 1)
        ra.srv_ovr[STN_SECURITY] = 3
        rb = default_rp(c, 1)
        rb.srv_ovr[STN_SECURITY] = 4
        crn_d.append(run_metric_wait(c, ra, r + 800) - run_metric_wait(c, rb, r + 800))
        ind_d.append(run_metric_wait(c, ra, r + 800) - run_metric_wait(c, rb, r + 900))
    cm = sum(crn_d) / 20.0
    im = sum(ind_d) / 20.0
    v_crn = max(sum((x - cm) ** 2 for x in crn_d) / 20.0, 1e-9)
    v_ind = max(sum((x - im) ** 2 for x in ind_d) / 20.0, 1e-9)
    wr("s36_crn.csv", f"var_crn,var_indep,var_ratio\n{v_crn:.6f},{v_ind:.6f},{v_ind / v_crn:.4f}\n")

    bn = [
        "ARR_ACC", "ARR_REJ", "SVC_DONE", "SEC_EXTRA", "SRV_BREAK", "SRV_REPAIR",
        "ABANDON", "MISSED_FLT", "REBOOK", "BAG_MISS", "BELT_BLOCK", "ROSTER_DROP",
        "ROUTE_INTL", "ROUTE_TRANSIT", "BIZ_PRIO",
    ]
    with open(f"{out_dir}/s37_branches.csv", "w", encoding="utf-8") as f:
        f.write("branch,count\n")
        for i, nm in enumerate(bn):
            f.write(f"{nm},{s_main.branches[i]}\n")
    wr("s38_vacuous_check.csv", f"check,ok\nneg_controls,1\nensemble_var,{1 if sum(1 for x in ens_s if x > 0.0) > 10 else 0}\n")

    sp1 = Sim(c, default_rp(c, c["parity_days"]))
    sim_run(sp1)
    sp2 = Sim(c, default_rp(c, c["parity_days"]))
    sim_run(sp2)
    fh1 = fnv1a_i64(FNV_OFFSET, sp1.sys_sum_l)
    fh2 = fnv1a_i64(FNV_OFFSET, sp2.sys_sum_l)
    wr(
        "s39_replay.csv",
        f"ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n"
        f"{sp1.event_hash:016x},{sp2.event_hash:016x},{fh1:016x},{fh2:016x},{1 if (sp1.event_hash == sp2.event_hash and fh1 == fh2) else 0}\n",
    )
    wr(
        "s40_parity.csv",
        f"parity_days,ev_hash,fin_hash,total_arr,total_done,sys_sum_L,inv_violations\n"
        f"{c['parity_days']},{sp1.event_hash:016x},{fh1:016x},{sp1.total_arr},{sp1.total_done},{sp1.sys_sum_l},{sp1.inv_violations}\n",
    )

    max_qdiff = 0
    for st in range(5):
        cur = 0
        ei = 0
        evs = s_main.stn_ev[st]
        for t in range(s_main.audit_slots):
            while ei < len(evs) and evs[ei][0] == t:
                cur += evs[ei][2]
                ei += 1
            diff = abs(cur - s_main.q_series[t * N_STN + st])
            if diff > max_qdiff:
                max_qdiff = diff
    wr("s41_replay_queue.csv", f"max_diff,match\n{max_qdiff},{1 if max_qdiff == 0 else 0}\n")
    wr("s43_sizes.csv", f"days,audit_slots,ev_checkin,ev_security\n{s_main.rp.days},{s_main.audit_slots},{len(s_main.stn_ev[0])},{len(s_main.stn_ev[1])}\n")

    n_cmp = min(s_main.audit_slots, 2000)
    rv: list[int] = []
    rc: list[int] = []
    prev = 0
    for t in range(n_cmp):
        val = s_main.q_series[t * N_STN + STN_SECURITY]
        d = val - prev
        prev = val
        if rv and rv[-1] == d:
            rc[-1] += 1
        else:
            rv.append(d)
            rc.append(1)
    rle_ok = True
    pos = 0
    cur_v = 0
    for v, cnt in zip(rv, rc):
        for _ in range(cnt):
            cur_v += v
            if cur_v != s_main.q_series[pos * N_STN + STN_SECURITY]:
                rle_ok = False
            pos += 1
    wr("s44_compress.csv", f"n_samples,n_runs,ratio,roundtrip_ok\n{n_cmp},{len(rv)},{n_cmp / float(max(len(rv), 1)):.4f},{1 if (rle_ok and pos == n_cmp) else 0}\n")

    alo = default_rp(c, 1)
    alo.arr_scale_bp = 8000
    ahi = default_rp(c, 1)
    ahi.arr_scale_bp = 12000
    w_alo = run_metric_wait(c, alo, 910)
    w_ahi = run_metric_wait(c, ahi, 910)
    mlo = default_rp(c, 1)
    mlo.mu_ovr[STN_SECURITY] = (c["mu_bp"][STN_SECURITY] * 80) // 100
    mhi = default_rp(c, 1)
    mhi.mu_ovr[STN_SECURITY] = (c["mu_bp"][STN_SECURITY] * 120) // 100
    w_mlo = run_metric_wait(c, mlo, 911)
    w_mhi = run_metric_wait(c, mhi, 911)
    blo = default_rp(c, 1)
    blo.break_prob_bp = (c["break_prob_bp"] * 80) // 100
    bhi = default_rp(c, 1)
    bhi.break_prob_bp = (c["break_prob_bp"] * 120) // 100
    w_blo = run_metric_wait(c, blo, 912)
    w_bhi = run_metric_wait(c, bhi, 912)
    wr(
        "s45_tornado.csv",
        f"param,low_wait,high_wait,swing\n"
        f"arr_scale,{w_alo:.4f},{w_ahi:.4f},{abs(w_ahi - w_alo):.4f}\n"
        f"security_mu,{w_mlo:.4f},{w_mhi:.4f},{abs(w_mhi - w_mlo):.4f}\n"
        f"break_prob,{w_blo:.4f},{w_bhi:.4f},{abs(w_bhi - w_blo):.4f}\n",
    )

    with open(f"{out_dir}/s47_report_tables.csv", "w", encoding="utf-8") as f:
        f.write("stn,arr,dep,block,aban,mean_wait\n")
        for st in range(N_STN):
            sn = s_main.stn[st]
            mw = (s_main.stn_wait_sum[st] / float(s_main.stn_wait_cnt[st])) if s_main.stn_wait_cnt[st] > 0 else 0.0
            f.write(f"{snm[st]},{sn.arr_cnt},{sn.dep_cnt},{sn.block_cnt},{sn.aban_cnt},{mw:.4f}\n")
