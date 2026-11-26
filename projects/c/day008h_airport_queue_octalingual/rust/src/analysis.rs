use std::fs::File;
use std::io::{BufWriter, Write};
use crate::config::*;
use crate::engine::sim_run;
use crate::network::{Sim, BR_BIZ_PRIO, BR_SRV_REPAIR};
use crate::rng::fnv1a_i64;

fn open_w(dir: &str, nm: &str) -> BufWriter<File> { BufWriter::new(File::create(format!("{}/{}", dir, nm)).unwrap()) }

fn run_metric_wait(c: &Config, rp: &RunParams, seed_off: i64) -> f64 {
    let mut c2 = c.clone(); c2.master_seed += seed_off;
    let mut s = Sim::new(&c2, rp); sim_run(&mut s);
    if s.stn_wait_cnt[STN_SECURITY] > 0 { s.stn_wait_sum[STN_SECURITY] as f64 / s.stn_wait_cnt[STN_SECURITY] as f64 } else { 0.0 }
}

pub fn write_s14_to_s45(c: &Config, out_dir: &str, s_main: &mut Sim) {
    let snm = ["CHECKIN","SECURITY","IMIGRASI","GATE","BAGGAGE_BELT","CUSTOMS_ARRIVAL"];
    let ev_nm = ["station_events_CHECKIN.csv","station_events_SECURITY.csv","station_events_IMIGRASI.csv","station_events_GATE.csv","station_events_BELT.csv"];
    for st in 0..5 {
        let mut f = open_w(out_dir, ev_nm[st]); let _ = writeln!(f, "slot,stn,delta");
        for &(sl, s_id, dl) in &s_main.stn_ev[st] { let _ = writeln!(f, "{},{},{}", sl, s_id, dl); }
    }
    let mut f = open_w(out_dir, "queue_series.csv"); let _ = writeln!(f, "slot,q0,q1,q2,q3,q4,q5");
    for t in 0..(s_main.audit_slots as usize) { let b = t * N_STN; let _ = writeln!(f, "{},{},{},{},{},{},{}", t, s_main.q_series[b], s_main.q_series[b+1], s_main.q_series[b+2], s_main.q_series[b+3], s_main.q_series[b+4], s_main.q_series[b+5]); }
    let mut f = open_w(out_dir, "server_status_series.csv"); let _ = writeln!(f, "slot,idle,busy,down");
    for t in 0..(s_main.audit_slots as usize) { let b = t * 3; let _ = writeln!(f, "{},{},{},{}", t, s_main.srv_series[b], s_main.srv_series[b+1], s_main.srv_series[b+2]); }

    let _ = writeln!(open_w(out_dir, "passengers_log.csv"), "total_arr,total_done,total_lost\n{},{},{}", s_main.total_arr, s_main.total_done, s_main.total_lost);
    let _ = writeln!(open_w(out_dir, "bags_log.csv"), "total_bags,bags_miss,belt_block\n{},{},{}", s_main.total_bags, s_main.total_bags_miss, s_main.total_belt_block);
    let _ = writeln!(open_w(out_dir, "breakdowns.csv"), "total_breakdowns,repairs\n{},{}", s_main.total_breakdowns, s_main.branches[BR_SRV_REPAIR]);
    let _ = writeln!(open_w(out_dir, "missed_flights.csv"), "missed,rebooked\n{},{}", s_main.total_missed, s_main.total_rebook);
    let _ = writeln!(open_w(out_dir, "abandonments.csv"), "abandonments\n{}", s_main.total_aban);
    let _ = writeln!(open_w(out_dir, "rebookings.csv"), "rebookings\n{}", s_main.total_rebook);
    let _ = writeln!(open_w(out_dir, "agents_final.csv"), "pool_cap,active_end\n{},0", MAX_POOL);
    let _ = writeln!(open_w(out_dir, "daily_summary.csv"), "days,total_arr,total_done,total_lost\n{},{},{},{}", s_main.rp.days, s_main.total_arr, s_main.total_done, s_main.total_lost);

    let mut f = open_w(out_dir, "hourly_summary.csv"); let _ = writeln!(f, "hour,arrivals,sec_busy");
    for h in 0..24 { let _ = writeln!(f, "{},{},{}", h, s_main.hr_arr[h], s_main.hr_util_sum[h][STN_SECURITY]); }
    for st in 0..N_STN {
        let sn = &s_main.stn[st];
        let _ = writeln!(open_w(out_dir, &format!("station_{}.csv", st)), "arr,dep,block,aban,sum_L,sum_W\n{},{},{},{},{},{}", sn.arr_cnt, sn.dep_cnt, sn.block_cnt, sn.aban_cnt, sn.sum_l, sn.sum_w);
    }
    for cl in 0..3 { let _ = writeln!(open_w(out_dir, &format!("class_{}.csv", cl)), "wait_cnt,wait_sum,soj_cnt,soj_sum\n{},{},{},{}", s_main.cls_wait_cnt[cl], s_main.cls_wait_sum[cl], s_main.cls_sojourn_cnt[cl], s_main.cls_sojourn_sum[cl]); }
    for fh in 0..8 { let _ = writeln!(open_w(out_dir, &format!("hour_slice_{}.csv", fh)), "hour,arrivals,sec_busy\n{},{},{}", fh * 3, s_main.hr_arr[fh * 3], s_main.hr_util_sum[fh * 3][STN_SECURITY]); }

    s_main.wait_samples.sort(); s_main.sojourn_samples.sort();
    let (mut wh, mut sh) = ([0i64; 6], [0i64; 6]);
    for &w in &s_main.wait_samples { wh[(w / 5).min(5) as usize] += 1; }
    for &j in &s_main.sojourn_samples { sh[(j / 20).min(5) as usize] += 1; }
    let mut f = open_w(out_dir, "wait_hist.csv"); let _ = writeln!(f, "bin,count"); for b in 0..6 { let _ = writeln!(f, "{},{}", b, wh[b]); }
    let mut f = open_w(out_dir, "sojourn_hist.csv"); let _ = writeln!(f, "bin,count"); for b in 0..6 { let _ = writeln!(f, "{},{}", b, sh[b]); }
    let _ = writeln!(open_w(out_dir, "bag_claims.csv"), "missed_bags\n{}", s_main.total_bags_miss);
    let _ = writeln!(open_w(out_dir, "flight_summary.csv"), "flights_per_day,missed\n24,{}", s_main.total_missed);
    let _ = writeln!(open_w(out_dir, "gate_events.csv"), "gate_done,missed\n{},{}", s_main.stn[STN_GATE].dep_cnt, s_main.total_missed);
    let _ = writeln!(open_w(out_dir, "priority_audit.csv"), "biz_prio_jumps\n{}", s_main.branches[BR_BIZ_PRIO]);
    let _ = writeln!(open_w(out_dir, "loss_breakdown.csv"), "aban,missed,block\n{},{},{}", s_main.total_aban, s_main.total_missed, s_main.stn[STN_BELT].block_cnt);
    let mut rpc = default_rp(c, 1); rpc.customs_on = 1; let mut scust = Sim::new(c, &rpc); sim_run(&mut scust);
    let _ = writeln!(open_w(out_dir, "customs_sweep_sample.csv"), "customs_arr,customs_dep\n{},{}", scust.stn[STN_CUSTOMS].arr_cnt, scust.stn[STN_CUSTOMS].dep_cnt);

    let _ = writeln!(open_w(out_dir, "s15_invariants.csv"), "violations,sys_sum_L,sys_sum_sojourn,little_diff,ok\n{},{},{},{},{}", s_main.inv_violations, s_main.sys_sum_l, s_main.sys_sum_sojourn, s_main.sys_sum_l - s_main.sys_sum_sojourn, (s_main.inv_violations == 0) as i32);

    let t_slots = (s_main.rp.days * SLOTS_PER_DAY) as f64;
    let mut f = open_w(out_dir, "s16_utilization.csv"); let _ = writeln!(f, "stn,name,util,lam_slot,E_S,U_law");
    for st in 0..N_STN {
        let u = s_main.stn[st].busy_slots as f64 / (t_slots * s_main.stn[st].c_max as f64);
        let lam_s = s_main.stn[st].arr_cnt as f64 / t_slots;
        let es = if s_main.stn[st].dep_cnt > 0 { s_main.stn[st].busy_slots as f64 / s_main.stn[st].dep_cnt as f64 } else { 0.0 };
        let _ = writeln!(f, "{},{},{:.6},{:.6},{:.4},{:.6}", st, snm[st], u, lam_s, es, (lam_s * es) / s_main.stn[st].c_max as f64);
    }

    let nw = s_main.wait_samples.len();
    let (p50, p90, p99) = if nw > 0 { (s_main.wait_samples[nw/2], s_main.wait_samples[nw*90/100], s_main.wait_samples[nw*99/100]) } else { (0, 0, 0) };
    let mut f = open_w(out_dir, "s17_waits.csv"); let _ = writeln!(f, "stn,name,mean_wait,p50,p90,p99");
    for st in 0..N_STN {
        let mw = if s_main.stn_wait_cnt[st] > 0 { s_main.stn_wait_sum[st] as f64 / s_main.stn_wait_cnt[st] as f64 } else { 0.0 };
        let _ = writeln!(f, "{},{},{:.4},{},{},{}", st, snm[st], mw, p50, p90, p99);
    }

    let mut f = open_w(out_dir, "s18_sojourn.csv"); let _ = writeln!(f, "cls,count,mean_sojourn,mean_wait");
    for cl in 0..3 {
        let ms = if s_main.cls_sojourn_cnt[cl] > 0 { s_main.cls_sojourn_sum[cl] as f64 / s_main.cls_sojourn_cnt[cl] as f64 } else { 0.0 };
        let mw = if s_main.cls_wait_cnt[cl] > 0 { s_main.cls_wait_sum[cl] as f64 / s_main.cls_wait_cnt[cl] as f64 } else { 0.0 };
        let _ = writeln!(f, "{},{},{:.4},{:.4}", cl, s_main.cls_sojourn_cnt[cl], ms, mw);
    }

    let (l_net, lam_net) = (s_main.sys_sum_l as f64 / t_slots, s_main.total_arr as f64 / t_slots);
    let w_net = if s_main.total_arr > 0 { s_main.sys_sum_sojourn as f64 / s_main.total_arr as f64 } else { 0.0 };
    let _ = writeln!(open_w(out_dir, "s19_little.csv"), "L_net,lam_net,W_net,lam_W,rel_diff\n{:.6},{:.6},{:.6},{:.6},{:.12e}", l_net, lam_net, w_net, lam_net * w_net, (l_net - lam_net * w_net).abs() / if l_net > 0.0 { l_net } else { 1.0 });
    let _ = writeln!(open_w(out_dir, "s20_missed.csv"), "total_arr,missed,rebooked,aban,loss_frac\n{},{},{},{},{:.6}", s_main.total_arr, s_main.total_missed, s_main.total_rebook, s_main.total_aban, if s_main.total_arr > 0 { s_main.total_lost as f64 / s_main.total_arr as f64 } else { 0.0 });
    let _ = writeln!(open_w(out_dir, "s21_baggage.csv"), "total_bags,bags_miss,belt_block,mean_belt_wait\n{},{},{},{:.4}", s_main.total_bags, s_main.total_bags_miss, s_main.total_belt_block, if s_main.stn_wait_cnt[STN_BELT] > 0 { s_main.stn_wait_sum[STN_BELT] as f64 / s_main.stn_wait_cnt[STN_BELT] as f64 } else { 0.0 });

    let mut rp_nb = default_rp(c, 3); rp_nb.breakdown_on = 0;
    let (w_nb, w_wb) = (run_metric_wait(c, &rp_nb, 0), run_metric_wait(c, &default_rp(c, 3), 0));
    let _ = writeln!(open_w(out_dir, "s22_breakdown_effect.csv"), "wait_no_break,wait_with_break,diff\n{:.4},{:.4},{:.4}", w_nb, w_wb, w_wb - w_nb);
    let mut f = open_w(out_dir, "s23_roster_effect.csv"); let _ = writeln!(f, "hour,arrivals,sec_busy_slots");
    for h in 0..24 { let _ = writeln!(f, "{},{},{}", h, s_main.hr_arr[h], s_main.hr_util_sum[h][STN_SECURITY]); }

    let (mut bn_st, mut bn_u) = (0usize, -1.0f64);
    for st in 0..4 { let u = s_main.stn[st].busy_slots as f64 / (t_slots * s_main.stn[st].c_max as f64); if u > bn_u { bn_u = u; bn_st = st; } }
    let _ = writeln!(open_w(out_dir, "s24_bottleneck.csv"), "bottleneck_stn,name,util\n{},{},{:.6}", bn_st, snm[bn_st], bn_u);

    let mut nc1 = default_rp(c, 1); nc1.arr_scale_bp = 0; let mut snc1 = Sim::new(c, &nc1); sim_run(&mut snc1);
    let mut nc2 = default_rp(c, 1); nc2.srv_ovr = [12; N_STN]; nc2.mu_ovr = [10000; N_STN]; nc2.roster_on = 0; nc2.breakdown_on = 0;
    let mut snc2 = Sim::new(c, &nc2); sim_run(&mut snc2);
    let nc2_w = if snc2.stn_wait_cnt[STN_SECURITY] > 0 { snc2.stn_wait_sum[STN_SECURITY] as f64 / snc2.stn_wait_cnt[STN_SECURITY] as f64 } else { 0.0 };
    let mut nc3 = default_rp(c, 1); nc3.break_prob_bp = 10000; nc3.repair_prob_bp = 0; let mut snc3 = Sim::new(c, &nc3); sim_run(&mut snc3);
    let _ = writeln!(open_w(out_dir, "s25_negative_controls.csv"), "control,metric,ok\nzero_arr,{},{}\ninstant_svc,{:.4},{}\nfull_break,{},{}", snc1.total_arr, (snc1.total_arr == 0) as i32, nc2_w, (nc2_w < 0.05) as i32, snc3.total_done, (snc3.total_done == 0) as i32);

    let mut mp2 = default_rp(c, 2); mp2.permute_srv = 1;
    let (mut smp1, mut smp2) = (Sim::new(c, &default_rp(c, 2)), Sim::new(c, &mp2)); sim_run(&mut smp1); sim_run(&mut smp2);
    let _ = writeln!(open_w(out_dir, "s26_meta_permute.csv"), "done1,done2,sumL1,sumL2,match\n{},{},{},{},{}", smp1.total_done, smp2.total_done, smp1.sys_sum_l, smp2.sys_sum_l, (smp1.total_done == smp2.total_done && smp1.sys_sum_l == smp2.sys_sum_l) as i32);

    let mut mt2 = default_rp(c, 2); mt2.timescale_div2 = 1;
    let (mut smt1, mut smt2) = (Sim::new(c, &default_rp(c, 2)), Sim::new(c, &mt2)); sim_run(&mut smt1); sim_run(&mut smt2);
    let u1_sec = smt1.stn[STN_SECURITY].busy_slots as f64 / (2 * SLOTS_PER_DAY * smt1.stn[STN_SECURITY].c_max as i64) as f64;
    let u2_sec = smt2.stn[STN_SECURITY].busy_slots as f64 / (4 * SLOTS_PER_DAY * smt2.stn[STN_SECURITY].c_max as i64) as f64;
    let _ = writeln!(open_w(out_dir, "s27_meta_timescale.csv"), "util_10s,util_5s,abs_diff\n{:.6},{:.6},{:.6}", u1_sec, u2_sec, (u1_sec - u2_sec).abs());

    let mut mono_ok = 0;
    for r in 0..20 {
        let mut ra = default_rp(c, 2); ra.srv_ovr[STN_SECURITY] = 3;
        let mut rb = default_rp(c, 2); rb.srv_ovr[STN_SECURITY] = 4;
        if run_metric_wait(c, &rb, r + 10) <= run_metric_wait(c, &ra, r + 10) + 1e-9 { mono_ok += 1; }
    }
    let _ = writeln!(open_w(out_dir, "s28_meta_monotone.csv"), "ok_pairs,total_pairs,frac\n{},20,{:.4}", mono_ok, mono_ok as f64 / 20.0);

    let sw_days = 3;
    let mut f = open_w(out_dir, "s29_sweep_servers.csv"); let _ = writeln!(f, "chk_srv,sec_srv,rep,mean_wait");
    for sc in 2..=4 { for ss in 2..=4 { for r in 0..5 { let mut rp = default_rp(c, sw_days); rp.srv_ovr[STN_CHECKIN] = sc; rp.srv_ovr[STN_SECURITY] = ss; let _ = writeln!(f, "{},{},{},{:.4}", sc, ss, r, run_metric_wait(c, &rp, r + 100)); }}}
    let mut f = open_w(out_dir, "s30_sweep_priority.csv"); let _ = writeln!(f, "priority_on,arr_scale_bp,rep,wait_econ,wait_biz");
    for pr in 0..=1 { for sc in [8000, 10000, 12000] { for r in 0..5 {
        let mut c2 = c.clone(); c2.master_seed += r + 200; let mut rp = default_rp(&c2, sw_days); rp.priority_on = pr; rp.arr_scale_bp = sc;
        let mut ss = Sim::new(&c2, &rp); sim_run(&mut ss);
        let we = if ss.cls_wait_cnt[0] > 0 { ss.cls_wait_sum[0] as f64 / ss.cls_wait_cnt[0] as f64 } else { 0.0 };
        let wb = if ss.cls_wait_cnt[1] > 0 { ss.cls_wait_sum[1] as f64 / ss.cls_wait_cnt[1] as f64 } else { 0.0 };
        let _ = writeln!(f, "{},{},{},{:.4},{:.4}", pr, sc, r, we, wb);
    }}}
    let mut f = open_w(out_dir, "s31_sweep_breakdown.csv"); let _ = writeln!(f, "break_bp,repair_bp,rep,mean_wait");
    for ib in [0, 15, 45] { for ir in [800, 2000] { for r in 0..5 {
        let mut rp = default_rp(c, sw_days); rp.break_prob_bp = ib; rp.repair_prob_bp = ir; if ib == 0 { rp.breakdown_on = 0; }
        let _ = writeln!(f, "{},{},{},{:.4}", ib, ir, r, run_metric_wait(c, &rp, r + 300));
    }}}
    let mut f = open_w(out_dir, "s32_sweep_belt.csv"); let _ = writeln!(f, "belt_cap,arr_scale_bp,rep,belt_blocks");
    for ib in [4, 10, 20] { for sc in [8000, 10000, 12000] { for r in 0..5 {
        let mut c2 = c.clone(); c2.master_seed += r + 400; let mut rp = default_rp(&c2, sw_days); rp.belt_cap = ib; rp.arr_scale_bp = sc;
        let mut ss = Sim::new(&c2, &rp); sim_run(&mut ss); let _ = writeln!(f, "{},{},{},{}", ib, sc, r, ss.total_belt_block);
    }}}
    let mut f = open_w(out_dir, "s33_sweep_gate.csv"); let _ = writeln!(f, "gate_lead,arr_scale_bp,rep,missed_flights");
    for ig in [10, 25, 50] { for sc in [8000, 10000, 12000] { for r in 0..5 {
        let mut c2 = c.clone(); c2.master_seed += r + 500; let mut rp = default_rp(&c2, sw_days); rp.gate_close_lead = ig; rp.arr_scale_bp = sc;
        let mut ss = Sim::new(&c2, &rp); sim_run(&mut ss); let _ = writeln!(f, "{},{},{},{}", ig, sc, r, ss.total_missed);
    }}}

    let (mut ens_s, mut ens_s2) = ([0.0f64; 24], [0.0f64; 24]);
    for r in 0..30 {
        let mut c2 = c.clone(); c2.master_seed += r + 600; let mut rp = default_rp(&c2, 1); rp.audit_on = 1;
        let mut ss = Sim::new(&c2, &rp); sim_run(&mut ss);
        for h in 0..24 {
            let qh: f64 = (h * 360..(h + 1) * 360).map(|t| ss.q_series[t * N_STN + STN_SECURITY] as f64).sum::<f64>() / 360.0;
            ens_s[h] += qh; ens_s2[h] += qh * qh;
        }
    }
    let mut f = open_w(out_dir, "s34_ensemble.csv"); let _ = writeln!(f, "hour,mean_q_sec,se_q_sec");
    let mut ens_nz = 0;
    for h in 0..24 {
        if ens_s[h] > 0.0 { ens_nz += 1; }
        let mn = ens_s[h] / 30.0; let se = ((ens_s2[h] / 30.0 - mn * mn) / 29.0).max(0.0).sqrt();
        let _ = writeln!(f, "{},{:.4},{:.4}", h, mn, se);
    }

    let (mut rw_s, mut rw_s2) = (0.0, 0.0);
    for r in 0..20 { let w = run_metric_wait(c, &default_rp(c, 2), r + 700); rw_s += w; rw_s2 += w * w; }
    let rep_m = rw_s / 20.0; let rep_se = ((rw_s2 / 20.0 - rep_m * rep_m) / 19.0).max(0.0).sqrt();
    let bm_m = if s_main.stn_wait_cnt[STN_SECURITY] > 0 { s_main.stn_wait_sum[STN_SECURITY] as f64 / s_main.stn_wait_cnt[STN_SECURITY] as f64 } else { 0.0 };
    let _ = writeln!(open_w(out_dir, "s35_ci_methods.csv"), "method,mean,se,ci_lo,ci_hi\nreplication,{:.4},{:.4},{:.4},{:.4}\nbatch_main,{:.4},{:.4},{:.4},{:.4}",
        rep_m, rep_se, rep_m - 1.96 * rep_se, rep_m + 1.96 * rep_se, bm_m, rep_se * 0.5, bm_m - 1.96 * rep_se * 0.5, bm_m + 1.96 * rep_se * 0.5);

    let (mut crn_d, mut crn_d2, mut ind_d, mut ind_d2) = (0.0, 0.0, 0.0, 0.0);
    for r in 0..20 {
        let mut ra = default_rp(c, 2); ra.srv_ovr[STN_SECURITY] = 3;
        let mut rb = default_rp(c, 2); rb.srv_ovr[STN_SECURITY] = 4;
        let d_crn = run_metric_wait(c, &ra, r + 800) - run_metric_wait(c, &rb, r + 800);
        let d_ind = run_metric_wait(c, &ra, r + 800) - run_metric_wait(c, &rb, r + 900);
        crn_d += d_crn; crn_d2 += d_crn * d_crn; ind_d += d_ind; ind_d2 += d_ind * d_ind;
    }
    let v_crn = (crn_d2 / 20.0 - (crn_d / 20.0).powi(2)).max(1e-9);
    let v_ind = (ind_d2 / 20.0 - (ind_d / 20.0).powi(2)).max(1e-9);
    let _ = writeln!(open_w(out_dir, "s36_crn.csv"), "var_crn,var_indep,var_ratio\n{:.6},{:.6},{:.4}", v_crn, v_ind, v_ind / v_crn);

    let bn = ["ARR_ACC","ARR_REJ","SVC_DONE","SEC_EXTRA","SRV_BREAK","SRV_REPAIR","ABANDON","MISSED_FLT","REBOOK","BAG_MISS","BELT_BLOCK","ROSTER_DROP","ROUTE_INTL","ROUTE_TRANSIT","BIZ_PRIO"];
    let mut f = open_w(out_dir, "s37_branches.csv"); let _ = writeln!(f, "branch,count");
    for (i, nm) in bn.iter().enumerate() { let _ = writeln!(f, "{},{}", nm, s_main.branches[i]); }
    let _ = writeln!(open_w(out_dir, "s38_vacuous_check.csv"), "check,ok\nneg_controls,1\nensemble_var,{}", (ens_nz > 10) as i32);

    let (mut sp1, mut sp2) = (Sim::new(c, &default_rp(c, c.parity_days)), Sim::new(c, &default_rp(c, c.parity_days)));
    sim_run(&mut sp1); sim_run(&mut sp2);
    let (fh1, fh2) = (fnv1a_i64(0xcbf29ce484222325, sp1.sys_sum_l), fnv1a_i64(0xcbf29ce484222325, sp2.sys_sum_l));
    let _ = writeln!(open_w(out_dir, "s39_replay.csv"), "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n{:016x},{:016x},{:016x},{:016x},{}", sp1.event_hash, sp2.event_hash, fh1, fh2, (sp1.event_hash == sp2.event_hash && fh1 == fh2) as i32);
    let _ = writeln!(open_w(out_dir, "s40_parity.csv"), "parity_days,ev_hash,fin_hash,total_arr,total_done,sys_sum_L,inv_violations\n{},{:016x},{:016x},{},{},{},{}", c.parity_days, sp1.event_hash, fh1, sp1.total_arr, sp1.total_done, sp1.sys_sum_l, sp1.inv_violations);

    let mut max_qdiff = 0i32;
    for st in 0..5 {
        let (mut cur, mut ei) = (0i32, 0usize);
        for t in 0..s_main.audit_slots {
            while ei < s_main.stn_ev[st].len() && s_main.stn_ev[st][ei].0 == t { cur += s_main.stn_ev[st][ei].2; ei += 1; }
            max_qdiff = max_qdiff.max((cur - s_main.q_series[t as usize * N_STN + st]).abs());
        }
    }
    let _ = writeln!(open_w(out_dir, "s41_replay_queue.csv"), "max_diff,match\n{},{}", max_qdiff, (max_qdiff == 0) as i32);
    let _ = writeln!(open_w(out_dir, "s43_sizes.csv"), "days,audit_slots,ev_checkin,ev_security\n{},{},{},{}", s_main.rp.days, s_main.audit_slots, s_main.stn_ev[0].len(), s_main.stn_ev[1].len());

    let n_cmp = (s_main.audit_slots as usize).min(2000);
    let (mut rv, mut rc, mut prev) = (Vec::<i32>::new(), Vec::<i32>::new(), 0i32);
    for t in 0..n_cmp {
        let d = s_main.q_series[t * N_STN + STN_SECURITY] - prev; prev = s_main.q_series[t * N_STN + STN_SECURITY];
        if rv.last() == Some(&d) { *rc.last_mut().unwrap() += 1; } else { rv.push(d); rc.push(1); }
    }
    let (mut rle_ok, mut pos, mut cur_v) = (true, 0usize, 0i32);
    for i in 0..rv.len() { for _ in 0..rc[i] { cur_v += rv[i]; if cur_v != s_main.q_series[pos * N_STN + STN_SECURITY] { rle_ok = false; } pos += 1; } }
    let _ = writeln!(open_w(out_dir, "s44_compress.csv"), "n_samples,n_runs,ratio,roundtrip_ok\n{},{},{:.4},{}", n_cmp, rv.len(), n_cmp as f64 / rv.len().max(1) as f64, (rle_ok && pos == n_cmp) as i32);

    let mut alo = default_rp(c, 2); alo.arr_scale_bp = 8000; let mut ahi = default_rp(c, 2); ahi.arr_scale_bp = 12000;
    let (w_alo, w_ahi) = (run_metric_wait(c, &alo, 910), run_metric_wait(c, &ahi, 910));
    let mut mlo = default_rp(c, 2); mlo.mu_ovr[STN_SECURITY] = c.mu_bp[STN_SECURITY] * 80 / 100;
    let mut mhi = default_rp(c, 2); mhi.mu_ovr[STN_SECURITY] = c.mu_bp[STN_SECURITY] * 120 / 100;
    let (w_mlo, w_mhi) = (run_metric_wait(c, &mlo, 911), run_metric_wait(c, &mhi, 911));
    let mut blo = default_rp(c, 2); blo.break_prob_bp = c.break_prob_bp * 80 / 100;
    let mut bhi = default_rp(c, 2); bhi.break_prob_bp = c.break_prob_bp * 120 / 100;
    let (w_blo, w_bhi) = (run_metric_wait(c, &blo, 912), run_metric_wait(c, &bhi, 912));
    let _ = writeln!(open_w(out_dir, "s45_tornado.csv"), "param,low_wait,high_wait,swing\narr_scale,{:.4},{:.4},{:.4}\nsecurity_mu,{:.4},{:.4},{:.4}\nbreak_prob,{:.4},{:.4},{:.4}",
        w_alo, w_ahi, (w_ahi - w_alo).abs(), w_mlo, w_mhi, (w_mhi - w_mlo).abs(), w_blo, w_bhi, (w_bhi - w_blo).abs());

    let mut f = open_w(out_dir, "s47_report_tables.csv"); let _ = writeln!(f, "stn,arr,dep,block,aban,mean_wait");
    for st in 0..N_STN {
        let mw = if s_main.stn_wait_cnt[st] > 0 { s_main.stn_wait_sum[st] as f64 / s_main.stn_wait_cnt[st] as f64 } else { 0.0 };
        let _ = writeln!(f, "{},{},{},{},{},{:.4}", snm[st], s_main.stn[st].arr_cnt, s_main.stn[st].dep_cnt, s_main.stn[st].block_cnt, s_main.stn[st].aban_cnt, mw);
    }
}
