use std::fs::File;
use std::io::{BufWriter, Write};
use crate::engine::*;
use crate::population::*;
use crate::rng::*;

fn open_bw(dir: &str, name: &str) -> BufWriter<File> {
    BufWriter::new(File::create(format!("{}/{}", dir, name)).unwrap())
}

fn make_rp(n: i64, days: i64, beta_bp: i64, interv_on: i32) -> RunParams {
    RunParams { n, days, beta_bp, interv_on, ..Default::default() }
}

fn net_degrees(s: &Sim) -> (f64, f64, f64) {
    let (mut shh, mut ssw, mut scom) = (0.0, 0.0, 0.0);
    for p in &s.pop {
        let h = p.hh_id as usize;
        shh += (s.hh_off[h + 1] - s.hh_off[h] - 1) as f64;
        if p.age_grp == 0 { ssw += 6.0; } else if p.age_grp == 1 { ssw += 5.0; }
        scom += if p.super_flag != 0 { 12.0 } else { 3.0 };
    }
    let nf = s.n as f64;
    (shh / nf, ssw / nf, scom / nf)
}

pub fn write_s01_to_s09(c: &Config, out_dir: &str) {
    let mut hrng = 0xcbf29ce484222325u64;
    let pairs = [(0u32, 0u32), (1, 2), (2, 1)];
    let exp_hv = [0xbe40b757fe0b6697u64, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb];
    let mut rng_ok = true;
    for i in 0..3 {
        let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(pairs[i].0, pairs[i].1);
        let sd = splitmix64(&mut x);
        let mut st = Xoshiro256::seed(sd);
        let mut buf = [0u8; 32];
        for j in 0..4 { buf[j * 8..(j + 1) * 8].copy_from_slice(&st.next_u64().to_le_bytes()); }
        let hv = fnv1a_bytes(&buf);
        if hv != exp_hv[i] { rng_ok = false; }
        hrng = fnv1a_i64(hrng, hv as i64);
    }
    let mut st0 = Sim::new(c, &make_rp(200, 10, 450, 0), false); sim_run(&mut st0);
    let inv_ok = st0.inv_violations == 0;
    let _ = writeln!(open_bw(out_dir, "s01_selftest.txt"), "ok,h_rng,inv_ok,tot_trans\n{},{:016x},{},{}", (rng_ok && inv_ok) as i32, hrng, inv_ok as i32, st0.total_trans);
    let _ = writeln!(open_bw(out_dir, "s02_config.csv"), "key,value\nmaster_seed,{}\nparity_N,{}\nparity_days,{}\nmain_days,{}\nbeta_bp_init,{}", c.master_seed, c.parity_n, c.parity_days, c.main_days, c.beta_bp_init);

    let sp = Sim::new(c, &make_rp(1000, 1, 450, 0), false);
    let _ = writeln!(open_bw(out_dir, "s03_population.csv"), "age_grp,count\n0,{}\n1,{}\n2,{}", sp.age_tot[0], sp.age_tot[1], sp.age_tot[2]);
    let mut f = open_bw(out_dir, "s03_households.csv"); let _ = writeln!(f, "hh_id,size");
    for h in 0..sp.n_hh.min(50) { let _ = writeln!(f, "{},{}", h, sp.hh_off[h + 1] - sp.hh_off[h]); }
    let mut f = open_bw(out_dir, "s03_schools_work.csv"); let _ = writeln!(f, "sw_id,size");
    for w in 0..sp.n_sw.min(50) { let _ = writeln!(f, "{},{}", w, sp.sw_off[w + 1] - sp.sw_off[w]); }
    let (chh, csw, ccom) = net_degrees(&sp);
    let _ = writeln!(open_bw(out_dir, "s04_network.csv"), "n_hh,n_sw,deg_hh,deg_sw,deg_comm\n{},{},{:.4},{:.4},{:.4}", sp.n_hh, sp.n_sw, chh, csw, ccom);

    let mut f = open_bw(out_dir, "s05_contact_script.csv"); let _ = writeln!(f, "case_id,vax,mask,quar,eff_bp");
    for k in 0..12 {
        let (vax, mask, quar) = (k % 3, (k / 3) % 2, (k / 6) % 2);
        let mut p = if quar != 0 { 0 } else { 450 };
        if mask != 0 { p = (p * (BP - c.mask_eff_bp)) / BP; }
        if vax == 1 { p = (p * (BP - c.vax_eff1_bp)) / BP; } else if vax == 2 { p = (p * (BP - c.vax_eff2_bp)) / BP; }
        let _ = writeln!(f, "{},{},{},{},{}", k, vax, mask, quar, p);
    }
    let mut f = open_bw(out_dir, "s06_progression_script.csv"); let _ = writeln!(f, "case_id,age_grp,hosp_rej,p_hosp_bp,p_die_bp");
    for k in 0..12 {
        let (ag, rej) = (k % 3, (k / 3) % 2);
        let ph = if ag == 0 { 120 } else if ag == 1 { 450 } else { 1800 };
        let mut pd = if ag == 0 { 100 } else if ag == 1 { 500 } else { 2000 };
        if rej != 0 { pd *= 2; }
        let _ = writeln!(f, "{},{},{},{},{}", k, ag, rej, ph, pd);
    }
    let _ = writeln!(open_bw(out_dir, "s07_interventions.csv"), "name,start_day,param_bp\nschool_close,{},10000\nlockdown,{},5000\nmask,{},{}\ntest,{},{}\nvax,{},{}\nvariant,{},{}",
        c.school_close_day, c.lockdown_day, c.mask_day, c.mask_eff_bp, c.test_day, c.test_prob_bp, c.vax_day, c.vax_daily_bp, c.variant_day, c.variant_mult_bp);

    let mut f = open_bw(out_dir, "s08_calibration.csv"); let _ = writeln!(f, "beta_bp,r0_analytic,attack_rate,calibrated");
    for bv in [300, 380, 450, 520, 600] {
        let mut sc = Sim::new(c, &make_rp(1000, 120, bv, 0), false); sim_run(&mut sc);
        let (ch, cs2, cc2) = net_degrees(&sc);
        let r0 = analytical_r0(c, bv, ch, cs2, cc2);
        let ar = sc.total_trans as f64 / sc.n as f64;
        let _ = writeln!(f, "{},{:.4},{:.4},{}", bv, r0, ar, ((1.5..=3.5).contains(&r0) && (0.20..=0.80).contains(&ar)) as i32);
    }
    let mut s9 = Sim::new(c, &make_rp(500, 30, 450, 0), false); sim_run(&mut s9);
    let _ = writeln!(open_bw(out_dir, "s09_baseline_small.csv"), "N,days,total_trans,deaths,inv_ok\n500,30,{},{},{}", s9.total_trans, s9.daily_cnt[29][ST_D as usize], (s9.inv_violations == 0) as i32);
}

fn write_daily_csv(s: &Sim, out_dir: &str, fname: &str) {
    let mut f = open_bw(out_dir, fname); let _ = writeln!(f, "day,S,E,Ia,Is,H,ICU,R,D");
    for d in 0..(s.rp.days as usize) {
        let _ = write!(f, "{}", d);
        for k in 0..8 { let _ = write!(f, ",{}", s.daily_cnt[d][k]); }
        let _ = writeln!(f);
    }
}

fn write_trans_csv(s: &Sim, out_dir: &str, fname: &str) {
    let mut f = open_bw(out_dir, fname); let _ = writeln!(f, "day,src_id,dst_id,layer,variant,gen_int");
    for r in &s.trans { let _ = writeln!(f, "{},{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4, r.5); }
}

pub fn write_s10_to_s37(s_no: &Sim, s_main: &Sim, c: &Config, out_dir: &str, beta_cal: i64) {
    let d_cnt = s_main.rp.days as usize;
    write_daily_csv(s_no, out_dir, "daily_counts_nointerv.csv");
    write_trans_csv(s_no, out_dir, "transmissions_log_nointerv.csv");
    write_daily_csv(s_main, out_dir, "daily_counts.csv");
    write_trans_csv(s_main, out_dir, "transmissions_log.csv");

    let mut f = open_bw(out_dir, "incidence_by_age.csv"); let _ = writeln!(f, "day,age_0_17,age_18_59,age_60p");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{},{}", d, s_main.daily_inc_age[d][0], s_main.daily_inc_age[d][1], s_main.daily_inc_age[d][2]); }
    let mut f = open_bw(out_dir, "hospital_daily.csv"); let _ = writeln!(f, "day,hosp_occ,icu_occ,hosp_rej_cum,icu_rej_cum");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{},{},{}", d, s_main.daily_hosp[d][0], s_main.daily_hosp[d][1], s_main.daily_hosp[d][2], s_main.daily_hosp[d][3]); }
    let mut f = open_bw(out_dir, "deaths_by_age.csv"); let _ = writeln!(f, "day,d_0_17,d_18_59,d_60p");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{},{}", d, s_main.daily_deaths_age[d][0], s_main.daily_deaths_age[d][1], s_main.daily_deaths_age[d][2]); }
    let mut f = open_bw(out_dir, "tests.csv"); let _ = writeln!(f, "day,tests,positives");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{}", d, s_main.daily_tests[d][0], s_main.daily_tests[d][1]); }
    let mut f = open_bw(out_dir, "quarantine.csv"); let _ = writeln!(f, "day,quarantined");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{}", d, s_main.daily_quar[d]); }
    let mut f = open_bw(out_dir, "vaccination.csv"); let _ = writeln!(f, "day,dose1,dose2");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{}", d, s_main.daily_vax[d][0], s_main.daily_vax[d][1]); }
    let mut f = open_bw(out_dir, "variant_share.csv"); let _ = writeln!(f, "day,v0_inf,v1_inf");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{}", d, s_main.daily_var[d][0], s_main.daily_var[d][1]); }
    let mut f = open_bw(out_dir, "generation_interval.csv"); let _ = writeln!(f, "idx,gen_int");
    for (i, r) in s_main.trans.iter().take(2000).enumerate() { let _ = writeln!(f, "{},{}", i, r.5); }
    let mut f = open_bw(out_dir, "agents_final.csv"); let _ = writeln!(f, "id,age_grp,status,vax_doses,offspring,ever_inf");
    for (i, p) in s_main.pop.iter().take(2000).enumerate() { let _ = writeln!(f, "{},{},{},{},{},{}", i, p.age_grp, p.status, p.vax_doses, p.offspring, p.ever_inf); }
    let mut f = open_bw(out_dir, "progression_log.csv"); let _ = writeln!(f, "day,pid,from_st,to_st");
    for r in &s_main.progs { let _ = writeln!(f, "{},{},{},{}", r.0, r.1, r.2, r.3); }

    for a in 0..3 {
        let mut f = open_bw(out_dir, &format!("daily_age_{}.csv", a)); let _ = writeln!(f, "day,inc,deaths");
        for d in 0..d_cnt { let _ = writeln!(f, "{},{},{}", d, s_main.daily_inc_age[d][a], s_main.daily_deaths_age[d][a]); }
    }
    for (l, nm) in ["layer_hh.csv", "layer_sw.csv", "layer_comm.csv"].iter().enumerate() { let _ = writeln!(open_bw(out_dir, nm), "layer,infections\n{},{}", l, s_main.layer_inf[l]); }
    let mut f = open_bw(out_dir, "hospital_cap_nointerv.csv"); let _ = writeln!(f, "day,hosp_occ,icu_occ");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{},{}", d, s_no.daily_hosp[d][0], s_no.daily_hosp[d][1]); }
    let mut f = open_bw(out_dir, "deaths_nointerv.csv"); let _ = writeln!(f, "day,deaths");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{}", d, s_no.daily_cnt[d][ST_D as usize]); }
    let mut oh_bins = [0i64; 6];
    for p in &s_no.pop { if p.ever_inf != 0 { oh_bins[(p.offspring as usize).min(5)] += 1; } }
    let mut f = open_bw(out_dir, "offspring_hist.csv"); let _ = writeln!(f, "bin,count");
    for b in 0..6 { let _ = writeln!(f, "{},{}", b, oh_bins[b]); }
    let mut f = open_bw(out_dir, "rt_nointerv.csv"); let _ = writeln!(f, "day,inc");
    for d in 0..d_cnt { let _ = writeln!(f, "{},{}", d, s_no.daily_inc_age[d].iter().sum::<i64>()); }

    let tot_viol = s_no.inv_violations + s_main.inv_violations;
    let _ = writeln!(open_bw(out_dir, "s12_invariants.csv"), "violations,ok\n{},{}", tot_viol, (tot_viol == 0) as i32);

    let mut f = open_bw(out_dir, "s13_rt.csv"); let _ = writeln!(f, "day,rt");
    for d in 0..d_cnt {
        let num: i64 = (0..5).filter(|&k| d >= k).map(|k| s_main.daily_inc_age[d - k].iter().sum::<i64>()).sum();
        let den: i64 = (4..9).filter(|&k| d >= k).map(|k| s_main.daily_inc_age[d - k].iter().sum::<i64>()).sum();
        let _ = writeln!(f, "{},{:.4}", d, if den > 0 { num as f64 / den as f64 } else { 0.0 });
    }

    let (chh, csw, ccom) = net_degrees(s_no);
    let r0_ana = analytical_r0(c, beta_cal, chh, csw, ccom);
    let i_early = (s_no.daily_inc_age[3].iter().sum::<i64>() + s_no.daily_inc_age[4].iter().sum::<i64>() + 1) as f64 / 2.0;
    let i_mid = (10..16.min(d_cnt)).map(|d| s_no.daily_inc_age[d].iter().sum::<i64>() as f64).sum::<f64>() / 6.0;
    let r_grow = if i_mid > i_early { (i_mid / i_early).ln() / 9.0 } else { 0.0 };
    let _ = writeln!(open_bw(out_dir, "s14_growth.csv"), "growth_r,doubling_days\n{:.6},{:.4}", r_grow, if r_grow > 0.0 { 2.0f64.ln() / r_grow } else { 0.0 });
    let r0_emp = 1.0 + r_grow * 4.5;
    let _ = writeln!(open_bw(out_dir, "s15_r0.csv"), "r0_analytic,r0_empirical,ratio\n{:.4},{:.4},{:.4}", r0_ana, r0_emp, if r0_ana > 0.0 { r0_emp / r0_ana } else { 0.0 });
    let z_ana = solve_final_size(r0_ana); let herd = if r0_ana > 1.0 { 1.0 - 1.0 / r0_ana } else { 0.0 };
    let _ = writeln!(open_bw(out_dir, "s16_final_size.csv"), "z_analytic,herd_thresh,attack_nointerv,attack_interv\n{:.4},{:.4},{:.4},{:.4}", z_ana, herd, s_no.total_trans as f64 / s_no.n as f64, s_main.total_trans as f64 / s_main.n as f64);

    let mut f = open_bw(out_dir, "s17_attack.csv");
    let _ = writeln!(f, "category,count,frac\nHH,{},{:.4}\nSW,{},{:.4}\nCOMM,{},{:.4}\nAGE0,{},{:.4}\nAGE1,{},{:.4}\nAGE2,{},{:.4}",
        s_main.layer_inf[0], s_main.layer_inf[0] as f64 / s_main.n as f64, s_main.layer_inf[1], s_main.layer_inf[1] as f64 / s_main.n as f64, s_main.layer_inf[2], s_main.layer_inf[2] as f64 / s_main.n as f64,
        s_main.age_inf[0], s_main.age_inf[0] as f64 / s_main.age_tot[0].max(1) as f64, s_main.age_inf[1], s_main.age_inf[1] as f64 / s_main.age_tot[1].max(1) as f64, s_main.age_inf[2], s_main.age_inf[2] as f64 / s_main.age_tot[2].max(1) as f64);

    let (mut m1, mut m2, mut n_inf) = (0.0, 0.0, 0i64);
    for p in &s_no.pop { if p.ever_inf != 0 { let o = p.offspring as f64; m1 += o; m2 += o * o; n_inf += 1; } }
    if n_inf > 0 { m1 /= n_inf as f64; m2 /= n_inf as f64; }
    let var_o = m2 - m1 * m1;
    let _ = writeln!(open_bw(out_dir, "s18_offspring.csv"), "n_inf,mean_offspring,var_offspring,k_disp\n{},{:.4},{:.4},{:.4}", n_inf, m1, var_o, if var_o > m1 && m1 > 0.0 { (m1 * m1) / (var_o - m1) } else { 10.0 });

    let pk_h_no = (0..d_cnt).map(|d| s_no.daily_hosp[d][0]).max().unwrap_or(0);
    let pk_i_no = (0..d_cnt).map(|d| s_no.daily_hosp[d][1]).max().unwrap_or(0);
    let pk_h_m = (0..d_cnt).map(|d| s_main.daily_hosp[d][0]).max().unwrap_or(0);
    let pk_i_m = (0..d_cnt).map(|d| s_main.daily_hosp[d][1]).max().unwrap_or(0);
    let _ = writeln!(open_bw(out_dir, "s19_hospital.csv"), "hosp_cap,icu_cap,peak_h_no,peak_i_no,rej_no,peak_h_main,peak_i_main,rej_main\n{},{},{},{},{},{},{},{}",
        s_main.hosp_cap, s_main.icu_cap, pk_h_no, pk_i_no, s_no.branches[9] + s_no.branches[11], pk_h_m, pk_i_m, s_main.branches[9] + s_main.branches[11]);

    let mut f = open_bw(out_dir, "s20_ifr.csv"); let _ = writeln!(f, "age_grp,infections,deaths,ifr");
    for a in 0..3 { let dt: i64 = (0..d_cnt).map(|d| s_main.daily_deaths_age[d][a]).sum(); let _ = writeln!(f, "{},{},{},{:.6}", a, s_main.age_inf[a], dt, if s_main.age_inf[a] > 0 { dt as f64 / s_main.age_inf[a] as f64 } else { 0.0 }); }
    let t_tot: i64 = (0..d_cnt).map(|d| s_main.daily_tests[d][0]).sum(); let t_pos: i64 = (0..d_cnt).map(|d| s_main.daily_tests[d][1]).sum();
    let _ = writeln!(open_bw(out_dir, "s21_testing.csv"), "tests,positives,quar_traces\n{},{},{}", t_tot, t_pos, s_main.branches[17]);

    let (mut vn, mut vi) = ([0i64; 3], [0i64; 3]);
    for p in &s_main.pop { let vd = (p.vax_doses as usize).min(2); vn[vd] += 1; if p.ever_inf != 0 { vi[vd] += 1; } }
    let mut f = open_bw(out_dir, "s22_vaccine.csv"); let _ = writeln!(f, "doses,n_people,n_inf,attack_rate");
    for d in 0..3 { let _ = writeln!(f, "{},{},{},{:.6}", d, vn[d], vi[d], if vn[d] > 0 { vi[d] as f64 / vn[d] as f64 } else { 0.0 }); }

    let v0: i64 = (0..d_cnt).map(|d| s_main.daily_var[d][0]).sum(); let v1: i64 = (0..d_cnt).map(|d| s_main.daily_var[d][1]).sum();
    let _ = writeln!(open_bw(out_dir, "s23_variant.csv"), "v0_total,v1_total,v1_share\n{},{},{:.6}", v0, v1, if v0 + v1 > 0 { v1 as f64 / (v0 + v1) as f64 } else { 0.0 });
    let _ = writeln!(open_bw(out_dir, "s24_waning.csv"), "waning_events,reinfections\n{},{}", s_main.branches[14], s_main.total_reinf);

    let mut r1 = make_rp(500, 30, beta_cal, 0); r1.contact_scale_bp = 0; let mut nc1 = Sim::new(c, &r1, false); sim_run(&mut nc1);
    let mut r2 = make_rp(500, 30, 0, 0); r2.contact_scale_bp = -1; let mut nc2 = Sim::new(c, &r2, false); sim_run(&mut nc2);
    let mut r3 = make_rp(500, 30, beta_cal, 1); r3.vax_perfect = 1; let mut nc3 = Sim::new(c, &r3, false); sim_run(&mut nc3);
    let _ = writeln!(open_bw(out_dir, "s25_negative_controls.csv"), "control,infections,ok\nzero_contact,{},{}\nzero_seed,{},{}\nperfect_vax,{},{}",
        nc1.total_trans, (nc1.total_trans == 0) as i32, nc2.total_trans, (nc2.total_trans == 0) as i32, nc3.total_trans, (nc3.total_trans == 0) as i32);

    let (nsw, dsw) = (1000i64, 90i64);
    let mut f = open_bw(out_dir, "s26_sweep_comp_test.csv"); let _ = writeln!(f, "comp_bp,test_bp,rep,attack_rate");
    for cv in [2000, 6000, 9500] { for tv in [1000, 5000, 9000] { for rp in 0..3 { let mut c2 = c.clone(); c2.master_seed += rp + 1; let mut r = make_rp(nsw, dsw, beta_cal, 1); r.compliance_ovr = cv; r.test_ovr = tv; let mut ss = Sim::new(&c2, &r, false); sim_run(&mut ss); let _ = writeln!(f, "{},{},{},{:.6}", cv, tv, rp, ss.total_trans as f64 / nsw as f64); } } }
    let mut f = open_bw(out_dir, "s27_sweep_vax.csv"); let _ = writeln!(f, "vax_day,vax_cov_bp,rep,attack_rate");
    for vd in [10, 45, 80] { for vc in [2000, 5500, 9000] { for rp in 0..3 { let mut c2 = c.clone(); c2.master_seed += rp + 10; let mut r = make_rp(nsw, dsw, beta_cal, 1); r.vax_day_ovr = vd; r.vax_cov_ovr = vc; let mut ss = Sim::new(&c2, &r, false); sim_run(&mut ss); let _ = writeln!(f, "{},{},{},{:.6}", vd, vc, rp, ss.total_trans as f64 / nsw as f64); } } }
    let mut f = open_bw(out_dir, "s28_sweep_school_lock.csv"); let _ = writeln!(f, "school_day,lock_day,rep,attack_rate");
    for sc in [15, 200] { for lk in [15, 45, 200] { for rp in 0..3 { let mut c2 = c.clone(); c2.master_seed += rp + 20; let mut r = make_rp(nsw, dsw, beta_cal, 1); r.school_ovr = sc; r.lock_ovr = lk; let mut ss = Sim::new(&c2, &r, false); sim_run(&mut ss); let _ = writeln!(f, "{},{},{},{:.6}", sc, lk, rp, ss.total_trans as f64 / nsw as f64); } } }
    let mut f = open_bw(out_dir, "s29_sweep_variant_mask.csv"); let _ = writeln!(f, "var_day,mask_day,rep,attack_rate");
    for vr in [15, 45, 200] { for mk in [10, 40, 200] { for rp in 0..3 { let mut c2 = c.clone(); c2.master_seed += rp + 30; let mut r = make_rp(nsw, dsw, beta_cal, 1); r.var_day_ovr = vr; r.mask_ovr = mk; let mut ss = Sim::new(&c2, &r, false); sim_run(&mut ss); let _ = writeln!(f, "{},{},{},{:.6}", vr, mk, rp, ss.total_trans as f64 / nsw as f64); } } }

    let (mut ens_s, mut ens_s2) = ([0.0f64; 90], [0.0f64; 90]);
    for rp in 0..20 { let mut c2 = c.clone(); c2.master_seed += rp + 100; let mut ss = Sim::new(&c2, &make_rp(nsw, dsw, beta_cal, 1), false); sim_run(&mut ss); for d in 0..90 { let inc = ss.daily_inc_age[d].iter().sum::<i64>() as f64; ens_s[d] += inc; ens_s2[d] += inc * inc; } }
    let mut f = open_bw(out_dir, "s30_ensemble.csv"); let _ = writeln!(f, "day,mean_inc,se_inc");
    for d in 0..90 { let mn = ens_s[d] / 20.0; let _ = writeln!(f, "{},{:.4},{:.4}", d, mn, ((ens_s2[d] / 20.0 - mn * mn) / 19.0).max(0.0).sqrt()); }

    let bn = ["INF_HH","INF_SW","INF_COMM","FAIL_VAX","FAIL_MASK","FAIL_QUAR","PROG_IA","PROG_IS","HOSP_ADM","HOSP_REJ","ICU_ADM","ICU_REJ","DEATH","RECOVER","WANING","TEST_POS","TEST_NEG","QUAR_TRACE","VAX_DOSE","VAX_OUT","VAR_IMPORT"];
    let mut f = open_bw(out_dir, "s31_branches.csv"); let _ = writeln!(f, "branch,count");
    for (i, nm) in bn.iter().enumerate() { let _ = writeln!(f, "{},{}", nm, s_main.branches[i]); }
    let _ = writeln!(open_bw(out_dir, "s32_vacuous_check.csv"), "check,ok\nneg_controls,{}\nensemble_var,{}", (nc1.total_trans==0 && nc2.total_trans==0 && nc3.total_trans==0) as i32, (ens_s.iter().filter(|&&x| x > 0.0).count() > 5) as i32);

    let mut sp1 = Sim::new(c, &make_rp(c.parity_n, c.parity_days, beta_cal, 1), false); sim_run(&mut sp1);
    let mut sp2 = Sim::new(c, &make_rp(c.parity_n, c.parity_days, beta_cal, 1), false); sim_run(&mut sp2);
    let (mut h1, mut h2) = (0xcbf29ce484222325u64, 0xcbf29ce484222325u64);
    for i in 0..(sp1.n as usize) { h1 = fnv1a_i64(h1, sp1.pop[i].status as i64); h2 = fnv1a_i64(h2, sp2.pop[i].status as i64); h1 = fnv1a_i64(h1, sp1.pop[i].offspring as i64); h2 = fnv1a_i64(h2, sp2.pop[i].offspring as i64); }
    let _ = writeln!(open_bw(out_dir, "s33_replay.csv"), "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n{:016x},{:016x},{:016x},{:016x},{}", sp1.event_hash, sp2.event_hash, h1, h2, (sp1.event_hash == sp2.event_hash && h1 == h2) as i32);
    let _ = writeln!(open_bw(out_dir, "s34_parity.csv"), "parity_N,parity_days,ev_hash,fin_hash,total_trans,deaths,inv_violations\n{},{},{:016x},{:016x},{},{},{}", c.parity_n, c.parity_days, sp1.event_hash, h1, sp1.total_trans, sp1.daily_cnt[c.parity_days as usize - 1][ST_D as usize], sp1.inv_violations);
    let _ = writeln!(open_bw(out_dir, "s36_sizes.csv"), "N,days,n_trans,n_progs\n{},{},{},{}", s_main.n, s_main.rp.days, s_main.trans.len(), s_main.progs.len());

    let mut cnt_r = [s_main.n, 0, 0, 0, 0, 0, 0, 0]; let mut pi = 0usize; let mut max_diff = 0i64;
    for d in 0..d_cnt {
        while pi < s_main.progs.len() && s_main.progs[pi].0 == d as i32 { cnt_r[s_main.progs[pi].2 as usize] -= 1; cnt_r[s_main.progs[pi].3 as usize] += 1; pi += 1; }
        for k in 0..8 { max_diff = max_diff.max((cnt_r[k] - s_main.daily_cnt[d][k]).abs()); }
    }
    let _ = writeln!(open_bw(out_dir, "s37_replay_counts.csv"), "max_diff,match\n{},{}", max_diff, (max_diff == 0) as i32);
}
