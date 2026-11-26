use std::fs::File;
use std::io::{BufWriter, Write};
use crate::config::*;
use crate::dtmc::*;
use crate::engine::sim_run;
use crate::network::Sim;
use crate::rng::*;

fn open_w(dir: &str, nm: &str) -> BufWriter<File> { BufWriter::new(File::create(format!("{}/{}", dir, nm)).unwrap()) }

fn run_single_sim(rng: &mut Xoshiro256, c: usize, k_cap: usize, lam_bp: i64, mu_bp: i64, slots: i64) -> (f64, f64, f64) {
    let mut q = 0usize; let (mut arr, mut blk, mut sum_l, mut sum_b) = (0i64, 0i64, 0i64, 0i64);
    for t in 0..slots {
        let a = rng.bernoulli(lam_bp) as usize;
        let m = q.min(c); let mut d = 0usize;
        for _ in 0..m { if rng.bernoulli(mu_bp) != 0 { d += 1; } }
        if q < k_cap { q += a; if a != 0 { arr += 1; } } else if a != 0 { blk += 1; }
        q -= d;
        if t >= slots / 5 { let mb = q.min(c); sum_l += q as i64; sum_b += mb as i64; }
    }
    let ns = (slots - slots / 5) as f64;
    let pb = blk as f64 / if arr + blk > 0 { (arr + blk) as f64 } else { 1.0 };
    (sum_l as f64 / ns, sum_b as f64 / (ns * c as f64), pb)
}

fn run_tandem_sim(rng: &mut Xoshiro256, k1: usize, k2: usize, lam_bp: i64, mu1_bp: i64, mu2_bp: i64, slots: i64) -> (f64, f64, f64, f64) {
    let (mut q1, mut q2) = (0usize, 0usize); let (mut sl1, mut sl2, mut su1, mut su2) = (0i64, 0i64, 0i64, 0i64);
    for t in 0..slots {
        let a = rng.bernoulli(lam_bp) as usize;
        let d1 = if q1 > 0 { rng.bernoulli(mu1_bp) as usize } else { 0 };
        let d2 = if q2 > 0 { rng.bernoulli(mu2_bp) as usize } else { 0 };
        if q1 < k1 { q1 += a; }
        q1 -= d1; q2 -= d2;
        if d1 != 0 && q2 < k2 { q2 += 1; }
        if t >= slots / 5 { sl1 += q1 as i64; sl2 += q2 as i64; if q1 > 0 { su1 += 1; } if q2 > 0 { su2 += 1; } }
    }
    let ns = (slots - slots / 5) as f64;
    (sl1 as f64 / ns, sl2 as f64 / ns, su1 as f64 / ns, su2 as f64 / ns)
}

pub fn write_s01_to_s13(c: &Config, out_dir: &str) {
    let mut hrng = 0xcbf29ce484222325u64;
    let pairs = [(0u32, 0u32), (1, 2), (2, 1)];
    let exp_hv = [0xbe40b757fe0b6697u64, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb];
    let mut rng_ok = true;
    for i in 0..3 {
        let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(pairs[i].0, pairs[i].1);
        let mut st = Xoshiro256::seed(splitmix64(&mut x));
        let mut buf = [0u8; 32];
        for j in 0..4 { buf[j * 8..(j + 1) * 8].copy_from_slice(&st.next_u64().to_le_bytes()); }
        let hv = fnv1a_bytes(&buf);
        if hv != exp_hv[i] { rng_ok = false; }
        hrng = fnv1a_i64(hrng, hv as i64);
    }
    let mut s0 = Sim::new(c, &default_rp(c, 1)); sim_run(&mut s0);
    let inv_ok = s0.inv_violations == 0;
    let _ = writeln!(open_w(out_dir, "s01_selftest.txt"), "ok,h_rng,inv_ok,tot_arr,tot_done\n{},{:016x},{},{},{}", (rng_ok && inv_ok) as i32, hrng, inv_ok as i32, s0.total_arr, s0.total_done);
    let _ = writeln!(open_w(out_dir, "s02_config.csv"), "key,value\nmaster_seed,{}\nslots_per_day,{}\nparity_days,{}\narr_scale_bp,{}\nbelt_cap,{}", c.master_seed, c.slots_per_day, c.parity_days, c.arr_scale_bp, c.belt_cap);

    let mut f = open_w(out_dir, "s03_flights.csv"); let _ = writeln!(f, "flight_id,dep_slot,gate_close_slot,seats");
    for h in 0..24 { let _ = writeln!(f, "{},{},{},180", h, h * 360 + 300, h * 360 + 300 - c.gate_close_lead as i32); }
    let mut f = open_w(out_dir, "s03_arrival_profile.csv"); let _ = writeln!(f, "hour,arr_bp");
    for h in 0..24 { let _ = writeln!(f, "{},{}", h, s0.arr_prof[h]); }

    let snm = ["CHECKIN","SECURITY","IMIGRASI","GATE","BAGGAGE_BELT","CUSTOMS_ARRIVAL"];
    let mut f = open_w(out_dir, "s04_stations.csv"); let _ = writeln!(f, "stn_id,name,c_max,mu_bp,cap,has_prio");
    for i in 0..N_STN { let _ = writeln!(f, "{},{},{},{},{},{}", i, snm[i], s0.stn[i].c_max, c.mu_bp[i], s0.stn[i].cap, s0.stn[i].has_prio); }
    let mut f = open_w(out_dir, "s04_roster.csv"); let _ = writeln!(f, "hour,checkin,security,imigrasi,gate,belt,customs");
    for h in 0..24 { let _ = writeln!(f, "{},{},{},{},{},{},{}", h, s0.roster[h][0], s0.roster[h][1], s0.roster[h][2], s0.roster[h][3], s0.roster[h][4], s0.roster[h][5]); }

    let mut f = open_w(out_dir, "s05_queue_script.csv"); let _ = writeln!(f, "case_id,mode,q_init,cap,prio,breakdown,gate_closed,accepted,served,lost");
    for k in 0..12 {
        let (q0, cap, brk, gc) = ((k % 4) + 1, if k >= 6 { 3 } else { 10 }, (k % 3 == 2) as i32, (k == 11) as i32);
        let acc = q0.min(cap); let srv = if gc != 0 { 0 } else if brk != 0 { acc - 1 } else { acc }; let lost = if gc != 0 { acc } else { q0 - acc };
        let _ = writeln!(f, "{},{},{},{},{},{},{},{},{},{}", k, k % 3, q0, cap, k % 2, brk, gc, acc, srv, lost);
    }
    let mut f = open_w(out_dir, "s06_bag_script.csv"); let _ = writeln!(f, "case_id,n_bags,belt_len,belt_cap,miss_bp,accepted,blocked,missed");
    for k in 0..12 {
        let (nb, bl, cap) = (k % 3, (k / 3) * 4, 10);
        let acc = if bl + nb <= cap { nb } else { (cap - bl).max(0) }; let blk = nb - acc; let ms = (k % 4 == 3 && acc > 0) as i32;
        let _ = writeln!(f, "{},{},{},{},{},{},{},{}", k, nb, bl, cap, c.bag_miss_bp, acc, blk, ms);
    }

    let mut f = open_w(out_dir, "s07_dtmc_single.csv"); let _ = writeln!(f, "c,K,lam,mu,p_block,L,L_q,U,W,W_q");
    let (lam, mu) = (0.08, 0.12);
    for cc in [1usize, 2, 3] { for kk in [5usize, 10, 60] {
        let dr = solve_dtmc_single(cc, kk, lam, mu);
        let _ = writeln!(f, "{},{},{:.4},{:.4},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12}", cc, kk, lam, mu, dr.p_block, dr.l, dr.l_q, dr.u, dr.w, dr.w_q);
    }}
    let (t1, t2, pbt) = solve_dtmc_tandem(6, 6, 0.08, 0.14, 0.12);
    let _ = writeln!(open_w(out_dir, "s08_dtmc_tandem.csv"), "K1,K2,lam,mu1,mu2,p_blk1,L1,U1,p_blk2,L2,U2,p_blk_tnd\n6,6,0.0800,0.1400,0.1200,{:.12},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12}",
        t1.p_block, t1.l, t1.u, t2.p_block, t2.l, t2.u, pbt);
    let (cf, d60) = (solve_geo_closed_form(lam, mu), solve_dtmc_single(1, 60, lam, mu));
    let _ = writeln!(open_w(out_dir, "s09_geo_closed_form.csv"), "lam,mu,L_closed,W_closed,U_closed,L_dtmc60,rel_err_L\n{:.4},{:.4},{:.12},{:.12},{:.12},{:.12},{:.12e}",
        lam, mu, cf.l, cf.w, cf.u, d60.l, (cf.l - d60.l).abs() / cf.l);

    let mut f = open_w(out_dir, "s10_single_validation.csv"); let _ = writeln!(f, "c,K,metric,dtmc,sim_mean,sim_se,z_score,within_3se");
    let val_slots = 20000i64;
    for (ic, &cc) in [1usize, 2, 3].iter().enumerate() {
        let kk = if cc == 1 { 10 } else { 5 };
        let dr = solve_dtmc_single(cc, kk, lam, mu);
        let (mut sl, mut sl2, mut su, mut su2, mut sp, mut sp2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for r in 0..50u32 {
            let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(500 + ic as u32, r);
            let mut rng = Xoshiro256::seed(splitmix64(&mut x));
            let (vl, vu, vp) = run_single_sim(&mut rng, cc, kk, 800, 1200, val_slots);
            sl += vl; sl2 += vl * vl; su += vu; su2 += vu * vu; sp += vp; sp2 += vp * vp;
        }
        for (nm, dt, s1, s2) in [("L", dr.l, sl, sl2), ("U", dr.u, su, su2), ("p_block", dr.p_block, sp, sp2)] {
            let mn = s1 / 50.0; let se = ((s2 / 50.0 - mn * mn) / 49.0).max(1e-18).sqrt(); let z = (mn - dt) / se;
            let _ = writeln!(f, "{},{},{},{:.6},{:.6},{:.6},{:.3},{}", cc, kk, nm, dt, mn, se, z, (z.abs() <= 3.0) as i32);
        }
    }

    let mut f = open_w(out_dir, "s11_tandem_validation.csv"); let _ = writeln!(f, "metric,dtmc,sim_mean,sim_se,z_score,within_3se");
    let (mut sl1, mut sl12, mut sl2, mut sl22, mut su1, mut su12, mut su2, mut su22) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for r in 0..50u32 {
        let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(600, r);
        let mut rng = Xoshiro256::seed(splitmix64(&mut x));
        let (l1, l2, u1, u2) = run_tandem_sim(&mut rng, 6, 6, 800, 1400, 1200, val_slots);
        sl1 += l1; sl12 += l1*l1; sl2 += l2; sl22 += l2*l2; su1 += u1; su12 += u1*u1; su2 += u2; su22 += u2*u2;
    }
    for (nm, dt, s1, s2) in [("L1", t1.l, sl1, sl12), ("L2", t2.l, sl2, sl22), ("U1", t1.u, su1, su12), ("U2", t2.u, su2, su22)] {
        let mn = s1 / 50.0; let se = ((s2 / 50.0 - mn * mn) / 49.0).max(1e-18).sqrt(); let z = (mn - dt) / se;
        let _ = writeln!(f, "{},{:.6},{:.6},{:.6},{:.3},{}", nm, dt, mn, se, z, (z.abs() <= 3.0) as i32);
    }

    let mut rpw = default_rp(c, 2); rpw.audit_on = 1;
    let mut sw = Sim::new(c, &rpw); sim_run(&mut sw);
    let n_blk = sw.audit_slots as usize / 5;
    let blk_y: Vec<f64> = (0..n_blk).map(|b| (0..5).map(|k| sw.q_series[(b * 5 + k) * N_STN + STN_SECURITY] as f64).sum::<f64>() / 5.0).collect();
    let (mut best_d, mut best_mser) = (0usize, 1e30f64);
    for d in (0..n_blk / 2).step_by(10) {
        let m = (n_blk - d) as f64; let mn = blk_y[d..].iter().sum::<f64>() / m;
        let s2: f64 = blk_y[d..].iter().map(|&y| (y - mn) * (y - mn)).sum();
        let mser = s2 / (m * m);
        if mser < best_mser { best_mser = mser; best_d = d * 5; }
    }
    let _ = writeln!(open_w(out_dir, "s12_warmup.csv"), "warmup_slots,mser_min\n{},{:.6}", best_d, best_mser);

    let mut f = open_w(out_dir, "s13_calibration.csv"); let _ = writeln!(f, "arr_scale_bp,max_util,loss_frac,calibrated");
    for scale in [15000, 18000, 22000, 26000, 30000] {
        let mut rpc = default_rp(c, 2); rpc.arr_scale_bp = scale;
        let mut sc = Sim::new(c, &rpc); sim_run(&mut sc);
        let max_u = (0..4).map(|st| sc.stn[st].busy_slots as f64 / (2 * SLOTS_PER_DAY * sc.stn[st].c_max as i64) as f64).fold(0.0f64, f64::max);
        let lf = if sc.total_arr > 0 { sc.total_lost as f64 / sc.total_arr as f64 } else { 0.0 };
        let _ = writeln!(f, "{},{:.4},{:.4},{}", scale, max_u, lf, ((0.55..=0.92).contains(&max_u) && lf < 0.15) as i32);
    }
}
