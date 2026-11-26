use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::time::Instant;

const MASTER_SEED: u64 = 20261126;
const W_WALKS: usize = 200_000;
const W_RUIN: usize = 20_000;
const BENCH_STEPS: usize = 10_000_000;
const T_EXACT_MAX: usize = 5000;
const TS_MSD: [usize; 13] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];
const TS_RET: [usize; 4] = [10, 100, 1000, 10000];

#[inline]
fn sm64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

struct Rng { s: [u64; 4] }
impl Rng {
    fn new(seed: u64) -> Self {
        let mut st = seed;
        let mut s = [0u64; 4];
        for v in &mut s { st = st.wrapping_add(0x9E3779B97F4A7C15); *v = sm64(st); }
        Self { s }
    }
    #[inline]
    fn u64(&mut self) -> u64 {
        let res = (self.s[0].wrapping_add(self.s[3])).rotate_left(23).wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0]; self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2]; self.s[0] ^= self.s[3];
        self.s[2] ^= t; self.s[3] = self.s[3].rotate_left(45);
        res
    }
    #[inline]
    fn u01(&mut self) -> f64 { ((self.u64() >> 11) as f64) * 2.0f64.powi(-53) }
    #[inline]
    fn rint(&mut self, n: usize) -> usize {
        let k = (self.u01() * (n as f64)).floor() as usize;
        if k >= n { n - 1 } else { k }
    }
}

fn task_seed(d: usize, st: usize, idx: usize) -> u64 {
    let mut h = 0x517CC1B727220A95u64;
    h = sm64(h ^ (d as u64));
    h = sm64(h ^ (st as u64));
    h = sm64(h ^ (idx as u64));
    sm64(MASTER_SEED ^ h)
}

fn fnv1a_i32(v: &[i32]) -> u64 {
    let mut h = 0xCBF29CE484222325u64;
    for &val in v {
        for b in val.to_le_bytes() { h ^= b as u64; h = h.wrapping_mul(0x100000001B3); }
    }
    h
}

struct State {
    msd: [[f64; 13]; 3],
    msd_se: [[f64; 13]; 3],
    ret_sim: [[f64; 4]; 3],
    ret_se: [[f64; 4]; 3],
    ret_ex: [[f64; 4]; 3],
    ruin_p: [f64; 6],
    ruin_pe: [f64; 6],
    s03_1024: Vec<Vec<f64>>,
    s04_10k: Vec<Vec<u8>>,
    stage_sec: [f64; 12],
}
impl State {
    fn new() -> Self {
        Self {
            msd: [[0.0; 13]; 3], msd_se: [[0.0; 13]; 3],
            ret_sim: [[0.0; 4]; 3], ret_se: [[0.0; 4]; 3], ret_ex: [[0.0; 4]; 3],
            ruin_p: [0.0; 6], ruin_pe: [0.0; 6],
            s03_1024: vec![vec![0.0; W_WALKS]; 3],
            s04_10k: vec![vec![0; W_WALKS]; 3],
            stage_sec: [0.0; 12],
        }
    }
}

fn s01(_st: &mut State) {
    let mut r = Rng::new(MASTER_SEED);
    let (x0, x1, x2) = (r.u64(), r.u64(), r.u64());
    let ok = x0 == 0xA9920E0E7DAB3EF0 && x1 == 0x033B9CEB4B3BA7A0 && x2 == 0xE7473A6A39AFCEA6;
    let mut w = [1.0, 1.0 / 3.0, 0.0, 0.0, 0.0];
    for n in 1..4 {
        let nf = n as f64;
        w[n + 1] = ((10.0 * nf * nf + 10.0 * nf + 3.0) * w[n] - nf * nf * w[n - 1]) / (9.0 * (nf + 1.0) * (nf + 1.0));
    }
    let enum_ok = (w[2] - 15.0 / 81.0).abs() < 1e-14 && (w[3] - 93.0 / 729.0).abs() < 1e-14 && (w[4] - 639.0 / 6561.0).abs() < 1e-14;
    let mut f = File::create("out/rust/s01_selftest.txt").unwrap();
    writeln!(f, "ok,x0,x1,x2,enum_ok\n{},{:016x},{:016x},{:016x},{}",
             if ok && enum_ok { 1 } else { 0 }, x0, x1, x2, if enum_ok { 1 } else { 0 }).unwrap();
    if !(ok && enum_ok) { std::process::exit(1); }
}

fn s02(_st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s02_parity.csv").unwrap());
    writeln!(f, "d,walk,x0,x1,x2,fnv1a_hex").unwrap();
    for d in 1..=3 {
        let mut coords = Vec::with_capacity(8 * d);
        let mut pos_all = [[0i32; 3]; 8];
        for w in 0..8 {
            let mut r = Rng::new(task_seed(d, 2, w));
            let mut p = [0i32; 3];
            for _ in 0..200 { let k = r.rint(2 * d); p[k >> 1] += if (k & 1) != 0 { -1 } else { 1 }; }
            pos_all[w] = p;
            for a in 0..d { coords.push(p[a]); }
        }
        let h = fnv1a_i32(&coords);
        for w in 0..8 {
            writeln!(f, "{},{},{},{},{},{:016x}", d, w, pos_all[w][0], pos_all[w][1], pos_all[w][2], h).unwrap();
        }
    }
}

fn s03(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s03_msd.csv").unwrap());
    writeln!(f, "d,t,msd,se,exact").unwrap();
    for d in 1..=3 {
        let mut sum1 = [0.0f64; 13];
        let mut sum2 = [0.0f64; 13];
        for w in 0..W_WALKS {
            let mut r = Rng::new(task_seed(d, 3, w));
            let (mut p, mut ti) = ([0i32; 3], 0usize);
            for t in 1..=4096 {
                let k = r.rint(2 * d); p[k >> 1] += if (k & 1) != 0 { -1 } else { 1 };
                if t == TS_MSD[ti] {
                    let r2 = (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) as f64;
                    sum1[ti] += r2; sum2[ti] += r2 * r2;
                    if ti == 10 { st.s03_1024[d - 1][w] = r2; }
                    ti += 1;
                }
            }
        }
        for ti in 0..13 {
            let m = sum1[ti] / (W_WALKS as f64);
            let se = ((sum2[ti] / (W_WALKS as f64) - m * m).max(0.0) / (W_WALKS as f64)).sqrt();
            st.msd[d - 1][ti] = m; st.msd_se[d - 1][ti] = se;
            writeln!(f, "{},{},{:.8},{:.8},{}", d, TS_MSD[ti], m, se, TS_MSD[ti]).unwrap();
        }
    }
}

fn s04(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s04_return_sim.csv").unwrap());
    writeln!(f, "d,T,p_hat,se").unwrap();
    for d in 1..=3 {
        let mut hits = [0usize; 4];
        for w in 0..W_WALKS {
            let mut r = Rng::new(task_seed(d, 4, w));
            let (mut p, mut ret_t) = ([0i32; 3], 0usize);
            for t in 1..=10000 {
                let k = r.rint(2 * d); p[k >> 1] += if (k & 1) != 0 { -1 } else { 1 };
                if (p[0] | p[1] | p[2]) == 0 { ret_t = t; break; }
            }
            st.s04_10k[d - 1][w] = if ret_t > 0 { 1 } else { 0 };
            if ret_t > 0 { for ti in 0..4 { if ret_t <= TS_RET[ti] { hits[ti] += 1; } } }
        }
        for ti in 0..4 {
            let ph = (hits[ti] as f64) / (W_WALKS as f64);
            let se = ((ph * (1.0 - ph)).max(0.0) / (W_WALKS as f64)).sqrt();
            st.ret_sim[d - 1][ti] = ph; st.ret_se[d - 1][ti] = se;
            writeln!(f, "{},{},{:.8},{:.8}", d, TS_RET[ti], ph, se).unwrap();
        }
    }
}

fn s05(st: &mut State) {
    let mut u1 = vec![1.0f64; T_EXACT_MAX + 1];
    let mut u2 = vec![1.0f64; T_EXACT_MAX + 1];
    let mut w3 = vec![1.0f64; T_EXACT_MAX + 1];
    let mut u3 = vec![1.0f64; T_EXACT_MAX + 1];
    let mut fv = vec![0.0f64; T_EXACT_MAX + 1];
    w3[1] = 1.0 / 3.0;
    for n in 1..T_EXACT_MAX {
        let nf = n as f64;
        w3[n + 1] = ((10.0 * nf * nf + 10.0 * nf + 3.0) * w3[n] - nf * nf * w3[n - 1]) / (9.0 * (nf + 1.0) * (nf + 1.0));
    }
    for n in 1..=T_EXACT_MAX {
        let nf = n as f64;
        u1[n] = u1[n - 1] * (2.0 * nf - 1.0) / (2.0 * nf);
        u2[n] = u1[n] * u1[n]; u3[n] = u1[n] * w3[n];
    }
    let us = [&u1, &u2, &u3];
    let mut f = BufWriter::new(File::create("out/rust/s05_return_exact.csv").unwrap());
    writeln!(f, "d,T,F_exact").unwrap();
    for d in 1..=3 {
        let u = us[d - 1]; let mut acc = 0.0;
        for n in 1..=T_EXACT_MAX {
            let mut s = u[n];
            for k in 1..n { s -= fv[k] * u[n - k]; }
            fv[n] = s; acc += s;
            for ti in 0..4 { if 2 * n == TS_RET[ti] { st.ret_ex[d - 1][ti] = acc; writeln!(f, "{},{},{:.10}", d, TS_RET[ti], acc).unwrap(); } }
        }
    }
}

fn s06(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s06_return_compare.csv").unwrap());
    writeln!(f, "d,T,p_hat,se,F_exact,z_score,p_inf_est").unwrap();
    let p3_inf = st.ret_sim[2][3] + (st.ret_sim[2][3] - st.ret_sim[2][2]) / (10.0f64.sqrt() - 1.0);
    for d in 1..=3 {
        for ti in 0..4 {
            let z = if st.ret_se[d - 1][ti] > 1e-12 { (st.ret_sim[d - 1][ti] - st.ret_ex[d - 1][ti]) / st.ret_se[d - 1][ti] } else { 0.0 };
            writeln!(f, "{},{},{:.8},{:.8},{:.10},{:.4},{:.8}", d, TS_RET[ti], st.ret_sim[d - 1][ti], st.ret_se[d - 1][ti],
                     st.ret_ex[d - 1][ti], z, if d == 3 && ti == 3 { p3_inf } else { 0.0 }).unwrap();
        }
    }
}

fn s07(st: &mut State) {
    let ns = [10, 10, 10, 50, 50, 50];
    let xs = [2, 5, 8, 10, 25, 40];
    let mut f = BufWriter::new(File::create("out/rust/s07_ruin.csv").unwrap());
    writeln!(f, "N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z").unwrap();
    for ci in 0..6 {
        let (n, x0) = (ns[ci], xs[ci]);
        let (mut wins, mut st1, mut st2) = (0usize, 0.0f64, 0.0f64);
        for w in 0..W_RUIN {
            let mut r = Rng::new(task_seed(1, 7, ci * 100000 + w));
            let (mut x, mut steps) = (x0, 0usize);
            while x > 0 && x < n { x += if r.rint(2) == 0 { 1 } else { -1 }; steps += 1; }
            if x == n { wins += 1; }
            st1 += steps as f64; st2 += (steps as f64) * (steps as f64);
        }
        let ph = (wins as f64) / (W_RUIN as f64);
        let pse = ((ph * (1.0 - ph)).max(0.0) / (W_RUIN as f64)).sqrt();
        let pe = (x0 as f64) / (n as f64);
        let th = st1 / (W_RUIN as f64);
        let tse = ((st2 / (W_RUIN as f64) - th * th).max(0.0) / (W_RUIN as f64)).sqrt();
        let te = (x0 * (n - x0)) as f64;
        st.ruin_p[ci] = ph; st.ruin_pe[ci] = pe;
        writeln!(f, "{},{},{:.8},{:.8},{:.8},{:.4},{:.4},{:.4},{:.4},{:.4}",
                 n, x0, ph, pse, pe, (ph - pe) / pse, th, tse, te, (th - te) / tse).unwrap();
    }
}

fn s08(_st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s08_endpoint.csv").unwrap());
    writeln!(f, "d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count").unwrap();
    let mut hist = [0usize; 40];
    let n_ep = 20_000;
    for d in 1..=3 {
        let (mut m1, mut m2, mut m4) = (0.0, 0.0, 0.0);
        for w in 0..n_ep {
            let mut r = Rng::new(task_seed(d, 8, w));
            let mut p = [0i32; 3];
            for _ in 0..1024 { let k = r.rint(2 * d); p[k >> 1] += if (k & 1) != 0 { -1 } else { 1 }; }
            let x = p[0] as f64; m1 += x; m2 += x * x; m4 += x * x * x * x;
            if d == 1 { let b = ((x + 128.0) / 6.4).floor().clamp(0.0, 39.0) as usize; hist[b] += 1; }
        }
        m1 /= n_ep as f64; let v = m2 / (n_ep as f64) - m1 * m1;
        writeln!(f, "{},{:.6},{:.6},{:.6},{:.6},-1,0,0,{}", d, v, 1024.0 / (d as f64), (m4 / (n_ep as f64)) / (v * v), 3.0 + ((d as f64) - 3.0) / 1024.0, n_ep).unwrap();
    }
    for (b, cnt) in hist.iter().enumerate() {
        writeln!(f, "1,0,0,0,0,{},{:.2},{:.2},{}", b, -128.0 + (b as f64) * 6.4, -128.0 + ((b + 1) as f64) * 6.4, cnt).unwrap();
    }
}

fn s09(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s09_ci.csv").unwrap());
    writeln!(f, "metric,d,estimate,boot_se,ci_lo,ci_hi").unwrap();
    let mut boot = [0.0f64; 200];
    for d in 1..=3 {
        let mut r1 = Rng::new(task_seed(d, 9, 0));
        for b in &mut boot {
            let mut s = 0.0; for _ in 0..W_WALKS { s += st.s03_1024[d - 1][r1.rint(W_WALKS)]; }
            *b = s / (W_WALKS as f64);
        }
        boot.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let bm = boot.iter().sum::<f64>() / 200.0;
        let bv = boot.iter().map(|v| (v - bm) * (v - bm)).sum::<f64>() / 199.0;
        writeln!(f, "msd_1024,{},{:.6},{:.6},{:.6},{:.6}", d, st.msd[d - 1][10], bv.sqrt(), boot[5], boot[194]).unwrap();

        let mut r2 = Rng::new(task_seed(d, 9, 1));
        for b in &mut boot {
            let mut s = 0usize; for _ in 0..W_WALKS { s += st.s04_10k[d - 1][r2.rint(W_WALKS)] as usize; }
            *b = (s as f64) / (W_WALKS as f64);
        }
        boot.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let bm = boot.iter().sum::<f64>() / 200.0;
        let bv = boot.iter().map(|v| (v - bm) * (v - bm)).sum::<f64>() / 199.0;
        writeln!(f, "return_10000,{},{:.6},{:.6},{:.6},{:.6}", d, st.ret_sim[d - 1][3], bv.sqrt(), boot[5], boot[194]).unwrap();
    }
}

fn s10(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s10_fit.csv").unwrap());
    writeln!(f, "d,slope,slope_se,intercept").unwrap();
    for d in 1..=3 {
        let mx = TS_MSD.iter().map(|&t| (t as f64).ln()).sum::<f64>() / 13.0;
        let my = st.msd[d - 1].iter().map(|&m| m.ln()).sum::<f64>() / 13.0;
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..13 { let dx = (TS_MSD[i] as f64).ln() - mx; num += dx * (st.msd[d - 1][i].ln() - my); den += dx * dx; }
        let sl = num / den; let ic = my - sl * mx;
        let rss = (0..13).map(|i| { let e = st.msd[d - 1][i].ln() - (sl * (TS_MSD[i] as f64).ln() + ic); e * e }).sum::<f64>();
        writeln!(f, "{},{:.6},{:.6},{:.6}", d, sl, ((rss / 11.0) / den).sqrt(), ic).unwrap();
    }
}

fn s11(_st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s11_bench.csv").unwrap());
    writeln!(f, "d,steps_per_rep,reps,ns_med,ns_mad").unwrap();
    let mut sink: i64 = 0;
    for d in 1..=3 {
        let mut vals = [0.0f64; 7];
        for rep in 0..7 {
            let mut r = Rng::new(task_seed(d, 11, rep));
            let mut p = [0i32; 3];
            let t0 = Instant::now();
            for _ in 0..BENCH_STEPS { let k = r.rint(2 * d); p[k >> 1] += if (k & 1) != 0 { -1 } else { 1 }; }
            vals[rep] = t0.elapsed().as_secs_f64() * 1e9 / (BENCH_STEPS as f64);
            sink += (p[0] + p[1] + p[2]) as i64;
        }
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = vals[3];
        let mut dev = vals.map(|v| (v - med).abs());
        dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        writeln!(f, "{},{},7,{:.4},{:.4}", d, BENCH_STEPS, med, dev[3]).unwrap();
    }
    std::hint::black_box(sink);
}

fn s12(st: &mut State) {
    let cols = ["#1f77b4", "#ff7f0e", "#2ca02c"];
    let mut s1 = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>");
    for d in 0..3 {
        s1.push_str(&format!("<polyline fill=\"none\" stroke=\"{}\" stroke-width=\"1.6\" points=\"", cols[d]));
        for i in 0..13 { s1.push_str(&format!("{:.1},{:.1} ", 60.0 + 500.0 * (TS_MSD[i] as f64).log10() / 4096.0f64.log10(), 310.0 - 280.0 * st.msd[d][i].max(1.0).log10() / 4096.0f64.log10())); }
        s1.push_str("\"/>");
    }
    s1.push_str("<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Rust MSD Log-Log (y=t)</text></svg>\n");
    fs::write("out/rust/figures/msd.svg", s1).unwrap();

    let mut s2 = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for d in 0..3 {
        s2.push_str(&format!("<polyline fill=\"none\" stroke=\"{}\" stroke-width=\"1.6\" points=\"", cols[d]));
        for i in 0..4 { s2.push_str(&format!("{:.1},{:.1} ", 60.0 + 500.0 * (i as f64 / 3.0), 310.0 - 280.0 * st.ret_ex[d][i])); }
        s2.push_str("\"/>");
        for i in 0..4 { s2.push_str(&format!("<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"4\" fill=\"{}\"/>", 60.0 + 500.0 * (i as f64 / 3.0), 310.0 - 280.0 * st.ret_sim[d][i], cols[d])); }
    }
    s2.push_str("<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Rust Return Prob F_d(T) Exact vs Sim</text></svg>\n");
    fs::write("out/rust/figures/return_prob.svg", s2).unwrap();

    let mut s3 = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>");
    for i in 0..6 { s3.push_str(&format!("<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"5\" fill=\"#1f77b4\"/>", 60.0 + 500.0 * st.ruin_pe[i], 310.0 - 280.0 * st.ruin_p[i])); }
    s3.push_str("<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Rust Gambler's Ruin Exact vs Sim</text></svg>\n");
    fs::write("out/rust/figures/ruin.svg", s3).unwrap();

    let tot: f64 = st.stage_sec.iter().sum();
    fs::write("out/rust/s12_summary.json", format!("{{\"lang\":\"rust\",\"W\":{},\"ruin_W\":{},\"bench_steps\":{},\"toolchain\":\"rustc 1.93.1\",\"total_sec\":{:.3}}}\n", W_WALKS, W_RUIN, BENCH_STEPS, tot)).unwrap();
}

fn main() {
    fs::create_dir_all("out/rust/figures").unwrap();
    let cmd = std::env::args().nth(1).unwrap_or_else(|| "all".to_string());
    let mut st = State::new();
    let fns: [(&str, fn(&mut State)); 12] = [
        ("s01", s01), ("s02", s02), ("s03", s03), ("s04", s04), ("s05", s05), ("s06", s06),
        ("s07", s07), ("s08", s08), ("s09", s09), ("s10", s10), ("s11", s11), ("s12", s12),
    ];
    if cmd == "all" {
        for (i, (nm, f)) in fns.iter().enumerate() {
            println!("[rust] tahap {}/12 ({})", i + 1, nm);
            let t0 = Instant::now(); f(&mut st); st.stage_sec[i] = t0.elapsed().as_secs_f64();
        }
        return;
    }
    for (nm, f) in fns { if cmd == nm { f(&mut st); return; } }
    std::process::exit(1);
}
