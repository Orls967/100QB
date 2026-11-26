use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::time::Instant;

const PI_VAL: f64 = 3.14159265358979323846;
const NS_ALL: [usize; 6] = [15, 31, 63, 127, 255, 511];
const MNAMES: [&str; 7] = ["jacobi", "gs", "rbgs", "sor", "cg", "pcg_ssor", "mg"];

fn sin_spec(mut x: f64) -> f64 {
    let mut sgn = 1.0;
    if x < 0.0 { x = -x; sgn = -sgn; }
    while x > PI_VAL { x -= PI_VAL; sgn = -sgn; }
    if x > 0.5 * PI_VAL { x = PI_VAL - x; }
    let x2 = x * x;
    let c = [
        1.0, -1.0/6.0, 1.0/120.0, -1.0/5040.0, 1.0/362880.0, -1.0/39916800.0,
        1.0/6227020800.0, -1.0/1307674368000.0, 1.0/355687428096000.0, -1.0/121645100408832000.0,
    ];
    let mut p = c[9];
    for k in (0..9).rev() { p = p * x2 + c[k]; }
    sgn * (x * p)
}
fn cos_spec(x: f64) -> f64 { sin_spec(0.5 * PI_VAL - x) }
fn omega_opt(n: usize) -> f64 {
    match n {
        15 => 1.673513677715992, 31 => 1.8214651907890225, 63 => 1.906454701582762,
        127 => 1.952093233850055, 255 => 1.9757544535797149, _ => 1.9878030696593354,
    }
}
fn fnv1a_f64(v: &[f64]) -> u64 {
    let mut h = 0xCBF29CE484222325u64;
    for &x in v {
        for b in x.to_bits().to_le_bytes() { h ^= b as u64; h = h.wrapping_mul(0x100000001B3); }
    }
    h
}
fn apply_a(n: usize, u: &[f64], out: &mut [f64]) {
    for i in 0..n {
        let r = i * n;
        for j in 0..n {
            let k = r + j;
            let w = if j > 0 { u[k - 1] } else { 0.0 };
            let e = if j + 1 < n { u[k + 1] } else { 0.0 };
            let s = if i > 0 { u[k - n] } else { 0.0 };
            let nn = if i + 1 < n { u[k + n] } else { 0.0 };
            out[k] = (((4.0 * u[k] - w) - e) - s) - nn;
        }
    }
}
fn dot_seq(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0; for i in 0..a.len() { s += a[i] * b[i]; } s
}
fn norm2_seq(a: &[f64]) -> f64 { dot_seq(a, a).sqrt() }
fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    let mut m = 0.0f64; for i in 0..a.len() { m = m.max((a[i] - b[i]).abs()); } m
}
fn build_prob(n: usize) -> (Vec<f64>, Vec<f64>) {
    let h = 1.0 / ((n + 1) as f64); let h2 = h * h;
    let sx: Vec<f64> = (0..n).map(|i| sin_spec(PI_VAL * ((i + 1) as f64) * h)).collect();
    let sy: Vec<f64> = (0..n).map(|i| sin_spec(2.0 * PI_VAL * ((i + 1) as f64) * h)).collect();
    let mut b = vec![0.0; n * n]; let mut ue = vec![0.0; n * n];
    for i in 0..n {
        let x = ((i + 1) as f64) * h; let xx = x * (1.0 - x);
        for j in 0..n {
            let y = ((j + 1) as f64) * h; let yy = y * (1.0 - y); let sxy = sx[i] * sy[j];
            ue[i * n + j] = sxy + xx * yy;
            b[i * n + j] = h2 * (5.0 * PI_VAL * PI_VAL * sxy + 2.0 * xx + 2.0 * yy);
        }
    }
    (b, ue)
}
fn compute_res(n: usize, u: &[f64], b: &[f64], r: &mut [f64]) -> f64 {
    apply_a(n, u, r); for i in 0..n * n { r[i] = b[i] - r[i]; } norm2_seq(r)
}
fn rbgs_sweep(n: usize, u: &mut [f64], b: &[f64]) {
    for col in 0..2 {
        for i in 0..n {
            let row = i * n;
            for j in 0..n {
                if ((i + j) & 1) != col { continue; }
                let k = row + j;
                let w = if j > 0 { u[k - 1] } else { 0.0 };
                let e = if j + 1 < n { u[k + 1] } else { 0.0 };
                let s = if i > 0 { u[k - n] } else { 0.0 };
                let nn = if i + 1 < n { u[k + n] } else { 0.0 };
                u[k] = 0.25 * ((((b[k] + w) + e) + s) + nn);
            }
        }
    }
}
fn restrict_fw(nf: usize, rf: &[f64]) -> Vec<f64> {
    let nc = (nf - 1) / 2; let mut rc = vec![0.0; nc * nc];
    for ic in 0..nc {
        let fi = 2 * ic + 1;
        for jc in 0..nc {
            let fj = 2 * jc + 1; let k = fi * nf + fj;
            let edge = ((rf[k - 1] + rf[k + 1]) + rf[k - nf]) + rf[k + nf];
            let corn = ((rf[k - nf - 1] + rf[k - nf + 1]) + rf[k + nf - 1]) + rf[k + nf + 1];
            rc[ic * nc + jc] = 4.0 * (0.25 * rf[k] + 0.125 * edge + 0.0625 * corn);
        }
    }
    rc
}
fn prolong_add(nf: usize, uf: &mut [f64], ec: &[f64]) {
    let nc = (nf - 1) / 2;
    for ic in 0..nc {
        let fi = 2 * ic + 1;
        for jc in 0..nc {
            let fj = 2 * jc + 1; let k = fi * nf + fj; let v = ec[ic * nc + jc];
            uf[k] += v;
            uf[k - 1] += 0.5 * v; uf[k + 1] += 0.5 * v; uf[k - nf] += 0.5 * v; uf[k + nf] += 0.5 * v;
            uf[k - nf - 1] += 0.25 * v; uf[k - nf + 1] += 0.25 * v; uf[k + nf - 1] += 0.25 * v; uf[k + nf + 1] += 0.25 * v;
        }
    }
}
fn solve_dense3(b: &[f64], x: &mut [f64]) {
    let mut m = [[0.0f64; 10]; 9];
    for i in 0..3 { for j in 0..3 {
        let r = i * 3 + j; m[r][r] = 4.0;
        if j > 0 { m[r][r - 1] = -1.0; } if j < 2 { m[r][r + 1] = -1.0; }
        if i > 0 { m[r][r - 3] = -1.0; } if i < 2 { m[r][r + 3] = -1.0; }
        m[r][9] = b[r];
    } }
    for c in 0..9 {
        let mut piv = c;
        for r in (c + 1)..9 { if m[r][c].abs() > m[piv][c].abs() { piv = r; } }
        if piv != c { m.swap(c, piv); }
        for r in (c + 1)..9 { let f = m[r][c] / m[c][c]; for k in c..10 { m[r][k] -= f * m[c][k]; } }
    }
    for r in (0..9).rev() {
        let mut s = m[r][9]; for k in (r + 1)..9 { s -= m[r][k] * x[k]; }
        x[r] = s / m[r][r];
    }
}
fn vcycle(n: usize, u: &mut [f64], b: &[f64]) {
    if n == 3 { solve_dense3(b, u); return; }
    rbgs_sweep(n, u, b); rbgs_sweep(n, u, b);
    let mut r = vec![0.0; n * n]; compute_res(n, u, b, &mut r);
    let rc = restrict_fw(n, &r);
    let nc = (n - 1) / 2; let mut ec = vec![0.0; nc * nc];
    vcycle(nc, &mut ec, &rc); prolong_add(n, u, &ec); rbgs_sweep(n, u, b);
}
fn fmg_solve(n: usize, u: &mut [f64], b: &[f64]) {
    if n == 3 { solve_dense3(b, u); return; }
    let bc = restrict_fw(n, b); let nc = (n - 1) / 2; let mut uc = vec![0.0; nc * nc];
    fmg_solve(nc, &mut uc, &bc); u.fill(0.0); prolong_add(n, u, &uc); vcycle(n, u, b);
}

#[derive(Clone, Default)]
struct SolveRes { iters: usize, rel_res: f64, err_inf: f64, rho_emp: f64, hist_hash: u64, sol_hash: u64 }
struct State { res: Vec<Vec<SolveRes>>, hist63: Vec<Vec<f64>>, disc_err: [f64; 6], stage_s: [f64; 19] }
impl State {
    fn new() -> Self { Self { res: vec![vec![SolveRes::default(); 6]; 7], hist63: vec![vec![]; 7], disc_err: [0.0; 6], stage_s: [0.0; 19] } }
}

fn run_solver(meth: usize, n: usize, st: &mut State, save63: bool) -> SolveRes {
    let n2 = n * n;
    let (b, ue) = build_prob(n);
    let (mut u, mut r, mut p, mut ap, mut z) = (vec![0.0; n2], vec![0.0; n2], vec![0.0; n2], vec![0.0; n2], vec![0.0; n2]);
    let nb = norm2_seq(&b); let om = omega_opt(n);
    let max_it = match meth { 0 => 100000, 1 | 2 => 50000, 3 => 10000, 4 => 5000, 5 => 3000, _ => 100 };
    let mut hist = Vec::with_capacity(256);
    let (mut rr, mut rz) = (0.0, 0.0);
    if meth == 4 || meth == 5 { r.copy_from_slice(&b); rr = dot_seq(&r, &r); }
    if meth == 4 { p.copy_from_slice(&r); }
    if meth == 5 {
        for i in 0..n { for j in 0..n { let k = i * n + j; z[k] = 0.25 * ((r[k] + if j > 0 { z[k - 1] } else { 0.0 }) + if i > 0 { z[k - n] } else { 0.0 }); } }
        for i in (0..n).rev() { for j in (0..n).rev() { let k = i * n + j; p[k] = z[k] + 0.25 * (if j + 1 < n { p[k + 1] } else { 0.0 } + if i + 1 < n { p[k + n] } else { 0.0 }); } }
        z.copy_from_slice(&p); rz = dot_seq(&r, &z);
    }
    for _ in 1..=max_it {
        let rel;
        if meth == 0 {
            for i in 0..n { for j in 0..n {
                let k = i * n + j;
                let w = if j > 0 { u[k - 1] } else { 0.0 }; let e = if j + 1 < n { u[k + 1] } else { 0.0 };
                let s = if i > 0 { u[k - n] } else { 0.0 }; let nn = if i + 1 < n { u[k + n] } else { 0.0 };
                z[k] = 0.25 * ((((b[k] + w) + e) + s) + nn);
            } }
            u.copy_from_slice(&z); rel = compute_res(n, &u, &b, &mut r) / nb;
        } else if meth == 1 || meth == 3 {
            for i in 0..n { for j in 0..n {
                let k = i * n + j;
                let w = if j > 0 { u[k - 1] } else { 0.0 }; let e = if j + 1 < n { u[k + 1] } else { 0.0 };
                let s = if i > 0 { u[k - n] } else { 0.0 }; let nn = if i + 1 < n { u[k + n] } else { 0.0 };
                let gs = 0.25 * ((((b[k] + w) + e) + s) + nn);
                u[k] = if meth == 1 { gs } else { (1.0 - om) * u[k] + om * gs };
            } }
            rel = compute_res(n, &u, &b, &mut r) / nb;
        } else if meth == 2 {
            rbgs_sweep(n, &mut u, &b); rel = compute_res(n, &u, &b, &mut r) / nb;
        } else if meth == 4 {
            apply_a(n, &p, &mut ap);
            let alpha = rr / dot_seq(&p, &ap);
            for k in 0..n2 { u[k] += alpha * p[k]; r[k] -= alpha * ap[k]; }
            let rr_new = dot_seq(&r, &r); rel = rr_new.sqrt() / nb; hist.push(rel);
            if rel <= 1e-8 { break; }
            let beta = rr_new / rr; for k in 0..n2 { p[k] = r[k] + beta * p[k]; }
            rr = rr_new; continue;
        } else if meth == 5 {
            apply_a(n, &p, &mut ap);
            let alpha = rz / dot_seq(&p, &ap);
            for k in 0..n2 { u[k] += alpha * p[k]; r[k] -= alpha * ap[k]; }
            rel = norm2_seq(&r) / nb; hist.push(rel);
            if rel <= 1e-8 { break; }
            for i in 0..n { for j in 0..n { let k = i * n + j; ap[k] = 0.25 * ((r[k] + if j > 0 { ap[k - 1] } else { 0.0 }) + if i > 0 { ap[k - n] } else { 0.0 }); } }
            for i in (0..n).rev() { for j in (0..n).rev() { let k = i * n + j; z[k] = ap[k] + 0.25 * (if j + 1 < n { z[k + 1] } else { 0.0 } + if i + 1 < n { z[k + n] } else { 0.0 }); } }
            let rz_new = dot_seq(&r, &z); let beta = rz_new / rz;
            for k in 0..n2 { p[k] = z[k] + beta * p[k]; }
            rz = rz_new; continue;
        } else {
            vcycle(n, &mut u, &b); rel = compute_res(n, &u, &b, &mut r) / nb;
        }
        hist.push(rel);
        if rel <= 1e-8 { break; }
    }
    let used = hist.len();
    let rho = if used >= 51 { (hist[used - 1] / hist[used - 51]).powf(0.02) } else if used >= 2 { (hist[used - 1] / hist[0]).powf(1.0 / ((used - 1) as f64)) } else { 0.0 };
    if save63 && n == 63 { st.hist63[meth] = hist.clone(); }
    SolveRes { iters: used, rel_res: *hist.last().unwrap(), err_inf: max_abs_diff(&u, &ue), rho_emp: rho, hist_hash: fnv1a_f64(&hist), sol_hash: fnv1a_f64(&u) }
}

fn s01(_st: &mut State) {
    let u = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]; let mut au = [0.0; 9];
    apply_a(3, &u, &mut au);
    let (h1, h2, h3) = (fnv1a_f64(&u), fnv1a_f64(&au), fnv1a_f64(&[1.0, -2.5, 3.25]));
    let ok = max_abs_diff(&au, &[-2.0,-1.0,4.0,3.0,0.0,7.0,16.0,11.0,22.0]) == 0.0 && h1 == 0xc56de10c6671e6c2 && h2 == 0x59ed72ec82ee959c && h3 == 0xafc5dc1ca48a2366;
    fs::write("out/rust/s01_selftest.txt", format!("ok,h_u3,h_Au3,h_test\n{},{:016x},{:016x},{:016x}\n", if ok {1} else {0}, h1, h2, h3)).unwrap();
    if !ok { std::process::exit(1); }
}
fn s02(_st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s02_problem.csv").unwrap());
    writeln!(f, "N,norm_b,h_b,h_ue").unwrap();
    for &n in &NS_ALL {
        let (b, ue) = build_prob(n);
        writeln!(f, "{},{:.17},{:016x},{:016x}", n, norm2_seq(&b), fnv1a_f64(&b), fnv1a_f64(&ue)).unwrap();
    }
}
fn s03(_st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s03_eigen.csv").unwrap());
    writeln!(f, "N,j,k,lambda_exact,rel_err_inf").unwrap();
    for &n in &NS_ALL[..3] {
        let h = 1.0 / ((n + 1) as f64); let mut v = vec![0.0; n * n]; let mut av = vec![0.0; n * n];
        for (j, k) in [(1,1), (1,2), (3,5)] {
            let lam = 4.0 - 2.0 * cos_spec((j as f64) * PI_VAL * h) - 2.0 * cos_spec((k as f64) * PI_VAL * h);
            for r in 0..n { for c in 0..n { v[r * n + c] = sin_spec((j as f64) * PI_VAL * ((r + 1) as f64) * h) * sin_spec((k as f64) * PI_VAL * ((c + 1) as f64) * h); } }
            apply_a(n, &v, &mut av);
            let me = (0..n * n).map(|q| (av[q] - lam * v[q]).abs()).fold(0.0f64, f64::max);
            writeln!(f, "{},{},{},{:.17},{:.3e}", n, j, k, lam, me / lam).unwrap();
        }
    }
}
fn run_m(st: &mut State, m: usize, fnm: &str, nmax: usize) {
    let mut f = BufWriter::new(File::create(fnm).unwrap());
    writeln!(f, "N,iters,rel_res,err_inf,rho_emp").unwrap();
    for i in 0..nmax {
        let r = run_solver(m, NS_ALL[i], st, NS_ALL[i] == 63);
        writeln!(f, "{},{},{:.10e},{:.10e},{:.10}", NS_ALL[i], r.iters, r.rel_res, r.err_inf, r.rho_emp).unwrap();
        st.res[m][i] = r;
    }
}
fn s04(st: &mut State) { run_m(st, 0, "out/rust/s04_jacobi.csv", 4); }
fn s05(st: &mut State) { run_m(st, 1, "out/rust/s05_gs.csv", 4); }
fn s06(st: &mut State) { run_m(st, 2, "out/rust/s06_rbgs.csv", 4); }
fn s07(st: &mut State) { run_m(st, 3, "out/rust/s07_sor.csv", 6); }
fn s08(st: &mut State) { run_m(st, 4, "out/rust/s08_cg.csv", 6); }
fn s09(st: &mut State) { run_m(st, 5, "out/rust/s09_pcg_ssor.csv", 6); }
fn s10(st: &mut State) { run_m(st, 6, "out/rust/s10_mg.csv", 6); }
fn s11(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s11_fmg.csv").unwrap());
    writeln!(f, "N,fmg_err_inf,disc_err_inf,ratio").unwrap();
    for (i, &n) in NS_ALL.iter().enumerate() {
        let (b, ue) = build_prob(n); let mut u = vec![0.0; n * n]; fmg_solve(n, &mut u, &b);
        let fe = max_abs_diff(&u, &ue); st.disc_err[i] = st.res[6][i].err_inf;
        writeln!(f, "{},{:.10e},{:.10e},{:.6}", n, fe, st.disc_err[i], fe / st.disc_err[i]).unwrap();
    }
}
fn s12(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s12_order.csv").unwrap());
    writeln!(f, "N,h,err_inf,slope").unwrap();
    let xs: Vec<f64> = NS_ALL.iter().map(|&n| (1.0 / ((n + 1) as f64)).ln()).collect();
    let ys: Vec<f64> = st.disc_err.iter().map(|&e| e.ln()).collect();
    let mx = xs.iter().sum::<f64>() / 6.0; let my = ys.iter().sum::<f64>() / 6.0;
    let sl = (0..6).map(|i| (xs[i] - mx) * (ys[i] - my)).sum::<f64>() / xs.iter().map(|&x| (x - mx).powi(2)).sum::<f64>();
    for i in 0..6 { writeln!(f, "{},{:.10e},{:.10e},{:.6}", NS_ALL[i], 1.0 / ((NS_ALL[i] + 1) as f64), st.disc_err[i], sl).unwrap(); }
}
fn s13(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s13_scaling.csv").unwrap());
    writeln!(f, "method,N,iters,slope,slope_se").unwrap();
    for m in 0..7 {
        let cnt = if m <= 2 { 4 } else { 6 };
        let xs: Vec<f64> = (0..cnt).map(|i| (NS_ALL[i] as f64).ln()).collect();
        let ys: Vec<f64> = (0..cnt).map(|i| (st.res[m][i].iters as f64).ln()).collect();
        let mx = xs.iter().sum::<f64>() / (cnt as f64); let my = ys.iter().sum::<f64>() / (cnt as f64);
        let den = xs.iter().map(|&x| (x - mx).powi(2)).sum::<f64>();
        let sl = (0..cnt).map(|i| (xs[i] - mx) * (ys[i] - my)).sum::<f64>() / den;
        let ic = my - sl * mx;
        let rss = (0..cnt).map(|i| (ys[i] - (sl * xs[i] + ic)).powi(2)).sum::<f64>();
        let se = ((rss / ((cnt - 2) as f64)) / den).sqrt();
        for i in 0..cnt { writeln!(f, "{},{},{},{:.6},{:.6}", MNAMES[m], NS_ALL[i], st.res[m][i].iters, sl, se).unwrap(); }
    }
}
fn s14(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s14_history.csv").unwrap());
    writeln!(f, "method,iter,rel_res").unwrap();
    for m in 0..7 {
        let len = st.hist63[m].len(); let step = if len > 200 { len / 200 } else { 1 };
        for k in (0..len).step_by(step) { writeln!(f, "{},{},{:.10e}", MNAMES[m], k + 1, st.hist63[m][k]).unwrap(); }
        if (len - 1) % step != 0 { writeln!(f, "{},{},{:.10e}", MNAMES[m], len, st.hist63[m][len - 1]).unwrap(); }
    }
}
fn s15(_st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s15_mg_factor.csv").unwrap());
    writeln!(f, "N,iters,mg_factor_last5").unwrap();
    for &n in &[63usize, 127] {
        let (b, _) = build_prob(n); let mut u = vec![0.0; n * n]; let mut r = vec![0.0; n * n];
        let nb = norm2_seq(&b); let mut h = Vec::new();
        for _ in 1..=20 { vcycle(n, &mut u, &b); let rel = compute_res(n, &u, &b, &mut r) / nb; h.push(rel); if rel <= 1e-8 { break; } }
        let it = h.len(); writeln!(f, "{},{},{:.8}", n, it, (h[it - 1] / h[it - 6]).powf(0.2)).unwrap();
    }
}
fn s16(_st: &mut State) {
    let flops = [9, 9, 9, 12, 21, 35, 36];
    let mut f = BufWriter::new(File::create("out/rust/s16_cost.csv").unwrap());
    writeln!(f, "method,N,flops_per_unk_iter,ns_med,ns_mad").unwrap();
    let n = 63usize; let n2 = n * n; let (b, _) = build_prob(n); let (mut u, mut z) = (vec![0.0; n2], vec![0.0; n2]);
    for m in 0..7 {
        let iters = if m == 6 { 10 } else { 80 }; let mut vals = [0.0; 7];
        for rp in 0..7 {
            u.fill(0.0); let t0 = Instant::now();
            for _ in 0..iters {
                if m == 6 { vcycle(n, &mut u, &b); } else if m == 2 { rbgs_sweep(n, &mut u, &b); }
                else { apply_a(n, &u, &mut z); for k in 0..n2 { u[k] += 0.25 * (b[k] - z[k]); } }
            }
            vals[rp] = t0.elapsed().as_secs_f64() * 1e9 / ((iters * n2) as f64);
        }
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap()); let med = vals[3];
        let mut dev = vals.map(|v| (v - med).abs()); dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        writeln!(f, "{},{},{},{:.4},{:.4}", MNAMES[m], n, flops[m], med, dev[3]).unwrap();
    }
}
fn s17(st: &mut State) {
    let mut f = BufWriter::new(File::create("out/rust/s17_parity.csv").unwrap());
    writeln!(f, "method,N,iters,hist_fnv1a,sol_fnv1a").unwrap();
    for m in 0..7 { for ni in 0..2 { writeln!(f, "{},{},{},{:016x},{:016x}", MNAMES[m], NS_ALL[ni], st.res[m][ni].iters, st.res[m][ni].hist_hash, st.res[m][ni].sol_hash).unwrap(); } }
}
fn s18(st: &mut State) {
    let cols = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2"];
    let mut s1 = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for m in 0..7 {
        s1.push_str(&format!("<polyline fill=\"none\" stroke=\"{}\" stroke-width=\"1.4\" points=\"", cols[m]));
        let len = st.hist63[m].len(); let step = if len > 150 { len / 150 } else { 1 };
        for k in (0..len).step_by(step) {
            s1.push_str(&format!("{:.1},{:.1} ", 60.0 + 500.0 * ((k + 1) as f64).log10() / (st.hist63[0].len() as f64).log10(), 310.0 - 280.0 * (-st.hist63[m][k].max(1e-12).log10() / 9.0).clamp(0.0, 1.0)));
        }
        s1.push_str("\"/>");
    }
    fs::write("out/rust/figures/history.svg", s1 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Rust Residual History (N=63)</text></svg>\n").unwrap();
    let mut s2 = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for m in 0..7 {
        let cnt = if m <= 2 { 4 } else { 6 };
        s2.push_str(&format!("<polyline fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\" points=\"", cols[m]));
        for i in 0..cnt { s2.push_str(&format!("{:.1},{:.1} ", 60.0 + 500.0 * ((NS_ALL[i] as f64).log10() - 15.0f64.log10()) / (511.0f64.log10() - 15.0f64.log10()), 310.0 - 280.0 * (st.res[m][i].iters as f64).log10() / 5.0)); }
        s2.push_str("\"/>");
    }
    fs::write("out/rust/figures/scaling.svg", s2 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Rust Iterations vs N</text></svg>\n").unwrap();
    let mut s3 = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"1.8\" points=\"");
    for i in 0..6 { s3.push_str(&format!("{:.1},{:.1} ", 60.0 + 500.0 * (i as f64 / 5.0), 310.0 - 280.0 * (st.disc_err[i].log10() + 5.0) / 4.0)); }
    fs::write("out/rust/figures/order.svg", s3 + "\"/><text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Rust Discretization Error vs h (Order ~2)</text></svg>\n").unwrap();
}
fn s19(st: &mut State) {
    let tot: f64 = st.stage_s.iter().sum();
    fs::write("out/rust/s19_summary.json", format!("{{\"lang\":\"rust\",\"N_max\":511,\"N_stat_max\":127,\"toolchain\":\"rustc 1.93.1\",\"total_sec\":{:.3}}}\n", tot)).unwrap();
}
fn main() {
    fs::create_dir_all("out/rust/figures").unwrap();
    let cmd = std::env::args().nth(1).unwrap_or_else(|| "all".to_string());
    let mut st = State::new();
    let fns: [(&str, fn(&mut State)); 19] = [("s01",s01),("s02",s02),("s03",s03),("s04",s04),("s05",s05),("s06",s06),("s07",s07),("s08",s08),("s09",s09),("s10",s10),("s11",s11),("s12",s12),("s13",s13),("s14",s14),("s15",s15),("s16",s16),("s17",s17),("s18",s18),("s19",s19)];
    if cmd == "all" {
        for (i, (nm, f)) in fns.iter().enumerate() { println!("[rust] tahap {}/19 ({})", i+1, nm); let t0 = Instant::now(); f(&mut st); st.stage_s[i] = t0.elapsed().as_secs_f64(); }
        return;
    }
    for (nm, f) in fns { if cmd == nm { f(&mut st); return; } }
    std::process::exit(1);
}
