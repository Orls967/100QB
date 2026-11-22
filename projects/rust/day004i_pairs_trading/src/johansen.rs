//! Bivariate Johansen trace cointegration test (unrestricted constant, rank 0 vs >= 1):
//! Cholesky reduction of `S01 * S11^{-1} * S10 v = lambda * S00 v` to symmetric 2x2 eigenproblem,
//! H0 critical value simulation (R=20,000, T=500), and size/power on spurious vs cointegrated pairs.

use crate::linalg::{cholesky_2x2, eig_sym_2x2_desc, inv_2x2, solve_lu};
use crate::ou::{simulate_ou_path, OuParams};
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;

/// Computes bivariate Johansen trace statistic for `H0: rank(Pi) = 0` with unrestricted constant.
pub fn johansen_trace_bivariate(x: &[f64], y: &[f64]) -> f64 {
    let t_len = x.len().min(y.len());
    if t_len < 5 {
        return 0.0;
    }
    let n = t_len - 1;
    let nf = n as f64;

    let mut m_dx = 0.0_f64;
    let mut m_dy = 0.0_f64;
    let mut m_lx = 0.0_f64;
    let mut m_ly = 0.0_f64;
    for t in 0..n {
        m_dx += x[t + 1] - x[t];
        m_dy += y[t + 1] - y[t];
        m_lx += x[t];
        m_ly += y[t];
    }
    m_dx /= nf;
    m_dy /= nf;
    m_lx /= nf;
    m_ly /= nf;

    let mut s00 = [[0.0_f64; 2]; 2];
    let mut s11 = [[0.0_f64; 2]; 2];
    let mut s01 = [[0.0_f64; 2]; 2];
    for t in 0..n {
        let r0 = [(x[t + 1] - x[t]) - m_dx, (y[t + 1] - y[t]) - m_dy];
        let r1 = [x[t] - m_lx, y[t] - m_ly];
        for i in 0..2 {
            for j in 0..2 {
                s00[i][j] += r0[i] * r0[j];
                s11[i][j] += r1[i] * r1[j];
                s01[i][j] += r0[i] * r1[j];
            }
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            s00[i][j] /= nf;
            s11[i][j] /= nf;
            s01[i][j] /= nf;
        }
    }

    let s11_inv = inv_2x2(&s11);
    // M = S01 * S11^{-1} * S01^T
    let mut tmp = [[0.0_f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            tmp[i][j] = s01[i][0] * s11_inv[0][j] + s01[i][1] * s11_inv[1][j];
        }
    }
    let mut m_mat = [[0.0_f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            m_mat[i][j] = tmp[i][0] * s01[j][0] + tmp[i][1] * s01[j][1];
        }
    }

    // Cholesky S00 = L L^T -> L^{-1} M L^{-T}
    let l = cholesky_2x2(&s00);
    let l_inv = [
        [1.0 / l[0][0], 0.0],
        [-l[1][0] / (l[0][0] * l[1][1]), 1.0 / l[1][1]],
    ];
    // sym = L^{-1} * M * L^{-T}
    let mut lm = [[0.0_f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            lm[i][j] = l_inv[i][0] * m_mat[0][j] + l_inv[i][1] * m_mat[1][j];
        }
    }
    let mut sym = [[0.0_f64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            sym[i][j] = lm[i][0] * l_inv[j][0] + lm[i][1] * l_inv[j][1];
        }
    }
    let eigs = eig_sym_2x2_desc(&sym);
    let mut trace_stat = 0.0_f64;
    for &lam in &eigs {
        let clamped = lam.clamp(0.0, 1.0 - 1e-15);
        trace_stat += -nf * (1.0 - clamped).ln();
    }
    trace_stat
}

#[derive(Clone, Debug)]
pub struct JohansenSummary {
    pub n_paths: usize,
    pub t_len: usize,
    pub q90: f64,
    pub q95: f64,
    pub se_q95: f64,
    pub q99: f64,
    pub q95_h0_drift: f64,
    pub se_q95_h0_drift: f64,
    pub spurious_size_at_q95: f64,
    pub spurious_size_at_q95_drift: f64,
    pub spurious_size_at_lit_1541: f64,
    pub coint_power_theta16_2yr: f64,
    pub coint_power_at_q95_drift: f64,
    pub coint_power_at_lit_1541: f64,
    pub brownian_q95_case1_no_const: f64,
    pub brownian_q95_case2_rest_const: f64,
    pub brownian_q95_case3_unrest_const: f64,
    pub brownian_q95_case3_nonzero_drift: f64,
    pub t_sweep_q95: Vec<(usize, f64, f64)>,
}

/// Simulates asymptotic Brownian functional `tr(int dB B^T (int B B^T dt)^{-1} int B dB^T)`
/// with `T=5000, R=20000` for Case 1 (raw B), Case 2 (augmented [B, 1]), Case 3 zero drift (demeaned B),
/// and Case 3 non-zero drift (`[B1 - m1, u - 0.5]`, Osterwald-Lenum 15.41).
pub fn johansen_brownian_functional_q95(
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let dt = 1.0 / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let res = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3600, r as u64));
        let mut b1 = 0.0_f64;
        let mut b2 = 0.0_f64;
        let mut sb1 = 0.0_f64;
        let mut sb2 = 0.0_f64;
        let mut b11 = 0.0_f64;
        let mut b12 = 0.0_f64;
        let mut b22 = 0.0_f64;
        let mut bd11 = 0.0_f64;
        let mut bd12 = 0.0_f64;
        let mut bd21 = 0.0_f64;
        let mut bd22 = 0.0_f64;
        let mut ub1 = 0.0_f64;
        let mut udb1 = 0.0_f64;
        let mut udb2 = 0.0_f64;
        for step_idx in 0..n_steps {
            let u = step_idx as f64 * dt - 0.5;
            let db1 = sqrt_dt * rg.normal();
            let db2 = sqrt_dt * rg.normal();
            sb1 += b1;
            sb2 += b2;
            b11 += b1 * b1;
            b12 += b1 * b2;
            b22 += b2 * b2;
            bd11 += b1 * db1;
            bd12 += b1 * db2;
            bd21 += b2 * db1;
            bd22 += b2 * db2;
            ub1 += u * b1;
            udb1 += u * db1;
            udb2 += u * db2;
            b1 += db1;
            b2 += db2;
        }
        let m1 = sb1 * dt;
        let m2 = sb2 * dt;
        let i11 = b11 * dt;
        let i12 = b12 * dt;
        let i22 = b22 * dt;

        // Case 1: raw B (2x2)
        let inv1 = inv_2x2(&[[i11, i12], [i12, i22]]);
        let c1 = bd11 * (inv1[0][0] * bd11 + inv1[0][1] * bd21)
            + bd21 * (inv1[1][0] * bd11 + inv1[1][1] * bd21)
            + bd12 * (inv1[0][0] * bd12 + inv1[0][1] * bd22)
            + bd22 * (inv1[1][0] * bd12 + inv1[1][1] * bd22);

        // Case 3: demeaned B (2x2), int B_dem dB^T = int B dB^T - [m1; m2]*[b1, b2]
        let inv3 = inv_2x2(&[
            [i11 - m1 * m1, i12 - m1 * m2],
            [i12 - m1 * m2, i22 - m2 * m2],
        ]);
        let d11 = bd11 - m1 * b1;
        let d12 = bd12 - m1 * b2;
        let d21 = bd21 - m2 * b1;
        let d22 = bd22 - m2 * b2;
        let c3 = d11 * (inv3[0][0] * d11 + inv3[0][1] * d21)
            + d21 * (inv3[1][0] * d11 + inv3[1][1] * d21)
            + d12 * (inv3[0][0] * d12 + inv3[0][1] * d22)
            + d22 * (inv3[1][0] * d12 + inv3[1][1] * d22);

        // Case 2: restricted constant => augmented 3x3 F = [B1, B2, 1]
        let f_mat = vec![vec![i11, i12, m1], vec![i12, i22, m2], vec![m1, m2, 1.0]];
        let col1 = solve_lu(&f_mat, &[bd11, bd21, b1]).unwrap_or(vec![0.0; 3]);
        let col2 = solve_lu(&f_mat, &[bd12, bd22, b2]).unwrap_or(vec![0.0; 3]);
        let c2 = bd11 * col1[0]
            + bd21 * col1[1]
            + b1 * col1[2]
            + bd12 * col2[0]
            + bd22 * col2[1]
            + b2 * col2[2];

        // Case 3 non-zero drift: F = [B1 - m1, u - 0.5]
        let ub1_int = ub1 * dt;
        let inv3d = inv_2x2(&[[i11 - m1 * m1, ub1_int], [ub1_int, 1.0 / 12.0]]);
        let e11 = bd11 - m1 * b1;
        let e12 = bd12 - m1 * b2;
        let e21 = udb1;
        let e22 = udb2;
        let c3d = e11 * (inv3d[0][0] * e11 + inv3d[0][1] * e21)
            + e21 * (inv3d[1][0] * e11 + inv3d[1][1] * e21)
            + e12 * (inv3d[0][0] * e12 + inv3d[0][1] * e22)
            + e22 * (inv3d[1][0] * e12 + inv3d[1][1] * e22);

        (c1, c2, c3, c3d)
    });
    let mut v1: Vec<f64> = res.iter().map(|x| x.0).collect();
    let mut v2: Vec<f64> = res.iter().map(|x| x.1).collect();
    let mut v3: Vec<f64> = res.iter().map(|x| x.2).collect();
    let mut v3d: Vec<f64> = res.iter().map(|x| x.3).collect();
    v1.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v2.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v3.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v3d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    (
        stats::quantile_and_se(&v1, 0.95).0,
        stats::quantile_and_se(&v2, 0.95).0,
        stats::quantile_and_se(&v3, 0.95).0,
        stats::quantile_and_se(&v3d, 0.95).0,
    )
}

pub fn run_johansen_study(
    t_len: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> JohansenSummary {
    let mut stats_h0 = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3500, r as u64));
        let mut x = Vec::with_capacity(t_len);
        let mut y = Vec::with_capacity(t_len);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..t_len {
            cx += rg.normal();
            cy += rg.normal();
            x.push(cx);
            y.push(cy);
        }
        johansen_trace_bivariate(&x, &y)
    });
    stats_h0.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (q90, _) = stats::quantile_and_se(&stats_h0, 0.90);
    let (q95, se_q95) = stats::quantile_and_se(&stats_h0, 0.95);
    let (q99, _) = stats::quantile_and_se(&stats_h0, 0.99);

    // Out-of-sample spurious size (independent seed) vs cointegrated power (2yr, theta=16)
    let dt = 1.0_f64 / 252.0;
    let sqrt_dt = dt.sqrt();
    let mut stats_h0_drift = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed ^ 0x3520, 0x3521, r as u64));
        let mut x = Vec::with_capacity(t_len);
        let mut y = Vec::with_capacity(t_len);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..t_len {
            cx += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            cy += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            x.push(cx);
            y.push(cy);
        }
        johansen_trace_bivariate(&x, &y)
    });
    stats_h0_drift.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (q95_h0_drift, se_q95_h0_drift) = stats::quantile_and_se(&stats_h0_drift, 0.95);

    let p_ou = OuParams::default_params();
    let oos = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed ^ 0x7777_0001, 0x3550, r as u64));
        let mut x_sp = Vec::with_capacity(t_len);
        let mut y_sp = Vec::with_capacity(t_len);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..t_len {
            cx += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            cy += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            x_sp.push(cx);
            y_sp.push(cy);
        }
        let tr_sp = johansen_trace_bivariate(&x_sp, &y_sp);
        let rej_sp = tr_sp > q95;
        let rej_sp_drift = tr_sp > q95_h0_drift;
        let rej_lit = tr_sp > 15.41;

        let n_2yr = 504;
        let spread = simulate_ou_path(p_ou, n_2yr - 1, true, &mut rg);
        let mut x_co = Vec::with_capacity(n_2yr);
        let mut y_co = Vec::with_capacity(n_2yr);
        let mut cxc = 0.0_f64;
        for &sp in &spread {
            cxc += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            x_co.push(cxc);
            y_co.push(0.1 + 1.2 * cxc + sp);
        }
        let tr_co = johansen_trace_bivariate(&x_co, &y_co);
        let rej_co = tr_co > q95;
        let rej_co_drift = tr_co > q95_h0_drift;
        let rej_co_lit = tr_co > 15.41;
        (
            rej_sp,
            rej_sp_drift,
            rej_lit,
            rej_co,
            rej_co_drift,
            rej_co_lit,
        )
    });
    let nf = oos.len() as f64;
    let sp_size = oos.iter().filter(|x| x.0).count() as f64 / nf;
    let sp_size_drift = oos.iter().filter(|x| x.1).count() as f64 / nf;
    let sp_lit = oos.iter().filter(|x| x.2).count() as f64 / nf;
    let co_pwr = oos.iter().filter(|x| x.3).count() as f64 / nf;
    let co_pwr_drift = oos.iter().filter(|x| x.4).count() as f64 / nf;
    let co_pwr_lit = oos.iter().filter(|x| x.5).count() as f64 / nf;

    let (bm1, bm2, bm3, bm3d) = johansen_brownian_functional_q95(5000, n_paths, n_workers, seed);
    let mut t_sweep_q95 = Vec::new();
    for &t_sw in &[250_usize, 500, 1000, 5000] {
        let r_sw = if t_sw >= 5000 { n_paths / 4 } else { n_paths };
        let mut v_sw = parallel::run_chunks(r_sw, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x3700 ^ (t_sw as u64), r as u64));
            let mut x = Vec::with_capacity(t_sw);
            let mut y = Vec::with_capacity(t_sw);
            let mut cx = 0.0_f64;
            let mut cy = 0.0_f64;
            for _ in 0..t_sw {
                cx += rg.normal();
                cy += rg.normal();
                x.push(cx);
                y.push(cy);
            }
            johansen_trace_bivariate(&x, &y)
        });
        v_sw.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let (qt, se_qt) = stats::quantile_and_se(&v_sw, 0.95);
        t_sweep_q95.push((t_sw, qt, se_qt));
    }

    JohansenSummary {
        n_paths,
        t_len,
        q90,
        q95,
        se_q95,
        q99,
        q95_h0_drift,
        se_q95_h0_drift,
        spurious_size_at_q95: sp_size,
        spurious_size_at_q95_drift: sp_size_drift,
        spurious_size_at_lit_1541: sp_lit,
        coint_power_theta16_2yr: co_pwr,
        coint_power_at_q95_drift: co_pwr_drift,
        coint_power_at_lit_1541: co_pwr_lit,
        brownian_q95_case1_no_const: bm1,
        brownian_q95_case2_rest_const: bm2,
        brownian_q95_case3_unrest_const: bm3,
        brownian_q95_case3_nonzero_drift: bm3d,
        t_sweep_q95,
    }
}
