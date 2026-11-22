//! Bivariate Johansen trace cointegration test (unrestricted constant, rank 0 vs >= 1):
//! Cholesky reduction of `S01 * S11^{-1} * S10 v = lambda * S00 v` to symmetric 2x2 eigenproblem,
//! H0 critical value simulation (R=20,000, T=500), and size/power on spurious vs cointegrated pairs.

use crate::linalg::{cholesky_2x2, eig_sym_2x2_desc, inv_2x2};
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
    pub spurious_size_at_q95: f64,
    pub coint_power_theta16_2yr: f64,
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

    // Out-of-sample spurious size vs cointegrated power (2yr, theta=16)
    let dt = 1.0_f64 / 252.0;
    let sqrt_dt = dt.sqrt();
    let p_ou = OuParams::default_params();
    let oos = parallel::run_chunks(n_paths / 2, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3550, r as u64));
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
        let rej_sp = johansen_trace_bivariate(&x_sp, &y_sp) > q95;

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
        let rej_co = johansen_trace_bivariate(&x_co, &y_co) > q95;
        (rej_sp, rej_co)
    });
    let nf = oos.len() as f64;
    let sp_size = oos.iter().filter(|x| x.0).count() as f64 / nf;
    let co_pwr = oos.iter().filter(|x| x.1).count() as f64 / nf;

    JohansenSummary {
        n_paths,
        t_len,
        q90,
        q95,
        se_q95,
        q99,
        spurious_size_at_q95: sp_size,
        coint_power_theta16_2yr: co_pwr,
    }
}
