//! Bagian 2: Dampak transien (Obizhaeva-Wang), solver tridiagonal O(N), kernel lain & manipulasi harga.

use std::time::Instant;

use crate::ac::AcParams;
use crate::linalg::{
    cholesky, dot, jacobi_extreme_eigenvalues, jacobi_min_eigenvalue, quad_form, solve_spd,
};

/// Parameter model Obizhaeva-Wang.
#[derive(Clone, Copy, Debug)]
pub struct OwParams {
    pub x_total: f64,
    pub t_horizon: f64,
    pub n_steps: usize,
    pub lambda_ow: f64,
    pub rho: f64,
}

impl Default for OwParams {
    fn default() -> Self {
        Self {
            x_total: 1_000_000.0,
            t_horizon: 1.0,
            n_steps: 50,
            lambda_ow: 5e-7,
            rho: 5.0,
        }
    }
}

/// Bangun matriks kernel ``G[i][j] = f(|t_i - t_j|)`` pada grid seragam ``t_i = (i+1)*Delta``.
pub fn build_kernel_matrix<F>(n: usize, t_horizon: f64, kernel_fn: F) -> Vec<Vec<f64>>
where
    F: Fn(f64) -> f64,
{
    let dt = t_horizon / (n as f64);
    let mut g = vec![vec![0.0; n]; n];
    for (i, row) in g.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let dist = ((i as f64) - (j as f64)).abs() * dt;
            *cell = kernel_fn(dist);
        }
    }
    g
}

/// Solusi umum eksak ``x* = X * G^{-1} 1 / (1^T G^{-1} 1)`` lewat Cholesky.
pub fn solve_transient_general(g: &[Vec<f64>], x_total: f64) -> Result<(Vec<f64>, f64), usize> {
    let n = g.len();
    let ones = vec![1.0; n];
    let g_inv_1 = solve_spd(g, &ones)?;
    let denom = dot(&ones, &g_inv_1);
    let x: Vec<f64> = g_inv_1.into_iter().map(|v| x_total * v / denom).collect();
    let cost = 0.5 * quad_form(g, &x);
    Ok((x, cost))
}

/// Solusi analitik O(N) dari ``G^{-1}`` tridiagonal untuk kernel eksponensial Obizhaeva-Wang.
pub fn solve_ow_exponential_analytic(p: &OwParams) -> (Vec<f64>, f64) {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let a = (-p.rho * dt).exp();
    let denom = 2.0 + ((n as f64) - 2.0) * (1.0 - a);
    let edge = p.x_total / denom;
    let interior = p.x_total * (1.0 - a) / denom;

    let mut x = vec![interior; n];
    if n >= 1 {
        x[0] = edge;
        x[n - 1] = edge;
    }
    let cost = 0.5 * p.lambda_ow * p.x_total * p.x_total * (1.0 + a) / denom;
    (x, cost)
}

/// Bandingkan waktu & akurasi Cholesky penuh O(N^3) vs solver analitik tridiagonal O(N).
pub fn benchmark_ow_cholesky_vs_analytic(p: &OwParams, n_large: usize) -> (f64, f64, f64) {
    let mut p_lg = *p;
    p_lg.n_steps = n_large;
    let g = build_kernel_matrix(n_large, p_lg.t_horizon, |d| {
        p_lg.lambda_ow * (-p_lg.rho * d).exp()
    });

    let t0 = Instant::now();
    let (x_chol, _) = solve_transient_general(&g, p_lg.x_total).expect("OW G harus SPD");
    let dt_chol = t0.elapsed().as_secs_f64();

    let t1 = Instant::now();
    let (x_ana, _) = solve_ow_exponential_analytic(&p_lg);
    let dt_ana = t1.elapsed().as_secs_f64();

    let max_rel_err = x_chol
        .iter()
        .zip(x_ana.iter())
        .map(|(a, b)| (a - b).abs() / p_lg.x_total)
        .fold(0.0f64, f64::max);

    (max_rel_err, dt_chol, dt_ana)
}

/// Uji kernel berosilasi ``G(t) = lambda_ow * exp(-rho_osc * t) * cos(omega * t)`` (+ diagonal kecil):
/// mengembalikan (is_pd, min_eig, x_star, min_x, cost).
pub fn analyze_oscillating_kernel(
    n: usize,
    t_horizon: f64,
    lambda_ow: f64,
    rho_osc: f64,
    omega: f64,
    diag_eps: f64,
    x_total: f64,
) -> (bool, f64, Vec<f64>, f64, f64) {
    let mut g = build_kernel_matrix(n, t_horizon, |d| {
        lambda_ow * (-rho_osc * d).exp() * (omega * d).cos()
    });
    for (i, row) in g.iter_mut().enumerate() {
        row[i] += diag_eps * lambda_ow;
    }
    let is_pd = cholesky(&g).is_ok();
    let (min_eig, _) = jacobi_min_eigenvalue(&g);
    let (x, cost) = solve_transient_general(&g, x_total).unwrap_or((vec![0.0; n], f64::NAN));
    let min_x = x.iter().copied().fold(f64::INFINITY, f64::min);
    (is_pd, min_eig, x, min_x, cost)
}

/// Uji kernel NON-positif-definit: Cholesky gagal dan ada arah penurunan ``1^T d = 0`` dengan ``d^T G d < 0``.
pub fn analyze_non_pd_kernel(
    n: usize,
    t_horizon: f64,
    lambda_ow: f64,
    x_total: f64,
) -> (bool, f64, f64, f64) {
    // Kernel yang melanggar positif-definit: G(d) = lambda_ow * (1 - 6 * d^2) (cekung di d=0 -> indefinit)
    let g = build_kernel_matrix(n, t_horizon, |d| lambda_ow * (1.0 - 6.0 * d * d));
    let chol_failed = cholesky(&g).is_err();
    let (min_eig, eigvec) = jacobi_min_eigenvalue(&g);

    // Proyeksikan eigvec ke subruang {1^T d = 0} agar x_twap + alpha * d tetap memenuhi kendala sum(x) = X
    let mean_v: f64 = eigvec.iter().sum::<f64>() / (n as f64);
    let d: Vec<f64> = eigvec.iter().map(|&v| v - mean_v).collect();
    let d_quad = quad_form(&g, &d);

    let x_base = vec![x_total / (n as f64); n];
    let alpha = 10.0 * x_total;
    let x_pert: Vec<f64> = x_base
        .iter()
        .zip(d.iter())
        .map(|(a, b)| a + alpha * b)
        .collect();
    let cost_pert = 0.5 * quad_form(&g, &x_pert);

    (chol_failed, min_eig, d_quad, cost_pert)
}

/// Gabungan kernel transien ``G`` + suku risiko Almgren-Chriss ``lambda * sigma^2 * dt * sum y_j^2``.
pub fn solve_transient_with_risk(
    g: &[Vec<f64>],
    p_ac: &AcParams,
) -> Result<(Vec<f64>, f64), usize> {
    let n = p_ac.n_steps;
    assert_eq!(g.len(), n);
    let dt = p_ac.t_horizon / (n as f64);
    let risk_coeff = 2.0 * p_ac.lambda * p_ac.sigma * p_ac.sigma * dt;

    let mut h = g.to_vec();
    for (a, row) in h.iter_mut().enumerate() {
        for (b, cell) in row.iter_mut().enumerate() {
            *cell += risk_coeff * (a.min(b) as f64);
        }
    }
    solve_transient_general(&h, p_ac.x_total)
}

/// Satu baris hasil scan kernel osilasi ``G(t) = lambda_ow e^{-rho t} cos(omega t)`` + ``eps * lambda_ow * I``.
#[derive(Clone, Debug)]
pub struct OscScanRow {
    pub rho: f64,
    pub omega: f64,
    pub eps: f64,
    pub is_pd: bool,
    pub min_eig_rel: f64,
    pub max_eig_rel: f64,
    pub cond: f64,
    pub min_x: f64,
}

/// Scan deterministik grid parameter (dinyatakan sebelum dijalankan di ``main.rs``): kembalikan semua baris
/// berlabel PD; pemanggil menyaring yang memiliki ``min_x < 0`` ("PD tetapi manipulasi").
pub fn scan_oscillating_kernels(
    n: usize,
    t_horizon: f64,
    lambda_ow: f64,
    x_total: f64,
    rhos: &[f64],
    omegas: &[f64],
    epss: &[f64],
) -> Vec<OscScanRow> {
    let mut rows = Vec::new();
    for &rho in rhos {
        for &omega in omegas {
            for &eps in epss {
                let mut g = build_kernel_matrix(n, t_horizon, |d| {
                    lambda_ow * (-rho * d).exp() * (omega * d).cos()
                });
                for (i, row) in g.iter_mut().enumerate() {
                    row[i] += eps * lambda_ow;
                }
                let is_pd = cholesky(&g).is_ok();
                let (lo, hi) = jacobi_extreme_eigenvalues(&g);
                let min_x = if is_pd {
                    solve_transient_general(&g, x_total)
                        .map(|(x, _)| x.iter().copied().fold(f64::INFINITY, f64::min))
                        .unwrap_or(f64::NAN)
                } else {
                    f64::NAN
                };
                rows.push(OscScanRow {
                    rho,
                    omega,
                    eps,
                    is_pd,
                    min_eig_rel: lo / lambda_ow,
                    max_eig_rel: hi / lambda_ow,
                    cond: hi / lo,
                    min_x,
                });
            }
        }
    }
    rows
}
