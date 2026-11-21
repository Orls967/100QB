//! PCA of weekly yield changes: theoretical covariance, Jacobi eigendecomposition,
//! T=5000 x R=100 Monte Carlo study, and aggregate/diagnostic acceptance evaluation.

use crate::afns::simulate_stationary_path;
use crate::kalman::ModelSpec;
use crate::linalg;
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;

pub fn theoretical_cov(spec: &ModelSpec, taus: &[f64], dt: f64) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let m = taus.len();
    let mut h_mat = linalg::zeros(m, 3);
    for (i, &tau) in taus.iter().enumerate() {
        let row = spec.afns.yield_loadings(tau);
        h_mat[i][0] = row[0];
        h_mat[i][1] = row[1];
        h_mat[i][2] = row[2];
    }
    let phi = spec.afns.transition_matrix_closed(dt);
    let q_dt = spec.afns.transition_cov_closed(dt);
    let q_inf = spec.afns.stationary_cov();

    let mut phi_minus_i = phi;
    for i in 0..3 {
        phi_minus_i[i][i] -= 1.0;
    }
    let mut tmp = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            tmp[i][j] = phi_minus_i[i][0] * q_inf[0][j]
                + phi_minus_i[i][1] * q_inf[1][j]
                + phi_minus_i[i][2] * q_inf[2][j];
        }
    }
    let mut sigma_dx = linalg::zeros(3, 3);
    let mut sigma_dx_approx = linalg::zeros(3, 3);
    for i in 0..3 {
        for j in 0..3 {
            sigma_dx_approx[i][j] = q_dt[i][j];
            sigma_dx[i][j] = q_dt[i][j]
                + tmp[i][0] * phi_minus_i[j][0]
                + tmp[i][1] * phi_minus_i[j][1]
                + tmp[i][2] * phi_minus_i[j][2];
        }
    }
    let h_t = linalg::transpose(&h_mat);
    let mut exact_cov = linalg::mat_mul(&h_mat, &linalg::mat_mul(&sigma_dx, &h_t));
    let mut approx_cov = linalg::mat_mul(&h_mat, &linalg::mat_mul(&sigma_dx_approx, &h_t));
    let noise_var = 2.0 * spec.sigma_eps * spec.sigma_eps;
    for i in 0..m {
        exact_cov[i][i] += noise_var;
        approx_cov[i][i] += noise_var;
    }
    (exact_cov, approx_cov)
}

pub fn count_sign_changes(v: &[f64]) -> usize {
    let mut cnt = 0usize;
    let mut prev_sign = 0i32;
    for &x in v {
        let s = if x > 1e-12 {
            1
        } else if x < -1e-12 {
            -1
        } else {
            0
        };
        if s == 0 {
            continue;
        }
        if prev_sign != 0 && s != prev_sign {
            cnt += 1;
        }
        prev_sign = s;
    }
    cnt
}

fn orient_eigenvectors(vecs: &mut [Vec<f64>]) {
    if !vecs.is_empty() && stats::mean(&vecs[0]) < 0.0 {
        for x in &mut vecs[0] {
            *x = -*x;
        }
    }
    if vecs.len() >= 2 && vecs[1].len() >= 2 && vecs[1][vecs[1].len() - 1] < vecs[1][0] {
        for x in &mut vecs[1] {
            *x = -*x;
        }
    }
    if vecs.len() >= 3 && !vecs[2].is_empty() && vecs[2][0] < 0.0 {
        for x in &mut vecs[2] {
            *x = -*x;
        }
    }
}

fn align_vector_sign(v: &[f64], reference: &[f64]) -> Vec<f64> {
    if linalg::dot(v, reference) < 0.0 {
        v.iter().map(|&x| -x).collect()
    } else {
        v.to_vec()
    }
}

pub fn sample_yield_diff_cov(yields: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let t_steps = yields.len();
    let m = yields[0].len();
    let n_diff = t_steps - 1;
    let mut means = vec![0.0; m];
    for t in 1..t_steps {
        for j in 0..m {
            means[j] += yields[t][j] - yields[t - 1][j];
        }
    }
    for mj in &mut means {
        *mj /= n_diff as f64;
    }
    let mut cov = linalg::zeros(m, m);
    for t in 1..t_steps {
        for i in 0..m {
            let di = (yields[t][i] - yields[t - 1][i]) - means[i];
            for j in i..m {
                let dj = (yields[t][j] - yields[t - 1][j]) - means[j];
                cov[i][j] += di * dj;
            }
        }
    }
    let denom = (n_diff - 1) as f64;
    for i in 0..m {
        for j in i..m {
            let val = cov[i][j] / denom;
            cov[i][j] = val;
            cov[j][i] = val;
        }
    }
    cov
}

#[derive(Clone, Debug)]
pub struct PcaReplicationSummary {
    pub tenors: Vec<f64>,
    pub t_steps: usize,
    pub num_reps: usize,
    pub theo_eigenvalues: Vec<f64>,
    pub approx_eigenvalues: Vec<f64>,
    pub theo_var_explained: Vec<f64>,
    pub theo_top3_cum_var: f64,
    pub theo_eigenvectors: Vec<Vec<f64>>,
    pub theo_sign_changes: [usize; 3],
    pub emp_mean_eigenvalues: Vec<f64>,
    pub emp_sd_eigenvalues: Vec<f64>,
    pub emp_mean_se: Vec<f64>,
    pub z_score_sample_mean: Vec<f64>,
    pub single_run_eigenvals: Vec<f64>,
    pub single_run_z_scores: Vec<f64>,
    pub noise_bulk_theo_mean: f64,
    pub noise_bulk_emp_mean: f64,
    pub noise_bulk_mean_se: f64,
    pub noise_bulk_z_score: f64,
    pub single_run_top3_cos: [f64; 3],
    pub emp_mean_top3_cos: [f64; 3],
    pub emp_min_top3_cos: [f64; 3],
    pub emp_mean_top3_vecs: Vec<Vec<f64>>,
    pub emp_sign_changes: [usize; 3],
    pub emp_var_explained: Vec<f64>,
    pub emp_top3_cum_var: f64,
}

#[derive(Clone, Debug)]
pub struct PcaAcceptance {
    pub top3_eigenvalues_within_3se: bool,
    pub noise_bulk_mean_within_3se: bool,
    pub diagnostic_all10_sorted_within_3se: bool,
    pub top3_cosine_gt_0p99: bool,
    pub sign_patterns_012: bool,
    pub aggregate_passed: bool,
}

pub fn run_pca_study(
    spec: &ModelSpec,
    taus: &[f64],
    dt: f64,
    t_steps: usize,
    n_reps: usize,
    n_workers: usize,
    seed: u64,
) -> PcaReplicationSummary {
    let m = taus.len();
    let (exact_cov, approx_cov) = theoretical_cov(spec, taus, dt);
    let (theo_vals, mut theo_vecs, _) =
        linalg::symmetric_eigen_jacobi(&exact_cov, 80, 1e-14).unwrap();
    let (approx_vals, _, _) = linalg::symmetric_eigen_jacobi(&approx_cov, 80, 1e-14).unwrap();
    orient_eigenvectors(&mut theo_vecs);

    let theo_sum: f64 = theo_vals.iter().sum();
    let theo_var_exp: Vec<f64> = theo_vals.iter().map(|&v| 100.0 * v / theo_sum).collect();
    let theo_top3_cum = theo_var_exp.iter().take(3).sum();
    let theo_signs = [
        count_sign_changes(&theo_vecs[0]),
        count_sign_changes(&theo_vecs[1]),
        count_sign_changes(&theo_vecs[2]),
    ];

    struct RepOut {
        vals: Vec<f64>,
        vecs3: [Vec<f64>; 3],
        cos3: [f64; 3],
        noise_mean: f64,
    }

    let rows = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x4A01, r as u64));
        let (_, y) =
            simulate_stationary_path(&spec.afns, taus, t_steps, dt, spec.sigma_eps, &mut rg);
        let cov = sample_yield_diff_cov(&y);
        let (vals, vecs, _) = linalg::symmetric_eigen_jacobi(&cov, 80, 1e-14).unwrap();
        let mut vecs3 = [Vec::new(), Vec::new(), Vec::new()];
        let mut cos3 = [0.0; 3];
        for k in 0..3 {
            let aligned = align_vector_sign(&vecs[k], &theo_vecs[k]);
            cos3[k] = linalg::dot(&aligned, &theo_vecs[k]).abs();
            vecs3[k] = aligned;
        }
        let noise_mean = if m > 3 { stats::mean(&vals[3..]) } else { 0.0 };
        RepOut {
            vals,
            vecs3,
            cos3,
            noise_mean,
        }
    });

    let mut emp_mean = vec![0.0; m];
    let mut emp_sd = vec![0.0; m];
    let mut emp_se = vec![0.0; m];
    let mut z_mean = vec![0.0; m];
    let mut single_z = vec![0.0; m];
    for j in 0..m {
        let col: Vec<f64> = rows.iter().map(|r| r.vals[j]).collect();
        emp_mean[j] = stats::mean(&col);
        emp_sd[j] = stats::sd(&col);
        emp_se[j] = stats::se(&col);
        z_mean[j] = (emp_mean[j] - theo_vals[j]).abs() / emp_se[j].max(1e-18);
        single_z[j] = (rows[0].vals[j] - theo_vals[j]).abs() / emp_sd[j].max(1e-18);
    }

    let emp_sum: f64 = emp_mean.iter().sum();
    let emp_var_exp: Vec<f64> = emp_mean.iter().map(|&v| 100.0 * v / emp_sum).collect();
    let emp_top3_cum = emp_var_exp.iter().take(3).sum();

    let noise_theo = if m > 3 {
        stats::mean(&theo_vals[3..])
    } else {
        0.0
    };
    let noise_cols: Vec<f64> = rows.iter().map(|r| r.noise_mean).collect();
    let noise_emp = stats::mean(&noise_cols);
    let noise_se = stats::se(&noise_cols);
    let noise_z = (noise_emp - noise_theo).abs() / noise_se.max(1e-18);

    let mut mean_top3_cos = [0.0; 3];
    let mut min_top3_cos = [1.0; 3];
    let mut mean_top3_vecs = vec![vec![0.0; m]; 3];
    for k in 0..3 {
        let cs: Vec<f64> = rows.iter().map(|r| r.cos3[k]).collect();
        mean_top3_cos[k] = stats::mean(&cs);
        min_top3_cos[k] = cs.iter().cloned().fold(f64::INFINITY, f64::min);
        for j in 0..m {
            mean_top3_vecs[k][j] = rows.iter().map(|r| r.vecs3[k][j]).sum::<f64>() / n_reps as f64;
        }
        let norm = linalg::dot(&mean_top3_vecs[k], &mean_top3_vecs[k]).sqrt();
        if norm > 0.0 {
            for x in &mut mean_top3_vecs[k] {
                *x /= norm;
            }
        }
    }
    let emp_signs = [
        count_sign_changes(&mean_top3_vecs[0]),
        count_sign_changes(&mean_top3_vecs[1]),
        count_sign_changes(&mean_top3_vecs[2]),
    ];

    PcaReplicationSummary {
        tenors: taus.to_vec(),
        t_steps,
        num_reps: n_reps,
        theo_eigenvalues: theo_vals,
        approx_eigenvalues: approx_vals,
        theo_var_explained: theo_var_exp,
        theo_top3_cum_var: theo_top3_cum,
        theo_eigenvectors: theo_vecs,
        theo_sign_changes: theo_signs,
        emp_mean_eigenvalues: emp_mean,
        emp_sd_eigenvalues: emp_sd,
        emp_mean_se: emp_se,
        z_score_sample_mean: z_mean,
        single_run_eigenvals: rows[0].vals.clone(),
        single_run_z_scores: single_z,
        noise_bulk_theo_mean: noise_theo,
        noise_bulk_emp_mean: noise_emp,
        noise_bulk_mean_se: noise_se,
        noise_bulk_z_score: noise_z,
        single_run_top3_cos: rows[0].cos3,
        emp_mean_top3_cos: mean_top3_cos,
        emp_min_top3_cos: min_top3_cos,
        emp_mean_top3_vecs: mean_top3_vecs,
        emp_sign_changes: emp_signs,
        emp_var_explained: emp_var_exp,
        emp_top3_cum_var: emp_top3_cum,
    }
}

pub fn evaluate_acceptance(sum: &PcaReplicationSummary) -> PcaAcceptance {
    let top3_eig = (0..3).all(|i| sum.z_score_sample_mean[i] <= 3.0);
    let noise_bulk = sum.noise_bulk_z_score <= 3.0;
    let all10 = sum.z_score_sample_mean.iter().all(|&z| z <= 3.0);
    let top3_cos = sum.emp_mean_top3_cos.iter().all(|&c| c > 0.99);
    let signs = sum.emp_sign_changes == [0, 1, 2];
    PcaAcceptance {
        top3_eigenvalues_within_3se: top3_eig,
        noise_bulk_mean_within_3se: noise_bulk,
        diagnostic_all10_sorted_within_3se: all10,
        top3_cosine_gt_0p99: top3_cos,
        sign_patterns_012: signs,
        aggregate_passed: top3_eig && noise_bulk && top3_cos && signs,
    }
}
