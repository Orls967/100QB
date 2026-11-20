use serde::Serialize;

use crate::kalman::kalman_log_likelihood;
use crate::models::model_a_log_prior;

#[derive(Debug, Clone, Serialize)]
pub struct MarginalSummary {
    pub mean: f64,
    pub sd: f64,
    pub q025: f64,
    pub q500: f64,
    pub q975: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct GridPosterior {
    pub a_grid: Vec<f64>,
    pub b_grid: Vec<f64>,
    pub marginal_a: Vec<f64>,
    pub marginal_b: Vec<f64>,
    pub summary_a: MarginalSummary,
    pub summary_b: MarginalSummary,
    pub edge_mass: f64,
}

fn trapezoidal_weights(n: usize, step: f64) -> Vec<f64> {
    let mut w = vec![step; n];
    w[0] *= 0.5;
    w[n - 1] *= 0.5;
    w
}

fn summarize_1d(grid: &[f64], pdf: &[f64], step: f64) -> MarginalSummary {
    let w = trapezoidal_weights(grid.len(), step);
    let norm: f64 = pdf.iter().zip(&w).map(|(p, wi)| p * wi).sum();
    let mean: f64 = grid
        .iter()
        .zip(pdf)
        .zip(&w)
        .map(|((&x, &p), &wi)| x * p * wi)
        .sum::<f64>()
        / norm;
    let var: f64 = grid
        .iter()
        .zip(pdf)
        .zip(&w)
        .map(|((&x, &p), &wi)| (x - mean) * (x - mean) * p * wi)
        .sum::<f64>()
        / norm;

    let mut cdf = Vec::with_capacity(grid.len());
    let mut acc = 0.0;
    cdf.push(0.0);
    for i in 1..grid.len() {
        acc += 0.5 * (pdf[i - 1] + pdf[i]) * step / norm;
        cdf.push(acc);
    }

    let interp_q = |q: f64| -> f64 {
        for i in 1..cdf.len() {
            if cdf[i] >= q {
                let frac = (q - cdf[i - 1]) / (cdf[i] - cdf[i - 1]).max(1e-15);
                return grid[i - 1] + frac * (grid[i] - grid[i - 1]);
            }
        }
        *grid.last().unwrap_or(&0.0)
    };

    MarginalSummary {
        mean,
        sd: var.sqrt(),
        q025: interp_q(0.025),
        q500: interp_q(0.50),
        q975: interp_q(0.975),
    }
}

/// 2D exact grid posterior for Model A over `a in [a_min, a_max]`, `b in [b_min, b_max]`.
pub fn compute_model_a_grid(
    obs: &[f64],
    a_range: (f64, f64),
    b_range: (f64, f64),
    n_pts: usize,
) -> GridPosterior {
    let da = (a_range.1 - a_range.0) / ((n_pts - 1) as f64);
    let db = (b_range.1 - b_range.0) / ((n_pts - 1) as f64);
    let a_grid: Vec<f64> = (0..n_pts).map(|i| a_range.0 + (i as f64) * da).collect();
    let b_grid: Vec<f64> = (0..n_pts).map(|j| b_range.0 + (j as f64) * db).collect();

    let mut log_post = vec![vec![0.0; n_pts]; n_pts];
    let mut max_lp = f64::NEG_INFINITY;
    for (i, &a) in a_grid.iter().enumerate() {
        let s_eta = a.exp();
        for (j, &b) in b_grid.iter().enumerate() {
            let lp = kalman_log_likelihood(obs, s_eta, b.exp()) + model_a_log_prior(a, b);
            log_post[i][j] = lp;
            if lp > max_lp {
                max_lp = lp;
            }
        }
    }

    let wa = trapezoidal_weights(n_pts, da);
    let wb = trapezoidal_weights(n_pts, db);
    let mut total_mass = 0.0;
    for i in 0..n_pts {
        for j in 0..n_pts {
            let val = (log_post[i][j] - max_lp).exp();
            log_post[i][j] = val;
            total_mass += val * wa[i] * wb[j];
        }
    }
    let inv_mass = 1.0 / total_mass;
    for row in &mut log_post {
        for v in row.iter_mut() {
            *v *= inv_mass;
        }
    }

    let mut edge_mass = 0.0;
    for i in 0..n_pts {
        for j in 0..n_pts {
            if i == 0 || i == n_pts - 1 || j == 0 || j == n_pts - 1 {
                edge_mass += log_post[i][j] * wa[i] * wb[j];
            }
        }
    }

    let marginal_a: Vec<f64> = (0..n_pts)
        .map(|i| (0..n_pts).map(|j| log_post[i][j] * wb[j]).sum())
        .collect();
    let marginal_b: Vec<f64> = (0..n_pts)
        .map(|j| (0..n_pts).map(|i| log_post[i][j] * wa[i]).sum())
        .collect();

    let summary_a = summarize_1d(&a_grid, &marginal_a, da);
    let summary_b = summarize_1d(&b_grid, &marginal_b, db);

    GridPosterior {
        a_grid,
        b_grid,
        marginal_a,
        marginal_b,
        summary_a,
        summary_b,
        edge_mass,
    }
}
