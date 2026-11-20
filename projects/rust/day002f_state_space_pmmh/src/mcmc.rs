use rayon::prelude::*;

use crate::diagnostics::{summarize_chains, ParamDiagnostics};
use crate::rng::{cholesky_lower, derive_seed, make_rng, sample_normal, sample_uniform};

pub struct McmcConfig {
    pub n_chains: usize,
    pub n_iter: usize,
    pub burn_in: usize,
    pub master_seed: u64,
    pub initial_step: Vec<f64>,
}

pub struct McmcResult {
    pub chains: Vec<Vec<Vec<f64>>>,
    pub acceptance_rates: Vec<f64>,
    pub mean_acceptance: f64,
    pub param_summaries: Vec<ParamDiagnostics>,
}

/// Online Haario covariance tracker used ONLY during burn-in, then frozen.
struct OnlineCov {
    dim: usize,
    count: usize,
    mean: Vec<f64>,
    m2: Vec<Vec<f64>>,
}

impl OnlineCov {
    fn new(dim: usize) -> Self {
        Self {
            dim,
            count: 0,
            mean: vec![0.0; dim],
            m2: vec![vec![0.0; dim]; dim],
        }
    }

    fn push(&mut self, x: &[f64]) {
        self.count += 1;
        let c = self.count as f64;
        let mut delta = vec![0.0; self.dim];
        for i in 0..self.dim {
            delta[i] = x[i] - self.mean[i];
            self.mean[i] += delta[i] / c;
        }
        for (i, row_m2) in self.m2.iter_mut().enumerate() {
            for (j, cell) in row_m2.iter_mut().enumerate() {
                *cell += delta[i] * (x[j] - self.mean[j]);
            }
        }
    }

    fn proposal_chol(&self, eps_reg: f64) -> Option<Vec<Vec<f64>>> {
        if self.count < self.dim + 10 {
            return None;
        }
        let scale = (2.38 * 2.38) / (self.dim as f64);
        let denom = (self.count - 1) as f64;
        let mut cov = vec![vec![0.0; self.dim]; self.dim];
        for (i, cov_row) in cov.iter_mut().enumerate() {
            for (j, cell) in cov_row.iter_mut().enumerate() {
                *cell = scale * (self.m2[i][j] / denom);
            }
            cov_row[i] += eps_reg;
        }
        cholesky_lower(&cov).ok()
    }
}

/// Run parallel RW-MH or PMMH across `config.n_chains` chains.
/// `log_target(theta, rng)` returns `Some(log_prior + log_lik)` or `None` if out of support.
pub fn run_mcmc_parallel<F>(
    config: &McmcConfig,
    init_points: &[Vec<f64>],
    log_target: F,
) -> Result<McmcResult, String>
where
    F: Fn(&[f64], &mut rand_chacha::ChaCha8Rng) -> Option<f64> + Sync + Send,
{
    let dim = config.initial_step.len();
    if init_points.len() != config.n_chains {
        return Err("init_points length must match n_chains".to_string());
    }

    let chain_outputs: Vec<(Vec<Vec<f64>>, f64)> = (0..config.n_chains)
        .into_par_iter()
        .map(|cid| {
            let mut rng = make_rng(derive_seed(config.master_seed, cid as u64));
            let mut curr = init_points[cid].clone();
            let mut curr_lp = log_target(&curr, &mut rng).unwrap_or(f64::NEG_INFINITY);

            let mut tracker = OnlineCov::new(dim);
            let mut chol: Vec<Vec<f64>> = (0..dim)
                .map(|i| {
                    let mut row = vec![0.0; dim];
                    row[i] = config.initial_step[i];
                    row
                })
                .collect();

            let mut kept = Vec::with_capacity(config.n_iter - config.burn_in);
            let mut accepted_post = 0usize;
            let mut z = vec![0.0; dim];
            let mut prop = vec![0.0; dim];

            for it in 0..config.n_iter {
                for zv in z.iter_mut() {
                    *zv = sample_normal(&mut rng);
                }
                for i in 0..dim {
                    let mut step = 0.0;
                    for j in 0..=i {
                        step += chol[i][j] * z[j];
                    }
                    prop[i] = curr[i] + step;
                }

                if let Some(prop_lp) = log_target(&prop, &mut rng) {
                    let log_acc = prop_lp - curr_lp;
                    if log_acc >= 0.0 || sample_uniform(&mut rng).ln() < log_acc {
                        curr.copy_from_slice(&prop);
                        curr_lp = prop_lp;
                        if it >= config.burn_in {
                            accepted_post += 1;
                        }
                    }
                }

                if it < config.burn_in {
                    tracker.push(&curr);
                    if it >= 100 && (it + 1) % 25 == 0 {
                        if let Some(new_chol) = tracker.proposal_chol(1e-6) {
                            chol = new_chol;
                        }
                    }
                } else {
                    kept.push(curr.clone());
                }
            }

            let acc_rate = (accepted_post as f64) / ((config.n_iter - config.burn_in) as f64);
            (kept, acc_rate)
        })
        .collect();

    let mut chains_by_param =
        vec![vec![Vec::with_capacity(config.n_iter - config.burn_in); config.n_chains]; dim];
    let mut acc_rates = Vec::with_capacity(config.n_chains);
    for (cid, (samples, acc)) in chain_outputs.iter().enumerate() {
        acc_rates.push(*acc);
        for s in samples {
            for p in 0..dim {
                chains_by_param[p][cid].push(s[p]);
            }
        }
    }

    let param_summaries: Vec<ParamDiagnostics> = chains_by_param
        .iter()
        .map(|c| summarize_chains(c))
        .collect();
    let mean_acceptance = acc_rates.iter().sum::<f64>() / (config.n_chains as f64);
    let chains = chain_outputs.into_iter().map(|(c, _)| c).collect();

    Ok(McmcResult {
        chains,
        acceptance_rates: acc_rates,
        mean_acceptance,
        param_summaries,
    })
}
