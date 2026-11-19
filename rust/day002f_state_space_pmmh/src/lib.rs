pub mod diagnostics;
pub mod grid;
pub mod kalman;
pub mod mcmc;
pub mod models;
pub mod pf;
pub mod plots;
pub mod rng;

#[cfg(test)]
mod tests {
    use super::diagnostics::{effective_sample_size, split_r_hat};
    use super::grid::compute_model_a_grid;
    use super::kalman::{brute_force_log_likelihood, kalman_log_likelihood};
    use super::mcmc::{run_mcmc_parallel, McmcConfig};
    use super::models::{
        model_a_log_prior, simulate_model_a, MODEL_A_SEED, MODEL_A_SIGMA_EPS, MODEL_A_SIGMA_ETA,
    };
    use super::pf::{pf_model_a, systematic_resample, PfWorkspace};
    use super::rng::{derive_seed, make_rng, sample_normal};

    #[test]
    fn kalman_matches_brute_force_on_t10() {
        let data = simulate_model_a(10, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS, 42);
        for &(se, sp) in &[(0.5, 1.0), (0.2, 0.8), (1.3, 0.4)] {
            let kf = kalman_log_likelihood(&data.obs, se, sp);
            let bf = brute_force_log_likelihood(&data.obs, se, sp).expect("SPD covariance");
            assert!((kf - bf).abs() < 1e-10, "kf={kf} bf={bf}");
        }
    }

    #[test]
    fn systematic_resampling_is_unbiased() {
        let weights = [0.05, 0.15, 0.50, 0.20, 0.10];
        let n_p = weights.len();
        let mut cdf = vec![0.0; n_p];
        let mut anc = vec![0usize; n_p];
        let mut counts = vec![0.0; n_p];
        let reps = 10_000usize;
        let mut rng = make_rng(77);

        for _ in 0..reps {
            systematic_resample(&weights, &mut cdf, &mut anc, &mut rng);
            for &a in &anc {
                counts[a] += 1.0;
            }
        }
        for (i, &w) in weights.iter().enumerate() {
            let emp_mean = counts[i] / (reps as f64);
            let expected = (n_p as f64) * w;
            assert!(
                (emp_mean - expected).abs() < 0.03,
                "i={i} emp={emp_mean} exp={expected}"
            );
        }
    }

    #[test]
    fn pf_estimator_unbiased_and_grid_normalized() {
        let data = simulate_model_a(60, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS, MODEL_A_SEED);
        let ll_exact = kalman_log_likelihood(&data.obs, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS);
        let mut ws = PfWorkspace::new(200);
        let mut rng = make_rng(123);
        let r_reps = 400usize;
        let mut ratios = Vec::with_capacity(r_reps);
        for _ in 0..r_reps {
            let ll_hat = pf_model_a(
                &data.obs,
                MODEL_A_SIGMA_ETA,
                MODEL_A_SIGMA_EPS,
                &mut ws,
                &mut rng,
            );
            ratios.push((ll_hat - ll_exact).exp());
        }
        let mean_r = ratios.iter().sum::<f64>() / (r_reps as f64);
        let var_r = ratios
            .iter()
            .map(|&r| (r - mean_r) * (r - mean_r))
            .sum::<f64>()
            / ((r_reps - 1) as f64);
        let se_r = (var_r / (r_reps as f64)).sqrt();
        assert!(
            (mean_r - 1.0).abs() < 3.0 * se_r,
            "mean_r={mean_r}, se_r={se_r}"
        );

        let data200 = simulate_model_a(200, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS, MODEL_A_SEED);
        let grid = compute_model_a_grid(&data200.obs, (-3.0, 1.5), (-1.5, 1.5), 100);
        assert!(grid.edge_mass < 1e-8, "edge_mass={}", grid.edge_mass);
    }

    #[test]
    fn diagnostics_ar1_ess_and_iid_rhat() {
        let mut rng = make_rng(999);
        let rho: f64 = 0.9;
        let sd = (1.0 - rho * rho).sqrt();
        let n = 10_000usize;
        let mut chains = Vec::new();
        for _ in 0..4 {
            let mut ch = Vec::with_capacity(n);
            let mut x = sample_normal(&mut rng);
            for _ in 0..n {
                x = rho * x + sd * sample_normal(&mut rng);
                ch.push(x);
            }
            chains.push(ch);
        }
        let ess = effective_sample_size(&chains);
        let theo = ((4 * n) as f64) * (1.0 - rho) / (1.0 + rho);
        assert!((ess - theo).abs() / theo < 0.15, "ess={ess}, theo={theo}");

        let iid: Vec<Vec<f64>> = (0..4)
            .map(|_| (0..2000).map(|_| sample_normal(&mut rng)).collect())
            .collect();
        let rhat = split_r_hat(&iid);
        assert!((rhat - 1.0).abs() < 0.02, "rhat={rhat}");
    }

    #[test]
    fn mcmc_reproducible_and_thread_invariant() {
        let data = simulate_model_a(50, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS, 42);
        let cfg = McmcConfig {
            n_chains: 4,
            n_iter: 400,
            burn_in: 100,
            master_seed: derive_seed(42, 1),
            initial_step: vec![0.1, 0.1],
        };
        let inits = vec![vec![-0.5, 0.0]; 4];
        let run1 = run_mcmc_parallel(&cfg, &inits, |th, _rng| {
            Some(
                kalman_log_likelihood(&data.obs, th[0].exp(), th[1].exp())
                    + model_a_log_prior(th[0], th[1]),
            )
        })
        .expect("mcmc");

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("pool");
        let run2 = pool
            .install(|| {
                run_mcmc_parallel(&cfg, &inits, |th, _rng| {
                    Some(
                        kalman_log_likelihood(&data.obs, th[0].exp(), th[1].exp())
                            + model_a_log_prior(th[0], th[1]),
                    )
                })
            })
            .expect("mcmc 1 thread");

        assert_eq!(run1.chains, run2.chains);
    }
}
