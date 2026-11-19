use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use serde::Serialize;

use day002f_state_space_pmmh::diagnostics::{mean_and_var, summarize_chains, ParamDiagnostics};
use day002f_state_space_pmmh::grid::{compute_model_a_grid, MarginalSummary};
use day002f_state_space_pmmh::kalman::{brute_force_log_likelihood, kalman_log_likelihood};
use day002f_state_space_pmmh::mcmc::{run_mcmc_parallel, McmcConfig};
use day002f_state_space_pmmh::models::{
    model_a_log_prior, model_b_log_prior, simulate_model_a, simulate_model_b, MODEL_A_SEED,
    MODEL_A_SIGMA_EPS, MODEL_A_SIGMA_ETA, MODEL_A_T, MODEL_B_MU, MODEL_B_PHI, MODEL_B_SEED,
    MODEL_B_SIGMA, MODEL_B_T,
};
use day002f_state_space_pmmh::pf::{pf_model_a, pf_model_b, pf_model_b_filtered, PfWorkspace};
use day002f_state_space_pmmh::plots::{
    plot_model_a_posteriors, plot_model_b_sv, plot_pf_unbiased_and_var,
};
use day002f_state_space_pmmh::rng::make_rng;

#[derive(Serialize)]
struct PfUnbiasedSummary {
    n_particles: usize,
    r_reps: usize,
    mean_ratio: f64,
    se_ratio: f64,
    z_score: f64,
    np_grid: Vec<usize>,
    var_ll: Vec<f64>,
    bias_ll: Vec<f64>,
    loglog_slope: f64,
    pf_calls_per_sec: f64,
    particle_steps_per_sec: f64,
}

#[derive(Serialize)]
struct McmcSummaryJson {
    acceptance: f64,
    n_particles: Option<usize>,
    params: Vec<ParamDiagnostics>,
}

#[derive(Serialize)]
struct ModelBComparison {
    np_calibrated: usize,
    pmmh_primary: McmcSummaryJson,
    pmmh_double_np: McmcSummaryJson,
    diff_in_combined_mcse: Vec<f64>,
}

#[derive(Serialize)]
struct TimingsSec {
    kalman_and_grid: f64,
    pf_unbiased_and_sweep: f64,
    model_a_mh: f64,
    model_a_pmmh: f64,
    model_b_pmmh: f64,
    plots_and_io: f64,
    total: f64,
}

#[derive(Serialize)]
struct ResultsPayload {
    kalman_vs_bruteforce_max_diff: f64,
    pf_unbiased: PfUnbiasedSummary,
    grid_a: MarginalSummary,
    grid_b: MarginalSummary,
    grid_edge_mass: f64,
    model_a_mh: McmcSummaryJson,
    model_a_pmmh: McmcSummaryJson,
    model_b: ModelBComparison,
    timings_sec: TimingsSec,
}

fn find_repo_root() -> Result<PathBuf, String> {
    let mut dir = std::env::current_dir().map_err(|e| e.to_string())?;
    loop {
        if dir.join("src").join("qb100").is_dir() {
            return Ok(dir);
        }
        if !dir.pop() {
            return Err("Could not find 100QB repo root (src/qb100)".to_string());
        }
    }
}

fn fit_loglog_slope(xs: &[usize], ys: &[f64]) -> f64 {
    let lx: Vec<f64> = xs.iter().map(|&n| (n as f64).ln()).collect();
    let ly: Vec<f64> = ys.iter().map(|&v| v.max(1e-16).ln()).collect();
    let mx = lx.iter().sum::<f64>() / (lx.len() as f64);
    let my = ly.iter().sum::<f64>() / (ly.len() as f64);
    let num: f64 = lx.iter().zip(&ly).map(|(&x, &y)| (x - mx) * (y - my)).sum();
    let den: f64 = lx.iter().map(|&x| (x - mx) * (x - mx)).sum();
    num / den
}

fn calibrate_np_model_a(obs: &[f64], a: f64, b: f64) -> usize {
    let s_eta = a.exp();
    let s_eps = b.exp();
    let mut rng = make_rng(505);
    for &np in &[100usize, 200, 500, 1000, 2000] {
        let mut ws = PfWorkspace::new(np);
        let lls: Vec<f64> = (0..100)
            .map(|_| pf_model_a(obs, s_eta, s_eps, &mut ws, &mut rng))
            .collect();
        let (_, var) = mean_and_var(&lls);
        if var.sqrt() <= 1.5 {
            return np;
        }
    }
    2000
}

fn calibrate_np_model_b(obs: &[f64], mu: f64, phi: f64, sigma: f64) -> usize {
    let mut rng = make_rng(606);
    for &np in &[100usize, 200, 500, 1000, 2000] {
        let mut ws = PfWorkspace::new(np);
        let lls: Vec<f64> = (0..100)
            .map(|_| pf_model_b(obs, mu, phi, sigma, &mut ws, &mut rng))
            .collect();
        let (_, var) = mean_and_var(&lls);
        if var.sqrt() <= 1.5 {
            return np;
        }
    }
    2000
}

fn transform_model_b_chains(
    chains: &[Vec<Vec<f64>>],
) -> (Vec<ParamDiagnostics>, Vec<f64>, Vec<f64>) {
    let n_chains = chains.len();
    let mut mu_ch = vec![Vec::new(); n_chains];
    let mut phi_ch = vec![Vec::new(); n_chains];
    let mut sig_ch = vec![Vec::new(); n_chains];
    let mut phi_flat = Vec::new();
    let mut sig_flat = Vec::new();

    for (cid, ch) in chains.iter().enumerate() {
        for s in ch {
            let mu = s[0];
            let phi = s[1].tanh();
            let sig = s[2].exp();
            mu_ch[cid].push(mu);
            phi_ch[cid].push(phi);
            sig_ch[cid].push(sig);
            phi_flat.push(phi);
            sig_flat.push(sig);
        }
    }
    (
        vec![
            summarize_chains(&mu_ch),
            summarize_chains(&phi_ch),
            summarize_chains(&sig_ch),
        ],
        phi_flat,
        sig_flat,
    )
}

fn run_all() -> Result<(), String> {
    let t_total = Instant::now();
    let repo_root = find_repo_root()?;
    let out_dir = repo_root.join("outputs").join("day002f");
    let fig_dir = repo_root.join("outputs").join("figures");
    fs::create_dir_all(&out_dir).map_err(|e| e.to_string())?;
    fs::create_dir_all(&fig_dir).map_err(|e| e.to_string())?;

    // 1. Kalman vs brute-force (T=10) & Model A 2D Grid
    let t0 = Instant::now();
    let data_a10 = simulate_model_a(10, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS, MODEL_A_SEED);
    let mut max_kf_diff = 0.0_f64;
    for &(se, sp) in &[(0.5, 1.0), (0.2, 0.8), (1.3, 0.4)] {
        let kf = kalman_log_likelihood(&data_a10.obs, se, sp);
        let bf = brute_force_log_likelihood(&data_a10.obs, se, sp)?;
        max_kf_diff = max_kf_diff.max((kf - bf).abs());
    }
    let data_a = simulate_model_a(
        MODEL_A_T,
        MODEL_A_SIGMA_ETA,
        MODEL_A_SIGMA_EPS,
        MODEL_A_SEED,
    );
    let grid = compute_model_a_grid(&data_a.obs, (-3.0, 1.5), (-1.5, 1.5), 200);
    let t_kalman_grid = t0.elapsed().as_secs_f64();

    // 2. PF Unbiasedness (Np=500, R=2000) & Np sweep
    let t1 = Instant::now();
    let ll_exact_a = kalman_log_likelihood(&data_a.obs, MODEL_A_SIGMA_ETA, MODEL_A_SIGMA_EPS);
    let r_reps = 2000usize;
    let mut ws_500 = PfWorkspace::new(500);
    let mut rng_pf = make_rng(777);
    let t_throughput = Instant::now();
    let mut ll_diffs = Vec::with_capacity(r_reps);
    let mut ratios = Vec::with_capacity(r_reps);
    for _ in 0..r_reps {
        let ll_hat = pf_model_a(
            &data_a.obs,
            MODEL_A_SIGMA_ETA,
            MODEL_A_SIGMA_EPS,
            &mut ws_500,
            &mut rng_pf,
        );
        ll_diffs.push(ll_hat - ll_exact_a);
        ratios.push((ll_hat - ll_exact_a).exp());
    }
    let dt_throughput = t_throughput.elapsed().as_secs_f64();
    let pf_calls_per_sec = (r_reps as f64) / dt_throughput;
    let particle_steps_per_sec = ((r_reps * 500 * MODEL_A_T) as f64) / dt_throughput;
    let (mean_ratio, var_ratio) = mean_and_var(&ratios);
    let se_ratio = (var_ratio / (r_reps as f64)).sqrt();
    let z_ratio = (mean_ratio - 1.0).abs() / se_ratio;

    let np_grid = vec![50usize, 100, 200, 500, 1000, 2000];
    let mut var_ll = Vec::with_capacity(np_grid.len());
    let mut bias_ll = Vec::with_capacity(np_grid.len());
    for &np in &np_grid {
        let mut ws = PfWorkspace::new(np);
        let lls: Vec<f64> = (0..500)
            .map(|_| {
                pf_model_a(
                    &data_a.obs,
                    MODEL_A_SIGMA_ETA,
                    MODEL_A_SIGMA_EPS,
                    &mut ws,
                    &mut rng_pf,
                )
            })
            .collect();
        let (m, v) = mean_and_var(&lls);
        var_ll.push(v);
        bias_ll.push(m - ll_exact_a);
    }
    let loglog_slope = fit_loglog_slope(&np_grid, &var_ll);
    let t_pf_sweep = t1.elapsed().as_secs_f64();

    // 3. Model A MCMC: Exact MH vs PMMH (4 chains, 50k iter, 10k burn-in)
    let inits_a = vec![
        vec![grid.summary_a.mean - 0.1, grid.summary_b.mean - 0.05],
        vec![grid.summary_a.mean + 0.1, grid.summary_b.mean + 0.05],
        vec![grid.summary_a.mean - 0.05, grid.summary_b.mean + 0.1],
        vec![grid.summary_a.mean + 0.05, grid.summary_b.mean - 0.1],
    ];
    let cfg_a_mh = McmcConfig {
        n_chains: 4,
        n_iter: 50_000,
        burn_in: 10_000,
        master_seed: 100,
        initial_step: vec![0.15, 0.10],
    };
    let t2 = Instant::now();
    let res_a_mh = run_mcmc_parallel(&cfg_a_mh, &inits_a, |th, _rng| {
        Some(
            kalman_log_likelihood(&data_a.obs, th[0].exp(), th[1].exp())
                + model_a_log_prior(th[0], th[1]),
        )
    })?;
    let t_a_mh = t2.elapsed().as_secs_f64();

    let np_a = calibrate_np_model_a(&data_a.obs, grid.summary_a.mean, grid.summary_b.mean);
    let cfg_a_pmmh = McmcConfig {
        n_chains: 4,
        n_iter: 50_000,
        burn_in: 10_000,
        master_seed: 200,
        initial_step: vec![0.15, 0.10],
    };
    let t3 = Instant::now();
    let res_a_pmmh = run_mcmc_parallel(&cfg_a_pmmh, &inits_a, |th, rng| {
        let mut ws = PfWorkspace::new(np_a);
        let ll = pf_model_a(&data_a.obs, th[0].exp(), th[1].exp(), &mut ws, rng);
        Some(ll + model_a_log_prior(th[0], th[1]))
    })?;
    let t_a_pmmh = t3.elapsed().as_secs_f64();

    // 4. Model B (SV) PMMH (4 chains, 20k iter, 5k burn-in) + exactness check with 2x Np
    let data_b = simulate_model_b(
        MODEL_B_T,
        MODEL_B_MU,
        MODEL_B_PHI,
        MODEL_B_SIGMA,
        MODEL_B_SEED,
    );
    let np_b = calibrate_np_model_b(&data_b.obs, MODEL_B_MU, MODEL_B_PHI, MODEL_B_SIGMA);
    let z_init = MODEL_B_PHI.atanh();
    let s_init = MODEL_B_SIGMA.ln();
    let inits_b = vec![
        vec![MODEL_B_MU - 0.1, z_init - 0.1, s_init - 0.1],
        vec![MODEL_B_MU + 0.1, z_init + 0.1, s_init + 0.1],
        vec![MODEL_B_MU - 0.05, z_init + 0.1, s_init - 0.05],
        vec![MODEL_B_MU + 0.05, z_init - 0.1, s_init + 0.05],
    ];
    let cfg_b1 = McmcConfig {
        n_chains: 4,
        n_iter: 20_000,
        burn_in: 5_000,
        master_seed: 300,
        initial_step: vec![0.15, 0.20, 0.15],
    };
    let t4 = Instant::now();
    let res_b1 = run_mcmc_parallel(&cfg_b1, &inits_b, |th, rng| {
        let phi = th[1].tanh();
        if phi.abs() >= 0.999 {
            return None;
        }
        let mut ws = PfWorkspace::new(np_b);
        let ll = pf_model_b(&data_b.obs, th[0], phi, th[2].exp(), &mut ws, rng);
        Some(ll + model_b_log_prior(th[0], th[1], th[2]))
    })?;
    let np_b2 = (np_b * 2).min(2000);
    let cfg_b2 = McmcConfig {
        n_chains: 4,
        n_iter: 20_000,
        burn_in: 5_000,
        master_seed: 400,
        initial_step: vec![0.15, 0.20, 0.15],
    };
    let res_b2 = run_mcmc_parallel(&cfg_b2, &inits_b, |th, rng| {
        let phi = th[1].tanh();
        if phi.abs() >= 0.999 {
            return None;
        }
        let mut ws = PfWorkspace::new(np_b2);
        let ll = pf_model_b(&data_b.obs, th[0], phi, th[2].exp(), &mut ws, rng);
        Some(ll + model_b_log_prior(th[0], th[1], th[2]))
    })?;
    let t_b_pmmh = t4.elapsed().as_secs_f64();

    let (b1_params, phi_flat, sig_flat) = transform_model_b_chains(&res_b1.chains);
    let (b2_params, _, _) = transform_model_b_chains(&res_b2.chains);
    let diff_in_mcse: Vec<f64> = b1_params
        .iter()
        .zip(&b2_params)
        .map(|(p1, p2)| (p1.mean - p2.mean).abs() / (p1.mcse * p1.mcse + p2.mcse * p2.mcse).sqrt())
        .collect();

    // 5. Plots & JSON output
    let t5 = Instant::now();
    let mh_a: Vec<f64> = res_a_mh
        .chains
        .iter()
        .flat_map(|c| c.iter().map(|s| s[0]))
        .collect();
    let mh_b: Vec<f64> = res_a_mh
        .chains
        .iter()
        .flat_map(|c| c.iter().map(|s| s[1]))
        .collect();
    let pmmh_a: Vec<f64> = res_a_pmmh
        .chains
        .iter()
        .flat_map(|c| c.iter().map(|s| s[0]))
        .collect();
    let pmmh_b: Vec<f64> = res_a_pmmh
        .chains
        .iter()
        .flat_map(|c| c.iter().map(|s| s[1]))
        .collect();

    plot_pf_unbiased_and_var(
        &fig_dir.join("day002f_pf-unbiased-and-var.png"),
        &ll_diffs,
        &np_grid,
        &var_ll,
        loglog_slope,
    )?;
    plot_model_a_posteriors(
        &fig_dir.join("day002f_model_a_posteriors.png"),
        &grid,
        &mh_a,
        &mh_b,
        &pmmh_a,
        &pmmh_b,
    )?;

    let mut rng_sv = make_rng(888);
    let sv_filt = pf_model_b_filtered(
        &data_b.obs,
        b1_params[0].mean,
        b1_params[1].mean,
        b1_params[2].mean,
        1000,
        &mut rng_sv,
    );
    plot_model_b_sv(
        &fig_dir.join("day002f_model_b_sv.png"),
        &data_b.states,
        &sv_filt.mean_h,
        &sv_filt.q025_h,
        &sv_filt.q975_h,
        &phi_flat,
        &sig_flat,
    )?;

    let t_plots = t5.elapsed().as_secs_f64();
    let payload = ResultsPayload {
        kalman_vs_bruteforce_max_diff: max_kf_diff,
        pf_unbiased: PfUnbiasedSummary {
            n_particles: 500,
            r_reps,
            mean_ratio,
            se_ratio,
            z_score: z_ratio,
            np_grid,
            var_ll,
            bias_ll,
            loglog_slope,
            pf_calls_per_sec,
            particle_steps_per_sec,
        },
        grid_a: grid.summary_a,
        grid_b: grid.summary_b,
        grid_edge_mass: grid.edge_mass,
        model_a_mh: McmcSummaryJson {
            acceptance: res_a_mh.mean_acceptance,
            n_particles: None,
            params: res_a_mh.param_summaries,
        },
        model_a_pmmh: McmcSummaryJson {
            acceptance: res_a_pmmh.mean_acceptance,
            n_particles: Some(np_a),
            params: res_a_pmmh.param_summaries,
        },
        model_b: ModelBComparison {
            np_calibrated: np_b,
            pmmh_primary: McmcSummaryJson {
                acceptance: res_b1.mean_acceptance,
                n_particles: Some(np_b),
                params: b1_params,
            },
            pmmh_double_np: McmcSummaryJson {
                acceptance: res_b2.mean_acceptance,
                n_particles: Some(np_b2),
                params: b2_params,
            },
            diff_in_combined_mcse: diff_in_mcse,
        },
        timings_sec: TimingsSec {
            kalman_and_grid: t_kalman_grid,
            pf_unbiased_and_sweep: t_pf_sweep,
            model_a_mh: t_a_mh,
            model_a_pmmh: t_a_pmmh,
            model_b_pmmh: t_b_pmmh,
            plots_and_io: t_plots,
            total: t_total.elapsed().as_secs_f64(),
        },
    };

    let json_path = out_dir.join("results.json");
    let json_str = serde_json::to_string_pretty(&payload).map_err(|e| e.to_string())?;
    fs::write(&json_path, format!("{json_str}\n")).map_err(|e| e.to_string())?;
    println!(
        "Done in {:.2}s -> {}",
        payload.timings_sec.total,
        Path::new(&json_path).display()
    );
    Ok(())
}

fn main() {
    if let Err(err) = run_all() {
        eprintln!("Error: {err}");
        std::process::exit(1);
    }
}
