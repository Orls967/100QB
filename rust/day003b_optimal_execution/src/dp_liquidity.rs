//! Bagian 3: Likuiditas stokastik dengan Dynamic Programming kuadratik eksak.

use crate::ac::{solve_ac_qp, AcParams};
use crate::parallel::parallel_fixed_chunks;
use crate::rng::Xoshiro256PlusPlus;
use crate::stats::{mean, mean_se, variance};

/// Parameter model DP likuiditas stokastik 2-rezim.
#[derive(Clone, Copy, Debug)]
pub struct LiquidityDpParams {
    pub x_total: f64,
    pub t_horizon: f64,
    pub n_steps: usize,
    pub sigma: f64,
    pub lambda: f64,
    pub eta_liquid: f64,
    pub eta_illiquid: f64,
    pub p_00: f64, // P(likuid -> likuid)
    pub p_11: f64, // P(ilikuid -> ilikuid)
}

impl Default for LiquidityDpParams {
    fn default() -> Self {
        Self {
            x_total: 1_000_000.0,
            t_horizon: 1.0,
            n_steps: 50,
            sigma: 0.95,
            lambda: 2e-5,
            eta_liquid: 1.5e-6,
            eta_illiquid: 6.0e-6,
            p_00: 0.9,
            p_11: 0.7,
        }
    }
}

/// Solusi rekursi mundur DP kuadratik: ``a[k][s]`` dan ``c[k][s]`` untuk k = 0..N-1.
#[derive(Clone, Debug)]
pub struct DpSolution {
    pub a: Vec<[f64; 2]>,
    pub c: Vec<[f64; 2]>,
    pub v0_liquid: f64,
    pub v0_illiquid: f64,
}

/// Turunkan rekursi mundur eksak ``V_k(y, s) = a_k(s) * y^2`` dan ``x_k = c_k(s) * y_{k-1}``.
pub fn solve_liquidity_dp(p: &LiquidityDpParams) -> DpSolution {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let etas = [p.eta_liquid / dt, p.eta_illiquid / dt];
    let risk = p.lambda * p.sigma * p.sigma * dt;
    let trans = [[p.p_00, 1.0 - p.p_00], [1.0 - p.p_11, p.p_11]];

    let mut a = vec![[0.0; 2]; n];
    let mut c = vec![[0.0; 2]; n];

    // Langkah terminal k = N-1: c_{N-1}(s) = 1, a_{N-1}(s) = eta_s / dt
    for s in 0..2 {
        c[n - 1][s] = 1.0;
        a[n - 1][s] = etas[s];
    }

    for k in (0..n - 1).rev() {
        for s in 0..2 {
            let cont = trans[s][0] * a[k + 1][0] + trans[s][1] * a[k + 1][1];
            let big_a = risk + cont;
            let ck = big_a / (etas[s] + big_a);
            c[k][s] = ck;
            a[k][s] = etas[s] * ck;
        }
    }

    let x2 = p.x_total * p.x_total;
    DpSolution {
        v0_liquid: a[0][0] * x2,
        v0_illiquid: a[0][1] * x2,
        a,
        c,
    }
}

/// Verifikasi independen brute-force dengan golden-section search pada grid y halus untuk N kecil (N=5).
pub fn verify_dp_bruteforce(p: &LiquidityDpParams, n_grid: usize) -> f64 {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let etas = [p.eta_liquid / dt, p.eta_illiquid / dt];
    let risk = p.lambda * p.sigma * p.sigma * dt;
    let trans = [[p.p_00, 1.0 - p.p_00], [1.0 - p.p_11, p.p_11]];

    let y_grid: Vec<f64> = (0..n_grid)
        .map(|i| p.x_total * (i as f64) / ((n_grid - 1) as f64))
        .collect();
    let dy = p.x_total / ((n_grid - 1) as f64);

    // Fungsi interpolasi kuadratik lokal pada grid seragam y_grid
    let interp_v = |vals: &[f64], y: f64| -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        let pos = (y / dy).min((n_grid - 1) as f64);
        let idx = (pos.floor() as usize).min(n_grid - 3);
        let y0 = y_grid[idx];
        let u = (y - y0) / dy;
        let f0 = vals[idx];
        let f1 = vals[idx + 1];
        let f2 = vals[idx + 2];
        0.5 * (u - 1.0) * (u - 2.0) * f0 - u * (u - 2.0) * f1 + 0.5 * u * (u - 1.0) * f2
    };

    let mut v_next = [vec![0.0; n_grid], vec![0.0; n_grid]];
    for (i, &y) in y_grid.iter().enumerate() {
        v_next[0][i] = etas[0] * y * y;
        v_next[1][i] = etas[1] * y * y;
    }

    let phi = 0.5 * (5.0f64.sqrt() - 1.0);
    for _k in (0..n - 1).rev() {
        let mut v_curr = [vec![0.0; n_grid], vec![0.0; n_grid]];
        for s in 0..2 {
            for (i, &y) in y_grid.iter().enumerate() {
                if i == 0 {
                    continue;
                }
                let eval_obj = |x: f64| -> f64 {
                    let rem = y - x;
                    let cont = trans[s][0] * interp_v(&v_next[0], rem)
                        + trans[s][1] * interp_v(&v_next[1], rem);
                    etas[s] * x * x + risk * rem * rem + cont
                };
                // Golden-section search pada [0, y]
                let mut lo = 0.0;
                let mut hi = y;
                let mut x1 = hi - phi * (hi - lo);
                let mut x2 = lo + phi * (hi - lo);
                let mut f1 = eval_obj(x1);
                let mut f2 = eval_obj(x2);
                for _ in 0..60 {
                    if f1 < f2 {
                        hi = x2;
                        x2 = x1;
                        f2 = f1;
                        x1 = hi - phi * (hi - lo);
                        f1 = eval_obj(x1);
                    } else {
                        lo = x1;
                        x1 = x2;
                        f1 = f2;
                        x2 = lo + phi * (hi - lo);
                        f2 = eval_obj(x2);
                    }
                }
                v_curr[s][i] = f1.min(f2);
            }
        }
        v_next = v_curr;
    }

    let dp_exact = solve_liquidity_dp(p);
    (v_next[0][n_grid - 1] - dp_exact.v0_liquid).abs() / dp_exact.v0_liquid
}

/// Hasil simulasi Monte Carlo DP likuiditas stokastik.
#[derive(Clone, Debug)]
pub struct DpSimResult {
    pub mean_j: f64,
    pub se_j: f64,
    pub v0_analytic: f64,
    pub z_score: f64,
    pub illiquid_frac: f64,
    pub y_k10_mean: f64,
    pub y_k10_std: f64,
    pub y_k25_mean: f64,
    pub y_k25_std: f64,
    pub y_k40_mean: f64,
    pub y_k40_std: f64,
    pub j_samples_head: Vec<f64>,
}

/// Simulasi 400.000 jalur rezim Markov secara paralel deterministik.
pub fn simulate_dp_paths(
    p: &LiquidityDpParams,
    sol: &DpSolution,
    n_paths: usize,
    n_threads: usize,
    seed: u64,
) -> DpSimResult {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let etas = [p.eta_liquid / dt, p.eta_illiquid / dt];
    let risk = p.lambda * p.sigma * p.sigma * dt;
    let stay_probs = [p.p_00, p.p_11];

    let chunks = parallel_fixed_chunks(n_paths, 100, n_threads, seed, &|_idx, count, cseed| {
        let mut rng = Xoshiro256PlusPlus::from_seed(cseed);
        let mut js = Vec::with_capacity(count);
        let mut y10 = Vec::with_capacity(count);
        let mut y25 = Vec::with_capacity(count);
        let mut y40 = Vec::with_capacity(count);
        let mut illiq_steps = 0usize;

        for _ in 0..count {
            let mut s = 0usize; // rezim awal: likuid
            let mut y = p.x_total;
            let mut total_j = 0.0;

            for k in 0..n {
                if s == 1 {
                    illiq_steps += 1;
                }
                let xk = sol.c[k][s] * y;
                y -= xk;
                total_j += etas[s] * xk * xk;
                if k + 1 < n {
                    total_j += risk * y * y;
                }
                if k == 9 {
                    y10.push(y);
                } else if k == 24 {
                    y25.push(y);
                } else if k == 39 {
                    y40.push(y);
                }
                let u = rng.next_f64();
                s = if u < stay_probs[s] { s } else { 1 - s };
            }
            js.push(total_j);
        }
        (js, illiq_steps, y10, y25, y40)
    });

    let mut all_j = Vec::with_capacity(n_paths);
    let mut all_y10 = Vec::with_capacity(n_paths);
    let mut all_y25 = Vec::with_capacity(n_paths);
    let mut all_y40 = Vec::with_capacity(n_paths);
    let mut total_illiq = 0usize;

    for (js, illiq, y10, y25, y40) in chunks {
        all_j.extend(js);
        total_illiq += illiq;
        all_y10.extend(y10);
        all_y25.extend(y25);
        all_y40.extend(y40);
    }

    let mean_j = mean(&all_j);
    let se_j = mean_se(&all_j);
    let z_score = (mean_j - sol.v0_liquid).abs() / se_j.max(1e-12);
    let illiquid_frac = (total_illiq as f64) / ((n_paths * n) as f64);

    DpSimResult {
        mean_j,
        se_j,
        v0_analytic: sol.v0_liquid,
        z_score,
        illiquid_frac,
        y_k10_mean: mean(&all_y10),
        y_k10_std: variance(&all_y10).sqrt(),
        y_k25_mean: mean(&all_y25),
        y_k25_std: variance(&all_y25).sqrt(),
        y_k40_mean: mean(&all_y40),
        y_k40_std: variance(&all_y40).sqrt(),
        j_samples_head: all_j.into_iter().take(2000).collect(),
    }
}

/// Evaluasi ekspektasi biaya sejati ``E[J(x)]`` untuk jadwal deterministik ``x`` di bawah rantai Markov rezim.
pub fn evaluate_static_schedule_under_markov(p: &LiquidityDpParams, x: &[f64]) -> f64 {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let risk = p.lambda * p.sigma * p.sigma * dt;
    let trans = [[p.p_00, 1.0 - p.p_00], [1.0 - p.p_11, p.p_11]];

    let mut pi = [1.0, 0.0]; // mulai di rezim likuid
    let mut y = p.x_total;
    let mut total = 0.0;
    for (k, &xk) in x.iter().enumerate() {
        let eta_k = pi[0] * p.eta_liquid + pi[1] * p.eta_illiquid;
        total += (eta_k / dt) * xk * xk;
        y -= xk;
        if k + 1 < n {
            total += risk * y * y;
        }
        let next0 = pi[0] * trans[0][0] + pi[1] * trans[1][0];
        let next1 = pi[0] * trans[0][1] + pi[1] * trans[1][1];
        pi = [next0, next1];
    }
    total
}

/// Hitung nilai adaptivitas DP vs strategi statis (stasioner, likuid-optimal, ilikuid-optimal).
#[derive(Clone, Debug)]
pub struct AdaptivityReport {
    pub j_dp: f64,
    pub j_static_stationary: f64,
    pub j_static_liquid: f64,
    pub j_static_illiquid: f64,
    pub adaptivity_pct_vs_stationary: f64,
    pub adaptivity_pct_vs_liquid: f64,
    pub adaptivity_pct_vs_illiquid: f64,
}

pub fn compute_adaptivity_value(p: &LiquidityDpParams) -> AdaptivityReport {
    let dp = solve_liquidity_dp(p);
    let p01 = 1.0 - p.p_00;
    let p10 = 1.0 - p.p_11;
    let pi_inf_0 = p10 / (p01 + p10);
    let pi_inf_1 = p01 / (p01 + p10);
    let eta_stat = pi_inf_0 * p.eta_liquid + pi_inf_1 * p.eta_illiquid;

    let make_static = |eta_val: f64| -> Vec<f64> {
        let ac_p = AcParams {
            x_total: p.x_total,
            t_horizon: p.t_horizon,
            n_steps: p.n_steps,
            s0: 100.0,
            sigma: p.sigma,
            eta: eta_val,
            gamma: 0.0,
            lambda: p.lambda,
        };
        solve_ac_qp(&ac_p).x
    };

    let j_stat = evaluate_static_schedule_under_markov(p, &make_static(eta_stat));
    let j_liq = evaluate_static_schedule_under_markov(p, &make_static(p.eta_liquid));
    let j_ill = evaluate_static_schedule_under_markov(p, &make_static(p.eta_illiquid));

    AdaptivityReport {
        j_dp: dp.v0_liquid,
        j_static_stationary: j_stat,
        j_static_liquid: j_liq,
        j_static_illiquid: j_ill,
        adaptivity_pct_vs_stationary: 100.0 * (j_stat - dp.v0_liquid) / j_stat,
        adaptivity_pct_vs_liquid: 100.0 * (j_liq - dp.v0_liquid) / j_liq,
        adaptivity_pct_vs_illiquid: 100.0 * (j_ill - dp.v0_liquid) / j_ill,
    }
}
