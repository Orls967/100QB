//! Bagian 1: Eksekusi optimal Almgren-Chriss (QP KKT, sinh diskrit, limit kontinu, frontier, MC).

use crate::linalg::{dot, solve_spd};
use crate::parallel::parallel_fixed_chunks;
use crate::rng::Xoshiro256PlusPlus;
use crate::stats::{mean, mean_se, variance, variance_se_batch_means};

/// Parameter model Almgren-Chriss.
#[derive(Clone, Copy, Debug)]
pub struct AcParams {
    pub x_total: f64,
    pub t_horizon: f64,
    pub n_steps: usize,
    pub s0: f64,
    pub sigma: f64,
    pub eta: f64,
    pub gamma: f64,
    pub lambda: f64,
}

impl Default for AcParams {
    fn default() -> Self {
        Self {
            x_total: 1_000_000.0,
            t_horizon: 1.0,
            n_steps: 50,
            s0: 100.0,
            sigma: 0.95,
            eta: 2.5e-6,
            gamma: 2.5e-7,
            lambda: 2e-5,
        }
    }
}

/// Hasil jadwal eksekusi Almgren-Chriss.
#[derive(Clone, Debug)]
pub struct AcSchedule {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub expected_c: f64,
    pub var_c: f64,
    pub std_c: f64,
}

/// Hitung E[C] dan Var[C] langsung dari jadwal ``x`` (panjang N).
pub fn evaluate_ac_moments(p: &AcParams, x: &[f64]) -> (f64, f64) {
    let n = p.n_steps;
    assert_eq!(x.len(), n);
    let dt = p.t_horizon / (n as f64);
    let eta_tilde = p.eta - 0.5 * p.gamma * dt;

    let sum_x2: f64 = x.iter().map(|&v| v * v).sum();
    let expected_c = 0.5 * p.gamma * p.x_total * p.x_total + (eta_tilde / dt) * sum_x2;

    let mut y_rem = p.x_total;
    let mut sum_y2 = 0.0;
    for &xk in x.iter().take(n - 1) {
        y_rem -= xk;
        sum_y2 += y_rem * y_rem;
    }
    let var_c = p.sigma * p.sigma * dt * sum_y2;
    (expected_c, var_c)
}

/// Bangun sisa holding ``y[0..=N]`` dari jadwal trade ``x[0..N]``.
pub fn holdings_from_trades(x_total: f64, x: &[f64]) -> Vec<f64> {
    let mut y = Vec::with_capacity(x.len() + 1);
    let mut rem = x_total;
    y.push(rem);
    for &xk in x {
        rem -= xk;
        y.push(rem);
    }
    y
}

/// Solusi QP eksak lewat sistem KKT + Cholesky buatan sendiri.
pub fn solve_ac_qp(p: &AcParams) -> AcSchedule {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let eta_tilde = p.eta - 0.5 * p.gamma * dt;
    assert!(eta_tilde > 0.0, "eta_tilde harus positif");

    let diag_term = 2.0 * eta_tilde / dt;
    let risk_coeff = 2.0 * p.lambda * p.sigma * p.sigma * dt;

    let mut h = vec![vec![0.0; n]; n];
    for (a, row) in h.iter_mut().enumerate() {
        for (b, cell) in row.iter_mut().enumerate() {
            let min_ab = a.min(b) as f64; // 0-indexed: min(a+1, b+1) - 1 == a.min(b)
            *cell = risk_coeff * min_ab;
            if a == b {
                *cell += diag_term;
            }
        }
    }

    let ones = vec![1.0; n];
    let h_inv_1 = solve_spd(&h, &ones).expect("Matriks H Almgren-Chriss harus SPD");
    let denom = dot(&ones, &h_inv_1);
    let x: Vec<f64> = h_inv_1.into_iter().map(|v| p.x_total * v / denom).collect();
    let y = holdings_from_trades(p.x_total, &x);
    let (expected_c, var_c) = evaluate_ac_moments(p, &x);
    AcSchedule {
        x,
        y,
        expected_c,
        var_c,
        std_c: var_c.sqrt(),
    }
}

/// Solusi tertutup sinh diskrit dari persamaan Euler-Lagrange diskrit.
pub fn solve_ac_sinh_discrete(p: &AcParams) -> AcSchedule {
    let n = p.n_steps;
    let dt = p.t_horizon / (n as f64);
    let eta_tilde = p.eta - 0.5 * p.gamma * dt;

    if p.lambda <= 0.0 {
        let x = vec![p.x_total / (n as f64); n];
        let y = holdings_from_trades(p.x_total, &x);
        let (expected_c, var_c) = evaluate_ac_moments(p, &x);
        return AcSchedule {
            x,
            y,
            expected_c,
            var_c,
            std_c: var_c.sqrt(),
        };
    }

    let kappa_disc = (p.lambda * p.sigma * p.sigma / eta_tilde).sqrt();
    let arg = 0.5 * kappa_disc * dt;
    let kappa_tilde = (2.0 / dt) * arg.asinh();

    let mut y = Vec::with_capacity(n + 1);
    let denom = (kappa_tilde * p.t_horizon).sinh();
    for j in 0..=n {
        let t_rem = p.t_horizon - (j as f64) * dt;
        let y_j = if denom.is_infinite() {
            if j == 0 {
                p.x_total
            } else {
                0.0
            }
        } else {
            p.x_total * (kappa_tilde * t_rem).sinh() / denom
        };
        y.push(y_j);
    }
    y[0] = p.x_total;
    y[n] = 0.0;

    let mut x = Vec::with_capacity(n);
    for j in 0..n {
        x.push(y[j] - y[j + 1]);
    }
    let (expected_c, var_c) = evaluate_ac_moments(p, &x);
    AcSchedule {
        x,
        y,
        expected_c,
        var_c,
        std_c: var_c.sqrt(),
    }
}

/// Limit kontinu ``y(t) = X * sinh(kappa * (T - t)) / sinh(kappa * T)`` dengan ``kappa^2 = lambda * sigma^2 / eta``.
pub fn ac_continuous_holding(p: &AcParams, t: f64) -> f64 {
    if p.lambda <= 0.0 {
        return p.x_total * (1.0 - t / p.t_horizon);
    }
    let kappa = (p.lambda * p.sigma * p.sigma / p.eta).sqrt();
    let denom = (kappa * p.t_horizon).sinh();
    if denom.is_infinite() {
        if t <= 0.0 {
            p.x_total
        } else {
            0.0
        }
    } else {
        p.x_total * (kappa * (p.t_horizon - t)).sinh() / denom
    }
}

/// Galat sup relatif ``max_j |y_j - y_cont(t_j)| / X`` untuk ``N`` tertentu.
pub fn ac_continuous_sup_error(p: &AcParams, n_steps: usize) -> f64 {
    let mut p_n = *p;
    p_n.n_steps = n_steps;
    let sched = solve_ac_sinh_discrete(&p_n);
    let dt = p.t_horizon / (n_steps as f64);
    let mut max_err = 0.0f64;
    for (j, &yj) in sched.y.iter().enumerate() {
        let yc = ac_continuous_holding(&p_n, (j as f64) * dt);
        max_err = max_err.max((yj - yc).abs() / p.x_total);
    }
    max_err
}

/// Bangkitkan frontier efisien (40 titik lambda pada skala log [1e-8, 1e-3]).
pub fn ac_efficient_frontier(p: &AcParams, n_points: usize) -> Vec<(f64, AcSchedule)> {
    let log_min = (1e-8f64).ln();
    let log_max = (1e-3f64).ln();
    (0..n_points)
        .map(|i| {
            let frac = (i as f64) / ((n_points - 1).max(1) as f64);
            let lam = (log_min + frac * (log_max - log_min)).exp();
            let mut pi = *p;
            pi.lambda = lam;
            (lam, solve_ac_qp(&pi))
        })
        .collect()
}

/// Hasil simulasi Monte Carlo harga Almgren-Chriss (antitetik, paralel deterministik).
#[derive(Clone, Debug)]
pub struct AcMcResult {
    pub mc_mean: f64,
    pub mc_mean_se: f64,
    pub mc_antithetic_mean: f64,
    pub mc_var: f64,
    pub mc_var_se: f64,
    pub exact_mean: f64,
    pub exact_var: f64,
    pub z_mean: f64,
    pub z_var: f64,
}

/// Simulasi Monte Carlo antitetik untuk jadwal ``x`` yang diberikan.
pub fn simulate_ac_mc(
    p: &AcParams,
    x: &[f64],
    n_paths: usize,
    n_threads: usize,
    seed: u64,
) -> AcMcResult {
    assert_eq!(n_paths % 2, 0, "n_paths harus genap untuk antitetik");
    let n_pairs = n_paths / 2;
    let n = x.len();
    let dt = p.t_horizon / (n as f64);
    let sqrt_dt = dt.sqrt();

    let (exact_mean, exact_var) = evaluate_ac_moments(p, x);

    // Simulasikan pasangan antitetik dalam 100 chunk tetap
    let chunk_vecs = parallel_fixed_chunks(n_pairs, 100, n_threads, seed, &|_idx, count, cseed| {
        let mut rng = Xoshiro256PlusPlus::from_seed(cseed);
        let mut pair_means = Vec::with_capacity(count);
        let mut raw_samples = Vec::with_capacity(count * 2);
        for _ in 0..count {
            let mut s_p = p.s0;
            let mut s_m = p.s0;
            let mut rev_p = 0.0;
            let mut rev_m = 0.0;
            for &xk in x {
                let p_exec_p = s_p - p.eta * xk / dt;
                let p_exec_m = s_m - p.eta * xk / dt;
                rev_p += xk * p_exec_p;
                rev_m += xk * p_exec_m;

                let z = rng.next_normal();
                let shock = p.sigma * sqrt_dt * z;
                s_p = s_p + shock - p.gamma * xk;
                s_m = s_m - shock - p.gamma * xk;
            }
            let c_p = p.x_total * p.s0 - rev_p;
            let c_m = p.x_total * p.s0 - rev_m;
            pair_means.push(0.5 * (c_p + c_m));
            // Rata-rata variansi per jalur tidak dipasangkan: simpan kedua jalur
            raw_samples.push(c_p);
            raw_samples.push(c_m);
        }
        (pair_means, raw_samples)
    });

    let mut all_pair_means = Vec::with_capacity(n_pairs);
    let mut all_raw = Vec::with_capacity(n_paths);
    for (pm, raw) in chunk_vecs {
        all_pair_means.extend(pm);
        all_raw.extend(raw);
    }

    let mc_antithetic_mean = mean(&all_pair_means);
    let mc_mean = mean(&all_raw);
    let mc_mean_se = mean_se(&all_raw);
    let mc_var = variance(&all_raw);
    let mc_var_se = variance_se_batch_means(&all_raw, 100);
    let z_mean =
        if exact_var == 0.0 && (mc_mean - exact_mean).abs() <= 1e-9 * exact_mean.abs().max(1.0) {
            0.0
        } else {
            (mc_mean - exact_mean).abs() / mc_mean_se.max(1e-12)
        };
    let z_var = if exact_var == 0.0 && mc_var <= 1e-12 {
        0.0
    } else {
        (mc_var - exact_var).abs() / mc_var_se.max(1e-12)
    };

    AcMcResult {
        mc_mean,
        mc_mean_se,
        mc_antithetic_mean,
        mc_var,
        mc_var_se,
        exact_mean,
        exact_var,
        z_mean,
        z_var,
    }
}

/// Hasil MC TANPA antitetik (jalur independen) untuk rata-rata dan variansi shortfall.
#[derive(Clone, Debug)]
pub struct AcPlainMcResult {
    pub mc_mean: f64,
    pub mc_mean_se: f64,
    pub exact_mean: f64,
    pub z_mean: f64,
    pub mc_var: f64,
    pub mc_var_se: f64,
    pub exact_var: f64,
    pub z_var: f64,
}

/// Simulasi MC jalur independen (tanpa antitetik). Karena shortfall linear di xi, mean dengan
/// antitetik identik dengan rumus eksak (uji vakum); di sini z berukuran O(1) dan dilaporkan apa adanya.
pub fn simulate_ac_mc_plain(
    p: &AcParams,
    x: &[f64],
    n_paths: usize,
    n_threads: usize,
    seed: u64,
) -> AcPlainMcResult {
    let n = x.len();
    let dt = p.t_horizon / (n as f64);
    let sqrt_dt = dt.sqrt();
    let (exact_mean, exact_var) = evaluate_ac_moments(p, x);

    let chunk_vecs = parallel_fixed_chunks(n_paths, 100, n_threads, seed, &|_idx, count, cseed| {
        let mut rng = Xoshiro256PlusPlus::from_seed(cseed);
        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            let mut s = p.s0;
            let mut rev = 0.0;
            for &xk in x {
                rev += xk * (s - p.eta * xk / dt);
                s = s + p.sigma * sqrt_dt * rng.next_normal() - p.gamma * xk;
            }
            out.push(p.x_total * p.s0 - rev);
        }
        out
    });
    let all: Vec<f64> = chunk_vecs.into_iter().flatten().collect();
    let mc_mean = mean(&all);
    let mc_mean_se = mean_se(&all);
    let mc_var = variance(&all);
    let mc_var_se = variance_se_batch_means(&all, 100);
    AcPlainMcResult {
        mc_mean,
        mc_mean_se,
        exact_mean,
        z_mean: (mc_mean - exact_mean).abs() / mc_mean_se.max(1e-12),
        mc_var,
        mc_var_se,
        exact_var,
        z_var: if exact_var == 0.0 {
            0.0
        } else {
            (mc_var - exact_var).abs() / mc_var_se.max(1e-12)
        },
    }
}

/// Holding kontinu dengan koefisien impact efektif ``eta_eff``: y(t) = X sinh(k (T-t)) / sinh(k T).
fn continuous_holding_eta(p: &AcParams, eta_eff: f64, t: f64) -> f64 {
    let kappa = (p.lambda * p.sigma * p.sigma / eta_eff).sqrt();
    p.x_total * (kappa * (p.t_horizon - t)).sinh() / (kappa * p.t_horizon).sinh()
}

/// Dekomposisi galat limit kontinu: total = d1 + d2 dengan
/// d1 = y_sinh-diskrit - y_kontinu(eta_tilde)  (murni diskretisasi),
/// d2 = y_kontinu(eta_tilde) - y_kontinu(eta)  (efek eta_tilde vs eta).
/// Mengembalikan (sup|d1|, sup|d2|, sup|total|, sup|d1+d2| terhadap tanda: apakah d1 dan d2 berlawanan tanda di
/// titik galat total maksimum) sebagai (e1, e2, e_total, opposite_sign_at_max).
pub fn ac_error_decomposition(p: &AcParams, n_steps: usize) -> (f64, f64, f64, bool) {
    let mut pn = *p;
    pn.n_steps = n_steps;
    let sched = solve_ac_sinh_discrete(&pn);
    let dt = p.t_horizon / (n_steps as f64);
    let eta_tilde = p.eta - 0.5 * p.gamma * dt;
    let (mut e1, mut e2, mut et) = (0.0f64, 0.0f64, 0.0f64);
    let mut opp = false;
    for (j, &yj) in sched.y.iter().enumerate() {
        let t = (j as f64) * dt;
        let y_et = continuous_holding_eta(&pn, eta_tilde, t);
        let y_e = ac_continuous_holding(&pn, t);
        let d1 = (yj - y_et) / p.x_total;
        let d2 = (y_et - y_e) / p.x_total;
        let tot = d1 + d2;
        e1 = e1.max(d1.abs());
        e2 = e2.max(d2.abs());
        if tot.abs() > et {
            et = tot.abs();
            opp = d1 * d2 < 0.0;
        }
    }
    (e1, e2, et, opp)
}
