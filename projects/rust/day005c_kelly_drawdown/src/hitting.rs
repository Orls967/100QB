use crate::gbm::GbmParams;
use crate::parallel::reduce_vec_fixed_blocks;
use crate::rng::{derive_seed, Rng};
use crate::special::norm_cdf;

#[derive(Clone, Debug)]
pub struct HittingGridRow {
    pub dt: f64,
    pub t: f64,
    pub x: f64,
    pub f: f64,
    pub exact_prob: f64,
    pub naive_mc: f64,
    pub naive_se: f64,
    pub naive_bias: f64,
    pub bridge_mc: f64,
    pub bridge_se: f64,
    pub bridge_z: f64,
}

#[derive(Clone, Debug)]
pub struct InfiniteHitRow {
    pub c: f64,
    pub f: f64,
    pub g: f64,
    pub s: f64,
    pub two_g_over_s2: f64,
    pub prob_hit_half: f64,
}

/// Infinite-horizon hitting probability `P(inf_{t>=0} W_t <= x * W_0)` for `x in (0, 1)`.
pub fn prob_hit_infinite(g: f64, s: f64, x: f64) -> f64 {
    if x >= 1.0 {
        return 1.0;
    }
    if x <= 0.0 || s <= 0.0 {
        return if g <= 0.0 { 1.0 } else { 0.0 };
    }
    if g <= 0.0 {
        return 1.0;
    }
    let exponent = 2.0 * g / (s * s);
    x.powf(exponent).clamp(0.0, 1.0)
}

/// Finite-horizon hitting probability `P(min_{0<=t<=T} W_t <= x * W_0)` for `x in (0, 1)`
/// derived via reflection principle + Girsanov.
pub fn prob_hit_finite(g: f64, s: f64, t: f64, x: f64) -> f64 {
    if x >= 1.0 {
        return 1.0;
    }
    if x <= 0.0 || t <= 0.0 || s <= 0.0 {
        return 0.0;
    }
    let a = x.ln();
    let st = s * t.sqrt();
    let t1 = norm_cdf((a - g * t) / st);
    let exp_arg = (2.0 * g * a) / (s * s);
    let t2 = if exp_arg < -700.0 {
        0.0
    } else {
        exp_arg.exp() * norm_cdf((a + g * t) / st)
    };
    (t1 + t2).clamp(0.0, 1.0)
}

#[inline]
pub fn bridge_cross_prob(y_i: f64, y_next: f64, a: f64, s2_dt: f64) -> f64 {
    if y_i <= a || y_next <= a {
        1.0
    } else {
        (-2.0 * (y_i - a) * (y_next - a) / s2_dt).exp()
    }
}

/// Shared-path Monte Carlo across `(T, x, f)` for each `dt` in `dts`.
/// Uses Rao-Blackwellized Brownian bridge survival `1 - prod(1 - p_cross)` for unbiased
/// minimum-variance estimation alongside naive discrete monitoring.
#[allow(clippy::too_many_arguments)]
pub fn run_hitting_grid(
    gbm: &GbmParams,
    dts: &[f64],
    ts: &[f64],
    xs: &[f64],
    fs: &[f64],
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> Vec<HittingGridRow> {
    let mut rows = Vec::with_capacity(dts.len() * ts.len() * xs.len() * fs.len());
    let max_t = ts.iter().copied().fold(0.0_f64, f64::max);

    for (d_idx, &dt) in dts.iter().enumerate() {
        let steps_per_year = (1.0 / dt).round() as usize;
        let total_steps = (max_t * steps_per_year as f64).round() as usize;
        let t_steps: Vec<usize> = ts
            .iter()
            .map(|&t| (t * steps_per_year as f64).round() as usize)
            .collect();
        let sqrt_dt = dt.sqrt();

        // Number of (T, x, f) combinations = n_combo; we store 4 accumulators per combo:
        // [sum_naive, sum_naive_sq, sum_bridge, sum_bridge_sq]
        let n_t = ts.len();
        let n_x = xs.len();
        let n_f = fs.len();
        let n_combo = n_t * n_x * n_f;
        let n_acc = n_combo * 4;

        let gs: Vec<f64> = fs.iter().map(|&f| gbm.drift_log(f)).collect();
        let ss: Vec<f64> = fs.iter().map(|&f| gbm.vol_log(f)).collect();
        let s2_dts: Vec<f64> = ss.iter().map(|&s| s * s * dt).collect();
        let a_xs: Vec<f64> = xs.iter().map(|&x| x.ln()).collect();

        let sums = reduce_vec_fixed_blocks(r_paths, n_workers, n_acc, |idx, acc| {
            let mut rng = Rng::new(derive_seed(main_seed, 200 + d_idx as u64, idx as u64));
            let mut y = vec![0.0_f64; n_f];
            let mut naive_hit = vec![false; n_f * n_x];
            let mut surv_prob = vec![1.0_f64; n_f * n_x];
            let mut t_ptr = 0_usize;

            for step in 1..=total_steps {
                let z = rng.normal();
                let dw = sqrt_dt * z;
                for fi in 0..n_f {
                    let y_prev = y[fi];
                    let y_next = y_prev + gs[fi] * dt + ss[fi] * dw;
                    y[fi] = y_next;
                    let s2_dt = s2_dts[fi];
                    for (xi, &a) in a_xs.iter().enumerate() {
                        let slot = fi * n_x + xi;
                        if y_next <= a {
                            naive_hit[slot] = true;
                            surv_prob[slot] = 0.0;
                        } else if surv_prob[slot] > 0.0 && y_prev > a {
                            let pc = (-2.0 * (y_prev - a) * (y_next - a) / s2_dt).exp();
                            surv_prob[slot] *= 1.0 - pc;
                        }
                    }
                }
                while t_ptr < n_t && step == t_steps[t_ptr] {
                    for fi in 0..n_f {
                        for xi in 0..n_x {
                            let slot = fi * n_x + xi;
                            let nv = if naive_hit[slot] { 1.0 } else { 0.0 };
                            let br = 1.0 - surv_prob[slot];
                            let c_idx = (t_ptr * n_x * n_f + xi * n_f + fi) * 4;
                            acc[c_idx] += nv;
                            acc[c_idx + 1] += nv * nv;
                            acc[c_idx + 2] += br;
                            acc[c_idx + 3] += br * br;
                        }
                    }
                    t_ptr += 1;
                }
            }
        });

        let rf = r_paths as f64;
        for (ti, &t) in ts.iter().enumerate() {
            for (xi, &x) in xs.iter().enumerate() {
                for (fi, &f) in fs.iter().enumerate() {
                    let c_idx = (ti * n_x * n_f + xi * n_f + fi) * 4;
                    let naive_mc = sums[c_idx] / rf;
                    let naive_var = ((sums[c_idx + 1] / rf) - naive_mc * naive_mc).max(0.0);
                    let naive_se = (naive_var / rf).sqrt();

                    let bridge_mc = sums[c_idx + 2] / rf;
                    let bridge_var = ((sums[c_idx + 3] / rf) - bridge_mc * bridge_mc).max(0.0);
                    let bridge_se = (bridge_var / rf).sqrt();

                    let exact_prob = prob_hit_finite(gs[fi], ss[fi], t, x);
                    let bridge_z = if bridge_se > 1e-15 {
                        (bridge_mc - exact_prob) / bridge_se
                    } else {
                        0.0
                    };
                    rows.push(HittingGridRow {
                        dt,
                        t,
                        x,
                        f,
                        exact_prob,
                        naive_mc,
                        naive_se,
                        naive_bias: naive_mc - exact_prob,
                        bridge_mc,
                        bridge_se,
                        bridge_z,
                    });
                }
            }
        }
    }
    rows
}

pub fn infinite_horizon_table(gbm: &GbmParams, cs: &[f64], x: f64) -> Vec<InfiniteHitRow> {
    let f_star = gbm.kelly_fraction();
    cs.iter()
        .map(|&c| {
            let f = c * f_star;
            let g = gbm.drift_log(f);
            let s = gbm.vol_log(f);
            let two_g_over_s2 = if s > 0.0 { 2.0 * g / (s * s) } else { 0.0 };
            let prob_hit_half = prob_hit_infinite(g, s, x);
            InfiniteHitRow {
                c,
                f,
                g,
                s,
                two_g_over_s2,
                prob_hit_half,
            }
        })
        .collect()
}
