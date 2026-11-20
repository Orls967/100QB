use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Distribution, StandardNormal};
use rayon::prelude::*;
use serde::Serialize;
use std::f64::consts::FRAC_1_SQRT_2;

const CHUNK_SIZE: usize = 4_000;
const CI_Z_95: f64 = 1.96;

#[derive(Debug, Clone, Copy, Serialize)]
pub struct OptionParams {
    pub s0: f64,
    pub k: f64,
    pub r: f64,
    pub sigma: f64,
    pub t: f64,
    pub n_steps: usize,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EstimatorResult {
    pub price: f64,
    pub se: f64,
    pub ci_low: f64,
    pub ci_high: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct ControlVariateResult {
    pub price: f64,
    pub se: f64,
    pub ci_low: f64,
    pub ci_high: f64,
    pub beta: f64,
    pub y_exact: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct DeltaResult {
    pub delta: f64,
    pub se: f64,
    pub ci_low: f64,
    pub ci_high: f64,
}

/// High-precision erf(x) using W. J. Cody's rational Chebyshev approximation (relative error < 1e-15).
pub fn erf_approx(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let abs_x = x.abs();
    let res = if abs_x < 0.46875 {
        const A: [f64; 5] = [
            3.161_123_743_870_565_5,
            113.864_154_151_050_16,
            377.485_237_685_302,
            3_209.377_589_138_469_4,
            0.185_777_706_184_603_15,
        ];
        const B: [f64; 4] = [
            23.601_290_952_344_12,
            244.024_637_934_444_17,
            1_282.616_526_077_372_3,
            2_844.236_833_439_171,
        ];
        let x2 = abs_x * abs_x;
        let num = (((A[4] * x2 + A[0]) * x2 + A[1]) * x2 + A[2]) * x2 + A[3];
        let den = (((x2 + B[0]) * x2 + B[1]) * x2 + B[2]) * x2 + B[3];
        abs_x * (num / den)
    } else if abs_x < 4.0 {
        const C: [f64; 9] = [
            0.564_188_496_988_670_1,
            8.883_149_794_388_376,
            66.119_190_637_141_63,
            298.635_138_197_400_1,
            881.952_221_241_769,
            1_712.047_612_634_070_6,
            2_051.078_377_826_071_5,
            1_230.339_354_797_997_2,
            2.153_115_354_744_038_5e-8,
        ];
        const D: [f64; 8] = [
            15.744_926_110_709_835,
            117.693_950_891_312_5,
            537.181_101_862_009_9,
            1_621.389_574_566_690_2,
            3_290.799_235_733_459_6,
            4_362.619_090_143_247,
            3_439.367_674_143_721_6,
            1_230.339_354_803_749_4,
        ];
        let mut num = C[8];
        let mut den = 1.0;
        for i in 0..8 {
            num = num * abs_x + C[i];
            den = den * abs_x + D[i];
        }
        let erfc_val = (-abs_x * abs_x).exp() * (num / den);
        1.0 - erfc_val
    } else if abs_x < 26.5 {
        const P: [f64; 6] = [
            0.305_326_634_961_232_34,
            0.360_344_899_949_804_43,
            0.125_781_726_111_229_24,
            0.016_083_785_148_742_276,
            6.587_491_615_298_378e-4,
            0.016_315_387_137_302_097,
        ];
        const Q: [f64; 5] = [
            2.568_520_192_289_822,
            1.872_952_849_923_460_4,
            0.527_905_102_951_428_4,
            0.060_518_341_312_441_32,
            0.002_335_204_976_268_691_8,
        ];
        let inv_x2 = 1.0 / (abs_x * abs_x);
        let mut num = P[5];
        let mut den = 1.0;
        for i in 0..5 {
            num = num * inv_x2 + P[i];
            den = den * inv_x2 + Q[i];
        }
        let erfc_val =
            ((-abs_x * abs_x).exp() / abs_x) * (0.564_189_583_547_756_3 - inv_x2 * (num / den));
        1.0 - erfc_val
    } else {
        1.0
    };
    if x < 0.0 {
        -res
    } else {
        res
    }
}

#[inline]
pub fn norm_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf_approx(x * FRAC_1_SQRT_2))
}

pub fn bs_call(p: &OptionParams) -> f64 {
    if p.t <= 0.0 || p.sigma <= 0.0 {
        return (p.s0 - p.k * (-p.r * p.t).exp()).max(0.0);
    }
    let vol_sqrt_t = p.sigma * p.t.sqrt();
    let d1 = ((p.s0 / p.k).ln() + (p.r + 0.5 * p.sigma * p.sigma) * p.t) / vol_sqrt_t;
    let d2 = d1 - vol_sqrt_t;
    p.s0 * norm_cdf(d1) - p.k * (-p.r * p.t).exp() * norm_cdf(d2)
}

pub fn bs_put(p: &OptionParams) -> f64 {
    if p.t <= 0.0 || p.sigma <= 0.0 {
        return (p.k * (-p.r * p.t).exp() - p.s0).max(0.0);
    }
    let vol_sqrt_t = p.sigma * p.t.sqrt();
    let d1 = ((p.s0 / p.k).ln() + (p.r + 0.5 * p.sigma * p.sigma) * p.t) / vol_sqrt_t;
    let d2 = d1 - vol_sqrt_t;
    p.k * (-p.r * p.t).exp() * norm_cdf(-d2) - p.s0 * norm_cdf(-d1)
}

pub fn geometric_asian_call(p: &OptionParams) -> Result<f64, String> {
    if p.n_steps < 1 {
        return Err("n_steps must be >= 1".to_string());
    }
    let n = p.n_steps as f64;
    let mu_ln_g = p.s0.ln() + (p.r - 0.5 * p.sigma * p.sigma) * p.t * (n + 1.0) / (2.0 * n);
    let var_ln_g = (p.sigma * p.sigma * p.t) * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    let sigma_ln_g = var_ln_g.sqrt();
    let expected_g = (mu_ln_g + 0.5 * var_ln_g).exp();
    let d1 = (mu_ln_g - p.k.ln() + var_ln_g) / sigma_ln_g;
    let d2 = d1 - sigma_ln_g;
    Ok((-p.r * p.t).exp() * (expected_g * norm_cdf(d1) - p.k * norm_cdf(d2)))
}

#[inline]
fn chunk_rng(seed: u64, chunk_idx: usize) -> ChaCha8Rng {
    let derived = seed.wrapping_add((chunk_idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    ChaCha8Rng::seed_from_u64(derived)
}

#[inline]
fn simulate_one_path_avgs(
    rng: &mut ChaCha8Rng,
    s0: f64,
    drift: f64,
    vol: f64,
    n_steps: usize,
) -> (f64, f64, f64) {
    let mut log_s = s0.ln();
    let mut sum_s = 0.0;
    let mut sum_log_s = 0.0;
    for _ in 0..n_steps {
        let z: f64 = StandardNormal.sample(rng);
        log_s += drift + vol * z;
        let s = log_s.exp();
        sum_s += s;
        sum_log_s += log_s;
    }
    let inv_n = 1.0 / (n_steps as f64);
    let s_t = log_s.exp();
    let arith_avg = sum_s * inv_n;
    let geom_avg = (sum_log_s * inv_n).exp();
    (s_t, arith_avg, geom_avg)
}

#[inline]
fn simulate_antithetic_pair_arith(
    rng: &mut ChaCha8Rng,
    s0: f64,
    drift: f64,
    vol: f64,
    n_steps: usize,
) -> (f64, f64) {
    let mut log_s_pos = s0.ln();
    let mut log_s_neg = log_s_pos;
    let mut sum_pos = 0.0;
    let mut sum_neg = 0.0;
    for _ in 0..n_steps {
        let z: f64 = StandardNormal.sample(rng);
        let inc = vol * z;
        log_s_pos += drift + inc;
        log_s_neg += drift - inc;
        sum_pos += log_s_pos.exp();
        sum_neg += log_s_neg.exp();
    }
    let inv_n = 1.0 / (n_steps as f64);
    (sum_pos * inv_n, sum_neg * inv_n)
}

fn summarize_moments(count: usize, sum: f64, sum_sq: f64) -> EstimatorResult {
    let n = count as f64;
    let mean = sum / n;
    let var = ((sum_sq - n * mean * mean) / (n - 1.0)).max(0.0);
    let se = (var / n).sqrt();
    EstimatorResult {
        price: mean,
        se,
        ci_low: mean - CI_Z_95 * se,
        ci_high: mean + CI_Z_95 * se,
    }
}

fn validate_mc_inputs(p: &OptionParams, n_paths: usize) -> Result<(), String> {
    if p.n_steps < 1 {
        return Err("n_steps must be >= 1".to_string());
    }
    if n_paths < 2 {
        return Err("n_paths must be >= 2".to_string());
    }
    Ok(())
}

pub fn price_european_mc(
    p: &OptionParams,
    n_paths: usize,
    seed: u64,
) -> Result<EstimatorResult, String> {
    validate_mc_inputs(p, n_paths)?;
    let dt = p.t / (p.n_steps as f64);
    let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    let vol = p.sigma * dt.sqrt();
    let disc = (-p.r * p.t).exp();
    let n_chunks = n_paths.div_ceil(CHUNK_SIZE);

    let (sum, sum_sq) = (0..n_chunks)
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK_SIZE;
            let len = (n_paths - start).min(CHUNK_SIZE);
            let mut rng = chunk_rng(seed, c);
            let mut local_sum = 0.0;
            let mut local_sq = 0.0;
            for _ in 0..len {
                let (s_t, _, _) = simulate_one_path_avgs(&mut rng, p.s0, drift, vol, p.n_steps);
                let payoff = (s_t - p.k).max(0.0) * disc;
                local_sum += payoff;
                local_sq += payoff * payoff;
            }
            (local_sum, local_sq)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    Ok(summarize_moments(n_paths, sum, sum_sq))
}

pub fn price_asian_naive(
    p: &OptionParams,
    n_paths: usize,
    seed: u64,
) -> Result<EstimatorResult, String> {
    validate_mc_inputs(p, n_paths)?;
    let dt = p.t / (p.n_steps as f64);
    let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    let vol = p.sigma * dt.sqrt();
    let disc = (-p.r * p.t).exp();
    let n_chunks = n_paths.div_ceil(CHUNK_SIZE);

    let (sum, sum_sq) = (0..n_chunks)
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK_SIZE;
            let len = (n_paths - start).min(CHUNK_SIZE);
            let mut rng = chunk_rng(seed, c);
            let mut local_sum = 0.0;
            let mut local_sq = 0.0;
            for _ in 0..len {
                let (_, arith_avg, _) =
                    simulate_one_path_avgs(&mut rng, p.s0, drift, vol, p.n_steps);
                let payoff = (arith_avg - p.k).max(0.0) * disc;
                local_sum += payoff;
                local_sq += payoff * payoff;
            }
            (local_sum, local_sq)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    Ok(summarize_moments(n_paths, sum, sum_sq))
}

pub fn price_asian_antithetic(
    p: &OptionParams,
    n_paths: usize,
    seed: u64,
) -> Result<EstimatorResult, String> {
    if p.n_steps < 1 {
        return Err("n_steps must be >= 1".to_string());
    }
    if n_paths < 4 || !n_paths.is_multiple_of(2) {
        return Err("n_paths must be even and >= 4 for antithetic variates".to_string());
    }
    let n_pairs = n_paths / 2;
    let dt = p.t / (p.n_steps as f64);
    let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    let vol = p.sigma * dt.sqrt();
    let disc = (-p.r * p.t).exp();
    let n_chunks = n_pairs.div_ceil(CHUNK_SIZE);

    let (sum, sum_sq) = (0..n_chunks)
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK_SIZE;
            let len = (n_pairs - start).min(CHUNK_SIZE);
            let mut rng = chunk_rng(seed, c);
            let mut local_sum = 0.0;
            let mut local_sq = 0.0;
            for _ in 0..len {
                let (avg_pos, avg_neg) =
                    simulate_antithetic_pair_arith(&mut rng, p.s0, drift, vol, p.n_steps);
                let pair_avg = 0.5 * ((avg_pos - p.k).max(0.0) + (avg_neg - p.k).max(0.0)) * disc;
                local_sum += pair_avg;
                local_sq += pair_avg * pair_avg;
            }
            (local_sum, local_sq)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    Ok(summarize_moments(n_pairs, sum, sum_sq))
}

pub fn price_asian_control_variate(
    p: &OptionParams,
    n_paths: usize,
    seed: u64,
) -> Result<ControlVariateResult, String> {
    validate_mc_inputs(p, n_paths)?;
    let y_exact = geometric_asian_call(p)?;
    let dt = p.t / (p.n_steps as f64);
    let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    let vol = p.sigma * dt.sqrt();
    let disc = (-p.r * p.t).exp();
    let n_chunks = n_paths.div_ceil(CHUNK_SIZE);

    let (sx, sy, sxx, syy, sxy) = (0..n_chunks)
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK_SIZE;
            let len = (n_paths - start).min(CHUNK_SIZE);
            let mut rng = chunk_rng(seed, c);
            let (mut lx, mut ly, mut lxx, mut lyy, mut lxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for _ in 0..len {
                let (_, arith_avg, geom_avg) =
                    simulate_one_path_avgs(&mut rng, p.s0, drift, vol, p.n_steps);
                let x = (arith_avg - p.k).max(0.0) * disc;
                let y = (geom_avg - p.k).max(0.0) * disc;
                lx += x;
                ly += y;
                lxx += x * x;
                lyy += y * y;
                lxy += x * y;
            }
            (lx, ly, lxx, lyy, lxy)
        })
        .reduce(
            || (0.0, 0.0, 0.0, 0.0, 0.0),
            |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3, a.4 + b.4),
        );

    let n = n_paths as f64;
    let mx = sx / n;
    let my = sy / n;
    let cov_xy = (sxy - n * mx * my) / (n - 1.0);
    let var_y = (syy - n * my * my) / (n - 1.0);
    if var_y <= 0.0 {
        return Err("var_y is zero; cannot compute control variate beta".to_string());
    }
    let beta = cov_xy / var_y;

    let z_mean = mx - beta * (my - y_exact);
    let var_x = (sxx - n * mx * mx) / (n - 1.0);
    let var_z = (var_x + beta * beta * var_y - 2.0 * beta * cov_xy).max(0.0);
    let se = (var_z / n).sqrt();

    Ok(ControlVariateResult {
        price: z_mean,
        se,
        ci_low: z_mean - CI_Z_95 * se,
        ci_high: z_mean + CI_Z_95 * se,
        beta,
        y_exact,
    })
}

pub fn delta_asian_pathwise(
    p: &OptionParams,
    n_paths: usize,
    seed: u64,
) -> Result<DeltaResult, String> {
    validate_mc_inputs(p, n_paths)?;
    let dt = p.t / (p.n_steps as f64);
    let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    let vol = p.sigma * dt.sqrt();
    let disc = (-p.r * p.t).exp();
    let n_chunks = n_paths.div_ceil(CHUNK_SIZE);

    let (sum, sum_sq) = (0..n_chunks)
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK_SIZE;
            let len = (n_paths - start).min(CHUNK_SIZE);
            let mut rng = chunk_rng(seed, c);
            let mut local_sum = 0.0;
            let mut local_sq = 0.0;
            for _ in 0..len {
                let (_, arith_avg, _) =
                    simulate_one_path_avgs(&mut rng, p.s0, drift, vol, p.n_steps);
                let sample = if arith_avg > p.k {
                    disc * (arith_avg / p.s0)
                } else {
                    0.0
                };
                local_sum += sample;
                local_sq += sample * sample;
            }
            (local_sum, local_sq)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    let res = summarize_moments(n_paths, sum, sum_sq);
    Ok(DeltaResult {
        delta: res.price,
        se: res.se,
        ci_low: res.ci_low,
        ci_high: res.ci_high,
    })
}

pub fn delta_asian_bump(
    p: &OptionParams,
    n_paths: usize,
    h_rel: f64,
    seed: u64,
) -> Result<DeltaResult, String> {
    validate_mc_inputs(p, n_paths)?;
    if h_rel <= 0.0 {
        return Err("h_rel must be > 0".to_string());
    }
    let h = h_rel * p.s0;
    let scale_up = (p.s0 + h) / p.s0;
    let scale_dn = (p.s0 - h) / p.s0;
    let dt = p.t / (p.n_steps as f64);
    let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    let vol = p.sigma * dt.sqrt();
    let disc = (-p.r * p.t).exp();
    let n_chunks = n_paths.div_ceil(CHUNK_SIZE);

    let (sum, sum_sq) = (0..n_chunks)
        .into_par_iter()
        .map(|c| {
            let start = c * CHUNK_SIZE;
            let len = (n_paths - start).min(CHUNK_SIZE);
            let mut rng = chunk_rng(seed, c);
            let mut local_sum = 0.0;
            let mut local_sq = 0.0;
            for _ in 0..len {
                let (_, arith_avg, _) =
                    simulate_one_path_avgs(&mut rng, p.s0, drift, vol, p.n_steps);
                let payoff_up = (arith_avg * scale_up - p.k).max(0.0);
                let payoff_dn = (arith_avg * scale_dn - p.k).max(0.0);
                let sample = disc * (payoff_up - payoff_dn) / (2.0 * h);
                local_sum += sample;
                local_sq += sample * sample;
            }
            (local_sum, local_sq)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    let res = summarize_moments(n_paths, sum, sum_sq);
    Ok(DeltaResult {
        delta: res.price,
        se: res.se,
        ci_low: res.ci_low,
        ci_high: res.ci_high,
    })
}

pub fn fit_loglog_slope(ns: &[usize], ses: &[f64]) -> Result<f64, String> {
    if ns.len() != ses.len() || ns.len() < 2 {
        return Err("ns and ses must have equal length >= 2".to_string());
    }
    let n = ns.len() as f64;
    let xs: Vec<f64> = ns.iter().map(|&v| (v as f64).log10()).collect();
    let ys: Vec<f64> = ses.iter().map(|&v| v.log10()).collect();
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let num: f64 = xs
        .iter()
        .zip(ys.iter())
        .map(|(x, y)| (x - mx) * (y - my))
        .sum();
    let den: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    if den <= 0.0 {
        return Err("zero variance in x".to_string());
    }
    Ok(num / den)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> OptionParams {
        OptionParams {
            s0: 100.0,
            k: 100.0,
            r: 0.05,
            sigma: 0.20,
            t: 1.0,
            n_steps: 252,
        }
    }

    #[test]
    fn test_bs_put_call_parity() {
        let p = OptionParams {
            k: 105.0,
            ..default_params()
        };
        let c = bs_call(&p);
        let put = bs_put(&p);
        let rhs = p.s0 - p.k * (-p.r * p.t).exp();
        assert!((c - put - rhs).abs() < 1e-10);
    }

    #[test]
    fn test_geometric_asian_matches_large_mc() {
        let p = default_params();
        let exact = geometric_asian_call(&p).expect("valid params");
        let n_paths = 50_000usize;
        let dt = p.t / (p.n_steps as f64);
        let drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
        let vol = p.sigma * dt.sqrt();
        let disc = (-p.r * p.t).exp();
        let n_chunks = n_paths.div_ceil(CHUNK_SIZE);

        let (sum, sum_sq) = (0..n_chunks)
            .into_par_iter()
            .map(|c| {
                let start = c * CHUNK_SIZE;
                let len = (n_paths - start).min(CHUNK_SIZE);
                let mut rng = chunk_rng(42, c);
                let mut s = 0.0;
                let mut sq = 0.0;
                for _ in 0..len {
                    let (_, _, geom_avg) =
                        simulate_one_path_avgs(&mut rng, p.s0, drift, vol, p.n_steps);
                    let g_pay = (geom_avg - p.k).max(0.0) * disc;
                    s += g_pay;
                    sq += g_pay * g_pay;
                }
                (s, sq)
            })
            .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        let mc = summarize_moments(n_paths, sum, sum_sq);
        assert!((mc.price - exact).abs() < 3.0 * mc.se);
    }

    #[test]
    fn test_control_variate_beats_naive_se() {
        let p = default_params();
        let naive = price_asian_naive(&p, 20_000, 42).expect("naive ok");
        let cv = price_asian_control_variate(&p, 20_000, 42).expect("cv ok");
        assert!(cv.se < naive.se / 10.0);
    }

    #[test]
    fn test_reproducible_same_seed() {
        let p = OptionParams {
            n_steps: 50,
            ..default_params()
        };
        let a = price_asian_control_variate(&p, 5_000, 7).expect("a ok");
        let b = price_asian_control_variate(&p, 5_000, 7).expect("b ok");
        assert_eq!(a.price, b.price);
        assert_eq!(a.se, b.se);
    }

    #[test]
    fn test_delta_pathwise_matches_bump() {
        let p = default_params();
        let pw = delta_asian_pathwise(&p, 50_000, 42).expect("pw ok");
        let bmp = delta_asian_bump(&p, 50_000, 0.01, 42).expect("bump ok");
        let se_comb = pw.se.hypot(bmp.se);
        assert!((pw.delta - bmp.delta).abs() < 3.0 * se_comb);
    }
}
