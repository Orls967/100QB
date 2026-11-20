use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct ParamDiagnostics {
    pub mean: f64,
    pub sd: f64,
    pub q025: f64,
    pub q500: f64,
    pub q975: f64,
    pub r_hat: f64,
    pub ess: f64,
    pub mcse: f64,
}

pub fn mean_and_var(xs: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let mean = xs.iter().sum::<f64>() / n;
    let var = xs.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0);
    (mean, var)
}

pub fn quantiles(samples: &[f64]) -> (f64, f64, f64) {
    let mut sorted = samples.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let n = sorted.len();
    let get_q = |p: f64| -> f64 {
        let pos = p * ((n - 1) as f64);
        let lo = pos.floor() as usize;
        let hi = pos.ceil() as usize;
        let frac = pos - (lo as f64);
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    };
    (get_q(0.025), get_q(0.50), get_q(0.975))
}

/// Split R-hat across `chains` (each chain split into two equal halves).
pub fn split_r_hat(chains: &[Vec<f64>]) -> f64 {
    let mut halves: Vec<&[f64]> = Vec::with_capacity(chains.len() * 2);
    for ch in chains {
        let half = ch.len() / 2;
        halves.push(&ch[..half]);
        halves.push(&ch[half..2 * half]);
    }
    let m = halves.len() as f64;
    let n = halves[0].len() as f64;

    let mut means = Vec::with_capacity(halves.len());
    let mut vars = Vec::with_capacity(halves.len());
    for h in &halves {
        let (mu, var) = mean_and_var(h);
        means.push(mu);
        vars.push(var);
    }

    let grand_mean = means.iter().sum::<f64>() / m;
    let b_var = n * means
        .iter()
        .map(|&mu| (mu - grand_mean) * (mu - grand_mean))
        .sum::<f64>()
        / (m - 1.0);
    let w_var = vars.iter().sum::<f64>() / m;
    let var_plus = ((n - 1.0) / n) * w_var + (1.0 / n) * b_var;
    (var_plus / w_var.max(1e-18)).sqrt()
}

/// Sample autocorrelation up to `max_lag` for a single chain.
pub fn autocorrelation(xs: &[f64], max_lag: usize) -> Vec<f64> {
    let n = xs.len();
    let (mean, _) = mean_and_var(xs);
    let denom: f64 = xs.iter().map(|&x| (x - mean) * (x - mean)).sum();
    if denom <= 1e-18 {
        return vec![1.0; max_lag + 1];
    }
    let mut acf = Vec::with_capacity(max_lag + 1);
    acf.push(1.0);
    for lag in 1..=max_lag.min(n - 1) {
        let mut num = 0.0;
        for t in 0..(n - lag) {
            num += (xs[t] - mean) * (xs[t + lag] - mean);
        }
        acf.push(num / denom);
    }
    acf
}

/// Effective Sample Size via Geyer's initial positive sequence across multiple chains.
pub fn effective_sample_size(chains: &[Vec<f64>]) -> f64 {
    let m = chains.len();
    let n = chains[0].len();
    let max_lag = (n / 2).min(1000);
    let mut mean_acf = vec![0.0; max_lag + 1];
    for ch in chains {
        let acf = autocorrelation(ch, max_lag);
        for (i, &v) in acf.iter().enumerate() {
            mean_acf[i] += v / (m as f64);
        }
    }

    // Geyer's initial positive sequence: sum consecutive pairs (rho_{2k-1} + rho_{2k}) while > 0
    let mut tau = 1.0;
    let mut k = 1usize;
    while k < max_lag {
        let pair = mean_acf[k] + mean_acf[k + 1];
        if pair <= 0.0 {
            break;
        }
        tau += 2.0 * pair;
        k += 2;
    }
    ((m * n) as f64) / tau.max(1.0)
}

pub fn summarize_chains(chains: &[Vec<f64>]) -> ParamDiagnostics {
    let pooled: Vec<f64> = chains.iter().flat_map(|c| c.iter().copied()).collect();
    let (mean, var) = mean_and_var(&pooled);
    let sd = var.sqrt();
    let (q025, q500, q975) = quantiles(&pooled);
    let r_hat = split_r_hat(chains);
    let ess = effective_sample_size(chains);
    let mcse = sd / ess.sqrt();
    ParamDiagnostics {
        mean,
        sd,
        q025,
        q500,
        q975,
        r_hat,
        ess,
        mcse,
    }
}
