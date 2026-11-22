pub fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().sum::<f64>() / xs.len() as f64
}

pub fn variance_sample(xs: &[f64]) -> f64 {
    if xs.len() < 2 {
        return 0.0;
    }
    let m = mean(xs);
    let ss: f64 = xs.iter().map(|&x| (x - m) * (x - m)).sum();
    ss / (xs.len() - 1) as f64
}

pub fn sd(xs: &[f64]) -> f64 {
    variance_sample(xs).sqrt()
}

pub fn se(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    sd(xs) / (xs.len() as f64).sqrt()
}

pub fn quantile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let idx = p.clamp(0.0, 1.0) * (sorted.len() - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let frac = idx - lo as f64;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

/// CVaR_alpha of loss `L = -pnl` (i.e. mean of worst `(1 - alpha)` losses).
pub fn cvar_loss(pnl: &[f64], alpha: f64) -> f64 {
    if pnl.is_empty() {
        return 0.0;
    }
    let mut losses: Vec<f64> = pnl.iter().map(|&x| -x).collect();
    losses.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let k = ((1.0 - alpha) * losses.len() as f64).ceil().max(1.0) as usize;
    losses[..k.min(losses.len())].iter().sum::<f64>() / k.min(losses.len()) as f64
}

/// Entropic risk measure `rho_lambda(X) = (1/lambda) * ln(E[exp(-lambda * X)])`
/// using log-sum-exp shift for numerical stability.
pub fn entropic_risk(pnl: &[f64], lambda: f64) -> f64 {
    if pnl.is_empty() {
        return 0.0;
    }
    let mut max_arg = f64::NEG_INFINITY;
    for &x in pnl {
        let a = -lambda * x;
        if a > max_arg {
            max_arg = a;
        }
    }
    let mut sum_exp = 0.0_f64;
    for &x in pnl {
        sum_exp += (-lambda * x - max_arg).exp();
    }
    (max_arg + (sum_exp / pnl.len() as f64).ln()) / lambda
}
