//! Summary statistics, quantiles, OLS regression with t-stats, R^2, and Durbin-Watson.

pub fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().sum::<f64>() / xs.len() as f64
}

pub fn sd(xs: &[f64]) -> f64 {
    if xs.len() < 2 {
        return 0.0;
    }
    let m = mean(xs);
    let ss: f64 = xs.iter().map(|&x| (x - m) * (x - m)).sum();
    (ss / (xs.len() - 1) as f64).sqrt()
}

pub fn se(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    sd(xs) / (xs.len() as f64).sqrt()
}

pub fn skewness(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 3 {
        return 0.0;
    }
    let m = mean(xs);
    let s = sd(xs);
    if s <= 1e-18 {
        return 0.0;
    }
    let m3: f64 = xs.iter().map(|&x| ((x - m) / s).powi(3)).sum::<f64>() / n as f64;
    m3
}

/// Computes quantile `q in [0, 1]` and its asymptotic SE via kernel density estimate at the quantile.
pub fn quantile_and_se(sorted: &[f64], q: f64) -> (f64, f64) {
    let n = sorted.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let pos = q * (n - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let w = pos - lo as f64;
    let val = sorted[lo] * (1.0 - w) + sorted[hi] * w;

    // Bandwidth h for density estimate at quantile
    let sd_est = (sorted[(0.75 * (n - 1) as f64) as usize]
        - sorted[(0.25 * (n - 1) as f64) as usize])
        / 1.349;
    let h = (0.9 * sd_est.max(1e-6) * (n as f64).powf(-0.2)).max(1e-4);
    let mut cnt = 0.0_f64;
    for &x in sorted {
        if (x - val).abs() <= h {
            cnt += 1.0;
        }
    }
    let f_hat = (cnt / (2.0 * h * n as f64)).max(1e-6);
    let se_q = (q * (1.0 - q) / (n as f64 * f_hat * f_hat)).sqrt();
    (val, se_q)
}

#[derive(Clone, Debug)]
pub struct OlsResult {
    pub alpha: f64,
    pub beta: f64,
    pub se_alpha: f64,
    pub se_beta: f64,
    pub t_beta: f64,
    pub r2: f64,
    pub durbin_watson: f64,
    pub residuals: Vec<f64>,
}

/// Simple OLS `y = alpha + beta * x + e`.
pub fn ols_simple(x: &[f64], y: &[f64]) -> OlsResult {
    let n = x.len().min(y.len());
    let mx = mean(&x[..n]);
    let my = mean(&y[..n]);
    let mut sxx = 0.0_f64;
    let mut sxy = 0.0_f64;
    let mut syy = 0.0_f64;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let beta = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    let alpha = my - beta * mx;
    let mut residuals = Vec::with_capacity(n);
    let mut rss = 0.0_f64;
    for i in 0..n {
        let e = y[i] - (alpha + beta * x[i]);
        rss += e * e;
        residuals.push(e);
    }
    let s2 = if n > 2 { rss / (n - 2) as f64 } else { 0.0 };
    let se_beta = if sxx > 0.0 { (s2 / sxx).sqrt() } else { 0.0 };
    let se_alpha = if sxx > 0.0 {
        (s2 * (1.0 / n as f64 + mx * mx / sxx)).sqrt()
    } else {
        0.0
    };
    let t_beta = if se_beta > 0.0 { beta / se_beta } else { 0.0 };
    let r2 = if syy > 0.0 { 1.0 - rss / syy } else { 1.0 };

    let mut diff_sq = 0.0_f64;
    for i in 1..n {
        let d = residuals[i] - residuals[i - 1];
        diff_sq += d * d;
    }
    let durbin_watson = if rss > 0.0 { diff_sq / rss } else { 2.0 };

    OlsResult {
        alpha,
        beta,
        se_alpha,
        se_beta,
        t_beta,
        r2,
        durbin_watson,
        residuals,
    }
}
