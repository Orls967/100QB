//! Summary statistics, OLS linear regression, and variance-reduction metrics.

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

/// Ordinary least squares linear regression `y = intercept + slope * x`.
/// Returns `(intercept, slope, se_slope, r2)`.
pub fn ols_line(x: &[f64], y: &[f64]) -> (f64, f64, f64, f64) {
    let n = x.len().min(y.len());
    if n < 2 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mx = mean(&x[..n]);
    let my = mean(&y[..n]);
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    let intercept = my - slope * mx;
    let mut rss = 0.0;
    for i in 0..n {
        let r = y[i] - (intercept + slope * x[i]);
        rss += r * r;
    }
    let se_slope = if n > 2 && sxx > 0.0 {
        ((rss / (n - 2) as f64) / sxx).sqrt()
    } else {
        0.0
    };
    let r2 = if syy > 0.0 { 1.0 - rss / syy } else { 1.0 };
    (intercept, slope, se_slope, r2)
}
