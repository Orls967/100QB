//! Descriptive statistics, quantiles, Wilson CI, OLS/WLS trends, VaR/ES, and paired differences.

pub const Z95: f64 = 1.959963984540054;

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

pub fn rmse(xs: &[f64], target: f64) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let ss: f64 = xs.iter().map(|&x| (x - target) * (x - target)).sum();
    (ss / xs.len() as f64).sqrt()
}

pub fn quantile(xs: &[f64], q: f64) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    let mut s = xs.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if q <= 0.0 {
        return s[0];
    }
    if q >= 1.0 {
        return s[s.len() - 1];
    }
    let pos = q * (s.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    s[lo] * (1.0 - frac) + s[hi] * frac
}

pub fn wilson_ci(k: usize, n: usize, z: f64) -> (f64, f64, f64) {
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    let nf = n as f64;
    let p = k as f64 / nf;
    let z2 = z * z;
    let denom = 1.0 + z2 / nf;
    let center = (p + z2 / (2.0 * nf)) / denom;
    let half = (z / denom) * ((p * (1.0 - p) / nf + z2 / (4.0 * nf * nf)).sqrt());
    ((p), (center - half).max(0.0), (center + half).min(1.0))
}

pub fn var_and_es(xs: &[f64], alpha: f64) -> (f64, f64) {
    if xs.is_empty() {
        return (0.0, 0.0);
    }
    let v = quantile(xs, alpha);
    let mut sum = 0.0;
    let mut cnt = 0usize;
    for &x in xs {
        if x <= v {
            sum += x;
            cnt += 1;
        }
    }
    let es = if cnt > 0 { sum / cnt as f64 } else { v };
    (v, es)
}

#[derive(Clone, Debug)]
pub struct PairedDiffSummary {
    pub n: usize,
    pub mean_diff: f64,
    pub sd_diff: f64,
    pub mcse: f64,
    pub ci95_low: f64,
    pub ci95_high: f64,
    pub positive_frac: f64,
    pub z_stat: f64,
}

pub fn paired_diff_ci(a: &[f64], b: &[f64]) -> PairedDiffSummary {
    let n = a.len().min(b.len());
    let mut diffs = Vec::with_capacity(n);
    let mut pos = 0usize;
    for i in 0..n {
        let d = a[i] - b[i];
        if d > 0.0 {
            pos += 1;
        }
        diffs.push(d);
    }
    let m = mean(&diffs);
    let s = sd(&diffs);
    let se_val = se(&diffs);
    let z = if se_val > 0.0 { m / se_val } else { 0.0 };
    PairedDiffSummary {
        n,
        mean_diff: m,
        sd_diff: s,
        mcse: se_val,
        ci95_low: m - Z95 * se_val,
        ci95_high: m + Z95 * se_val,
        positive_frac: if n > 0 { pos as f64 / n as f64 } else { 0.0 },
        z_stat: z,
    }
}

pub fn ols_slope_intercept(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len().min(y.len());
    if n < 2 {
        return (0.0, 0.0);
    }
    let mx = mean(&x[..n]);
    let my = mean(&y[..n]);
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        num += dx * (y[i] - my);
        den += dx * dx;
    }
    if den == 0.0 {
        return (0.0, my);
    }
    let slope = num / den;
    (slope, my - slope * mx)
}

pub fn wls_slope_intercept(x: &[f64], y: &[f64], w: &[f64]) -> (f64, f64) {
    let n = x.len().min(y.len()).min(w.len());
    if n < 2 {
        return (0.0, 0.0);
    }
    let sw: f64 = w[..n].iter().sum();
    if sw <= 0.0 {
        return ols_slope_intercept(x, y);
    }
    let mx = (0..n).map(|i| w[i] * x[i]).sum::<f64>() / sw;
    let my = (0..n).map(|i| w[i] * y[i]).sum::<f64>() / sw;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        num += w[i] * dx * (y[i] - my);
        den += w[i] * dx * dx;
    }
    if den == 0.0 {
        return (0.0, my);
    }
    let slope = num / den;
    (slope, my - slope * mx)
}

pub fn correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return 0.0;
    }
    let mx = mean(&x[..n]);
    let my = mean(&y[..n]);
    let mut num = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        num += dx * dy;
        sx += dx * dx;
        sy += dy * dy;
    }
    let den = (sx * sy).sqrt();
    if den > 0.0 {
        num / den
    } else {
        0.0
    }
}
