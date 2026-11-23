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

/// Linear-interpolated quantile on a slice (clones and sorts).
pub fn quantile(xs: &[f64], q: f64) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.total_cmp(b));
    let pos = q.clamp(0.0, 1.0) * (v.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    v[lo] * (1.0 - frac) + v[hi] * frac
}
