//! Statistik deskriptif dan estimasi standard error variansi via batch means.

/// Rata-rata sampel.
pub fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().sum::<f64>() / (xs.len() as f64)
}

/// Variansi sampel tak-bias (ddof = 1).
pub fn variance(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 2 {
        return 0.0;
    }
    let m = mean(xs);
    let ss: f64 = xs.iter().map(|&x| (x - m) * (x - m)).sum();
    ss / ((n - 1) as f64)
}

/// Standard error rata-rata.
pub fn mean_se(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    (variance(xs) / (xs.len() as f64)).sqrt()
}

/// Standard error dari variansi sampel menggunakan metode ``n_batches`` batch means.
pub fn variance_se_batch_means(xs: &[f64], n_batches: usize) -> f64 {
    let n = xs.len();
    let b = n_batches.min(n / 2).max(2);
    let batch_size = n / b;
    let mut batch_vars = Vec::with_capacity(b);
    for i in 0..b {
        let slice = &xs[i * batch_size..(i + 1) * batch_size];
        batch_vars.push(variance(slice));
    }
    (variance(&batch_vars) / (b as f64)).sqrt()
}
