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

pub fn sharpe(xs: &[f64]) -> f64 {
    let s = sd(xs);
    if s > 1e-15 {
        mean(xs) / s
    } else {
        0.0
    }
}

pub fn corr(xs: &[f64], ys: &[f64]) -> f64 {
    if xs.len() != ys.len() || xs.len() < 2 {
        return 0.0;
    }
    let mx = mean(xs);
    let my = mean(ys);
    let mut sxy = 0.0_f64;
    let mut sxx = 0.0_f64;
    let mut syy = 0.0_f64;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mx;
        let dy = y - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    let den = (sxx * syy).sqrt();
    if den > 1e-18 {
        sxy / den
    } else {
        0.0
    }
}
