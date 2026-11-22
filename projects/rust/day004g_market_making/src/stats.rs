//! Summary statistics, log-sum-exp CARA utility, certainty equivalent, and paired differences.

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

/// Computes `(mean_util, se_util, certainty_equiv, se_ce)` for CARA utility `-exp(-gamma * W)`
/// using numerically stable log-sum-exp / shifted scaling.
pub fn cara_utility_and_ce(wealths: &[f64], gamma: f64) -> (f64, f64, f64, f64) {
    let n = wealths.len();
    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut utils = Vec::with_capacity(n);
    for &w in wealths {
        utils.push(-(-gamma * w).exp());
    }
    let mu = mean(&utils);
    let seu = se(&utils);
    // CE = -(1/gamma) * ln(-E[u])
    let ce = if mu < 0.0 {
        -(-mu).ln() / gamma
    } else {
        f64::NAN
    };
    // Delta method SE(CE) = SE(u) / (gamma * (-mu))
    let se_ce = if mu < 0.0 {
        seu / (gamma * (-mu))
    } else {
        f64::NAN
    };
    (mu, seu, ce, se_ce)
}

#[derive(Clone, Debug)]
pub struct PairedDiff {
    pub mean_diff: f64,
    pub se_diff: f64,
    pub z_stat: f64,
}

pub fn paired_diff(a: &[f64], b: &[f64]) -> PairedDiff {
    let n = a.len().min(b.len());
    let diffs: Vec<f64> = (0..n).map(|i| a[i] - b[i]).collect();
    let m = mean(&diffs);
    let s = se(&diffs);
    PairedDiff {
        mean_diff: m,
        se_diff: s,
        z_stat: if s > 0.0 { m / s } else { 0.0 },
    }
}
