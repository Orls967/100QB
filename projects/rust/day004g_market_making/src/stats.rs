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

/// Computes Richardson extrapolation of order `p_order` from `ce2` (step 2*dt) and `ce3` (step dt),
/// together with honest CRN SE (delta-method influence covariance) and independent-sample SE.
pub fn richardson_ce_and_se(
    w2: &[f64],
    ce2: f64,
    se2: f64,
    w3: &[f64],
    ce3: f64,
    se3: f64,
    gamma: f64,
    p_order: f64,
) -> (f64, f64, f64) {
    let factor = 2.0_f64.powf(p_order);
    let denom = factor - 1.0;
    let c3 = factor / denom;
    let c2 = 1.0 / denom;
    let ce_rich = c3 * ce3 - c2 * ce2;
    let se_indep = (c3 * c3 * se3 * se3 + c2 * c2 * se2 * se2).sqrt();

    let n = w2.len().min(w3.len());
    if n < 2 {
        return (ce_rich, se_indep, se_indep);
    }
    let e2: Vec<f64> = (0..n).map(|i| (-gamma * w2[i]).exp()).collect();
    let e3: Vec<f64> = (0..n).map(|i| (-gamma * w3[i]).exp()).collect();
    let m2 = mean(&e2);
    let m3 = mean(&e3);
    let mut comb = Vec::with_capacity(n);
    for i in 0..n {
        let phi2 = -(e2[i] / m2 - 1.0) / gamma;
        let phi3 = -(e3[i] / m3 - 1.0) / gamma;
        comb.push(c3 * phi3 - c2 * phi2);
    }
    let se_crn = se(&comb);
    (ce_rich, se_crn, se_indep)
}

/// Standard error of sample variance using the fourth central moment.
pub fn se_of_variance(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 2 {
        return 0.0;
    }
    let m = mean(xs);
    let mut m2 = 0.0;
    let mut m4 = 0.0;
    for &x in xs {
        let d = x - m;
        let d2 = d * d;
        m2 += d2;
        m4 += d2 * d2;
    }
    let nf = n as f64;
    let s2 = m2 / (nf - 1.0);
    let mu4 = m4 / nf;
    let var_of_s2 = ((mu4 - s2 * s2 * (nf - 3.0) / (nf - 1.0)) / nf).max(0.0);
    var_of_s2.sqrt()
}
