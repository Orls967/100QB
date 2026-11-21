//! Semiannual par coupon generation and recursive discount factor bootstrapping.

use crate::ns::NsParams;

pub fn semi_annual_grid(max_year: f64) -> Vec<f64> {
    let n = (max_year * 2.0).round() as usize;
    (0..n).map(|i| 0.5 * (i + 1) as f64).collect()
}

pub fn par_rates_from_curve(
    curve: &NsParams,
    maturities: &[f64],
    delta: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = maturities.len();
    let mut par_rates = Vec::with_capacity(n);
    let mut true_df = Vec::with_capacity(n);
    let mut true_yields = Vec::with_capacity(n);
    let mut annuity = 0.0;
    for &tau in maturities {
        let y = curve.yield_at(tau);
        let df = (-tau * y).exp();
        annuity += delta * df;
        true_yields.push(y);
        true_df.push(df);
        par_rates.push((1.0 - df) / annuity);
    }
    (par_rates, true_df, true_yields)
}

pub fn bootstrap_from_par_rates(
    maturities: &[f64],
    par_rates: &[f64],
    delta: f64,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = maturities.len();
    if n == 0 || par_rates.len() != n {
        return Err("bootstrap_from_par_rates: invalid slice lengths".into());
    }
    let mut dfs = Vec::with_capacity(n);
    let mut zero_yields = Vec::with_capacity(n);
    let mut prev_annuity = 0.0;
    for i in 0..n {
        let c = par_rates[i];
        let denom = 1.0 + delta * c;
        if denom <= 0.0 {
            return Err(format!("non-positive denominator at k={i}"));
        }
        let p = (1.0 - c * prev_annuity) / denom;
        if p <= 0.0 {
            return Err(format!("non-positive discount factor {p} at k={i}"));
        }
        dfs.push(p);
        zero_yields.push(-p.ln() / maturities[i]);
        prev_annuity += delta * p;
    }
    Ok((dfs, zero_yields))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bootstrap_recovers_true_dfs() {
        let ns = NsParams::default_ns();
        let mats = semi_annual_grid(30.0);
        let (pars, true_df, true_y) = par_rates_from_curve(&ns, &mats, 0.5);
        let (dfs, ys) = bootstrap_from_par_rates(&mats, &pars, 0.5).unwrap();
        for i in 0..mats.len() {
            assert!((dfs[i] - true_df[i]).abs() < 1e-12);
            assert!((ys[i] - true_y[i]).abs() < 1e-12);
        }
    }
}
