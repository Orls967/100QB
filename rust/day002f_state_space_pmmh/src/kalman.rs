use crate::rng::{cholesky_lower, forward_solve};

/// Exact Kalman filter log-likelihood for Model A (`x_0 ~ N(0, 1)`).
pub fn kalman_log_likelihood(obs: &[f64], sigma_eta: f64, sigma_eps: f64) -> f64 {
    let q_var = sigma_eta * sigma_eta;
    let r_var = sigma_eps * sigma_eps;
    let mut m = 0.0;
    let mut p = 1.0;
    let mut ll = 0.0;
    let log_2pi = (2.0 * std::f64::consts::PI).ln();

    for &y in obs {
        let p_pred = p + q_var;
        let v = y - m;
        let s = p_pred + r_var;
        ll += -0.5 * (log_2pi + s.ln() + (v * v) / s);
        let k = p_pred / s;
        m += k * v;
        p = (1.0 - k) * p_pred;
    }
    ll
}

/// Independent brute-force multivariate Gaussian log-density for Model A on small `T`:
/// `Sigma_{ts} = 1 + min(t, s) * sigma_eta^2 + sigma_eps^2 * 1{t == s}` (1-based `t, s`).
pub fn brute_force_log_likelihood(
    obs: &[f64],
    sigma_eta: f64,
    sigma_eps: f64,
) -> Result<f64, String> {
    let n = obs.len();
    let q_var = sigma_eta * sigma_eta;
    let r_var = sigma_eps * sigma_eps;
    let mut cov = vec![vec![0.0; n]; n];
    for (i, row) in cov.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let min_ts = (i.min(j) + 1) as f64;
            *cell = 1.0 + min_ts * q_var + if i == j { r_var } else { 0.0 };
        }
    }
    let lower = cholesky_lower(&cov)?;
    let z = forward_solve(&lower, obs);
    let quad: f64 = z.iter().map(|v| v * v).sum();
    let log_det: f64 = 2.0
        * lower
            .iter()
            .enumerate()
            .map(|(i, row)| row[i].ln())
            .sum::<f64>();
    let log_2pi = (2.0 * std::f64::consts::PI).ln();
    Ok(-0.5 * ((n as f64) * log_2pi + log_det + quad))
}
