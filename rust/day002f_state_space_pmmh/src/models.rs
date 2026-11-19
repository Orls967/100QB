use crate::rng::{make_rng, normal_log_pdf, sample_normal};

pub const MODEL_A_T: usize = 200;
pub const MODEL_A_SIGMA_ETA: f64 = 0.5;
pub const MODEL_A_SIGMA_EPS: f64 = 1.0;
pub const MODEL_A_SEED: u64 = 42;

pub const MODEL_B_T: usize = 500;
pub const MODEL_B_MU: f64 = -1.0;
pub const MODEL_B_PHI: f64 = 0.95;
pub const MODEL_B_SIGMA: f64 = 0.25;
pub const MODEL_B_SEED: u64 = 43;

pub struct SyntheticSeries {
    pub states: Vec<f64>,
    pub obs: Vec<f64>,
}

/// Simulate Model A (local level): `x_0 ~ N(0, 1)`, `x_t = x_{t-1} + sigma_eta * eta_t`,
/// `y_t = x_t + sigma_eps * eps_t` for `t = 1..=T`.
pub fn simulate_model_a(
    t_len: usize,
    sigma_eta: f64,
    sigma_eps: f64,
    seed: u64,
) -> SyntheticSeries {
    let mut rng = make_rng(seed);
    let mut states = Vec::with_capacity(t_len);
    let mut obs = Vec::with_capacity(t_len);
    let mut x = sample_normal(&mut rng);
    for _ in 0..t_len {
        x += sigma_eta * sample_normal(&mut rng);
        let y = x + sigma_eps * sample_normal(&mut rng);
        states.push(x);
        obs.push(y);
    }
    SyntheticSeries { states, obs }
}

/// Simulate Model B (stochastic volatility): `h_0 ~ N(mu, sigma^2/(1-phi^2))`,
/// `h_t = mu + phi*(h_{t-1}-mu) + sigma*eta_t`, `y_t = exp(h_t/2)*eps_t` for `t = 1..=T`.
pub fn simulate_model_b(t_len: usize, mu: f64, phi: f64, sigma: f64, seed: u64) -> SyntheticSeries {
    let mut rng = make_rng(seed);
    let mut states = Vec::with_capacity(t_len);
    let mut obs = Vec::with_capacity(t_len);
    let sd0 = sigma / (1.0 - phi * phi).sqrt();
    let mut h = mu + sd0 * sample_normal(&mut rng);
    for _ in 0..t_len {
        h = mu + phi * (h - mu) + sigma * sample_normal(&mut rng);
        let y = (0.5 * h).exp() * sample_normal(&mut rng);
        states.push(h);
        obs.push(y);
    }
    SyntheticSeries { states, obs }
}

#[inline]
pub fn model_a_log_prior(a: f64, b: f64) -> f64 {
    normal_log_pdf(a, 0.0, 1.0) + normal_log_pdf(b, 0.0, 1.0)
}

#[inline]
pub fn model_b_log_prior(mu: f64, z_phi: f64, s_sigma: f64) -> f64 {
    normal_log_pdf(mu, 0.0, 1.0)
        + normal_log_pdf(z_phi, 1.5, 1.0)
        + normal_log_pdf(s_sigma, -1.5, 1.0)
}
