use rand_chacha::ChaCha8Rng;

use crate::rng::{normal_log_pdf, sample_normal, sample_uniform};

/// Reusable buffers for allocation-free bootstrap particle filtering.
pub struct PfWorkspace {
    pub particles: Vec<f64>,
    pub next_particles: Vec<f64>,
    pub log_weights: Vec<f64>,
    pub weights: Vec<f64>,
    pub cdf: Vec<f64>,
    pub ancestors: Vec<usize>,
}

impl PfWorkspace {
    pub fn new(n_particles: usize) -> Self {
        Self {
            particles: vec![0.0; n_particles],
            next_particles: vec![0.0; n_particles],
            log_weights: vec![0.0; n_particles],
            weights: vec![0.0; n_particles],
            cdf: vec![0.0; n_particles],
            ancestors: vec![0; n_particles],
        }
    }
}

/// Systematic resampling using pre-allocated `cdf` and `ancestors` buffers.
pub fn systematic_resample(
    weights: &[f64],
    cdf: &mut [f64],
    ancestors: &mut [usize],
    rng: &mut ChaCha8Rng,
) {
    let n = weights.len();
    let mut acc = 0.0;
    for i in 0..n {
        acc += weights[i];
        cdf[i] = acc;
    }
    cdf[n - 1] = 1.0;

    let inv_n = 1.0 / (n as f64);
    let u0 = sample_uniform(rng) * inv_n;
    let mut i = 0usize;
    for (j, anc) in ancestors.iter_mut().enumerate() {
        let u = u0 + (j as f64) * inv_n;
        while i + 1 < n && u > cdf[i] {
            i += 1;
        }
        *anc = i;
    }
}

/// Normalize log-weights and return `log_mean_exp` increment for the marginal likelihood.
#[inline]
fn normalize_log_weights(log_weights: &[f64], weights: &mut [f64]) -> f64 {
    let n = log_weights.len();
    let mut max_lw = f64::NEG_INFINITY;
    for &lw in log_weights {
        if lw > max_lw {
            max_lw = lw;
        }
    }
    let mut sum_w = 0.0;
    for i in 0..n {
        let w = (log_weights[i] - max_lw).exp();
        weights[i] = w;
        sum_w += w;
    }
    let inv_sum = 1.0 / sum_w;
    for w in weights.iter_mut() {
        *w *= inv_sum;
    }
    max_lw + (sum_w / (n as f64)).ln()
}

/// Bootstrap Particle Filter for Model A (local level). Zero heap allocations inside the loop.
pub fn pf_model_a(
    obs: &[f64],
    sigma_eta: f64,
    sigma_eps: f64,
    ws: &mut PfWorkspace,
    rng: &mut ChaCha8Rng,
) -> f64 {
    let n = ws.particles.len();
    for p in ws.particles.iter_mut() {
        *p = sample_normal(rng);
    }

    let mut ll = 0.0;
    for &y in obs {
        for i in 0..n {
            let x_prop = ws.particles[i] + sigma_eta * sample_normal(rng);
            ws.next_particles[i] = x_prop;
            ws.log_weights[i] = normal_log_pdf(y, x_prop, sigma_eps);
        }
        ll += normalize_log_weights(&ws.log_weights, &mut ws.weights);
        systematic_resample(&ws.weights, &mut ws.cdf, &mut ws.ancestors, rng);
        for i in 0..n {
            ws.particles[i] = ws.next_particles[ws.ancestors[i]];
        }
    }
    ll
}

pub struct SvFilteredState {
    pub log_lik: f64,
    pub mean_h: Vec<f64>,
    pub q025_h: Vec<f64>,
    pub q975_h: Vec<f64>,
}

/// Bootstrap Particle Filter for Model B (stochastic volatility). Zero heap allocations inside the loop.
pub fn pf_model_b(
    obs: &[f64],
    mu: f64,
    phi: f64,
    sigma: f64,
    ws: &mut PfWorkspace,
    rng: &mut ChaCha8Rng,
) -> f64 {
    let n = ws.particles.len();
    let sd0 = sigma / (1.0 - phi * phi).sqrt();
    for p in ws.particles.iter_mut() {
        *p = mu + sd0 * sample_normal(rng);
    }

    let log_2pi = (2.0 * std::f64::consts::PI).ln();
    let mut ll = 0.0;
    for &y in obs {
        let y_sq = y * y;
        for i in 0..n {
            let h_prop = mu + phi * (ws.particles[i] - mu) + sigma * sample_normal(rng);
            ws.next_particles[i] = h_prop;
            ws.log_weights[i] = -0.5 * (log_2pi + h_prop + y_sq * (-h_prop).exp());
        }
        ll += normalize_log_weights(&ws.log_weights, &mut ws.weights);
        systematic_resample(&ws.weights, &mut ws.cdf, &mut ws.ancestors, rng);
        for i in 0..n {
            ws.particles[i] = ws.next_particles[ws.ancestors[i]];
        }
    }
    ll
}

/// Run Model B PF and also record filtered `h_t` mean and 2.5%/97.5% bands at each step.
pub fn pf_model_b_filtered(
    obs: &[f64],
    mu: f64,
    phi: f64,
    sigma: f64,
    n_particles: usize,
    rng: &mut ChaCha8Rng,
) -> SvFilteredState {
    let mut ws = PfWorkspace::new(n_particles);
    let sd0 = sigma / (1.0 - phi * phi).sqrt();
    for p in ws.particles.iter_mut() {
        *p = mu + sd0 * sample_normal(rng);
    }

    let log_2pi = (2.0 * std::f64::consts::PI).ln();
    let mut ll = 0.0;
    let mut mean_h = Vec::with_capacity(obs.len());
    let mut q025_h = Vec::with_capacity(obs.len());
    let mut q975_h = Vec::with_capacity(obs.len());
    let mut sorted = vec![0.0; n_particles];

    for &y in obs {
        let y_sq = y * y;
        for i in 0..n_particles {
            let h_prop = mu + phi * (ws.particles[i] - mu) + sigma * sample_normal(rng);
            ws.next_particles[i] = h_prop;
            ws.log_weights[i] = -0.5 * (log_2pi + h_prop + y_sq * (-h_prop).exp());
        }
        ll += normalize_log_weights(&ws.log_weights, &mut ws.weights);
        systematic_resample(&ws.weights, &mut ws.cdf, &mut ws.ancestors, rng);
        let mut sum = 0.0;
        for (i, s_val) in sorted.iter_mut().enumerate() {
            let val = ws.next_particles[ws.ancestors[i]];
            ws.particles[i] = val;
            *s_val = val;
            sum += val;
        }
        sorted.sort_by(|a, b| a.total_cmp(b));
        mean_h.push(sum / (n_particles as f64));
        q025_h.push(sorted[((0.025 * ((n_particles - 1) as f64)).round()) as usize]);
        q975_h.push(sorted[((0.975 * ((n_particles - 1) as f64)).round()) as usize]);
    }

    SvFilteredState {
        log_lik: ll,
        mean_h,
        q025_h,
        q975_h,
    }
}
