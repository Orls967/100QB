use crate::rng::Rng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Activation {
    Tanh,
    Softplus,
    Relu,
}

impl Activation {
    #[inline]
    pub fn apply(self, z: f64) -> f64 {
        match self {
            Activation::Tanh => z.tanh(),
            Activation::Softplus => {
                if z > 30.0 {
                    z
                } else if z < -30.0 {
                    z.exp()
                } else {
                    (1.0 + z.exp()).ln()
                }
            }
            Activation::Relu => z.max(0.0),
        }
    }

    #[inline]
    pub fn deriv_from_z(self, z: f64, a: f64) -> f64 {
        match self {
            Activation::Tanh => 1.0 - a * a,
            Activation::Softplus => {
                if z >= 0.0 {
                    1.0 / (1.0 + (-z).exp())
                } else {
                    let ez = z.exp();
                    ez / (1.0 + ez)
                }
            }
            Activation::Relu => {
                if z > 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
        }
    }
}

/// 2-hidden-layer MLP policy: `d_in -> hidden -> hidden -> 1` (linear output `delta in R`).
/// Parameters are stored in a flat vector `[W1, b1, W2, b2, W3, b3]` for Adam & FD checks.
#[derive(Clone, Debug)]
pub struct MlpPolicy {
    pub d_in: usize,
    pub hidden: usize,
    pub act: Activation,
    pub params: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct StepCache {
    pub x: Vec<f64>,
    pub z1: Vec<f64>,
    pub h1: Vec<f64>,
    pub z2: Vec<f64>,
    pub h2: Vec<f64>,
    pub out: f64,
}

impl StepCache {
    pub fn new(d_in: usize, hidden: usize) -> Self {
        Self {
            x: vec![0.0; d_in],
            z1: vec![0.0; hidden],
            h1: vec![0.0; hidden],
            z2: vec![0.0; hidden],
            h2: vec![0.0; hidden],
            out: 0.0,
        }
    }
}

impl MlpPolicy {
    pub fn n_params_for(d_in: usize, hidden: usize) -> usize {
        hidden * d_in + hidden + hidden * hidden + hidden + hidden + 1
    }

    pub fn new(d_in: usize, hidden: usize, act: Activation, seed: u64) -> Self {
        let mut rg = Rng::new(seed);
        let n_p = Self::n_params_for(d_in, hidden);
        let mut params = vec![0.0_f64; n_p];

        let s1 = (1.0 / d_in as f64).sqrt();
        let s2 = (1.0 / hidden as f64).sqrt();
        let s3 = (1.0 / hidden as f64).sqrt();

        let off_w1 = 0;
        let off_b1 = off_w1 + hidden * d_in;
        let off_w2 = off_b1 + hidden;
        let off_b2 = off_w2 + hidden * hidden;
        let off_w3 = off_b2 + hidden;

        for p in &mut params[off_w1..off_b1] {
            *p = s1 * rg.normal();
        }
        for p in &mut params[off_w2..off_b2] {
            *p = s2 * rg.normal();
        }
        for p in &mut params[off_w3..off_w3 + hidden] {
            *p = s3 * rg.normal();
        }
        Self {
            d_in,
            hidden,
            act,
            params,
        }
    }

    #[inline]
    pub fn n_params(&self) -> usize {
        self.params.len()
    }

    pub fn set_output_bias(&mut self, bias: f64) {
        let idx = self.params.len() - 1;
        self.params[idx] = bias;
    }

    pub fn zero_output_layer(&mut self) {
        let (_, _, _, _, ow3, _) = self.offsets();
        for p in &mut self.params[ow3..] {
            *p = 0.0;
        }
    }

    #[inline]
    fn offsets(&self) -> (usize, usize, usize, usize, usize, usize) {
        let d = self.d_in;
        let h = self.hidden;
        let w1 = 0;
        let b1 = w1 + h * d;
        let w2 = b1 + h;
        let b2 = w2 + h * h;
        let w3 = b2 + h;
        let b3 = w3 + h;
        (w1, b1, w2, b2, w3, b3)
    }

    pub fn forward_step(&self, x: &[f64], cache: &mut StepCache) -> f64 {
        debug_assert_eq!(x.len(), self.d_in);
        cache.x.copy_from_slice(x);
        let d = self.d_in;
        let h = self.hidden;
        let (ow1, ob1, ow2, ob2, ow3, ob3) = self.offsets();

        for i in 0..h {
            let mut s = self.params[ob1 + i];
            let row = ow1 + i * d;
            for (j, &xj) in x.iter().enumerate().take(d) {
                s += self.params[row + j] * xj;
            }
            cache.z1[i] = s;
            cache.h1[i] = self.act.apply(s);
        }

        for i in 0..h {
            let mut s = self.params[ob2 + i];
            let row = ow2 + i * h;
            for j in 0..h {
                s += self.params[row + j] * cache.h1[j];
            }
            cache.z2[i] = s;
            cache.h2[i] = self.act.apply(s);
        }

        let mut out = self.params[ob3];
        for i in 0..h {
            out += self.params[ow3 + i] * cache.h2[i];
        }
        cache.out = out;
        out
    }

    /// Accumulates gradients into `grad` and writes `dL/dx` into `d_in_out` (overwritten).
    pub fn backward_step(
        &self,
        cache: &StepCache,
        d_out: f64,
        grad: &mut [f64],
        d_in_out: &mut [f64],
        dz2_buf: &mut [f64],
        dz1_buf: &mut [f64],
    ) {
        let d = self.d_in;
        let h = self.hidden;
        let (ow1, ob1, ow2, ob2, ow3, ob3) = self.offsets();

        // Output layer
        grad[ob3] += d_out;
        for i in 0..h {
            grad[ow3 + i] += d_out * cache.h2[i];
            let dh2 = d_out * self.params[ow3 + i];
            dz2_buf[i] = dh2 * self.act.deriv_from_z(cache.z2[i], cache.h2[i]);
        }

        // Hidden layer 2
        for i in 0..h {
            let dz2 = dz2_buf[i];
            grad[ob2 + i] += dz2;
            let row = ow2 + i * h;
            for j in 0..h {
                grad[row + j] += dz2 * cache.h1[j];
            }
        }
        for (j, dz1_j) in dz1_buf.iter_mut().enumerate().take(h) {
            let mut dh1 = 0.0_f64;
            for (i, &dz2_i) in dz2_buf.iter().enumerate().take(h) {
                dh1 += self.params[ow2 + i * h + j] * dz2_i;
            }
            *dz1_j = dh1 * self.act.deriv_from_z(cache.z1[j], cache.h1[j]);
        }

        // Hidden layer 1
        d_in_out.fill(0.0);
        for i in 0..h {
            let dz1 = dz1_buf[i];
            grad[ob1 + i] += dz1;
            let row = ow1 + i * d;
            for j in 0..d {
                grad[row + j] += dz1 * cache.x[j];
                d_in_out[j] += self.params[row + j] * dz1;
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum LossType {
    Mse,
    Entropic { lambda: f64 },
    Cvar { alpha: f64 },
}

/// Computes `(loss, d_loss / d_pnl_i)` for terminal net PnL vector `pnl`.
pub fn loss_and_d_pnl(pnl: &[f64], loss_type: LossType) -> (f64, Vec<f64>) {
    let n = pnl.len();
    let mut d_pnl = vec![0.0_f64; n];
    if n == 0 {
        return (0.0, d_pnl);
    }
    let nf = n as f64;
    match loss_type {
        LossType::Mse => {
            let mut sum_sq = 0.0_f64;
            for (i, &x) in pnl.iter().enumerate() {
                sum_sq += x * x;
                d_pnl[i] = 2.0 * x / nf;
            }
            (sum_sq / nf, d_pnl)
        }
        LossType::Entropic { lambda } => {
            let mut max_arg = f64::NEG_INFINITY;
            for &x in pnl {
                let a = -lambda * x;
                if a > max_arg {
                    max_arg = a;
                }
            }
            let mut sum_exp = 0.0_f64;
            for (i, &x) in pnl.iter().enumerate() {
                let e = (-lambda * x - max_arg).exp();
                d_pnl[i] = e;
                sum_exp += e;
            }
            let loss = (max_arg + (sum_exp / nf).ln()) / lambda;
            for d in &mut d_pnl {
                *d = -(*d) / sum_exp;
            }
            (loss, d_pnl)
        }
        LossType::Cvar { alpha } => {
            let mut idx_loss: Vec<(usize, f64)> =
                pnl.iter().enumerate().map(|(i, &x)| (i, -x)).collect();
            idx_loss.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            let k = ((1.0 - alpha) * nf).ceil().max(1.0) as usize;
            let k = k.min(n);
            let mut sum_l = 0.0_f64;
            let w = -1.0 / (k as f64);
            for &(idx, l) in idx_loss.iter().take(k) {
                sum_l += l;
                d_pnl[idx] = w;
            }
            (sum_l / (k as f64), d_pnl)
        }
    }
}

/// Runs a multi-step recurrent trajectory `pnl = p0 + sum_t delta_t * ds_t - c * sum_t |delta_t - delta_{t-1}| * s_t - payoff`
/// and returns the maximum relative error of analytical BPTT gradient vs central FD + Richardson extrapolation.
pub fn measure_bptt_grad_max_rel_err(loss_type: LossType, cost: f64, act: Activation) -> f64 {
    let n_steps = 10_usize;
    let n_paths = 40_usize;
    let d_in = 3_usize; // [ln(S_t/K), tau, delta_{t-1}]
    let hidden = 8_usize;
    let mut policy = MlpPolicy::new(d_in, hidden, act, 12345);
    let p0 = 0.05_f64;

    let mut rg = Rng::new(999);
    let mut s_paths = vec![vec![1.0_f64; n_steps + 1]; n_paths];
    for path in &mut s_paths {
        let mut s = 1.0_f64;
        for t in 0..n_steps {
            let dt = 0.25 / n_steps as f64;
            let z = rg.normal();
            s *= (-0.5 * 0.04 * dt + 0.2 * dt.sqrt() * z).exp();
            path[t + 1] = s;
        }
    }

    let eval_loss = |pol: &MlpPolicy| -> f64 {
        let mut cache = StepCache::new(d_in, hidden);
        let mut pnl = Vec::with_capacity(n_paths);
        for path in &s_paths {
            let mut prev_d = 0.0_f64;
            let mut wealth = p0;
            for t in 0..n_steps {
                let tau = 0.25 * (n_steps - t) as f64 / n_steps as f64;
                let x = [path[t].ln(), tau, prev_d];
                let d = pol.forward_step(&x, &mut cache);
                wealth += d * (path[t + 1] - path[t]) - cost * (d - prev_d).abs() * path[t];
                prev_d = d;
            }
            let payoff = (path[n_steps] - 1.0).max(0.0);
            pnl.push(wealth - payoff);
        }
        loss_and_d_pnl(&pnl, loss_type).0
    };

    let mut pnl = Vec::with_capacity(n_paths);
    let mut caches: Vec<Vec<StepCache>> = Vec::with_capacity(n_paths);
    for path in &s_paths {
        let mut prev_d = 0.0_f64;
        let mut wealth = p0;
        let mut path_caches = Vec::with_capacity(n_steps);
        for t in 0..n_steps {
            let tau = 0.25 * (n_steps - t) as f64 / n_steps as f64;
            let mut c = StepCache::new(d_in, hidden);
            let d = policy.forward_step(&[path[t].ln(), tau, prev_d], &mut c);
            wealth += d * (path[t + 1] - path[t]) - cost * (d - prev_d).abs() * path[t];
            prev_d = d;
            path_caches.push(c);
        }
        let payoff = (path[n_steps] - 1.0).max(0.0);
        pnl.push(wealth - payoff);
        caches.push(path_caches);
    }
    let (_loss0, d_pnl) = loss_and_d_pnl(&pnl, loss_type);
    let mut grad = vec![0.0_f64; policy.n_params()];
    let mut d_in_buf = vec![0.0_f64; d_in];
    let mut dz2 = vec![0.0_f64; hidden];
    let mut dz1 = vec![0.0_f64; hidden];
    for (p_idx, path) in s_paths.iter().enumerate() {
        let dp = d_pnl[p_idx];
        let mut d_prev_delta_next = 0.0_f64;
        for t in (0..n_steps).rev() {
            let c = &caches[p_idx][t];
            let d_cur = c.out;
            let d_prev = c.x[2];
            let sgn_cur = if d_cur >= d_prev { 1.0 } else { -1.0 };
            let d_direct = dp * ((path[t + 1] - path[t]) - cost * sgn_cur * path[t]);
            let d_out = d_direct + d_prev_delta_next;
            policy.backward_step(c, d_out, &mut grad, &mut d_in_buf, &mut dz2, &mut dz1);
            d_prev_delta_next = d_in_buf[2] + dp * cost * sgn_cur * path[t];
        }
    }

    let h = 1e-5_f64;
    let check_indices = [
        0,
        d_in * hidden - 1,
        d_in * hidden,
        policy.n_params() - hidden - 1,
        policy.n_params() - 2,
        policy.n_params() - 1,
    ];
    let mut max_rel = 0.0_f64;
    for &idx in &check_indices {
        let orig = policy.params[idx];
        policy.params[idx] = orig + h;
        let fp1 = eval_loss(&policy);
        policy.params[idx] = orig - h;
        let fm1 = eval_loss(&policy);
        let d1 = (fp1 - fm1) / (2.0 * h);

        policy.params[idx] = orig + 0.5 * h;
        let fp2 = eval_loss(&policy);
        policy.params[idx] = orig - 0.5 * h;
        let fm2 = eval_loss(&policy);
        let d2 = (fp2 - fm2) / h;

        policy.params[idx] = orig;
        let fd_rich = (4.0 * d2 - d1) / 3.0;
        let err = (grad[idx] - fd_rich).abs();
        let scale = grad[idx].abs().max(fd_rich.abs()).max(1e-8);
        max_rel = max_rel.max(err / scale);
    }
    max_rel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bptt_gradient_check_mse_entropic_cvar() {
        for (lt, c, act) in [
            (LossType::Mse, 0.0, Activation::Tanh),
            (LossType::Mse, 0.002, Activation::Softplus),
            (LossType::Entropic { lambda: 10.0 }, 0.001, Activation::Tanh),
            (LossType::Cvar { alpha: 0.95 }, 0.001, Activation::Tanh),
        ] {
            let rel = measure_bptt_grad_max_rel_err(lt, c, act);
            assert!(rel < 1e-6, "BPTT grad mismatch: rel={rel}");
        }
    }
}
