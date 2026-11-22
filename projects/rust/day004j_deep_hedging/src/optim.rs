/// Adam optimizer with standard hyperparameters `beta1 = 0.9, beta2 = 0.999, eps = 1e-8`.
#[derive(Clone, Debug)]
pub struct Adam {
    pub lr: f64,
    pub beta1: f64,
    pub beta2: f64,
    pub eps: f64,
    pub t: u64,
    pub m: Vec<f64>,
    pub v: Vec<f64>,
}

impl Adam {
    pub fn new(n_params: usize, lr: f64) -> Self {
        Self {
            lr,
            beta1: 0.9,
            beta2: 0.999,
            eps: 1e-8,
            t: 0,
            m: vec![0.0; n_params],
            v: vec![0.0; n_params],
        }
    }

    pub fn step(&mut self, params: &mut [f64], grads: &[f64]) {
        debug_assert_eq!(params.len(), grads.len());
        self.t += 1;
        let b1t = 1.0 - self.beta1.powi(self.t as i32);
        let b2t = 1.0 - self.beta2.powi(self.t as i32);
        for i in 0..params.len() {
            let g = grads[i];
            self.m[i] = self.beta1 * self.m[i] + (1.0 - self.beta1) * g;
            self.v[i] = self.beta2 * self.v[i] + (1.0 - self.beta2) * g * g;
            let m_hat = self.m[i] / b1t;
            let v_hat = self.v[i] / b2t;
            params[i] -= self.lr * m_hat / (v_hat.sqrt() + self.eps);
        }
    }
}

/// Runs Adam on a 2D quadratic and a 6-point OLS problem, returning `(quad_max_err, ols_max_err)`.
pub fn measure_adam_errors() -> (f64, f64) {
    let mut x = vec![0.0_f64, 0.0_f64];
    let mut opt = Adam::new(2, 0.05);
    for _ in 0..2000 {
        let g = [6.0 * (x[0] - 2.0), 1.0 * (x[1] + 4.0)];
        opt.step(&mut x, &g);
    }
    let quad_err = (x[0] - 2.0).abs().max((x[1] + 4.0).abs());

    let us = [-2.0_f64, -1.0, 0.0, 1.0, 2.0, 3.0];
    let ys: Vec<f64> = us.iter().map(|&u| 1.5 + 2.5 * u).collect();
    let mut beta = vec![0.0_f64, 0.0_f64];
    let mut opt_ols = Adam::new(2, 0.05);
    for _ in 0..2500 {
        let mut g = [0.0_f64, 0.0_f64];
        for (&u, &y) in us.iter().zip(ys.iter()) {
            let err = (beta[0] + beta[1] * u) - y;
            g[0] += 2.0 * err / us.len() as f64;
            g[1] += 2.0 * err * u / us.len() as f64;
        }
        opt_ols.step(&mut beta, &g);
    }
    let ols_err = (beta[0] - 1.5).abs().max((beta[1] - 2.5).abs());
    (quad_err, ols_err)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adam_converges_quadratic_and_ols() {
        let (quad_err, ols_err) = measure_adam_errors();
        assert!(quad_err < 1e-6, "quad_err={quad_err}");
        assert!(ols_err < 1e-6, "ols_err={ols_err}");
    }
}
