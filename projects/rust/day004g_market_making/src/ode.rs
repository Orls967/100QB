//! Pure Rust RK4 and Richardson extrapolation ODE integrator.

pub fn rk4_system<F>(t0: f64, t1: f64, y0: &[f64], n_steps: usize, deriv: F) -> Vec<f64>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let dim = y0.len();
    let mut y = y0.to_vec();
    if n_steps == 0 || t1 == t0 {
        return y;
    }
    let h = (t1 - t0) / n_steps as f64;
    let mut k1 = vec![0.0; dim];
    let mut k2 = vec![0.0; dim];
    let mut k3 = vec![0.0; dim];
    let mut k4 = vec![0.0; dim];
    let mut tmp = vec![0.0; dim];

    for step in 0..n_steps {
        let t = t0 + step as f64 * h;
        deriv(t, &y, &mut k1);
        for i in 0..dim {
            tmp[i] = y[i] + 0.5 * h * k1[i];
        }
        deriv(t + 0.5 * h, &tmp, &mut k2);
        for i in 0..dim {
            tmp[i] = y[i] + 0.5 * h * k2[i];
        }
        deriv(t + 0.5 * h, &tmp, &mut k3);
        for i in 0..dim {
            tmp[i] = y[i] + h * k3[i];
        }
        deriv(t + h, &tmp, &mut k4);
        for i in 0..dim {
            y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    y
}

pub fn rk4_richardson<F>(t0: f64, t1: f64, y0: &[f64], n_coarse: usize, deriv: F) -> Vec<f64>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let y1 = rk4_system(t0, t1, y0, n_coarse, &deriv);
    let y2 = rk4_system(t0, t1, y0, n_coarse * 2, &deriv);
    y1.iter()
        .zip(y2.iter())
        .map(|(&a, &b)| (16.0 * b - a) / 15.0)
        .collect()
}
