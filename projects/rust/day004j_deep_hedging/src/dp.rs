use crate::bs::norm_pdf;

/// Exact Dynamic Programming solver for entropic risk hedging under GBM
/// with proportional transaction costs `c` at `n_steps = 10`.
///
/// Bellman equation in certainty-equivalent log-value `V_t(x, d_prev)`:
/// `V_t(x, d_prev) = min_d { c * |d - d_prev| * S + (1/lambda) * ln E[ exp(-lambda * d * (S' - S) + lambda * V_{t+1}(x', d)) ] }`
/// At `t = 0, x = ln(S0/K), d_prev = 0`, `V_0` is the exact DP indifference price.
pub struct DpResult {
    pub indifference_price: f64,
    pub x_grid: Vec<f64>,
    pub d_grid: Vec<f64>,
    /// Optimal policy `d_star[t][ix * n_d + id_prev]`
    pub policy: Vec<Vec<f64>>,
    pub n_steps: usize,
    pub k: f64,
}

impl DpResult {
    /// Bilinear interpolation of `delta_t(x, d_prev)`.
    pub fn act(&self, t: usize, s: f64, d_prev: f64) -> f64 {
        let x = (s / self.k).ln();
        let nx = self.x_grid.len();
        let nd = self.d_grid.len();
        let x_min = self.x_grid[0];
        let x_max = self.x_grid[nx - 1];
        let d_min = self.d_grid[0];
        let d_max = self.d_grid[nd - 1];

        let tx = ((x.clamp(x_min, x_max) - x_min) / (x_max - x_min)) * (nx - 1) as f64;
        let ix0 = (tx.floor() as usize).min(nx - 2);
        let fx = tx - ix0 as f64;

        let td = ((d_prev.clamp(d_min, d_max) - d_min) / (d_max - d_min)) * (nd - 1) as f64;
        let id0 = (td.floor() as usize).min(nd - 2);
        let fd = td - id0 as f64;

        let pol = &self.policy[t];
        let v00 = pol[ix0 * nd + id0];
        let v01 = pol[ix0 * nd + id0 + 1];
        let v10 = pol[(ix0 + 1) * nd + id0];
        let v11 = pol[(ix0 + 1) * nd + id0 + 1];
        (1.0 - fx) * ((1.0 - fd) * v00 + fd * v01) + fx * ((1.0 - fd) * v10 + fd * v11)
    }
}

#[inline]
fn interp_1d(grid: &[f64], vals: &[f64], x: f64) -> f64 {
    let n = grid.len();
    let x_min = grid[0];
    let x_max = grid[n - 1];
    if x <= x_min {
        return vals[0];
    }
    if x >= x_max {
        return vals[n - 1];
    }
    let t = ((x - x_min) / (x_max - x_min)) * (n - 1) as f64;
    let i0 = (t.floor() as usize).min(n - 2);
    let f = t - i0 as f64;
    vals[i0] * (1.0 - f) + vals[i0 + 1] * f
}

#[allow(clippy::too_many_arguments)]
pub fn solve_entropic_dp_gbm(
    s0: f64,
    k: f64,
    r: f64,
    sigma: f64,
    t_mat: f64,
    n_steps: usize,
    cost: f64,
    lambda: f64,
) -> DpResult {
    let nx = 81_usize;
    let nd = 81_usize;
    let x_min = -0.55_f64;
    let x_max = 0.55_f64;
    let d_min = -0.05_f64;
    let d_max = 1.05_f64;

    let x_grid: Vec<f64> = (0..nx)
        .map(|i| x_min + (x_max - x_min) * i as f64 / (nx - 1) as f64)
        .collect();
    let d_grid: Vec<f64> = (0..nd)
        .map(|i| d_min + (d_max - d_min) * i as f64 / (nd - 1) as f64)
        .collect();

    // 25-point normal quadrature on [-5, 5]
    let nz = 25_usize;
    let z_max = 5.0_f64;
    let dz = (2.0 * z_max) / (nz - 1) as f64;
    let mut z_pts = Vec::with_capacity(nz);
    let mut z_wts = Vec::with_capacity(nz);
    let mut w_sum = 0.0_f64;
    for iz in 0..nz {
        let z = -z_max + iz as f64 * dz;
        let w = if iz == 0 || iz == nz - 1 { 0.5 } else { 1.0 } * dz * norm_pdf(z);
        z_pts.push(z);
        z_wts.push(w);
        w_sum += w;
    }
    for w in &mut z_wts {
        *w /= w_sum;
    }

    let dt = t_mat / n_steps as f64;
    let drift = (r - 0.5 * sigma * sigma) * dt;
    let vol_dt = sigma * dt.sqrt();

    // Terminal value V_n(x, d) = (K * exp(x) - K)_+
    let mut v_next = vec![0.0_f64; nx * nd];
    for (ix, &x) in x_grid.iter().enumerate() {
        let pay = (k * x.exp() - k).max(0.0);
        for id in 0..nd {
            v_next[ix * nd + id] = pay;
        }
    }

    let mut policy = vec![vec![0.0_f64; nx * nd]; n_steps];
    let mut v_cur = vec![0.0_f64; nx * nd];
    let mut cont = vec![0.0_f64; nx * nd];
    let mut v_slice = vec![0.0_f64; nx];

    for t in (0..n_steps).rev() {
        // Compute continuation value C(ix, id_new) = (1/lambda) * ln E[exp(-lambda * d * dS + lambda * V_{t+1}(x', d))]
        for (id, &d) in d_grid.iter().enumerate() {
            for ix in 0..nx {
                v_slice[ix] = v_next[ix * nd + id];
            }
            for (ix, &x) in x_grid.iter().enumerate() {
                let s = k * x.exp();
                let mut max_arg = f64::NEG_INFINITY;
                let mut args = [0.0_f64; 25];
                for iz in 0..nz {
                    let dx = drift + vol_dt * z_pts[iz];
                    let s_next = s * dx.exp();
                    let v_n = interp_1d(&x_grid, &v_slice, x + dx);
                    let arg = lambda * (v_n - d * (s_next - s));
                    args[iz] = arg;
                    if arg > max_arg {
                        max_arg = arg;
                    }
                }
                let mut sum_e = 0.0_f64;
                for iz in 0..nz {
                    sum_e += z_wts[iz] * (args[iz] - max_arg).exp();
                }
                cont[ix * nd + id] = (max_arg + sum_e.ln()) / lambda;
            }
        }

        // Minimize over d_new for each (ix, id_prev)
        for (ix, &x) in x_grid.iter().enumerate() {
            let s = k * x.exp();
            for (id_prev, &d_prev) in d_grid.iter().enumerate() {
                let mut best_val = f64::INFINITY;
                let mut best_d = d_prev;
                for (id_new, &d_new) in d_grid.iter().enumerate() {
                    let val = cost * (d_new - d_prev).abs() * s + cont[ix * nd + id_new];
                    if val < best_val {
                        best_val = val;
                        best_d = d_new;
                    }
                }
                v_cur[ix * nd + id_prev] = best_val;
                policy[t][ix * nd + id_prev] = best_d;
            }
        }
        v_next.copy_from_slice(&v_cur);
    }

    // Read V_0(x = ln(s0/k), d_prev = 0)
    let x0 = (s0 / k).ln();
    let d0_idx = ((0.0 - d_min) / (d_max - d_min) * (nd - 1) as f64).round() as usize;
    let v0_slice: Vec<f64> = (0..nx).map(|ix| v_next[ix * nd + d0_idx]).collect();
    let indifference_price = interp_1d(&x_grid, &v0_slice, x0);

    DpResult {
        indifference_price,
        x_grid,
        d_grid,
        policy,
        n_steps,
        k,
    }
}
