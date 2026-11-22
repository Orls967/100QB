//! Forward Kolmogorov ODE solver for inventory distribution, expected executions, and inventory variance.

use crate::glt::QuoteTable;

#[derive(Clone, Debug)]
pub struct KolmogorovResult {
    pub expected_executions: f64,
    pub mean_inventory: f64,
    pub var_inventory: f64,
    pub boundary_prob_terminal: f64,
    pub terminal_dist: Vec<f64>,
}

pub fn solve_forward_kolmogorov(table: &QuoteTable, n_steps: usize) -> KolmogorovResult {
    let p = table.params;
    let dim = p.dim();
    let dt = p.t_horizon / n_steps as f64;
    let mut prob = vec![0.0; dim];
    prob[p.idx(0)] = 1.0;

    let rhs_and_rate = |t: f64, state: &[f64], dstate: &mut [f64]| -> f64 {
        dstate.fill(0.0);
        let mut total_rate = 0.0;
        for i in 0..dim {
            let q = i as isize - p.q_max;
            let (db, da) = table.interp_quotes(t, q);
            let lb = if q < p.q_max {
                p.a_intens * (-p.k * db).exp()
            } else {
                0.0
            };
            let la = if q > -p.q_max {
                p.a_intens * (-p.k * da).exp()
            } else {
                0.0
            };
            let pi = state[i];
            total_rate += pi * (lb + la);
            if i + 1 < dim {
                dstate[i + 1] += lb * pi;
                dstate[i] -= lb * pi;
            }
            if i > 0 {
                dstate[i - 1] += la * pi;
                dstate[i] -= la * pi;
            }
        }
        total_rate
    };

    let mut k1 = vec![0.0; dim];
    let mut k2 = vec![0.0; dim];
    let mut k3 = vec![0.0; dim];
    let mut k4 = vec![0.0; dim];
    let mut tmp = vec![0.0; dim];
    let mut exp_exec = 0.0;

    for step in 0..n_steps {
        let t = step as f64 * dt;
        let r1 = rhs_and_rate(t, &prob, &mut k1);
        for i in 0..dim {
            tmp[i] = prob[i] + 0.5 * dt * k1[i];
        }
        let r2 = rhs_and_rate(t + 0.5 * dt, &tmp, &mut k2);
        for i in 0..dim {
            tmp[i] = prob[i] + 0.5 * dt * k2[i];
        }
        let r3 = rhs_and_rate(t + 0.5 * dt, &tmp, &mut k3);
        for i in 0..dim {
            tmp[i] = prob[i] + dt * k3[i];
        }
        let r4 = rhs_and_rate(t + dt, &tmp, &mut k4);
        for i in 0..dim {
            prob[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        exp_exec += (dt / 6.0) * (r1 + 2.0 * r2 + 2.0 * r3 + r4);
    }

    let mut mean_q = 0.0;
    let mut mean_q2 = 0.0;
    for (i, &pi) in prob.iter().enumerate() {
        let q = (i as isize - p.q_max) as f64;
        mean_q += q * pi;
        mean_q2 += q * q * pi;
    }

    KolmogorovResult {
        expected_executions: exp_exec,
        mean_inventory: mean_q,
        var_inventory: mean_q2 - mean_q * mean_q,
        boundary_prob_terminal: prob[0] + prob[dim - 1],
        terminal_dist: prob,
    }
}
