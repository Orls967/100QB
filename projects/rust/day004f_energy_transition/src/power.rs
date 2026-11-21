//! Three-regime electricity spot price model with compound Poisson jumps and 60-day burn-in.

use crate::linalg;
use crate::parallel;
use crate::rng::{derive_seed, Rng};

pub const NUM_REGIMES: usize = 3;
pub const HOURS_PER_DAY: usize = 24;
pub const REGIME_NORMAL: usize = 0;
pub const REGIME_SCARCITY: usize = 1;
pub const REGIME_OVERSUPPLY: usize = 2;
pub const REGIME_NAMES: [&str; 3] = ["normal", "scarcity", "oversupply"];

#[derive(Clone, Copy, Debug)]
pub struct JumpParams {
    pub lambda: f64,
    pub mean_jump: f64,
    pub sign: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct PowerModelSpec {
    pub transition_matrix: [[f64; 3]; 3],
    pub mu: [f64; 3],
    pub phi: [f64; 3],
    pub sigma: [f64; 3],
    pub jumps: [JumpParams; 3],
    pub seasonal_amp: f64,
    pub seasonal_peak_hour: f64,
    pub burn_in_days: usize,
}

impl PowerModelSpec {
    pub fn default_spec() -> Self {
        Self {
            transition_matrix: [
                [0.985, 0.010, 0.005],
                [0.250, 0.740, 0.010],
                [0.150, 0.050, 0.800],
            ],
            mu: [45.0, 90.0, -25.0],
            phi: [0.82, 0.60, 0.75],
            sigma: [18.0, 100.0, 40.0],
            jumps: [
                JumpParams {
                    lambda: 0.002,
                    mean_jump: 0.0,
                    sign: 0.0,
                },
                JumpParams {
                    lambda: 0.08,
                    mean_jump: 120.0,
                    sign: 1.0,
                },
                JumpParams {
                    lambda: 0.05,
                    mean_jump: 60.0,
                    sign: -1.0,
                },
            ],
            seasonal_amp: 15.0,
            seasonal_peak_hour: 19.0,
            burn_in_days: 60,
        }
    }

    #[inline]
    pub fn seasonal_mean(&self, r: usize, hour: usize) -> f64 {
        let h = (hour % HOURS_PER_DAY) as f64;
        self.mu[r]
            + self.seasonal_amp
                * (2.0 * std::f64::consts::PI * (h - self.seasonal_peak_hour) / 24.0).cos()
    }

    #[inline]
    pub fn step(
        &self,
        r_now: usize,
        h_now: usize,
        p_now: f64,
        rg: &mut Rng,
    ) -> (usize, f64, usize, f64) {
        let row = self.transition_matrix[r_now];
        let u = rg.uniform();
        let r_next = if u < row[0] {
            0
        } else if u < row[0] + row[1] {
            1
        } else {
            2
        };

        let jp = self.jumps[r_next];
        let n_jumps = rg.poisson(jp.lambda);
        let j_val = if jp.sign == 0.0 || n_jumps == 0 {
            0.0
        } else {
            let mut sum = 0.0;
            for _ in 0..n_jumps {
                sum += rg.exp(jp.mean_jump);
            }
            jp.sign * sum
        };

        let m_now = self.seasonal_mean(r_next, h_now);
        let m_next = self.seasonal_mean(r_next, h_now + 1);
        let z = rg.normal();
        let p_next = m_next + self.phi[r_next] * (p_now - m_now) + self.sigma[r_next] * z + j_val;
        (r_next, p_next, n_jumps, j_val)
    }

    pub fn one_step_conditional_moments(
        &self,
        r_now: usize,
        h_now: usize,
        p_now: f64,
    ) -> (f64, f64) {
        let row = self.transition_matrix[r_now];
        let mut cond_means = [0.0; 3];
        let mut cond_vars = [0.0; 3];
        let mut mean_total = 0.0;
        for r_next in 0..3 {
            let (jm, jv) = jump_theoretical_moments(self.jumps[r_next]);
            let m_now = self.seasonal_mean(r_next, h_now);
            let m_next = self.seasonal_mean(r_next, h_now + 1);
            cond_means[r_next] = m_next + self.phi[r_next] * (p_now - m_now) + jm;
            cond_vars[r_next] = self.sigma[r_next] * self.sigma[r_next] + jv;
            mean_total += row[r_next] * cond_means[r_next];
        }
        let mut var_total = 0.0;
        for r_next in 0..3 {
            let diff = cond_means[r_next] - mean_total;
            var_total += row[r_next] * (cond_vars[r_next] + diff * diff);
        }
        (mean_total, var_total)
    }
}

pub fn jump_theoretical_moments(jp: JumpParams) -> (f64, f64) {
    if jp.sign == 0.0 || jp.mean_jump <= 0.0 {
        return (0.0, 0.0);
    }
    let m = jp.sign * jp.lambda * jp.mean_jump;
    let v = 2.0 * jp.lambda * jp.mean_jump * jp.mean_jump;
    (m, v)
}

pub fn stationary_distribution(p: &[[f64; 3]; 3]) -> Result<[f64; 3], String> {
    let a = vec![
        vec![p[0][0] - 1.0, p[1][0], p[2][0]],
        vec![p[0][1], p[1][1] - 1.0, p[2][1]],
        vec![1.0, 1.0, 1.0],
    ];
    let rhs = vec![0.0, 0.0, 1.0];
    let sol = linalg::solve_lu(&a, &rhs)?;
    Ok([sol[0], sol[1], sol[2]])
}

pub fn stationary_residual(p: &[[f64; 3]; 3], pi: [f64; 3]) -> f64 {
    let mut max_err = 0.0f64;
    for j in 0..3 {
        let mut s = 0.0;
        for i in 0..3 {
            s += pi[i] * p[i][j];
        }
        max_err = max_err.max((s - pi[j]).abs());
    }
    max_err.max((pi.iter().sum::<f64>() - 1.0).abs())
}

pub fn sample_stationary_regime(pi: [f64; 3], rg: &mut Rng) -> usize {
    let u = rg.uniform();
    if u < pi[0] {
        0
    } else if u < pi[0] + pi[1] {
        1
    } else {
        2
    }
}

#[derive(Clone, Debug)]
pub struct DayPath {
    pub day_idx: usize,
    pub prices: [f64; HOURS_PER_DAY],
    pub regimes: [usize; HOURS_PER_DAY],
    pub num_jumps: [usize; HOURS_PER_DAY],
    pub jump_vals: [f64; HOURS_PER_DAY],
}

pub fn simulate_eval_days(
    n_days: usize,
    n_workers: usize,
    spec: &PowerModelSpec,
    pi: [f64; 3],
    eval_seed: u64,
    salt: u64,
) -> Vec<DayPath> {
    parallel::run_chunks(n_days, n_workers, |d| {
        let mut rg = Rng::new(derive_seed(eval_seed, salt, d as u64));
        let mut r = sample_stationary_regime(pi, &mut rg);
        let mut p = spec.seasonal_mean(r, 0);
        for step in 0..(spec.burn_in_days * HOURS_PER_DAY) {
            (r, p, _, _) = spec.step(r, step, p, &mut rg);
        }
        let mut day = DayPath {
            day_idx: d,
            prices: [0.0; HOURS_PER_DAY],
            regimes: [0; HOURS_PER_DAY],
            num_jumps: [0; HOURS_PER_DAY],
            jump_vals: [0.0; HOURS_PER_DAY],
        };
        for h in 0..HOURS_PER_DAY {
            let (rn, pn, nj, jv) = spec.step(r, h, p, &mut rg);
            r = rn;
            p = pn;
            day.prices[h] = pn;
            day.regimes[h] = rn;
            day.num_jumps[h] = nj;
            day.jump_vals[h] = jv;
        }
        day
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SEED_DAY004E_VERIFY;

    #[test]
    fn test_stationary_and_moments() {
        let spec = PowerModelSpec::default_spec();
        let pi = stationary_distribution(&spec.transition_matrix).unwrap();
        assert!(stationary_residual(&spec.transition_matrix, pi) < 1e-13);

        let mut rg = Rng::new(derive_seed(SEED_DAY004E_VERIFY, 0x9901, 0));
        let n = 100_000;
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        for _ in 0..n {
            let (_, p_next, _, _) = spec.step(REGIME_SCARCITY, 23, 210.0, &mut rg);
            sum += p_next;
            sum2 += p_next * p_next;
        }
        let emp_m = sum / n as f64;
        let emp_v = (sum2 - n as f64 * emp_m * emp_m) / (n - 1) as f64;
        let (tm, tv) = spec.one_step_conditional_moments(REGIME_SCARCITY, 23, 210.0);
        assert!((emp_m - tm).abs() < 4.0 * (tv / n as f64).sqrt());
        assert!((emp_v - tv).abs() / tv < 0.08);
    }
}
