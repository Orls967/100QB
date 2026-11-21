//! Log-linear discount factor interpolation, natural cubic spline on yields, and par-noise stability study.

use crate::bootstrap;
use crate::linalg;
use crate::ns::NsParams;
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;

pub struct LogLinearP {
    taus: Vec<f64>,
    ln_p: Vec<f64>,
}

impl LogLinearP {
    pub fn new(taus: &[f64], dfs: &[f64]) -> Result<Self, String> {
        if taus.is_empty() || taus.len() != dfs.len() {
            return Err("LogLinearP: invalid dimensions".into());
        }
        let mut t = Vec::with_capacity(taus.len() + 1);
        let mut lp = Vec::with_capacity(dfs.len() + 1);
        t.push(0.0);
        lp.push(0.0);
        for i in 0..taus.len() {
            if taus[i] <= t[i] || dfs[i] <= 0.0 {
                return Err(format!("LogLinearP: invalid knot/df at {i}"));
            }
            t.push(taus[i]);
            lp.push(dfs[i].ln());
        }
        Ok(Self { taus: t, ln_p: lp })
    }

    pub fn eval(&self, tau: f64) -> (f64, f64) {
        let n = self.taus.len();
        if tau <= 0.0 {
            let f0 = -(self.ln_p[1] - self.ln_p[0]) / (self.taus[1] - self.taus[0]);
            return (f0, f0);
        }
        let pos = self.taus.partition_point(|&x| x < tau);
        let idx = pos.saturating_sub(1).min(n - 2);
        let dt = self.taus[idx + 1] - self.taus[idx];
        let slope = (self.ln_p[idx + 1] - self.ln_p[idx]) / dt;
        let ln_p_tau = self.ln_p[idx] + slope * (tau - self.taus[idx]);
        (-ln_p_tau / tau, -slope)
    }
}

pub struct NaturalCubicSplineY {
    taus: Vec<f64>,
    yields: Vec<f64>,
    m2: Vec<f64>,
}

impl NaturalCubicSplineY {
    pub fn new(taus: &[f64], yields: &[f64]) -> Result<Self, String> {
        let n = taus.len();
        if n < 3 || yields.len() != n {
            return Err("NaturalCubicSplineY: need at least 3 knots".into());
        }
        let mut h = Vec::with_capacity(n - 1);
        for i in 0..(n - 1) {
            let hi = taus[i + 1] - taus[i];
            if hi <= 0.0 {
                return Err("NaturalCubicSplineY: knots must be strictly increasing".into());
            }
            h.push(hi);
        }
        let m = n - 2;
        let mut diag = vec![0.0; m];
        let mut rhs = vec![0.0; m];
        let mut sub = vec![0.0; m.saturating_sub(1)];
        let mut sup = vec![0.0; m.saturating_sub(1)];
        for row in 0..m {
            let i = row + 1;
            diag[row] = 2.0 * (h[i - 1] + h[i]);
            rhs[row] =
                6.0 * ((yields[i + 1] - yields[i]) / h[i] - (yields[i] - yields[i - 1]) / h[i - 1]);
            if row > 0 {
                sub[row - 1] = h[i - 1];
            }
            if row + 1 < m {
                sup[row] = h[i];
            }
        }
        let inner_m2 = linalg::solve_tridiagonal(&sub, &diag, &sup, &rhs)?;
        let mut m2 = vec![0.0; n];
        m2[1..(n - 1)].copy_from_slice(&inner_m2);
        Ok(Self {
            taus: taus.to_vec(),
            yields: yields.to_vec(),
            m2,
        })
    }

    pub fn eval(&self, tau: f64) -> (f64, f64) {
        let n = self.taus.len();
        let pos = self.taus.partition_point(|&x| x < tau);
        let idx = pos.saturating_sub(1).min(n - 2);
        let hi = self.taus[idx + 1] - self.taus[idx];
        let a = (self.taus[idx + 1] - tau) / hi;
        let b = (tau - self.taus[idx]) / hi;
        let y = a * self.yields[idx]
            + b * self.yields[idx + 1]
            + (hi * hi / 6.0)
                * ((a * a * a - a) * self.m2[idx] + (b * b * b - b) * self.m2[idx + 1]);
        let dy = (self.yields[idx + 1] - self.yields[idx]) / hi
            + (hi / 6.0)
                * (-(3.0 * a * a - 1.0) * self.m2[idx] + (3.0 * b * b - 1.0) * self.m2[idx + 1]);
        (y, y + tau * dy)
    }
}

#[derive(Clone, Debug)]
pub struct MethodMetrics {
    pub neg_forward_mean_count: f64,
    pub neg_forward_max_count: usize,
    pub neg_forward_frac_reps: f64,
    pub forward_tv_mean: f64,
    pub forward_tv_sd: f64,
    pub forward_jump_rms_mean: f64,
    pub forward_jump_rms_sd: f64,
    pub forward_d2_rms_mean: f64,
    pub forward_rmse_mean: f64,
    pub forward_rmse_sd: f64,
    pub yield_rmse_mean: f64,
    pub yield_rmse_sd: f64,
}

#[derive(Clone, Debug)]
pub struct NoiseExperimentSummary {
    pub num_reps: usize,
    pub sigma_par: f64,
    pub num_knots: usize,
    pub num_eval_grid: usize,
    pub log_linear_p: MethodMetrics,
    pub cubic_spline: MethodMetrics,
    pub sample_eval_taus: Vec<f64>,
    pub sample_true_y: Vec<f64>,
    pub sample_true_f: Vec<f64>,
    pub sample_loglin_y: Vec<f64>,
    pub sample_loglin_f: Vec<f64>,
    pub sample_spline_y: Vec<f64>,
    pub sample_spline_f: Vec<f64>,
}

fn eval_metrics_on_grid<F>(
    eval_fn: F,
    eval_taus: &[f64],
    true_y: &[f64],
    true_f: &[f64],
) -> (usize, f64, f64, f64, f64, f64, Vec<f64>, Vec<f64>)
where
    F: Fn(f64) -> (f64, f64),
{
    let m = eval_taus.len();
    let mut ys = Vec::with_capacity(m);
    let mut fs = Vec::with_capacity(m);
    let mut neg_cnt = 0usize;
    for &tau in eval_taus {
        let (y, f) = eval_fn(tau);
        if f < 0.0 {
            neg_cnt += 1;
        }
        ys.push(y);
        fs.push(f);
    }
    let mut tv = 0.0;
    let mut j_sq = 0.0;
    for j in 1..m {
        let df = fs[j] - fs[j - 1];
        tv += df.abs();
        j_sq += df * df;
    }
    let jump_rms = (j_sq / (m - 1) as f64).sqrt();
    let mut d2_sq = 0.0;
    for j in 1..(m - 1) {
        let d2 = fs[j + 1] - 2.0 * fs[j] + fs[j - 1];
        d2_sq += d2 * d2;
    }
    let d2_rms = (d2_sq / (m - 2) as f64).sqrt();
    let mut f_err2 = 0.0;
    let mut y_err2 = 0.0;
    for j in 0..m {
        let df = fs[j] - true_f[j];
        let dy = ys[j] - true_y[j];
        f_err2 += df * df;
        y_err2 += dy * dy;
    }
    let f_rmse = (f_err2 / m as f64).sqrt();
    let y_rmse = (y_err2 / m as f64).sqrt();
    (neg_cnt, tv, jump_rms, d2_rms, f_rmse, y_rmse, ys, fs)
}

pub fn run_noise_experiment(
    curve: &NsParams,
    n_reps: usize,
    sigma_par: f64,
    n_workers: usize,
    base_seed: u64,
) -> NoiseExperimentSummary {
    let knots = bootstrap::semi_annual_grid(30.0);
    let (true_pars, _, _) = bootstrap::par_rates_from_curve(curve, &knots, 0.5);
    let mut eval_taus = Vec::new();
    let mut t = 0.5;
    while t <= 30.0001 {
        eval_taus.push(t);
        t += 0.1;
    }
    let true_y: Vec<f64> = eval_taus.iter().map(|&x| curve.yield_at(x)).collect();
    let true_f: Vec<f64> = eval_taus.iter().map(|&x| curve.forward(x)).collect();

    struct Row {
        ll_neg: usize,
        ll_tv: f64,
        ll_j: f64,
        ll_d2: f64,
        ll_fr: f64,
        ll_yr: f64,
        sp_neg: usize,
        sp_tv: f64,
        sp_j: f64,
        sp_d2: f64,
        sp_fr: f64,
        sp_yr: f64,
        ll_y: Vec<f64>,
        ll_f: Vec<f64>,
        sp_y: Vec<f64>,
        sp_f: Vec<f64>,
    }

    let rows = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(base_seed, 0x1A01, r as u64));
        let noisy_pars: Vec<f64> = true_pars
            .iter()
            .map(|&p| p + sigma_par * rg.normal())
            .collect();
        let (dfs, zys) = bootstrap::bootstrap_from_par_rates(&knots, &noisy_pars, 0.5)
            .expect("bootstrap failed on 0.5bp noise");
        let ll = LogLinearP::new(&knots, &dfs).unwrap();
        let sp = NaturalCubicSplineY::new(&knots, &zys).unwrap();
        let (ll_neg, ll_tv, ll_j, ll_d2, ll_fr, ll_yr, ll_y, ll_f) =
            eval_metrics_on_grid(|tau| ll.eval(tau), &eval_taus, &true_y, &true_f);
        let (sp_neg, sp_tv, sp_j, sp_d2, sp_fr, sp_yr, sp_y, sp_f) =
            eval_metrics_on_grid(|tau| sp.eval(tau), &eval_taus, &true_y, &true_f);
        Row {
            ll_neg,
            ll_tv,
            ll_j,
            ll_d2,
            ll_fr,
            ll_yr,
            sp_neg,
            sp_tv,
            sp_j,
            sp_d2,
            sp_fr,
            sp_yr,
            ll_y: if r == 0 { ll_y } else { Vec::new() },
            ll_f: if r == 0 { ll_f } else { Vec::new() },
            sp_y: if r == 0 { sp_y } else { Vec::new() },
            sp_f: if r == 0 { sp_f } else { Vec::new() },
        }
    });

    let summarize =
        |neg: Vec<usize>, tv: Vec<f64>, j: Vec<f64>, d2: Vec<f64>, fr: Vec<f64>, yr: Vec<f64>| {
            let max_neg = *neg.iter().max().unwrap_or(&0);
            let frac_reps = neg.iter().filter(|&&c| c > 0).count() as f64 / n_reps.max(1) as f64;
            let mean_neg = neg.iter().sum::<usize>() as f64 / n_reps.max(1) as f64;
            MethodMetrics {
                neg_forward_mean_count: mean_neg,
                neg_forward_max_count: max_neg,
                neg_forward_frac_reps: frac_reps,
                forward_tv_mean: stats::mean(&tv),
                forward_tv_sd: stats::sd(&tv),
                forward_jump_rms_mean: stats::mean(&j),
                forward_jump_rms_sd: stats::sd(&j),
                forward_d2_rms_mean: stats::mean(&d2),
                forward_rmse_mean: stats::mean(&fr),
                forward_rmse_sd: stats::sd(&fr),
                yield_rmse_mean: stats::mean(&yr),
                yield_rmse_sd: stats::sd(&yr),
            }
        };

    let ll_metrics = summarize(
        rows.iter().map(|r| r.ll_neg).collect(),
        rows.iter().map(|r| r.ll_tv).collect(),
        rows.iter().map(|r| r.ll_j).collect(),
        rows.iter().map(|r| r.ll_d2).collect(),
        rows.iter().map(|r| r.ll_fr).collect(),
        rows.iter().map(|r| r.ll_yr).collect(),
    );
    let sp_metrics = summarize(
        rows.iter().map(|r| r.sp_neg).collect(),
        rows.iter().map(|r| r.sp_tv).collect(),
        rows.iter().map(|r| r.sp_j).collect(),
        rows.iter().map(|r| r.sp_d2).collect(),
        rows.iter().map(|r| r.sp_fr).collect(),
        rows.iter().map(|r| r.sp_yr).collect(),
    );

    NoiseExperimentSummary {
        num_reps: n_reps,
        sigma_par,
        num_knots: knots.len(),
        num_eval_grid: eval_taus.len(),
        log_linear_p: ll_metrics,
        cubic_spline: sp_metrics,
        sample_eval_taus: eval_taus,
        sample_true_y: true_y,
        sample_true_f: true_f,
        sample_loglin_y: rows[0].ll_y.clone(),
        sample_loglin_f: rows[0].ll_f.clone(),
        sample_spline_y: rows[0].sp_y.clone(),
        sample_spline_f: rows[0].sp_f.clone(),
    }
}
