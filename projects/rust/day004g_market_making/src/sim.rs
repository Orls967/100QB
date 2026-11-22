//! High-speed deterministic CRN forward simulator for market making evaluation.

use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;
use crate::strategies::Strategy;

#[derive(Clone, Copy, Debug)]
pub struct MarketSpec {
    pub t_horizon: f64,
    pub sigma: f64,
    pub mu: f64,
    pub a_intens: f64,
    pub k: f64,
    pub gamma: f64,
    pub q_max: isize,
    pub eps_adverse: f64,
}

#[derive(Clone, Debug)]
pub struct SimSummary {
    pub n_paths: usize,
    pub dt: f64,
    pub mean_wealth: f64,
    pub sd_wealth: f64,
    pub sharpe: f64,
    pub mean_inventory: f64,
    pub sd_inventory: f64,
    pub var_inventory: f64,
    pub mean_util: f64,
    pub se_util: f64,
    pub certainty_equiv: f64,
    pub se_ce: f64,
    pub mean_executions: f64,
    pub boundary_touch_frac: f64,
    pub wealths: Vec<f64>,
}

pub fn run_mc_sim(
    mkt: MarketSpec,
    strat: &Strategy<'_>,
    n_paths: usize,
    n_steps: usize,
    n_workers: usize,
    seed: u64,
    salt: u64,
) -> SimSummary {
    let dt = mkt.t_horizon / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let drift_step = mkt.mu * dt;
    let vol_step = mkt.sigma * sqrt_dt;

    struct PathOut {
        wealth: f64,
        q_final: isize,
        execs: usize,
        touched_bound: bool,
    }

    let rows = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, salt, r as u64));
        let mut s = 0.0_f64;
        let mut x = 0.0_f64;
        let mut q = 0_isize;
        let mut execs = 0_usize;
        let mut touched = false;

        for step in 0..n_steps {
            let t = step as f64 * dt;
            let (db, da) = strat.quotes(t, q);

            let z = rg.normal();
            let u_b = rg.uniform();
            let u_a = rg.uniform();

            if q < mkt.q_max && db.is_finite() {
                let prob_b = (mkt.a_intens * (-mkt.k * db).exp() * dt).min(1.0);
                if u_b < prob_b {
                    x -= s - db;
                    q += 1;
                    execs += 1;
                    s -= mkt.eps_adverse;
                }
            }
            if q > -mkt.q_max && da.is_finite() {
                let prob_a = (mkt.a_intens * (-mkt.k * da).exp() * dt).min(1.0);
                if u_a < prob_a {
                    x += s + da;
                    q -= 1;
                    execs += 1;
                    s += mkt.eps_adverse;
                }
            }
            if q == mkt.q_max || q == -mkt.q_max {
                touched = true;
            }
            s += drift_step + vol_step * z;
        }

        PathOut {
            wealth: x + (q as f64) * s,
            q_final: q,
            execs,
            touched_bound: touched,
        }
    });

    let wealths: Vec<f64> = rows.iter().map(|r| r.wealth).collect();
    let qs: Vec<f64> = rows.iter().map(|r| r.q_final as f64).collect();
    let exs: Vec<f64> = rows.iter().map(|r| r.execs as f64).collect();
    let bnd_cnt = rows.iter().filter(|r| r.touched_bound).count();

    let mw = stats::mean(&wealths);
    let sw = stats::sd(&wealths);
    let mq = stats::mean(&qs);
    let sq = stats::sd(&qs);
    let (mu, seu, ce, se_ce) = stats::cara_utility_and_ce(&wealths, mkt.gamma);

    SimSummary {
        n_paths,
        dt,
        mean_wealth: mw,
        sd_wealth: sw,
        sharpe: if sw > 0.0 { mw / sw } else { 0.0 },
        mean_inventory: mq,
        sd_inventory: sq,
        var_inventory: sq * sq,
        mean_util: mu,
        se_util: seu,
        certainty_equiv: ce,
        se_ce,
        mean_executions: stats::mean(&exs),
        boundary_touch_frac: bnd_cnt as f64 / n_paths.max(1) as f64,
        wealths,
    }
}

#[derive(Clone, Debug)]
pub struct ExamplePath {
    pub t: Vec<f64>,
    pub s: Vec<f64>,
    pub q: Vec<f64>,
    pub mtm_wealth: Vec<f64>,
}

pub fn simulate_single_example_path(
    mkt: MarketSpec,
    strat: &Strategy<'_>,
    n_steps: usize,
    seed: u64,
) -> ExamplePath {
    let dt = mkt.t_horizon / n_steps as f64;
    let vol_step = mkt.sigma * dt.sqrt();
    let mut rg = Rng::new(derive_seed(seed, 0xE101, 0));
    let mut s = 0.0_f64;
    let mut x = 0.0_f64;
    let mut q = 0_isize;

    let mut ts = Vec::with_capacity(n_steps + 1);
    let mut ss = Vec::with_capacity(n_steps + 1);
    let mut qs = Vec::with_capacity(n_steps + 1);
    let mut ws = Vec::with_capacity(n_steps + 1);
    ts.push(0.0);
    ss.push(s);
    qs.push(0.0);
    ws.push(0.0);

    for step in 0..n_steps {
        let t = step as f64 * dt;
        let (db, da) = strat.quotes(t, q);
        let z = rg.normal();
        let u_b = rg.uniform();
        let u_a = rg.uniform();
        if q < mkt.q_max
            && db.is_finite()
            && u_b < (mkt.a_intens * (-mkt.k * db).exp() * dt).min(1.0)
        {
            x -= s - db;
            q += 1;
            s -= mkt.eps_adverse;
        }
        if q > -mkt.q_max
            && da.is_finite()
            && u_a < (mkt.a_intens * (-mkt.k * da).exp() * dt).min(1.0)
        {
            x += s + da;
            q -= 1;
            s += mkt.eps_adverse;
        }
        s += mkt.mu * dt + vol_step * z;
        ts.push((step + 1) as f64 * dt);
        ss.push(s);
        qs.push(q as f64);
        ws.push(x + (q as f64) * s);
    }
    ExamplePath {
        t: ts,
        s: ss,
        q: qs,
        mtm_wealth: ws,
    }
}
