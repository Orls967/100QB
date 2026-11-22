//! Bagian 3 misspecification studies:
//! (1) Adverse selection (eps in {0, 0.02, 0.05, 0.1}) + shifted vs exact AS-HJB recovery,
//! (2) Wrong volatility (0.5*sigma, sigma, 1.5*sigma),
//! (3) Wrong intensity (A or k +/-30%, comparing GLT vs AS classic),
//! (4) Non-zero drift (mu in {-0.2, -0.05, 0, 0.05, 0.2}).

use crate::glt::{GltParams, QuoteTable};
use crate::hjb_fd;
use crate::sim::{run_mc_sim, MarketSpec, SimSummary};
use crate::strategies::Strategy;

#[derive(Clone, Debug)]
pub struct AdverseRow {
    pub eps: f64,
    pub unadjusted_glt: SimSummary,
    pub const_shifted_glt: SimSummary,
    pub exact_as_hjb: SimSummary,
}

#[derive(Clone, Debug)]
pub struct SigmaMisspecRow {
    pub sigma_mult: f64,
    pub sigma_model: f64,
    pub summary: SimSummary,
}

#[derive(Clone, Debug)]
pub struct IntensityMisspecRow {
    pub label: &'static str,
    pub a_true: f64,
    pub k_true: f64,
    pub glt_summary: SimSummary,
    pub as_summary: SimSummary,
}

#[derive(Clone, Debug)]
pub struct DriftMisspecRow {
    pub mu_true: f64,
    pub summary: SimSummary,
}

#[derive(Clone, Debug)]
pub struct MisspecResults {
    pub adverse: Vec<AdverseRow>,
    pub wrong_sigma: Vec<SigmaMisspecRow>,
    pub wrong_intensity: Vec<IntensityMisspecRow>,
    pub wrong_drift: Vec<DriftMisspecRow>,
}

pub fn run_misspec_studies(
    base_p: GltParams,
    base_table: &QuoteTable,
    n_paths: usize,
    n_steps: usize,
    n_workers: usize,
    seed: u64,
) -> MisspecResults {
    let base_mkt = MarketSpec {
        t_horizon: base_p.t_horizon,
        sigma: base_p.sigma,
        mu: 0.0,
        a_intens: base_p.a_intens,
        k: base_p.k,
        gamma: base_p.gamma,
        q_max: base_p.q_max,
        eps_adverse: 0.0,
    };

    // 1. Adverse selection
    let mut adverse = Vec::new();
    for &eps in &[0.0, 0.02, 0.05, 0.1] {
        let mut mkt = base_mkt;
        mkt.eps_adverse = eps;
        let unadj = run_mc_sim(
            mkt,
            &Strategy::GltOptimal { table: base_table },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3100,
        );
        let shifted = run_mc_sim(
            mkt,
            &Strategy::GltConstantShift {
                table: base_table,
                shift: eps,
            },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3100,
        );
        let (_, bid_as, ask_as) = hjb_fd::solve_nonlinear_hjb_theta(&base_p, 5e-4, eps);
        let exact_as = run_mc_sim(
            mkt,
            &Strategy::GltAsExact {
                p: base_p,
                bid_table: bid_as,
                ask_table: ask_as,
            },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3100,
        );
        adverse.push(AdverseRow {
            eps,
            unadjusted_glt: unadj,
            const_shifted_glt: shifted,
            exact_as_hjb: exact_as,
        });
    }

    // 2. Wrong sigma
    let mut wrong_sigma = Vec::new();
    for &mult in &[0.5, 1.0, 1.5] {
        let mut p_mod = base_p;
        p_mod.sigma = base_p.sigma * mult;
        let tab_mod = QuoteTable::build(p_mod, 1000);
        let sum = run_mc_sim(
            base_mkt,
            &Strategy::GltOptimal { table: &tab_mod },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3200,
        );
        wrong_sigma.push(SigmaMisspecRow {
            sigma_mult: mult,
            sigma_model: p_mod.sigma,
            summary: sum,
        });
    }

    // 3. Wrong intensity (A or k +/-30%, GLT vs AS)
    let int_cases = [
        ("A_-30pct", 0.7 * base_p.a_intens, base_p.k),
        ("A_+30pct", 1.3 * base_p.a_intens, base_p.k),
        ("k_-30pct", base_p.a_intens, 0.7 * base_p.k),
        ("k_+30pct", base_p.a_intens, 1.3 * base_p.k),
    ];
    let mut wrong_intensity = Vec::new();
    for (lbl, a_true, k_true) in int_cases {
        let mut mkt = base_mkt;
        mkt.a_intens = a_true;
        mkt.k = k_true;
        let glt_sum = run_mc_sim(
            mkt,
            &Strategy::GltOptimal { table: base_table },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3300,
        );
        let as_sum = run_mc_sim(
            mkt,
            &Strategy::AvellanedaStoikov { p: base_p },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3300,
        );
        wrong_intensity.push(IntensityMisspecRow {
            label: lbl,
            a_true,
            k_true,
            glt_summary: glt_sum,
            as_summary: as_sum,
        });
    }

    // 4. Non-zero drift
    let mut wrong_drift = Vec::new();
    for &mu in &[-0.2, -0.05, 0.0, 0.05, 0.2] {
        let mut mkt = base_mkt;
        mkt.mu = mu;
        let sum = run_mc_sim(
            mkt,
            &Strategy::GltOptimal { table: base_table },
            n_paths,
            n_steps,
            n_workers,
            seed,
            0x3400,
        );
        wrong_drift.push(DriftMisspecRow {
            mu_true: mu,
            summary: sum,
        });
    }

    MisspecResults {
        adverse,
        wrong_sigma,
        wrong_intensity,
        wrong_drift,
    }
}
