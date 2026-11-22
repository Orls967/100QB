use crate::conc::ConcRange;
use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};
use crate::stats;

#[derive(Clone, Debug)]
pub struct StrategyOutcome {
    pub name: String,
    pub mean_pnl: f64,
    pub sd_pnl: f64,
    pub se_pnl: f64,
    pub sharpe: f64,
    pub frac_in_range: f64,
    pub mean_gas: f64,
    pub mean_rebalances: f64,
    pub pnl: Vec<f64>,
}

/// Compares 5 LP/hedging strategies over `t_days = 30` with CRN:
/// (1) Hodl, (2) Full-range LP, (3) Passive Conc +-10%, (4) Rebalanced Conc +-10% (every `rebal_days`, gas `g`),
/// (5) Delta-hedged Conc +-10%.
#[allow(clippy::too_many_arguments)]
pub fn compare_lp_strategies(
    sigma: f64,
    phi: f64,
    daily_vol_tvl: f64,
    gas_cost: f64,
    rebal_days: usize,
    t_days: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<StrategyOutcome> {
    let steps_per_day = 24_usize; // hourly
    let n_steps = t_days * steps_per_day;
    let dt = (t_days as f64 / 365.0) / n_steps as f64;
    let drift = -0.5 * sigma * sigma * dt;
    let vol_dt = sigma * dt.sqrt();

    let full = ConcRange::full_range();
    let conc_base = ConcRange::normalized(0.9, 1.1);
    let ce = conc_base.capital_efficiency();
    let rebal_every = rebal_days * steps_per_day;

    // Rows: [hodl, full, passive_conc, rebal_conc, hedged_conc, in_range_pas, in_range_reb, gas, n_reb]
    let rows = run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x5001, r as u64));
        let mut p = 1.0_f64;

        let mut fee_full = 0.0_f64;
        let mut fee_pas = 0.0_f64;
        let mut fee_reb = 0.0_f64;
        let mut in_pas = 0_usize;
        let mut in_reb = 0_usize;

        let mut cur_range = conc_base;
        let mut wealth_reb_adj = 0.0_f64;
        let mut gas_tot = 0.0_f64;
        let mut n_reb = 0.0_f64;

        let mut hedge_pnl = 0.0_f64;

        let fee_per_step_base = phi * daily_vol_tvl / steps_per_day as f64;

        for step in 0..n_steps {
            let p_next = p * (drift + vol_dt * rg.normal()).exp();

            // Fees accrue proportionally to capital efficiency when in range
            fee_full += fee_per_step_base;
            if p > conc_base.pa && p < conc_base.pb {
                in_pas += 1;
                fee_pas += fee_per_step_base * ce;
            }
            if p > cur_range.pa && p < cur_range.pb {
                in_reb += 1;
                fee_reb += fee_per_step_base * ce;
            }

            // Delta hedge on passive concentrated LP
            let d_lp = conc_base.reserves(p).0;
            hedge_pnl -= d_lp * (p_next - p);

            p = p_next;

            if (step + 1) % rebal_every == 0
                && (step + 1) < n_steps
                && (p <= cur_range.pa || p >= cur_range.pb)
            {
                let v_cur = cur_range.value(p);
                let new_r = ConcRange {
                    pa: 0.9 * p,
                    pb: 1.1 * p,
                    l: (v_cur / (2.0 * p.sqrt() - (0.9 * p).sqrt() - p / (1.1 * p).sqrt())),
                };
                wealth_reb_adj -= gas_cost;
                gas_tot += gas_cost;
                n_reb += 1.0;
                cur_range = new_r;
            }
        }

        let pnl_hodl = 0.5 * p + 0.5 - 1.0;
        let pnl_full = full.value(p) + fee_full - 1.0;
        let pnl_pas = conc_base.value(p) + fee_pas - 1.0;
        let pnl_reb = cur_range.value(p) + fee_reb + wealth_reb_adj - 1.0;
        let pnl_hedged = conc_base.value(p) + fee_pas + hedge_pnl - 1.0;

        (
            [pnl_hodl, pnl_full, pnl_pas, pnl_reb, pnl_hedged],
            in_pas as f64 / n_steps as f64,
            in_reb as f64 / n_steps as f64,
            gas_tot,
            n_reb,
        )
    });

    let names = [
        "hodl",
        "full_range_lp",
        "passive_conc_10pct",
        "rebalanced_conc_10pct",
        "delta_hedged_conc_10pct",
    ];
    let mut out = Vec::with_capacity(5);
    for (idx, &nm) in names.iter().enumerate() {
        let pnl: Vec<f64> = rows.iter().map(|x| x.0[idx]).collect();
        let m = stats::mean(&pnl);
        let sd = stats::sd(&pnl);
        let se = stats::se(&pnl);
        let sh = stats::sharpe(&pnl);
        let (fin, mgas, mreb) = match idx {
            0 => (1.0, 0.0, 0.0),
            1 => (1.0, 0.0, 0.0),
            2 => (
                stats::mean(&rows.iter().map(|x| x.1).collect::<Vec<_>>()),
                0.0,
                0.0,
            ),
            3 => (
                stats::mean(&rows.iter().map(|x| x.2).collect::<Vec<_>>()),
                stats::mean(&rows.iter().map(|x| x.3).collect::<Vec<_>>()),
                stats::mean(&rows.iter().map(|x| x.4).collect::<Vec<_>>()),
            ),
            _ => (
                stats::mean(&rows.iter().map(|x| x.1).collect::<Vec<_>>()),
                0.0,
                0.0,
            ),
        };
        out.push(StrategyOutcome {
            name: nm.to_string(),
            mean_pnl: m,
            sd_pnl: sd,
            se_pnl: se,
            sharpe: sh,
            frac_in_range: fin,
            mean_gas: mgas,
            mean_rebalances: mreb,
            pnl,
        });
    }
    out
}
