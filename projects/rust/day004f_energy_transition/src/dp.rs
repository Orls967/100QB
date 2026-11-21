//! Stochastic 24-hour Bellman DP, price grids (17, 33, expanded 38), transition kernel builder,
//! deterministic perfect-foresight DP, threshold policy, and policy agreement.

use crate::battery::{
    BatterySpec, ACTIONS, ACTION_CHARGE_1MW, ACTION_CHARGE_2MW, ACTION_DISCHARGE_09MW,
    ACTION_DISCHARGE_18MW, ACTION_IDLE, NUM_ACTIONS, NUM_SOC_STATES,
};
use crate::parallel;
use crate::power::{DayPath, PowerModelSpec, HOURS_PER_DAY, NUM_REGIMES};
use crate::rng::{derive_state_seed, Rng};
use crate::stats;

#[derive(Clone, Debug)]
pub struct PriceGrid {
    pub name: &'static str,
    pub num_bins: usize,
    pub edges: Vec<f64>,
    pub rep_prices: Vec<f64>,
    pub lower_tail_at: f64,
    pub upper_tail_at: f64,
}

impl PriceGrid {
    fn build(name: &'static str, edges: Vec<f64>, lower_rep: f64, upper_rep: f64) -> Self {
        let num_bins = edges.len() + 1;
        let mut rep_prices = Vec::with_capacity(num_bins);
        rep_prices.push(lower_rep);
        for b in 1..(num_bins - 1) {
            rep_prices.push(0.5 * (edges[b - 1] + edges[b]));
        }
        rep_prices.push(upper_rep);
        let lower_tail_at = edges[0];
        let upper_tail_at = edges[edges.len() - 1];
        Self {
            name,
            num_bins,
            edges,
            rep_prices,
            lower_tail_at,
            upper_tail_at,
        }
    }

    pub fn coarse_17() -> Self {
        Self::build(
            "coarse_17",
            vec![
                -500.0, -200.0, -100.0, -50.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 120.0, 200.0,
                350.0, 600.0, 1000.0, 2000.0,
            ],
            -600.0,
            2500.0,
        )
    }

    pub fn fine_33() -> Self {
        Self::build(
            "fine_33",
            vec![
                -500.0, -350.0, -200.0, -150.0, -100.0, -75.0, -50.0, -35.0, -20.0, -10.0, 0.0,
                10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 140.0, 160.0, 200.0,
                275.0, 350.0, 475.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0,
            ],
            -600.0,
            2500.0,
        )
    }

    pub fn expanded_tail_38() -> Self {
        Self::build(
            "expanded_tail_38",
            vec![
                -1000.0, -750.0, -500.0, -350.0, -200.0, -150.0, -100.0, -75.0, -50.0, -35.0,
                -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0,
                140.0, 160.0, 200.0, 275.0, 350.0, 475.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0,
                2500.0, 3500.0, 5000.0,
            ],
            -1150.0,
            5500.0,
        )
    }

    #[inline]
    pub fn bin_index(&self, price: f64) -> usize {
        if price < self.edges[0] {
            return 0;
        }
        let last = self.edges.len() - 1;
        if price >= self.edges[last] {
            return self.num_bins - 1;
        }
        self.edges.partition_point(|&e| e <= price)
    }
}

#[derive(Clone, Debug)]
pub struct TransitionKernel {
    pub grid: PriceGrid,
    pub draws_per_state: usize,
    pub num_states: usize,
    pub target_stride: usize,
    pub total_draws: usize,
    pub probs: Vec<f64>,
    pub lower_tail_mass: f64,
    pub upper_tail_mass: f64,
    pub lower_tail_count: usize,
    pub upper_tail_count: usize,
    pub lower_tail_mean_p: f64,
    pub upper_tail_mean_p: f64,
    pub build_seconds: f64,
}

impl TransitionKernel {
    #[inline]
    pub fn state_index(&self, h: usize, r: usize, b: usize) -> usize {
        (h * NUM_REGIMES + r) * self.grid.num_bins + b
    }
}

pub fn build_kernel(
    spec: &PowerModelSpec,
    grid: PriceGrid,
    draws_per_state: usize,
    n_workers: usize,
    kernel_seed: u64,
    salt: u64,
    crn_interior_to_33: bool,
) -> TransitionKernel {
    let t0 = std::time::Instant::now();
    let n_bins = grid.num_bins;
    let n_states = HOURS_PER_DAY * NUM_REGIMES * n_bins;
    let stride = NUM_REGIMES * n_bins;

    struct ChunkOut {
        probs: Vec<f64>,
        lo_cnt: usize,
        hi_cnt: usize,
        lo_sum: f64,
        hi_sum: f64,
    }

    let chunks = parallel::run_chunks(n_states, n_workers, |st_idx| {
        let h = st_idx / (NUM_REGIMES * n_bins);
        let rem = st_idx % (NUM_REGIMES * n_bins);
        let r = rem / n_bins;
        let b = rem % n_bins;

        let (use_salt, use_bin) = if crn_interior_to_33 && (3..=33).contains(&b) {
            (0x4E33, b - 2)
        } else {
            (salt, b)
        };
        let mut rg = Rng::new(derive_state_seed(
            kernel_seed,
            use_salt,
            h as u64,
            r as u64,
            use_bin as u64,
        ));
        let p_now = grid.rep_prices[b];
        let mut counts = vec![0usize; stride];
        let mut lo_cnt = 0usize;
        let mut hi_cnt = 0usize;
        let mut lo_sum = 0.0;
        let mut hi_sum = 0.0;

        for _ in 0..draws_per_state {
            let (r_next, p_next, _, _) = spec.step(r, h, p_now, &mut rg);
            let b_next = grid.bin_index(p_next);
            counts[r_next * n_bins + b_next] += 1;
            if b_next == 0 {
                lo_cnt += 1;
                lo_sum += p_next;
            } else if b_next + 1 == n_bins {
                hi_cnt += 1;
                hi_sum += p_next;
            }
        }
        let inv_n = 1.0 / draws_per_state as f64;
        let probs = counts.into_iter().map(|c| c as f64 * inv_n).collect();
        ChunkOut {
            probs,
            lo_cnt,
            hi_cnt,
            lo_sum,
            hi_sum,
        }
    });

    let mut all_probs = Vec::with_capacity(n_states * stride);
    let mut tot_lo = 0usize;
    let mut tot_hi = 0usize;
    let mut sum_lo = 0.0;
    let mut sum_hi = 0.0;
    for ch in chunks {
        all_probs.extend(ch.probs);
        tot_lo += ch.lo_cnt;
        tot_hi += ch.hi_cnt;
        sum_lo += ch.lo_sum;
        sum_hi += ch.hi_sum;
    }
    let tot_draws = n_states * draws_per_state;
    let lo_mean = if tot_lo > 0 {
        sum_lo / tot_lo as f64
    } else {
        grid.rep_prices[0]
    };
    let hi_mean = if tot_hi > 0 {
        sum_hi / tot_hi as f64
    } else {
        grid.rep_prices[n_bins - 1]
    };

    TransitionKernel {
        grid,
        draws_per_state,
        num_states: n_states,
        target_stride: stride,
        total_draws: tot_draws,
        probs: all_probs,
        lower_tail_mass: tot_lo as f64 / tot_draws as f64,
        upper_tail_mass: tot_hi as f64 / tot_draws as f64,
        lower_tail_count: tot_lo,
        upper_tail_count: tot_hi,
        lower_tail_mean_p: lo_mean,
        upper_tail_mean_p: hi_mean,
        build_seconds: t0.elapsed().as_secs_f64(),
    }
}

#[derive(Clone, Debug)]
pub struct DpPolicy {
    pub grid: PriceGrid,
    pub battery: BatterySpec,
    pub num_bins: usize,
    pub values: Vec<f64>,
    pub actions: Vec<usize>,
    pub max_bellman_residual: f64,
    pub terminal_max_error: f64,
    pub solve_seconds: f64,
}

impl DpPolicy {
    #[inline]
    pub fn value_at(&self, h: usize, r: usize, b: usize, s: usize) -> f64 {
        self.values[((h * NUM_REGIMES + r) * self.num_bins + b) * NUM_SOC_STATES + s]
    }

    #[inline]
    pub fn action_at(&self, h: usize, r: usize, b: usize, s: usize) -> usize {
        self.actions[((h * NUM_REGIMES + r) * self.num_bins + b) * NUM_SOC_STATES + s]
    }
}

pub fn solve_stochastic_dp(kernel: &TransitionKernel, bat: BatterySpec) -> DpPolicy {
    let t0 = std::time::Instant::now();
    let n_bins = kernel.grid.num_bins;
    let n_soc = bat.num_soc_states;
    let stride = kernel.target_stride;

    let mut values = vec![0.0; (HOURS_PER_DAY + 1) * NUM_REGIMES * n_bins * n_soc];
    let mut actions = vec![ACTION_IDLE; HOURS_PER_DAY * NUM_REGIMES * n_bins * n_soc];

    let term_vals: [f64; NUM_SOC_STATES] = std::array::from_fn(|s| bat.terminal_value(s));
    for r in 0..NUM_REGIMES {
        for b in 0..n_bins {
            let off = ((HOURS_PER_DAY * NUM_REGIMES + r) * n_bins + b) * n_soc;
            values[off..(off + n_soc)].copy_from_slice(&term_vals);
        }
    }

    let mut ev_cont = [0.0; NUM_SOC_STATES];
    for h in (0..HOURS_PER_DAY).rev() {
        let next_base = (h + 1) * stride * n_soc;
        let cur_base = h * stride * n_soc;
        for r in 0..NUM_REGIMES {
            for b in 0..n_bins {
                let st_idx = (h * NUM_REGIMES + r) * n_bins + b;
                let probs = &kernel.probs[st_idx * stride..(st_idx + 1) * stride];
                let rep_price = kernel.grid.rep_prices[b];

                ev_cont.fill(0.0);
                for (target_idx, &pr) in probs.iter().enumerate() {
                    if pr == 0.0 {
                        continue;
                    }
                    let v_slice = &values
                        [next_base + target_idx * n_soc..next_base + (target_idx + 1) * n_soc];
                    for s_next in 0..n_soc {
                        ev_cont[s_next] += pr * v_slice[s_next];
                    }
                }

                let act_cf: [f64; NUM_ACTIONS] =
                    std::array::from_fn(|a| bat.hourly_cashflow(rep_price, ACTIONS[a]));
                let tb_off = (r * n_bins + b) * n_soc;
                for s in 0..n_soc {
                    let mut best_val = f64::NEG_INFINITY;
                    let mut best_act = ACTION_IDLE;
                    for (a_idx, act) in ACTIONS.iter().enumerate() {
                        if let Some(s_next) = bat.next_soc_idx(s, *act) {
                            let cand = act_cf[a_idx] + ev_cont[s_next];
                            if cand > best_val {
                                best_val = cand;
                                best_act = act.id;
                            }
                        }
                    }
                    values[cur_base + tb_off + s] = best_val;
                    actions[cur_base + tb_off + s] = best_act;
                }
            }
        }
    }

    let mut pol = DpPolicy {
        grid: kernel.grid.clone(),
        battery: bat,
        num_bins: n_bins,
        values,
        actions,
        max_bellman_residual: 0.0,
        terminal_max_error: 0.0,
        solve_seconds: 0.0,
    };

    let mut max_term = 0.0f64;
    for r in 0..NUM_REGIMES {
        for b in 0..n_bins {
            for s in 0..n_soc {
                max_term = max_term
                    .max((pol.value_at(HOURS_PER_DAY, r, b, s) - bat.terminal_value(s)).abs());
            }
        }
    }
    let mut max_res = 0.0f64;
    for h in 0..HOURS_PER_DAY {
        for r in 0..NUM_REGIMES {
            for b in 0..n_bins {
                let st_idx = (h * NUM_REGIMES + r) * n_bins + b;
                let probs = &kernel.probs[st_idx * stride..(st_idx + 1) * stride];
                let rep_price = kernel.grid.rep_prices[b];
                for s in 0..n_soc {
                    let mut best_opt = f64::NEG_INFINITY;
                    for act in ACTIONS {
                        if let Some(s_next) = bat.next_soc_idx(s, act) {
                            let mut ev = 0.0;
                            for (target_idx, &pr) in probs.iter().enumerate() {
                                if pr != 0.0 {
                                    ev += pr
                                        * pol.value_at(
                                            h + 1,
                                            target_idx / n_bins,
                                            target_idx % n_bins,
                                            s_next,
                                        );
                                }
                            }
                            best_opt = best_opt.max(bat.hourly_cashflow(rep_price, act) + ev);
                        }
                    }
                    max_res = max_res.max((pol.value_at(h, r, b, s) - best_opt).abs());
                }
            }
        }
    }

    pol.terminal_max_error = max_term;
    pol.max_bellman_residual = max_res;
    pol.solve_seconds = t0.elapsed().as_secs_f64();
    pol
}

pub fn solve_deterministic_dp(
    price_path: &[f64],
    bat: BatterySpec,
) -> (f64, f64, Vec<usize>, Vec<usize>) {
    let h_len = price_path.len();
    let n_soc = bat.num_soc_states;
    let mut v = vec![vec![0.0; n_soc]; h_len + 1];
    let mut act_tab = vec![vec![ACTION_IDLE; n_soc]; h_len];
    for s in 0..n_soc {
        v[h_len][s] = bat.terminal_value(s);
    }
    for h in (0..h_len).rev() {
        let p = price_path[h];
        for s in 0..n_soc {
            let mut best_val = f64::NEG_INFINITY;
            let mut best_act = ACTION_IDLE;
            for act in ACTIONS {
                if let Some(s_next) = bat.next_soc_idx(s, act) {
                    let cand = bat.hourly_cashflow(p, act) + v[h + 1][s_next];
                    if cand > best_val {
                        best_val = cand;
                        best_act = act.id;
                    }
                }
            }
            v[h][s] = best_val;
            act_tab[h][s] = best_act;
        }
    }
    let mut actions = Vec::with_capacity(h_len);
    let mut soc_path = Vec::with_capacity(h_len + 1);
    let mut soc = bat.initial_soc_idx;
    soc_path.push(soc);
    for h in 0..h_len {
        let aid = act_tab[h][soc];
        actions.push(aid);
        soc = bat.next_soc_idx(soc, ACTIONS[aid]).unwrap();
        soc_path.push(soc);
    }
    let gross_total = v[0][bat.initial_soc_idx];
    let net_profit = gross_total - bat.terminal_value(bat.initial_soc_idx);
    (net_profit, gross_total, actions, soc_path)
}

#[derive(Clone, Debug)]
pub struct StrategyResult {
    pub name: String,
    pub daily_profits: Vec<f64>,
    pub daily_gross_values: Vec<f64>,
    pub daily_cashflows: Vec<f64>,
    pub daily_grid_throughput: Vec<f64>,
    pub daily_stored_throughput: Vec<f64>,
    pub daily_efc: Vec<f64>,
    pub daily_regret_vs_pf: Vec<f64>,
    pub mean_profit: f64,
    pub median_profit: f64,
    pub sd_profit: f64,
    pub mcse_profit: f64,
    pub var5_profit: f64,
    pub es5_profit: f64,
    pub mean_gross_value: f64,
    pub mean_cashflow: f64,
    pub mean_throughput: f64,
    pub mean_efc: f64,
    pub mean_regret_vs_pf: f64,
    pub mcse_regret_vs_pf: f64,
    pub max_regret_vs_pf: f64,
    pub action_counts: [usize; NUM_ACTIONS],
    pub total_profit_neg_price: f64,
    pub mean_profit_neg_price_per_day: f64,
    pub total_profit_scarcity_or_jump: f64,
    pub mean_profit_scarcity_or_jump_per_day: f64,
    pub soc_violation_count: usize,
    pub power_violation_count: usize,
    pub pf_dom_violations: usize,
    pub max_dom_shortfall: f64,
}

#[derive(Clone, Debug)]
pub struct DayTrace {
    pub power_mw: [f64; HOURS_PER_DAY],
    pub soc_mwh: [f64; HOURS_PER_DAY + 1],
    pub net_profit: f64,
}

fn run_single_day<F>(
    day: &DayPath,
    bat: BatterySpec,
    mut decide: F,
) -> (
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    usize,
    usize,
    [usize; NUM_ACTIONS],
    DayTrace,
)
where
    F: FnMut(usize, usize) -> usize,
{
    let mut soc = bat.initial_soc_idx;
    let mut trace = DayTrace {
        power_mw: [0.0; HOURS_PER_DAY],
        soc_mwh: [0.0; HOURS_PER_DAY + 1],
        net_profit: 0.0,
    };
    trace.soc_mwh[0] = BatterySpec::soc_mwh(soc);
    let mut cash_sum = 0.0;
    let mut grid_th = 0.0;
    let mut stored_th = 0.0;
    let mut neg_cf = 0.0;
    let mut sj_cf = 0.0;
    let mut soc_viol = 0;
    let mut pow_viol = 0;
    let mut cnts = [0usize; NUM_ACTIONS];

    for h in 0..HOURS_PER_DAY {
        let mut aid = decide(h, soc);
        let mut act = ACTIONS[aid];
        if !bat.is_feasible(soc, act) {
            if act.power_mw > bat.max_charge_mw + 1e-12
                || (-act.power_mw) > bat.max_discharge_mw + 1e-12
            {
                pow_viol += 1;
            } else {
                soc_viol += 1;
            }
            aid = ACTION_IDLE;
            act = ACTIONS[aid];
        }
        let p = day.prices[h];
        let cf = bat.hourly_cashflow(p, act);
        cash_sum += cf;
        grid_th += act.grid_through_mw;
        stored_th += act.stored_delta.abs();
        if p < 0.0 {
            neg_cf += cf;
        }
        if day.regimes[h] == crate::power::REGIME_SCARCITY || day.num_jumps[h] > 0 {
            sj_cf += cf;
        }
        cnts[aid] += 1;
        soc = (soc as isize + act.delta_soc_idx) as usize;
        trace.power_mw[h] = act.power_mw;
        trace.soc_mwh[h + 1] = BatterySpec::soc_mwh(soc);
    }
    let gross = cash_sum + bat.terminal_value(soc);
    let net = gross - bat.terminal_value(bat.initial_soc_idx);
    trace.net_profit = net;
    (
        net, gross, cash_sum, grid_th, stored_th, neg_cf, sj_cf, soc_viol, pow_viol, cnts, trace,
    )
}

pub fn eval_strategy<F>(name: &str, days: &[DayPath], bat: BatterySpec, decide: F) -> StrategyResult
where
    F: Fn(usize, &DayPath, usize, usize) -> usize,
{
    let n = days.len();
    let mut daily_profits = Vec::with_capacity(n);
    let mut daily_gross = Vec::with_capacity(n);
    let mut daily_cf = Vec::with_capacity(n);
    let mut daily_gth = Vec::with_capacity(n);
    let mut daily_sth = Vec::with_capacity(n);
    let mut daily_efc = Vec::with_capacity(n);
    let mut tot_neg = 0.0;
    let mut tot_sj = 0.0;
    let mut soc_v = 0;
    let mut pow_v = 0;
    let mut act_cnts = [0usize; NUM_ACTIONS];

    for (i, day) in days.iter().enumerate() {
        let (net, gross, cf, gth, sth, ncf, sjcf, sv, pv, cnts, _) =
            run_single_day(day, bat, |h, s| decide(i, day, h, s));
        daily_profits.push(net);
        daily_gross.push(gross);
        daily_cf.push(cf);
        daily_gth.push(gth);
        daily_sth.push(sth);
        daily_efc.push(bat.equivalent_full_cycles(sth));
        tot_neg += ncf;
        tot_sj += sjcf;
        soc_v += sv;
        pow_v += pv;
        for a in 0..NUM_ACTIONS {
            act_cnts[a] += cnts[a];
        }
    }

    let (var5, es5) = stats::var_and_es(&daily_profits, 0.05);
    let fn_d = n.max(1) as f64;
    StrategyResult {
        name: name.into(),
        mean_profit: stats::mean(&daily_profits),
        median_profit: stats::quantile(&daily_profits, 0.5),
        sd_profit: stats::sd(&daily_profits),
        mcse_profit: stats::se(&daily_profits),
        var5_profit: var5,
        es5_profit: es5,
        mean_gross_value: stats::mean(&daily_gross),
        mean_cashflow: stats::mean(&daily_cf),
        mean_throughput: stats::mean(&daily_gth),
        mean_efc: stats::mean(&daily_efc),
        mean_regret_vs_pf: 0.0,
        mcse_regret_vs_pf: 0.0,
        max_regret_vs_pf: 0.0,
        action_counts: act_cnts,
        total_profit_neg_price: tot_neg,
        mean_profit_neg_price_per_day: tot_neg / fn_d,
        total_profit_scarcity_or_jump: tot_sj,
        mean_profit_scarcity_or_jump_per_day: tot_sj / fn_d,
        soc_violation_count: soc_v,
        power_violation_count: pow_v,
        pf_dom_violations: 0,
        max_dom_shortfall: 0.0,
        daily_profits,
        daily_gross_values: daily_gross,
        daily_cashflows: daily_cf,
        daily_grid_throughput: daily_gth,
        daily_stored_throughput: daily_sth,
        daily_efc,
        daily_regret_vs_pf: vec![0.0; n],
    }
}

pub fn eval_stochastic_dp_helper(
    name: &str,
    days: &[DayPath],
    policy: &DpPolicy,
    bat: BatterySpec,
) -> StrategyResult {
    eval_strategy(name, days, bat, |_, d, h, s| {
        let b = policy.grid.bin_index(d.prices[h]);
        policy.action_at(h, d.regimes[h], b, s)
    })
}

pub fn threshold_action(price: f64, soc_idx: usize, bat: BatterySpec) -> usize {
    if price <= 0.0 {
        if bat.is_feasible(soc_idx, ACTIONS[ACTION_CHARGE_2MW]) {
            return ACTION_CHARGE_2MW;
        }
        if bat.is_feasible(soc_idx, ACTIONS[ACTION_CHARGE_1MW]) {
            return ACTION_CHARGE_1MW;
        }
        return ACTION_IDLE;
    }
    if price >= 200.0 {
        if bat.is_feasible(soc_idx, ACTIONS[ACTION_DISCHARGE_18MW]) {
            return ACTION_DISCHARGE_18MW;
        }
        if bat.is_feasible(soc_idx, ACTIONS[ACTION_DISCHARGE_09MW]) {
            return ACTION_DISCHARGE_09MW;
        }
        return ACTION_IDLE;
    }
    ACTION_IDLE
}

pub fn eval_perfect_foresight(days: &[DayPath], bat: BatterySpec) -> StrategyResult {
    let pf_acts: Vec<Vec<usize>> = days
        .iter()
        .map(|d| solve_deterministic_dp(&d.prices, bat).2)
        .collect();
    eval_strategy("perfect_foresight_dp", days, bat, |i, _, h, _| {
        pf_acts[i][h]
    })
}

pub fn eval_no_battery(days: &[DayPath]) -> StrategyResult {
    let n = days.len();
    let mut act_cnts = [0usize; NUM_ACTIONS];
    act_cnts[ACTION_IDLE] = n * HOURS_PER_DAY;
    StrategyResult {
        name: "no_battery".into(),
        daily_profits: vec![0.0; n],
        daily_gross_values: vec![0.0; n],
        daily_cashflows: vec![0.0; n],
        daily_grid_throughput: vec![0.0; n],
        daily_stored_throughput: vec![0.0; n],
        daily_efc: vec![0.0; n],
        daily_regret_vs_pf: vec![0.0; n],
        mean_profit: 0.0,
        median_profit: 0.0,
        sd_profit: 0.0,
        mcse_profit: 0.0,
        var5_profit: 0.0,
        es5_profit: 0.0,
        mean_gross_value: 0.0,
        mean_cashflow: 0.0,
        mean_throughput: 0.0,
        mean_efc: 0.0,
        mean_regret_vs_pf: 0.0,
        mcse_regret_vs_pf: 0.0,
        max_regret_vs_pf: 0.0,
        action_counts: act_cnts,
        total_profit_neg_price: 0.0,
        mean_profit_neg_price_per_day: 0.0,
        total_profit_scarcity_or_jump: 0.0,
        mean_profit_scarcity_or_jump_per_day: 0.0,
        soc_violation_count: 0,
        power_violation_count: 0,
        pf_dom_violations: 0,
        max_dom_shortfall: 0.0,
    }
}

pub fn compute_regret_vs_pf(res: &mut StrategyResult, pf: &StrategyResult, tol: f64) {
    let n = res.daily_profits.len();
    let mut max_reg = f64::NEG_INFINITY;
    let mut max_short = 0.0f64;
    let mut viol = 0usize;
    for i in 0..n {
        let reg = pf.daily_profits[i] - res.daily_profits[i];
        res.daily_regret_vs_pf[i] = reg;
        max_reg = max_reg.max(reg);
        if reg < -tol {
            viol += 1;
            max_short = max_short.max(-reg);
        }
    }
    res.mean_regret_vs_pf = stats::mean(&res.daily_regret_vs_pf);
    res.mcse_regret_vs_pf = stats::se(&res.daily_regret_vs_pf);
    res.max_regret_vs_pf = if n == 0 { 0.0 } else { max_reg };
    res.pf_dom_violations = viol;
    res.max_dom_shortfall = max_short;
}

pub fn policy_action_agreement(
    days: &[DayPath],
    pol_ref: &DpPolicy,
    pol_other: &DpPolicy,
    bat: BatterySpec,
) -> (f64, usize, usize, f64, usize, usize) {
    let mut diff_steps = 0usize;
    let mut total_steps = 0usize;
    for day in days {
        let mut s_ref = bat.initial_soc_idx;
        let mut s_oth = bat.initial_soc_idx;
        for h in 0..HOURS_PER_DAY {
            let b_ref = pol_ref.grid.bin_index(day.prices[h]);
            let b_oth = pol_other.grid.bin_index(day.prices[h]);
            let a_ref = pol_ref.action_at(h, day.regimes[h], b_ref, s_ref);
            let a_oth = pol_other.action_at(h, day.regimes[h], b_oth, s_oth);
            if a_ref != a_oth {
                diff_steps += 1;
            }
            total_steps += 1;
            if let Some(sn) = bat.next_soc_idx(s_ref, ACTIONS[a_ref]) {
                s_ref = sn;
            }
            if let Some(sn) = bat.next_soc_idx(s_oth, ACTIONS[a_oth]) {
                s_oth = sn;
            }
        }
    }
    let eval_agree = (total_steps - diff_steps) as f64 / total_steps.max(1) as f64;

    let mut int_agree = 0usize;
    let mut int_total = 0usize;
    for b_ref in 1..(pol_ref.num_bins - 1) {
        let rep = pol_ref.grid.rep_prices[b_ref];
        let b_oth = pol_other.grid.bin_index(rep);
        if (pol_other.grid.rep_prices[b_oth] - rep).abs() > 1e-12 {
            continue;
        }
        for h in 0..HOURS_PER_DAY {
            for r in 0..NUM_REGIMES {
                for s in 0..bat.num_soc_states {
                    int_total += 1;
                    if pol_ref.action_at(h, r, b_ref, s) == pol_other.action_at(h, r, b_oth, s) {
                        int_agree += 1;
                    }
                }
            }
        }
    }
    let int_frac = int_agree as f64 / int_total.max(1) as f64;
    (
        eval_agree,
        diff_steps,
        total_steps,
        int_frac,
        int_agree,
        int_total,
    )
}

pub fn trace_example_day(
    day: &DayPath,
    bat: BatterySpec,
    policy: Option<&DpPolicy>,
    use_threshold: bool,
    use_pf: bool,
) -> DayTrace {
    let pf_acts = if use_pf {
        solve_deterministic_dp(&day.prices, bat).2
    } else {
        Vec::new()
    };
    let (_, _, _, _, _, _, _, _, _, _, tr) = run_single_day(day, bat, |h, soc| {
        if use_pf {
            pf_acts[h]
        } else if use_threshold {
            threshold_action(day.prices[h], soc, bat)
        } else if let Some(pol) = policy {
            let b = pol.grid.bin_index(day.prices[h]);
            pol.action_at(h, day.regimes[h], b, soc)
        } else {
            ACTION_IDLE
        }
    });
    tr
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SEED_DAY004E_KERNEL;

    fn brute_force_best(prices: &[f64], bat: BatterySpec) -> (f64, Vec<usize>) {
        let h_len = prices.len();
        let mut best_gross = f64::NEG_INFINITY;
        let mut seq = vec![ACTION_IDLE; h_len];
        let mut best_seq = vec![ACTION_IDLE; h_len];
        fn dfs(
            h: usize,
            soc: usize,
            cash: f64,
            prices: &[f64],
            bat: BatterySpec,
            seq: &mut [usize],
            best_gross: &mut f64,
            best_seq: &mut [usize],
        ) {
            if h == prices.len() {
                let tot = cash + bat.terminal_value(soc);
                if tot > *best_gross {
                    *best_gross = tot;
                    best_seq.copy_from_slice(seq);
                }
                return;
            }
            for act in ACTIONS {
                if let Some(sn) = bat.next_soc_idx(soc, act) {
                    seq[h] = act.id;
                    dfs(
                        h + 1,
                        sn,
                        cash + bat.hourly_cashflow(prices[h], act),
                        prices,
                        bat,
                        seq,
                        best_gross,
                        best_seq,
                    );
                }
            }
        }
        dfs(
            0,
            bat.initial_soc_idx,
            0.0,
            prices,
            bat,
            &mut seq,
            &mut best_gross,
            &mut best_seq,
        );
        (best_gross, best_seq)
    }

    #[test]
    fn test_bellman_matches_bruteforce_and_worker_determinism() {
        let bat = BatterySpec::default_spec();
        let prices = [-40.0, 35.0, 240.0, 10.0];
        let (bf_gross, bf_acts) = brute_force_best(&prices, bat);
        let (_, dp_gross, dp_acts, _) = solve_deterministic_dp(&prices, bat);
        assert!((dp_gross - bf_gross).abs() < 1e-12);
        assert_eq!(dp_acts, bf_acts);

        let spec = PowerModelSpec::default_spec();
        let g_exp = PriceGrid::expanded_tail_38();
        assert!(g_exp.lower_tail_at <= -1000.0 && g_exp.upper_tail_at >= 5000.0);
        assert_eq!(g_exp.bin_index(-5000.0), 0);
        assert_eq!(g_exp.bin_index(9000.0), 37);

        let k1 = build_kernel(
            &spec,
            g_exp.clone(),
            100,
            1,
            SEED_DAY004E_KERNEL,
            0x4E38,
            true,
        );
        let k4 = build_kernel(&spec, g_exp, 100, 4, SEED_DAY004E_KERNEL, 0x4E38, true);
        assert_eq!(k1.probs.len(), k4.probs.len());
        for (a, b) in k1.probs.iter().zip(k4.probs.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}
