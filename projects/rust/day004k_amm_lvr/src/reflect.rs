use crate::blocks::block_arb_outcome;
use crate::parallel::reduce_tuple4_fixed_blocks;
use crate::rng::{derive_seed, Rng};

/// Theoretical continuous-time reflected Brownian motion limit (Delta -> 0, phi > 0):
/// - half-width delta = -ln(1 - phi) approx phi
/// - total local time rate E[dL^+ + dL^-]/dt = sigma^2 / (2 * delta)
/// - annualized arbitrage volume / V = sigma^2 / (8 * delta)
/// - annualized LP fee income / V = phi * sigma^2 / (8 * delta)
#[inline]
pub fn reflected_limit_rates(sigma: f64, phi: f64) -> (f64, f64, f64) {
    if phi <= 1e-15 {
        return (0.0, 0.0, 0.0);
    }
    let delta = -(1.0 - phi).ln();
    let local_time_rate = sigma * sigma / (2.0 * delta);
    let arb_vol_rate = 0.25 * local_time_rate;
    let fee_rate = phi * arb_vol_rate;
    (delta, arb_vol_rate, fee_rate)
}

/// Multi-step discrete-block simulation on `[0, t_horizon]` with step `dt_sec` seconds
/// to observe convergence of `(arb_profit_rate, lp_fee_rate, arb_vol_rate)` as `Delta -> 0`.
pub fn mc_reflected_block_sim(
    sigma: f64,
    phi: f64,
    dt_sec: f64,
    t_horizon_days: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let sec_yr = 365.0 * 86_400.0;
    let dt = dt_sec / sec_yr;
    let t_yr = t_horizon_days / 365.0;
    let n_steps = (t_yr / dt).round().max(1.0) as usize;
    let mean_z = -0.5 * sigma * sigma * dt;
    let sd_z = sigma * dt.sqrt();
    let l = 0.5_f64;

    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3001 ^ (dt_sec as u64), r as u64));
        let mut p_ext = 1.0_f64;
        let mut p_pool = 1.0_f64;
        let mut tot_prof = 0.0_f64;
        let mut tot_fee = 0.0_f64;
        let mut tot_vol = 0.0_f64;
        for _ in 0..n_steps {
            p_ext *= (mean_z + sd_z * rg.normal()).exp();
            let z = (p_ext / p_pool).ln();
            let (pr, fe, vo, p_new) = block_arb_outcome(z, phi, l, p_pool);
            p_pool = p_new;
            tot_prof += pr;
            tot_fee += fe;
            tot_vol += vo;
        }
        let ann = 1.0 / (n_steps as f64 * dt);
        let p_a = tot_prof * ann;
        [p_a, tot_fee * ann, tot_vol * ann, p_a * p_a]
    });
    let nf = n_paths as f64;
    let m_pr = sums[0] / nf;
    let m_fe = sums[1] / nf;
    let m_vo = sums[2] / nf;
    let se_pr = ((sums[3] / nf - m_pr * m_pr).max(0.0) / nf).sqrt();
    (m_pr, se_pr, m_fe, m_vo)
}

/// Simulates Poisson noise-trader flow with daily volume/TVL `daily_vol_tvl` over `t_days`
/// alongside arbitrageur blocks (`dt = 12s`) to verify breakeven `phi * total_vol / V >= sigma^2 / 8`.
pub fn mc_noise_trader_breakeven(
    sigma: f64,
    phi: f64,
    daily_vol_tvl: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64) {
    let dt_day = 1.0 / 365.0;
    let n_blocks_day = 7200_usize; // 12s blocks per day
    let dt_blk = dt_day / n_blocks_day as f64;
    let mean_z = -0.5 * sigma * sigma * dt_blk;
    let sd_z = sigma * dt_blk.sqrt();
    let l = 0.5_f64;
    let avg_trades_per_day = 24.0_f64; // Poisson arrival rate
    let trade_size = if daily_vol_tvl > 0.0 {
        daily_vol_tvl / avg_trades_per_day
    } else {
        0.0
    };

    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(
            seed,
            0x3100 ^ ((daily_vol_tvl * 100.0) as u64),
            r as u64,
        ));
        let mut p_ext = 1.0_f64;
        let mut p_pool = 1.0_f64;
        let mut tot_arb_prof = 0.0_f64;
        let mut tot_fee = 0.0_f64;
        for _ in 0..n_blocks_day {
            p_ext *= (mean_z + sd_z * rg.normal()).exp();
            let z = (p_ext / p_pool).ln();
            let (pr, fe, _vo, p_new) = block_arb_outcome(z, phi, l, p_pool);
            p_pool = p_new;
            tot_arb_prof += pr;
            tot_fee += fe;
        }
        if daily_vol_tvl > 0.0 {
            let n_noise = rg.poisson(avg_trades_per_day);
            tot_fee += (n_noise as f64) * trade_size * phi;
        }
        let net_lp = (tot_fee - tot_arb_prof) * 365.0;
        [
            net_lp,
            tot_fee * 365.0,
            tot_arb_prof * 365.0,
            net_lp * net_lp,
        ]
    });
    let nf = n_paths as f64;
    let m_net = sums[0] / nf;
    let m_fee = sums[1] / nf;
    let se_net = ((sums[3] / nf - m_net * m_net).max(0.0) / nf).sqrt();
    (m_net, se_net, m_fee)
}
