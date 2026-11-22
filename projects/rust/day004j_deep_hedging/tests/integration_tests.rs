use day004j_deep_hedging::bs::bs_call;
use day004j_deep_hedging::dp::solve_entropic_dp_gbm;
use day004j_deep_hedging::fourier::heston_call_fourier;
use day004j_deep_hedging::heston::{generate_heston_qe_paths, HestonParams};
use day004j_deep_hedging::nn::{loss_and_d_pnl, Activation, LossType, MlpPolicy, StepCache};
use day004j_deep_hedging::optim::Adam;
use day004j_deep_hedging::parallel::{reduce_vec_fixed_blocks, run_chunks};
use day004j_deep_hedging::rng::{derive_seed, Rng};
use day004j_deep_hedging::stats;

fn train_20_iters(n_workers: usize) -> Vec<f64> {
    let n_paths = 256_usize;
    let n_steps = 10_usize;
    let d_in = 3_usize;
    let hidden = 16_usize;
    let seed = 20261122_u64;

    let paths: Vec<Vec<f64>> = run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x1001, r as u64));
        let mut s = 1.0_f64;
        let mut out = Vec::with_capacity(n_steps + 1);
        out.push(s);
        let dt = 0.25 / n_steps as f64;
        for _ in 0..n_steps {
            s *= (-0.5 * 0.04 * dt + 0.2 * dt.sqrt() * rg.normal()).exp();
            out.push(s);
        }
        out
    });

    let mut pol = MlpPolicy::new(d_in, hidden, Activation::Tanh, seed ^ 0x99);
    let mut opt = Adam::new(pol.n_params(), 1e-3);
    let p0 = 0.04_f64;
    let cost = 0.001_f64;

    for _ in 0..20 {
        let pnl = run_chunks(n_paths, n_workers, |r| {
            let path = &paths[r];
            let mut cache = StepCache::new(d_in, hidden);
            let mut prev_d = 0.0_f64;
            let mut wealth = p0;
            for t in 0..n_steps {
                let tau = 0.25 * (n_steps - t) as f64 / n_steps as f64;
                let d = pol.forward_step(&[path[t].ln(), tau, prev_d], &mut cache);
                wealth += d * (path[t + 1] - path[t]) - cost * (d - prev_d).abs() * path[t];
                prev_d = d;
            }
            wealth - (path[n_steps] - 1.0).max(0.0)
        });
        let (_loss, d_pnl) = loss_and_d_pnl(&pnl, LossType::Entropic { lambda: 10.0 });
        let (_s, grad) =
            reduce_vec_fixed_blocks(n_paths, pol.n_params(), n_workers, |r, local_g| {
                let path = &paths[r];
                let mut caches = Vec::with_capacity(n_steps);
                let mut prev_d = 0.0_f64;
                for (t, &s_t) in path.iter().enumerate().take(n_steps) {
                    let tau = 0.25 * (n_steps - t) as f64 / n_steps as f64;
                    let mut c = StepCache::new(d_in, hidden);
                    let d = pol.forward_step(&[s_t.ln(), tau, prev_d], &mut c);
                    prev_d = d;
                    caches.push(c);
                }
                let dp = d_pnl[r];
                let mut d_in_buf = vec![0.0_f64; d_in];
                let mut dz2 = vec![0.0_f64; hidden];
                let mut dz1 = vec![0.0_f64; hidden];
                let mut d_prev_delta_next = 0.0_f64;
                for t in (0..n_steps).rev() {
                    let c = &caches[t];
                    let d_cur = c.out;
                    let d_prev = c.x[2];
                    let sgn = if d_cur >= d_prev { 1.0 } else { -1.0 };
                    let d_direct = dp * ((path[t + 1] - path[t]) - cost * sgn * path[t]);
                    pol.backward_step(
                        c,
                        d_direct + d_prev_delta_next,
                        local_g,
                        &mut d_in_buf,
                        &mut dz2,
                        &mut dz1,
                    );
                    d_prev_delta_next = d_in_buf[2] + dp * cost * sgn * path[t];
                }
                0.0
            });
        opt.step(&mut pol.params, &grad);
    }
    pol.params
}

#[test]
fn test_worker_invariance_1_vs_8_bit_identical() {
    let p1 = train_20_iters(1);
    let p8 = train_20_iters(8);
    assert_eq!(p1.len(), p8.len());
    for (i, (&a, &b)) in p1.iter().zip(p8.iter()).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "param[{i}] differs between 1 and 8 workers: {a} vs {b}"
        );
    }
}

#[test]
fn test_dp_zero_cost_near_bs_and_cost_monotonicity() {
    let dp0 = solve_entropic_dp_gbm(1.0, 1.0, 0.0, 0.2, 0.25, 10, 0.0, 10.0);
    let bs = bs_call(1.0, 1.0, 0.0, 0.2, 0.25);
    assert!(
        (dp0.indifference_price - bs).abs() < 0.004,
        "dp0={} bs={bs}",
        dp0.indifference_price
    );
    let dp_c = solve_entropic_dp_gbm(1.0, 1.0, 0.0, 0.2, 0.25, 10, 0.002, 10.0);
    assert!(
        dp_c.indifference_price > dp0.indifference_price,
        "dp_c={} dp0={}",
        dp_c.indifference_price,
        dp0.indifference_price
    );
}

#[test]
fn test_heston_qe_martingale_and_fourier_3se() {
    let hp = HestonParams {
        s0: 1.0,
        k: 1.0,
        r: 0.0,
        kappa: 2.0,
        theta: 0.04,
        xi: 0.5,
        rho: -0.7,
        v0: 0.04,
        t_mat: 0.25,
        n_steps: 30,
    };
    let paths = generate_heston_qe_paths(&hp, 40_000, 8, 20261122, 0x9999);
    let st: Vec<f64> = paths.iter().map(|p| p.s[30]).collect();
    let pay: Vec<f64> = paths.iter().map(|p| (p.s[30] - 1.0).max(0.0)).collect();
    let m_st = stats::mean(&st);
    let se_st = stats::se(&st);
    assert!(
        (m_st - 1.0).abs() <= 3.0 * se_st,
        "E[S_T]={m_st} se={se_st}"
    );
    let mc_c = stats::mean(&pay);
    let se_c = stats::se(&pay);
    let f_c = heston_call_fourier(1.0, 1.0, 0.0, 0.04, 2.0, 0.04, 0.5, -0.7, 0.25);
    assert!(
        (mc_c - f_c).abs() <= 3.0 * se_c,
        "mc={mc_c} fourier={f_c} se={se_c}"
    );
}
