use day005d_merton_jump_diffusion::merton::{simulate_merton_mc, MertonParams};
use day005d_merton_jump_diffusion::pide::price_call_pide_cn;
use day005d_merton_jump_diffusion::tail_risk::simulate_short_put_tail_risk;

#[test]
fn integration_worker_invariance_and_seed_determinism() {
    let p = MertonParams::default();
    let ks = [80.0, 100.0, 120.0];

    let mc1 = simulate_merton_mc(&p, &ks, 1.0, 4096, 20261123, 1);
    let mc4 = simulate_merton_mc(&p, &ks, 1.0, 4096, 20261123, 4);
    let mc_diff_seed = simulate_merton_mc(&p, &ks, 1.0, 4096, 99999999, 4);

    assert_eq!(mc1.martingale_mean.to_bits(), mc4.martingale_mean.to_bits());
    assert_ne!(
        mc1.martingale_mean.to_bits(),
        mc_diff_seed.martingale_mean.to_bits()
    );
    for (r1, r4) in mc1.rows.iter().zip(mc4.rows.iter()) {
        assert_eq!(r1.call_mc.to_bits(), r4.call_mc.to_bits());
        assert_eq!(r1.put_mc.to_bits(), r4.put_mc.to_bits());
    }

    let tr1 = simulate_short_put_tail_risk(&p, 90.0, 1.0 / 12.0, 2048, 20261123, 1);
    let tr4 = simulate_short_put_tail_risk(&p, 90.0, 1.0 / 12.0, 2048, 20261123, 4);
    assert_eq!(tr1.var99_merton.to_bits(), tr4.var99_merton.to_bits());
    assert_eq!(tr1.es99_merton.to_bits(), tr4.es99_merton.to_bits());
}

#[test]
fn integration_pide_vs_mc_and_series_3se() {
    let p = MertonParams::default();
    let mc = simulate_merton_mc(&p, &[90.0, 100.0, 110.0], 1.0, 50_000, 20261123, 4);
    assert!(mc.martingale_z.abs() < 3.0);
    for row in &mc.rows {
        assert!(row.call_z.abs() < 3.0, "k={} call_z={}", row.k, row.call_z);
        let pide = price_call_pide_cn(&p, row.k, 1.0, 401, 200);
        assert!((pide - row.call_exact).abs() < 1e-2);
    }
}
