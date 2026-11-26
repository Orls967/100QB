mod solver;

use solver::{
    apply_h, build_dense_h, exp_s2_and_corr, fnv1a, jacobi_eigs, jw_xx_ground, lanczos, Sector,
};
use std::env;
use std::fs;
use std::time::Instant;

fn selftest() -> bool {
    // Hand N=4, Delta=1 periodic, S^z=0: symmetric {|N>=(|0101>+|1010>)/sqrt(2),
    // |D>=(|0011>+|0110>+|1100>+|1001>)/2} -> H_sym=[[-1, sqrt(2)],[sqrt(2), 0]], E0=-2.
    let s4 = Sector::new(4, 2);
    let (e4, _, _, _, _) = lanczos(&s4, 1.0, true);
    let mut h4 = build_dense_h(&s4, 1.0);
    let (ev4, _) = jacobi_eigs(&mut h4, s4.states.len(), false);
    // Single-bond 2-spin singlet: P(1/2) bond H_1 = [[-0.25, 0.5], [0.5, -0.25]] -> E0 = -0.75.
    let mut h_bond = [-0.25, 0.5, 0.5, -0.25];
    let (ev_b, _) = jacobi_eigs(&mut h_bond, 2, false);
    let ok = (e4 - (-2.0)).abs() < 1e-11
        && (ev4[0] - (-2.0)).abs() < 1e-11
        && (ev_b[0] - (-0.75)).abs() < 1e-12;
    println!("selftest ok={ok} E0(4,1)={e4:.12} E0_bond2={:.12}", ev_b[0]);
    ok
}

fn run_ed(out_dir: &str) -> (f64, f64, f64, f64, f64, usize) {
    let deltas = [0.0, 0.5, 1.0, 2.0];
    let (mut max_jac, mut max_jw, mut max_s2, mut max_nn, mut max_sumr, mut peak_vecs) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0usize);
    let mut ed_csv = String::from("N,dim,delta,E0_lanczos,E0_jacobi,diff_jac,iters,sec\n");
    let mut jw_csv = String::from("N,E0_lanczos,E0_jw,abs_diff\n");
    let mut id_csv = String::from("N,S2,nn_corr,E0_over_N,sum_S0_Sr\n");
    let mut prev_e = [0.0f64; 11];
    let mut mono_ok = 1;
    for (ni, n) in (4..=24).step_by(2).enumerate() {
        let sec = Sector::new(n, n / 2);
        let mut last_e = f64::INFINITY;
        for &del in &deltas {
            let need_v = (del == 1.0) || (del == 0.0 && n == 24);
            let t0 = Instant::now();
            let (e0, it, psi, rhist, pk) = lanczos(&sec, del, need_v);
            let dt = t0.elapsed().as_secs_f64();
            if n == 24 && pk > peak_vecs {
                peak_vecs = pk;
            }
            if e0 > last_e {
                mono_ok = 0;
            }
            last_e = e0;
            let mut ej = f64::NAN;
            if n <= 12 {
                let mut hm = build_dense_h(&sec, del);
                ej = jacobi_eigs(&mut hm, sec.states.len(), false).0[0];
                max_jac = max_jac.max((e0 - ej).abs());
            }
            ed_csv.push_str(&format!(
                "{n},{},{del:.1},{e0:.12},{ej:.12},{:.3e},{it},{dt:.4}\n",
                sec.states.len(),
                (e0 - ej).abs()
            ));
            if del == 0.0 {
                let ejw = jw_xx_ground(n);
                max_jw = max_jw.max((e0 - ejw).abs());
                jw_csv.push_str(&format!(
                    "{n},{e0:.12},{ejw:.12},{:.3e}\n",
                    (e0 - ejw).abs()
                ));
            }
            if del == 1.0 {
                prev_e[ni] = e0;
                let (s2, corr) = exp_s2_and_corr(&sec, &psi);
                let sum_r: f64 = corr.iter().sum();
                max_s2 = max_s2.max(s2.abs());
                max_nn = max_nn.max((corr[1] - e0 / (n as f64)).abs());
                max_sumr = max_sumr.max(sum_r.abs());
                id_csv.push_str(&format!(
                    "{n},{s2:.3e},{:.12},{:.12},{sum_r:.3e}\n",
                    corr[1],
                    e0 / (n as f64)
                ));
                if n == 24 {
                    let mut lc = String::from("iter,ritz_res\n");
                    for (k, r) in rhist.iter().enumerate() {
                        lc.push_str(&format!("{},{r:.6e}\n", k + 1));
                    }
                    fs::write(format!("{out_dir}/lanczos_conv_N24.csv"), lc).unwrap();
                    let mut w = vec![0.0; psi.len()];
                    apply_h(&sec, 1.0, &psi, &mut w);
                }
            }
        }
    }
    let _ = (prev_e[0], mono_ok);
    let h = fnv1a(0xcbf2_9ce4_8422_2325, ed_csv.as_bytes());
    fs::write(format!("{out_dir}/ed_energies.csv"), ed_csv).unwrap();
    fs::write(format!("{out_dir}/jw_check.csv"), jw_csv).unwrap();
    fs::write(format!("{out_dir}/identities.csv"), id_csv).unwrap();
    println!("ed done (hash={h:016x}): max_jac={max_jac:.2e} max_jw={max_jw:.2e} max_S2={max_s2:.2e} max_nn={max_nn:.2e} max_sumr={max_sumr:.2e} mono={mono_ok} peak_v={peak_vecs}");
    (max_jac, max_jw, max_s2, max_nn, max_sumr, peak_vecs)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("selftest");
    fs::create_dir_all("out").unwrap();
    match cmd {
        "selftest" => std::process::exit(if selftest() { 0 } else { 1 }),
        "ed" => {
            run_ed("out");
        }
        _ => {
            eprintln!("pending cmd: {cmd}");
            std::process::exit(1);
        }
    }
}
