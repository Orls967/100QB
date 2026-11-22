//! Unified CLI runner for day004g: GLT exact solvers, nonlinear HJB FD,
//! MC simulation (R=200,000 CRN), forward Kolmogorov, strategy comparison,
//! misspecification studies, JSON export, and 7 SVG figures.

use std::path::PathBuf;
use std::time::Instant;

use day004g_market_making::glt::{
    approx_bid_quote, eigendecompose_m, quotes_from_v, solve_v_jacobi, solve_v_mat_exp,
    solve_v_rk4_richardson, theta_from_v, GltParams, QuoteTable,
};
use day004g_market_making::hjb_fd::solve_nonlinear_hjb_theta;
use day004g_market_making::jsonout::{write_json_file, JsonVal};
use day004g_market_making::kolmogorov::solve_forward_kolmogorov;
use day004g_market_making::misspec::run_misspec_studies;
use day004g_market_making::parallel;
use day004g_market_making::rng::MAIN_SEED;
use day004g_market_making::sim::{run_mc_sim, simulate_single_example_path, MarketSpec};
use day004g_market_making::stats;
use day004g_market_making::strategies::Strategy;
use day004g_market_making::svg::{self, LineSeries, PALETTE};

fn find_repo_root() -> PathBuf {
    let mut dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    for _ in 0..8 {
        if dir.join("outputs").is_dir() && dir.join("projects").join("rust").is_dir() {
            return dir;
        }
        if !dir.pop() {
            break;
        }
    }
    PathBuf::from(".")
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let quick = args.iter().any(|a| a == "--quick" || a == "-quick");
    let n_workers = parallel::default_workers();
    let root = find_repo_root();
    let out_json = if quick {
        std::env::temp_dir().join("day004g_quick_results.json")
    } else {
        root.join("outputs").join("day004g").join("results.json")
    };
    let fig_dir = if quick {
        std::env::temp_dir().join("day004g_quick_figures")
    } else {
        root.join("outputs").join("figures")
    };

    let t_total = Instant::now();
    let p = GltParams::default_params();
    let n_paths = if quick { 5_000 } else { 200_000 };

    // ==================== BAGIAN 1: Exact Solvers & HJB ====================
    let t1 = Instant::now();
    let eig = eigendecompose_m(&p);
    let mut max_rel_je = 0.0f64;
    let mut max_rel_jr = 0.0f64;
    for &tau in &[0.1, 0.5, 1.0] {
        let vj = solve_v_jacobi(&p, &eig, tau);
        let ve = solve_v_mat_exp(&p, tau);
        let vr = solve_v_rk4_richardson(&p, tau, 6000);
        for i in 0..p.dim() {
            max_rel_je = max_rel_je.max((vj[i] - ve[i]).abs() / vj[i].abs().max(1e-15));
            max_rel_jr = max_rel_jr.max((vj[i] - vr[i]).abs() / vj[i].abs().max(1e-15));
        }
    }

    let v_t0 = solve_v_jacobi(&p, &eig, p.t_horizon);
    let theta_exact = theta_from_v(&p, &v_t0);
    let exact_ce = theta_exact[p.idx(0)];

    // Nonlinear HJB dt convergence
    let dts = [1e-3, 5e-4, 2.5e-4, 1.25e-4];
    let mut hjb_errs = Vec::new();
    for &dt in &dts {
        let (th_fd, _, _) = solve_nonlinear_hjb_theta(&p, dt, 0.0);
        let err = th_fd
            .iter()
            .zip(theta_exact.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        hjb_errs.push(err);
    }
    let hjb_orders: Vec<f64> = (1..dts.len())
        .map(|i| (hjb_errs[i - 1] / hjb_errs[i]).log2())
        .collect();

    // Symmetry & monotonicity check + sample quotes at t=0, 0.5, 0.9 for q in {-10, 0, 10}
    let mut sym_err = 0.0f64;
    let mut mono_ok = true;
    let mut quote_samples = Vec::new();
    for &t in &[0.0, 0.5, 0.9] {
        let vt = solve_v_jacobi(&p, &eig, p.t_horizon - t);
        let (db, da) = quotes_from_v(&p, &vt);
        for q in -p.q_max..p.q_max {
            sym_err = sym_err.max((db[p.idx(q)] - da[p.idx(-q)]).abs());
            if q + 1 < p.q_max && db[p.idx(q + 1)] <= db[p.idx(q)] {
                mono_ok = false;
            }
        }
        for &q in &[-10_isize, 0, 10] {
            let b = db[p.idx(q)];
            let a = da[p.idx(q)];
            quote_samples.push(JsonVal::Obj(vec![
                ("t", t.into()),
                ("q", q.into()),
                ("bid_delta", b.into()),
                ("ask_delta", a.into()),
                ("spread", (b + a).into()),
            ]));
        }
    }

    // Stationary limit (T=50 & ground state) + spectral gap
    let spectral_gap = eig.eigenvalues[1] - eig.eigenvalues[0];
    let v_t50 = solve_v_jacobi(&p, &eig, 50.0);
    let (db_t1, _) = quotes_from_v(&p, &v_t0);
    let (db_t50, _) = quotes_from_v(&p, &v_t50);
    let (db_ground, _) = quotes_from_v(&p, &eig.eigenvectors[0]);
    let max_t50_vs_ground = (1..(p.dim() - 1))
        .map(|i| (db_t50[i] - db_ground[i]).abs())
        .fold(0.0f64, f64::max);
    let max_t1_vs_ground = (1..(p.dim() - 1))
        .map(|i| (db_t1[i] - db_ground[i]).abs())
        .fold(0.0f64, f64::max);

    // Closed-form approximation error at q in {-10, -5, 0, 5, 10}
    let approx_errs: Vec<JsonVal> = [-10_isize, -5, 0, 5, 10]
        .iter()
        .map(|&q| {
            let ex = db_t1[p.idx(q)];
            let ap = approx_bid_quote(&p, q);
            JsonVal::Obj(vec![
                ("q", q.into()),
                ("exact_bid", ex.into()),
                ("approx_bid", ap.into()),
                ("abs_err", (ex - ap).abs().into()),
            ])
        })
        .collect();

    // Q sensitivity: Q in {10, 20, 40}
    let q_sens: Vec<JsonVal> = [10_isize, 20, 40]
        .iter()
        .map(|&qm| {
            let mut pq = p;
            pq.q_max = qm;
            let eq = eigendecompose_m(&pq);
            let vq = solve_v_jacobi(&pq, &eq, pq.t_horizon);
            let thq = theta_from_v(&pq, &vq);
            let (dbq, _) = quotes_from_v(&pq, &vq);
            JsonVal::Obj(vec![
                ("q_max", qm.into()),
                ("theta_0_0", thq[pq.idx(0)].into()),
                ("bid_q0", dbq[pq.idx(0)].into()),
                ("bid_q5", dbq[pq.idx(5)].into()),
            ])
        })
        .collect();
    let part1_secs = t1.elapsed().as_secs_f64();

    // ==================== BAGIAN 2: MC, Kolmogorov & Strategies ====================
    let t2 = Instant::now();
    let table_1000 = QuoteTable::build(p, 1000);
    let table_2000 = QuoteTable::build(p, 2000);
    let mkt = MarketSpec {
        t_horizon: p.t_horizon,
        sigma: p.sigma,
        mu: 0.0,
        a_intens: p.a_intens,
        k: p.k,
        gamma: p.gamma,
        q_max: p.q_max,
        eps_adverse: 0.0,
    };

    // Utility dt anchor: dt in {1e-3, 5e-4, 2.5e-4}
    let sim_dt_1e3 = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table_1000 },
        n_paths,
        1000,
        n_workers,
        MAIN_SEED,
        0x2001,
    );
    let sim_dt_5e4 = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table_2000 },
        n_paths,
        2000,
        n_workers,
        MAIN_SEED,
        0x2001,
    );
    let table_4000 = QuoteTable::build(p, 4000);
    let sim_dt_25e4 = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table_4000 },
        n_paths,
        4000,
        n_workers,
        MAIN_SEED,
        0x2001,
    );
    let ce_rich = (4.0 * sim_dt_25e4.certainty_equiv - sim_dt_5e4.certainty_equiv) / 3.0; // O(dt) Richardson
    let mc_vs_exact_z =
        (sim_dt_25e4.certainty_equiv - exact_ce).abs() / sim_dt_25e4.se_ce.max(1e-12);

    // Forward Kolmogorov anchor at dt = 5e-4 (2000 steps)
    let kol = solve_forward_kolmogorov(&table_2000, 2000);

    // Strategy comparison (CRN at dt = 5e-4)
    let strat_glt = &sim_dt_5e4;
    let strat_as = run_mc_sim(
        mkt,
        &Strategy::AvellanedaStoikov { p },
        n_paths,
        2000,
        n_workers,
        MAIN_SEED,
        0x2001,
    );
    let strat_sym = run_mc_sim(
        mkt,
        &Strategy::SymmetricFixed {
            delta: p.half_spread_zero_inv(),
            q_max: p.q_max,
        },
        n_paths,
        2000,
        n_workers,
        MAIN_SEED,
        0x2001,
    );
    let strat_tight = run_mc_sim(
        mkt,
        &Strategy::SymmetricFixed {
            delta: 0.1,
            q_max: p.q_max,
        },
        n_paths,
        2000,
        n_workers,
        MAIN_SEED,
        0x2001,
    );
    let diff_glt_as = stats::paired_diff(&strat_glt.wealths, &strat_as.wealths);
    let diff_glt_sym = stats::paired_diff(&strat_glt.wealths, &strat_sym.wealths);
    let diff_glt_tight = stats::paired_diff(&strat_glt.wealths, &strat_tight.wealths);
    let part2_secs = t2.elapsed().as_secs_f64();

    // ==================== BAGIAN 3: Misspecification ====================
    let t3 = Instant::now();
    let misspec = run_misspec_studies(p, &table_2000, n_paths, 2000, n_workers, MAIN_SEED);
    let part3_secs = t3.elapsed().as_secs_f64();

    // ==================== BAGIAN 4: Figures & Criteria ====================
    let t4 = Instant::now();
    let ex_path = simulate_single_example_path(
        mkt,
        &Strategy::GltOptimal { table: &table_2000 },
        2000,
        MAIN_SEED,
    );

    // 7 SVG Figures
    let qs_vec: Vec<f64> = (-p.q_max..=p.q_max).map(|q| q as f64).collect();
    let bid_t0: Vec<f64> = (-p.q_max..p.q_max).map(|q| db_t1[p.idx(q)]).collect();
    let bid_approx: Vec<f64> = (-p.q_max..p.q_max)
        .map(|q| approx_bid_quote(&p, q))
        .collect();
    let qs_no_edge: Vec<f64> = (-p.q_max..p.q_max).map(|q| q as f64).collect();
    svg::write_line_chart(
        &fig_dir.join("day004g_quotes-and-spread.svg"),
        "Kuotes Bid Optimal GLT vs Aproksimasi Tertutup (t=0, Q=20)",
        "Inventori q",
        "Jarak Kuote Bid delta^b",
        &[
            LineSeries {
                label: "GLT Eksak (t=0)",
                x: &qs_no_edge,
                y: &bid_t0,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Aproksimasi Tertutup",
                x: &qs_no_edge,
                y: &bid_approx,
                color: PALETTE[1],
                dashed: true,
            },
        ],
    )
    .unwrap();

    svg::write_line_chart(
        &fig_dir.join("day004g_hjb-dt-convergence.svg"),
        "Konvergensi HJB Non-Linear (theta-form) terhadap Langkah dt",
        "Langkah Waktu dt",
        "Galat Maksimum |theta_FD - theta_eksak|",
        &[LineSeries {
            label: "Max Error theta(0,q)",
            x: &dts,
            y: &hjb_errs,
            color: PALETTE[0],
            dashed: false,
        }],
    )
    .unwrap();

    let mc_dts = [1e-3, 5e-4, 2.5e-4];
    let mc_ces = [
        sim_dt_1e3.certainty_equiv,
        sim_dt_5e4.certainty_equiv,
        sim_dt_25e4.certainty_equiv,
    ];
    let exact_line = [exact_ce, exact_ce, exact_ce];
    svg::write_line_chart(
        &fig_dir.join("day004g_mc-utility-dt-extrapolation.svg"),
        "Anchor Utilitas CARA: Ekivalen Kepastian MC vs Nilai Eksak theta(0,0)",
        "Langkah Simulasi dt",
        "Certainty Equivalent",
        &[
            LineSeries {
                label: "MC Certainty Equiv",
                x: &mc_dts,
                y: &mc_ces,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Eksak theta(0,0)",
                x: &mc_dts,
                y: &exact_line,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    svg::write_bar_chart(
        &fig_dir.join("day004g_strategy-comparison.svg"),
        "Perbandingan Ekivalen Kepastian Antar-Strategi (R=200.000 CRN)",
        "Certainty Equivalent",
        &[
            "GLT Optimal",
            "Avellaneda-Stoikov",
            "Simetris Tetap",
            "Ketat (0.1)",
        ],
        &[
            strat_glt.certainty_equiv,
            strat_as.certainty_equiv,
            strat_sym.certainty_equiv,
            strat_tight.certainty_equiv,
        ],
        PALETTE[0],
    )
    .unwrap();

    let eps_vec: Vec<f64> = misspec.adverse.iter().map(|r| r.eps).collect();
    let ce_unadj: Vec<f64> = misspec
        .adverse
        .iter()
        .map(|r| r.unadjusted_glt.certainty_equiv)
        .collect();
    let ce_shift: Vec<f64> = misspec
        .adverse
        .iter()
        .map(|r| r.const_shifted_glt.certainty_equiv)
        .collect();
    let ce_exact_as: Vec<f64> = misspec
        .adverse
        .iter()
        .map(|r| r.exact_as_hjb.certainty_equiv)
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004g_adverse-selection-misspec.svg"),
        "Dampak Seleksi Merugikan (eps) pada Ekivalen Kepastian",
        "Adverse Selection eps",
        "Certainty Equivalent",
        &[
            LineSeries {
                label: "GLT Tanpa Penyesuaian",
                x: &eps_vec,
                y: &ce_unadj,
                color: PALETTE[1],
                dashed: false,
            },
            LineSeries {
                label: "GLT Geser +eps",
                x: &eps_vec,
                y: &ce_shift,
                color: PALETTE[3],
                dashed: true,
            },
            LineSeries {
                label: "HJB AS Eksak",
                x: &eps_vec,
                y: &ce_exact_as,
                color: PALETTE[2],
                dashed: false,
            },
        ],
    )
    .unwrap();

    let sig_mults: Vec<f64> = misspec.wrong_sigma.iter().map(|r| r.sigma_mult).collect();
    let sig_ces: Vec<f64> = misspec
        .wrong_sigma
        .iter()
        .map(|r| r.summary.certainty_equiv)
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004g_sigma-and-drift-misspec.svg"),
        "Sensitivitas Ekivalen Kepastian terhadap Salah Spesifikasi Volatilitas (sigma_model / sigma)",
        "Pengali Volatilitas Model",
        "Certainty Equivalent",
        &[LineSeries { label: "CE(sigma_model)", x: &sig_mults, y: &sig_ces, color: PALETTE[4], dashed: false }],
    ).unwrap();

    svg::write_line_chart(
        &fig_dir.join("day004g_example-path-s-q-wealth.svg"),
        "Lintasan Contoh GLT (S_t, Inventori q_t/10, Kekayaan MtM)",
        "Waktu t (Hari)",
        "Nilai",
        &[
            LineSeries {
                label: "Harga S_t",
                x: &ex_path.t,
                y: &ex_path.s,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Inventori q_t / 10",
                x: &ex_path.t,
                y: &ex_path.q.iter().map(|x| x * 0.1).collect::<Vec<_>>(),
                color: PALETTE[1],
                dashed: true,
            },
            LineSeries {
                label: "Kekayaan MtM / 10",
                x: &ex_path.t,
                y: &ex_path
                    .mtm_wealth
                    .iter()
                    .map(|x| x * 0.1)
                    .collect::<Vec<_>>(),
                color: PALETTE[2],
                dashed: false,
            },
        ],
    )
    .unwrap();
    let _ = qs_vec;
    let part4_secs = t4.elapsed().as_secs_f64();
    let total_secs = t_total.elapsed().as_secs_f64();

    // Criteria evaluation (separate from technical completion)
    let c1_solvers = max_rel_je < 1e-10 && max_rel_jr < 1e-10;
    let c2_hjb_conv = hjb_errs[3] < 2e-3 && hjb_orders.iter().all(|&o| o > 0.85 && o < 1.15);
    let c3_sym_mono = sym_err < 1e-12 && mono_ok;
    let c4_mc_utility =
        (sim_dt_25e4.certainty_equiv - exact_ce).abs() <= 3.0 * sim_dt_25e4.se_ce + 0.05;
    let c5_kolmogorov = (sim_dt_5e4.mean_executions - kol.expected_executions).abs() < 0.5
        && (sim_dt_5e4.var_inventory - kol.var_inventory).abs() < 0.5;
    let c6_glt_dom = strat_glt.certainty_equiv >= strat_as.certainty_equiv - 3.0 * strat_as.se_ce
        && strat_glt.certainty_equiv >= strat_sym.certainty_equiv - 3.0 * strat_sym.se_ce
        && strat_glt.certainty_equiv >= strat_tight.certainty_equiv - 3.0 * strat_tight.se_ce;
    let c7_as_recovery = misspec.adverse[3].exact_as_hjb.certainty_equiv
        > misspec.adverse[3].unadjusted_glt.certainty_equiv;
    let c8_budget = total_secs < 300.0 || quick;

    let crit_str = |ok: bool| if ok { "tercapai" } else { "belum tercapai" };

    let strat_json = |name: &'static str, s: &day004g_market_making::sim::SimSummary| {
        JsonVal::Obj(vec![
            ("strategy", name.into()),
            ("mean_wealth", s.mean_wealth.into()),
            ("sd_wealth", s.sd_wealth.into()),
            ("sharpe", s.sharpe.into()),
            ("mean_inventory", s.mean_inventory.into()),
            ("sd_inventory", s.sd_inventory.into()),
            ("var_inventory", s.var_inventory.into()),
            ("mean_util", s.mean_util.into()),
            ("se_util", s.se_util.into()),
            ("certainty_equiv", s.certainty_equiv.into()),
            ("se_ce", s.se_ce.into()),
            ("mean_executions", s.mean_executions.into()),
            ("boundary_touch_frac", s.boundary_touch_frac.into()),
        ])
    };

    let doc = JsonVal::Obj(vec![
        (
            "metadata",
            JsonVal::Obj(vec![
                ("session", "day004g".into()),
                ("date_wita", "2026-11-22".into()),
                ("language", "rust".into()),
                ("quick_mode", quick.into()),
                ("workers_used", n_workers.into()),
                ("main_seed", (MAIN_SEED as usize).into()),
                ("n_paths", n_paths.into()),
                (
                    "elapsed_seconds",
                    JsonVal::Obj(vec![
                        ("part1_exact_hjb", part1_secs.into()),
                        ("part2_mc_kolmogorov_strategies", part2_secs.into()),
                        ("part3_misspecification", part3_secs.into()),
                        ("part4_figures", part4_secs.into()),
                        ("total", total_secs.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part1_exact_and_hjb",
            JsonVal::Obj(vec![
                ("alpha", p.alpha().into()),
                ("eta", p.eta().into()),
                ("theta_0_0_exact", exact_ce.into()),
                ("max_rel_err_jacobi_vs_matexp", max_rel_je.into()),
                ("max_rel_err_jacobi_vs_rk4_rich", max_rel_jr.into()),
                ("hjb_dts", dts.to_vec().into()),
                ("hjb_max_errors", hjb_errs.clone().into()),
                ("hjb_convergence_orders", hjb_orders.clone().into()),
                ("quote_symmetry_max_err", sym_err.into()),
                ("quote_monotonicity_ok", mono_ok.into()),
                ("quote_samples", quote_samples.into()),
                ("spectral_gap", spectral_gap.into()),
                ("max_quote_diff_t50_vs_ground", max_t50_vs_ground.into()),
                ("max_quote_diff_t1_vs_ground", max_t1_vs_ground.into()),
                ("closed_form_approx_errors", approx_errs.into()),
                ("q_sensitivity", q_sens.into()),
            ]),
        ),
        (
            "part2_mc_and_strategies",
            JsonVal::Obj(vec![
                ("mc_dt_1e3_ce", sim_dt_1e3.certainty_equiv.into()),
                ("mc_dt_5e4_ce", sim_dt_5e4.certainty_equiv.into()),
                ("mc_dt_25e4_ce", sim_dt_25e4.certainty_equiv.into()),
                ("mc_ce_richardson", ce_rich.into()),
                ("mc_vs_exact_z", mc_vs_exact_z.into()),
                (
                    "kolmogorov_expected_executions",
                    kol.expected_executions.into(),
                ),
                (
                    "mc_expected_executions_dt5e4",
                    sim_dt_5e4.mean_executions.into(),
                ),
                ("kolmogorov_var_inventory", kol.var_inventory.into()),
                ("mc_var_inventory_dt5e4", sim_dt_5e4.var_inventory.into()),
                (
                    "strategies",
                    vec![
                        strat_json("glt_optimal", strat_glt),
                        strat_json("avellaneda_stoikov", &strat_as),
                        strat_json("symmetric_fixed", &strat_sym),
                        strat_json("tight_fixed_0.1", &strat_tight),
                    ]
                    .into(),
                ),
                (
                    "paired_wealth_glt_minus_as",
                    JsonVal::Obj(vec![
                        ("mean_diff", diff_glt_as.mean_diff.into()),
                        ("se_diff", diff_glt_as.se_diff.into()),
                        ("z_stat", diff_glt_as.z_stat.into()),
                    ]),
                ),
                (
                    "paired_wealth_glt_minus_sym",
                    JsonVal::Obj(vec![
                        ("mean_diff", diff_glt_sym.mean_diff.into()),
                        ("se_diff", diff_glt_sym.se_diff.into()),
                        ("z_stat", diff_glt_sym.z_stat.into()),
                    ]),
                ),
                (
                    "paired_wealth_glt_minus_tight",
                    JsonVal::Obj(vec![
                        ("mean_diff", diff_glt_tight.mean_diff.into()),
                        ("se_diff", diff_glt_tight.se_diff.into()),
                        ("z_stat", diff_glt_tight.z_stat.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part3_misspecification",
            JsonVal::Obj(vec![
                (
                    "adverse_selection",
                    misspec
                        .adverse
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("eps", r.eps.into()),
                                (
                                    "unadjusted_glt",
                                    strat_json("unadjusted_glt", &r.unadjusted_glt),
                                ),
                                (
                                    "const_shifted_glt",
                                    strat_json("const_shifted_glt", &r.const_shifted_glt),
                                ),
                                ("exact_as_hjb", strat_json("exact_as_hjb", &r.exact_as_hjb)),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "wrong_sigma",
                    misspec
                        .wrong_sigma
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("sigma_mult", r.sigma_mult.into()),
                                ("sigma_model", r.sigma_model.into()),
                                ("summary", strat_json("glt_wrong_sigma", &r.summary)),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "wrong_intensity",
                    misspec
                        .wrong_intensity
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("label", r.label.into()),
                                ("a_true", r.a_true.into()),
                                ("k_true", r.k_true.into()),
                                ("glt", strat_json("glt", &r.glt_summary)),
                                ("as_classic", strat_json("as_classic", &r.as_summary)),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "wrong_drift",
                    misspec
                        .wrong_drift
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("mu_true", r.mu_true.into()),
                                ("summary", strat_json("glt_wrong_drift", &r.summary)),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        ("run_completed_without_errors", true.into()),
        (
            "criteria",
            JsonVal::Obj(vec![
                (
                    "c1_three_solvers_within_1e10",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c1_solvers).into()),
                        ("max_rel_err_je", max_rel_je.into()),
                        ("max_rel_err_jr", max_rel_jr.into()),
                        ("threshold", 1e-10_f64.into()),
                    ]),
                ),
                (
                    "c2_hjb_nonlinear_order1",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c2_hjb_conv).into()),
                        ("min_err", hjb_errs[3].into()),
                        ("orders", hjb_orders.into()),
                    ]),
                ),
                (
                    "c3_quote_symmetry_and_monotonicity",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c3_sym_mono).into()),
                        ("sym_err", sym_err.into()),
                        ("monotonic", mono_ok.into()),
                    ]),
                ),
                (
                    "c4_mc_utility_matches_exact",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c4_mc_utility).into()),
                        ("mc_ce", sim_dt_25e4.certainty_equiv.into()),
                        ("exact_ce", exact_ce.into()),
                        ("z_score", mc_vs_exact_z.into()),
                    ]),
                ),
                (
                    "c5_kolmogorov_matches_mc",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c5_kolmogorov).into()),
                        (
                            "exec_diff",
                            (sim_dt_5e4.mean_executions - kol.expected_executions)
                                .abs()
                                .into(),
                        ),
                        (
                            "var_q_diff",
                            (sim_dt_5e4.var_inventory - kol.var_inventory).abs().into(),
                        ),
                    ]),
                ),
                (
                    "c6_glt_dominates_benchmarks_ce",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c6_glt_dom).into()),
                        ("glt_ce", strat_glt.certainty_equiv.into()),
                        ("as_ce", strat_as.certainty_equiv.into()),
                        ("sym_ce", strat_sym.certainty_equiv.into()),
                        ("tight_ce", strat_tight.certainty_equiv.into()),
                    ]),
                ),
                (
                    "c7_as_adjustment_improves_ce",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c7_as_recovery).into()),
                        (
                            "eps_0_1_unadj_ce",
                            misspec.adverse[3].unadjusted_glt.certainty_equiv.into(),
                        ),
                        (
                            "eps_0_1_exact_as_ce",
                            misspec.adverse[3].exact_as_hjb.certainty_equiv.into(),
                        ),
                    ]),
                ),
                (
                    "c8_runtime_under_300s",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c8_budget).into()),
                        ("total_seconds", total_secs.into()),
                        ("budget_seconds", 300.0_f64.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(&out_json, &doc).expect("failed to write results.json");
    println!(
        "day004g selesai dalam {:.2}s | workers={n_workers} | CE_exact={exact_ce:.4} | CE_MC={:.4} | c1..c8={}/{}/{}/{}/{}/{}/{}/{}",
        total_secs,
        sim_dt_25e4.certainty_equiv,
        c1_solvers,
        c2_hjb_conv,
        c3_sym_mono,
        c4_mc_utility,
        c5_kolmogorov,
        c6_glt_dom,
        c7_as_recovery,
        c8_budget
    );
}
