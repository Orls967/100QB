# Changelog

Semua perubahan penting pada project **100QB** dicatat di sini, per hari.
Format mengikuti [Keep a Changelog](https://keepachangelog.com/id/1.1.0/).

## [Hari 003b] — 2026-11-21 · Eksekusi optimal: Almgren–Chriss, dampak transien, likuiditas stokastik (Rust, Bonus)

### Ditambahkan
- Project Rust 2021 murni `rust/day003b_optimal_execution/` (`#![forbid(unsafe_code)]`, nol crate eksternal runtime): Cholesky, Thomas tridiagonal $O(N)$, Jacobi eigenvalue simetris, RNG `xoshiro256++` + Box–Muller, paralelisme `std::thread` deterministik (chunk tetap), penulis JSON (17 digit) & SVG buatan sendiri.
- Bagian 1–3: Almgren–Chriss (QP KKT vs `sinh` diskrit, limit kontinu, frontier 40 titik, MC `400.000` jalur), Obizhaeva–Wang (solusi umum, anchor eksak, benchmark $N=2000$ $O(N)$ vs $O(N^3)$, power-law, osilasi PD dengan manipulasi harga, non-PD, gabungan), dan DP likuiditas stokastik 2-rezim (rekursi kuadratik, brute-force $N=5$, simulasi `400.000` jalur, nilai adaptivitas `14.09%`, sensitivitas persistensi).
- Output `outputs/day003b/results.json`, 7 grafik SVG `outputs/figures/day003b_*.svg`, serta dokumen penjelasan `docs/project-explanations/day003b_quant_optimal-execution.md`.

### Diperbarui
- `README.md` (tabel progres baris `3b` Rust), `CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 003a] — 2026-11-21 · Volatilitas lokal Dupire dari permukaan SSVI (Python, Bonus)

### Ditambahkan
- Paket `src/qb100/day003a_local_vol/` (`bs.py`, `ssvi.py`, `dupire.py`, `grids.py`, `pde_backward.py`, `pde_forward.py`, `mc.py`, `barrier.py`, `greeks.py`, `calibration.py`, `plotting.py`, `runner.py`): permukaan SSVI bebas-arbitrase beserta turunan analitik & uji kontrafaktual, volatilitas lokal Dupire analitik + verifikasi beda hingga orde 4 & anchor $\eta \to 0$, dua solver PDE independen (backward di $x = \ln S$ dan forward Dupire di $K$, grid tak-seragam `sinh`, Crank–Nicolson + 4 langkah Rannacher), Monte Carlo log-Euler antitetik (`N = 400.000`, `200` langkah/tahun, PCG64) & uji *weak order*, *down-and-out call* (Reiner–Rubinstein, PDE absorbing, MC jembatan Brownian), Greeks (delta, gamma, vega-bucket vs pergeseran paralel), dan kalibrasi balik Levenberg–Marquardt buatan sendiri ($R = 50$ replikasi).
- Notebook `days/day003_quant_volatilitas-lokal-dupire-ssvi.ipynb`, unit test `tests/test_day003a_*.py`, `outputs/day003a/results.json` (presisi 17 digit), dan 7 grafik di `outputs/figures/day003a_*.png`.

### Diperbarui
- `README.md` (tabel progres baris `3`) dan `CHANGELOG.md`.

## [Hari 002j] — 2026-11-20 · Longevity risk: Lee–Carter, anuitas stokastik & modal solvabilitas (Ruby, Bonus)

### Ditambahkan
- Project Ruby 3.3+ murni `ruby/day002j_longevity_risk/` (tanpa gem runtime, YJIT aktif, paralelisme `Process.fork` deterministik): SVD Jacobi satu sisi, estimasi Lee–Carter (tanpa & dengan noise $R=200$), fit Gompertz–Makeham (Nelder–Mead), proyeksi stokastik $k_t$ ($N=20.000$) & bootstrap parametrik, anuitas kohor usia 65 (rekursi eksak vs simulasi 200.000 umur), Hukum Total Variansi portofolio, kurva $\mathrm{CV}(N)$, modal VaR 99.5% / ES 99% (simulasi penuh & aproksimasi normal bersyarat), serta stres deterministik 0.5%.
- Unit test `minitest` (`test/test_longevity.rb`) dan `rubocop` lulus bersih.
- Output `outputs/day002j/results.json` dan 6 grafik SVG (`day002j_mortality-surface.svg`, `day002j_life-expectancy-fan.svg`, `day002j_k-moments-qq.svg`, `day002j_cv-vs-n.svg`, `day002j_liability-hist-var.svg`, `day002j_capital-vs-n.svg`) di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2j` Ruby) dan `CHANGELOG.md`.

## [Hari 002i] — 2026-11-20 · Model short-rate Hull–White 1-faktor: pohon trinomial, Jamshidian, kalibrasi & Bermudan LSM (C++, Bonus)

### Ditambahkan
- Project C++20 `cpp/day002i_hull_white/` (hanya standard library): pohon trinomial Hull–White terkalibrasi ke kurva Nelson–Siegel via induksi maju Arrow–Debreu ($m \in \{12..192\}$, galat kurva $< 1.5 \times 10^{-15}$), pricing analitik & pohon ZBP/ZBC dan cap 5y, dekomposisi Jamshidian (Brent) vs pohon vs Monte Carlo antitetik OU + *control variate* martingale (`N = 200.000`, `dt = 1/52`), kalibrasi Nelder–Mead bebas-noise & 100 replikasi bernoise paralel (`std::thread`), swaption Bermudan 10y (pohon + Richardson vs *two-pass* LSM `N = 100.000`), serta DV01 pergeseran paralel $\pm 1\text{ bp}$.
- Unit test `doctest` (`tests/test_hw.cpp`) lulus di build Release dan UBSan (`-fsanitize=undefined`).
- Output `outputs/day002i/results.json` dan grafik `day002i_tree-convergence.svg`, `day002i_curve-and-swaptions.svg`, `day002i_calib-scatter.svg`, `day002i_bermudan-bounds-dv01.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2i` C++) dan `CHANGELOG.md`.

## [Hari 002h] — 2026-11-20 · Hawkes + limit order book: simulasi eksak, MLE, mesin pencocokan & fakta stilistik (C++, Bonus)

### Ditambahkan
- Project C++20 `cpp/day002h_hawkes_orderbook/` (hanya standard library): proses titik Hawkes (*Ogata thinning*, eksak *Dassios–Zhao*, uji penskalaan ulang waktu KS & Ljung–Box, MLE $O(n)$ + BFGS paralel `std::thread` $R=200$), mesin *limit order book* berkinerja tinggi (*pool allocator*, antrean FIFO intrusif, bitmap 2-level 64-bit, `0` alokasi *heap* setelah *warmup*, `23.9M ops/detik`) vs mesin referensi $O(n)$ pada 1.000.000 order, dan simulasi pasar gabungan ($T=20.000$) dengan fakta stilistik lengkap.
- Unit test `doctest` (`tests/test_hawkes_ob.cpp`) lulus di build Release dan UBSan (`-fsanitize=undefined`).
- Output `outputs/day002h/results.json` dan grafik `day002h_hawkes-fano-qq.svg`, `day002h_hawkes-mle-scatter.svg`, `day002h_ob-latency.svg`, `day002h_market-stylized-facts.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2h` C++) dan `CHANGELOG.md`.

## [Hari 002g] — 2026-11-20 · AAD untuk Greeks Monte Carlo: tape mode-balik, pathwise & keranjang 10 aset (C++, Bonus)

### Ditambahkan
- Project C++20 `cpp/day002g_aad_greeks/` tanpa library AD/numerik eksternal: *tape* mode-balik (`mark`/`rewind_to_mark`), `Dual<T>` mode maju, *forward-over-reverse* Gamma (`Dual<Number>`), faktorisasi Cholesky terdiferensiasi, dan mesin Monte Carlo paralel deterministik `std::thread` (`200.000` path) untuk GBM 1-langkah, Opsi Asia 252-langkah, *Digital call* (penghalusan $\Phi((S_T-K)/\varepsilon)$), keranjang aritmetik 10 aset (67 sensitivitas vs beda hingga CRN), dan keranjang geometrik 10 aset (vs gradien rumus tertutup).
- Unit test `doctest` (`tests/test_aad.cpp`) lulus di build Release dan UBSan (`-fsanitize=undefined`).
- Output `outputs/day002g/results.json` dan grafik `day002g_aad-vs-analytic.svg`, `day002g_digital-bias-variance.svg`, `day002g_basket-speedup.svg`, `day002g_tape-growth.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2g` C++) dan `CHANGELOG.md`.

## [Hari 002f] — 2026-11-20 · Inferensi Bayesian state-space: Kalman, particle filter & PMMH (Rust, Bonus)

### Ditambahkan
- Crate Rust `rust/day002f_state_space_pmmh/` (Kalman filter eksak + verifikasi *brute force* $T=10$ via Cholesky, bootstrap particle filter tanpa alokasi di loop panas + resampling sistematik, posterior grid 2D trapesium, Haario adaptive MH & PMMH paralel `rayon` deterministik, split $\hat{R}$, ESS *initial positive sequence*, MCSE, dan Model B Stochastic Volatility).
- 5 unit test `cargo test` (Kalman vs *brute force*, resampling sistematik tak bias, ketidakbiasan PF & normalisasi grid, validasi ESS AR(1) & $\hat{R}$ i.i.d., reproduktibilitas & invariansi jumlah thread).
- Output `outputs/day002f/results.json` dan grafik `day002f_pf-unbiased-and-var.png`, `day002f_model_a_posteriors.png`, `day002f_model_b_sv.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2f`) dan `CHANGELOG.md`.

## [Hari 002e] — 2026-11-20 · Model Heston: pricing Fourier, kalibrasi & Monte Carlo QE (Bonus)

### Ditambahkan
- Paket `src/qb100/heston/` (`cf.py`, `pricing.py`, `calibration.py`, `simulation.py`): fungsi karakteristik *Little Heston Trap* (Albrecher et al. 2007), *pricing* Fourier Gil-Pelaez & COS (Fang–Oosterlee 2008), *implied volatility* (Newton + Brent), kalibrasi *least-squares* dengan *multi-start* & analisis identifiabilitas Jacobian ($100$ replikasi noisy), serta Monte Carlo QE (Andersen 2008) vs Euler *full-truncation* pada kasus pelanggaran kondisi Feller, lengkap dengan *control variate* $S_T$ dan Greeks CRN.
- Unit test `tests/test_heston_pricing.py` & `tests/test_heston_calib_mc.py` (6 test).
- Notebook `days/day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb` beserta grafik `day002e_smile-and-fourier.png`, `day002e_calibration-noisy.png`, dan `day002e_qe-vs-euler.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2e`) dan `CHANGELOG.md`.

## [Hari 002d] — 2026-11-20 · Risiko kredit Vasicek: distribusi kerugian eksak, LHP & importance sampling (Bonus)

### Ditambahkan
- Modul `src/qb100/credit.py` (kopula Gaussian satu faktor Vasicek, DP konvolusi PMF bersyarat + Gauss–Hermite probabilis 200 node, verifikasi varians bivariat $\Phi_2$ dan `scipy.integrate.quad`, VaR/ES diskret, pendekatan LHP, dan 3 estimator Monte Carlo: Naif, Conditional MC via `CubicSpline` 2001-node DP, serta Importance Sampling pada $Z \sim \mathcal{N}(-\theta, 1)$).
- Unit test `tests/test_credit.py` (4 test: PMF/mean/varians bivariat, kesesuaian `quad`, konvergensi LHP pada $N=5000$, dan ketiga estimator tidak bias/reprodusibel/VRF).
- Notebook `days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb` beserta grafik `day002d_pmf-and-tail.png`, `day002d_lhp-vs-exact.png`, `day002d_is-variance-and-vrf.png`, dan `day002d_convergence-loglog.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2d`) dan `CHANGELOG.md`.

## [Hari 002c] — 2026-11-20 · HMM dua rezim: forward-backward, Baum–Welch & Viterbi (Bonus)

### Ditambahkan
- Modul `src/qb100/hmm.py` dengan NumPy murni (simulator HMM Gaussian, forward naif vs forward–backward ber-skala Rabiner, dekode Viterbi ruang log, enumerasi *brute force* untuk verifikasi, Baum–Welch EM *batched* `(S, K)` dengan inisialisasi kuantil $|r|$ deterministik & pengurutan state menurut $\sigma$, seleksi BIC, diagnostik Viterbi, dan diagram reliabilitas).
- Unit test `tests/test_hmm.py` (3 test: verifikasi *brute force* $T=10$ multi-seed, *overflow* versi naif vs stabilitas versi ber-skala $T=2000$, dan monotonisitas/reprodusibilitas EM).
- Notebook `days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb` beserta grafik `day002c_viterbi-and-em-ll.png`, `day002c_param-recovery.png`, dan `day002c_calibration-and-bic.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2c`) dan `CHANGELOG.md`.

## [Hari 002b] — 2026-11-20 · Estimasi kovarians & portofolio varians minimum: Wishart, Marchenko–Pastur, Ledoit–Wolf (Bonus)

### Ditambahkan
- Modul `src/qb100/covariance.py` (model faktor sintetis, kovarians sampel, shrinkage Ledoit–Wolf target identitas buatan sendiri, GMV bentuk tertutup & `SLSQP`, frontier *long-only*, spektrum Marchenko–Pastur, eksperimen Monte Carlo batched & sapuan rasio $q = N/T$).
- Unit test `tests/test_covariance.py` (5 test: SPD model faktor, GMV tertutup vs numerik, sifat & asimtotik Ledoit–Wolf, rasio Wishart & reprodusibilitas, batas eigenvalue Marchenko–Pastur).
- Notebook `days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb` beserta grafik `day002b_wishart-lw.png`, `day002b_marchenko-pastur.png`, `day002b_q-sweep.png`, dan `day002b_long-only-frontier.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2b`) dan `CHANGELOG.md`.

## [Hari 001i] — 2026-11-19 · Opsi Asia (Ruby): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Project Ruby murni `ruby/day001i_asian_options/` (Ruby 3.3+ / YJIT, nol dependensi runtime, `Random.new` + Box–Muller polar buatan sendiri, `Math.erfc`, Black-Scholes, Kemna–Vorst diskret, paralelisasi deterministik `Process.fork` + `IO.pipe` + `Marshal`, buffer pra-alokasi, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator SVG murni Ruby) beserta 7 unit test Minitest (40 assertions) dan RuboCop bersih (0 offense).
- Output hasil `outputs/day001i/results.json` dan grafik `day001i_asian-convergence-loglog.svg`, `day001i_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `scripts/update_readme.py` (menambah `ruby` ke `LANG_DIRS` beserta test di `tests/test_scripts.py`), `README.md` (tabel progres baris `1i`), dan `.gitignore` (`.bundle/`, `ruby/**/vendor/bundle/`).

## [Hari 001h] — 2026-11-19 · Opsi Asia (C): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Project C17 `c/day001h_asian_options/` (CMake, `xoshiro256**` + Polar Marsaglia, `erfc`, Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `pthread`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator SVG murni `fprintf`) beserta 7 unit test.
- Output hasil `outputs/day001h/results.json` dan grafik `day001h_asian-convergence-loglog.svg`, `day001h_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1h` otomatis terbaca tanpa mengubah `update_readme.py`) dan `.gitignore` (`c/**/build*/`).

## [Hari 001g] — 2026-11-19 · Opsi Asia (C++): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Project C++20 `cpp/day001g_asian_options/` (CMake, `Xoshiro256**` + Polar Marsaglia, `std::erfc`, Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `std::thread`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator SVG murni) beserta 6 unit test `doctest`.
- Output hasil `outputs/day001g/results.json` dan grafik `day001g_asian-convergence-loglog.svg`, `day001g_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1g` otomatis terbaca tanpa mengubah `update_readme.py`) dan `.gitignore` (`cpp/**/build*/`, `cmake-build-*/`).

## [Hari 001f] — 2026-11-19 · Opsi Asia (TypeScript): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Package TypeScript `typescript/day001f_asian_options/` (PRNG `Xoshiro128++` + Polar Marsaglia, `erf` Cody, Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `worker_threads`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator grafik SVG murni) beserta 6 unit test Vitest.
- Output hasil `outputs/day001f/results.json` dan grafik `day001f_asian-convergence-loglog.svg`, `day001f_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `scripts/update_readme.py` (`typescript` di `LANG_DIRS`), `README.md` (tabel progres baris `1f`), dan `.gitignore` (`node_modules/`, `typescript/**/dist/`).

## [Hari 001e] — 2026-11-19 · Opsi Asia (Go): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Module Go `go/day001e_asian_options/` (Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `sync.WaitGroup` + `math/rand/v2` PCG, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope konvergensi log-log) beserta 5 unit test.
- Output hasil `outputs/day001e/results.json` dan grafik `day001e_asian-convergence-loglog.png`, `day001e_asian-estimators-ci.png` di `outputs/figures/`.

### Diperbarui
- Generalisasi `scripts/update_readme.py` untuk membaca `*/meta.json` lintas bahasa (`rust`, `go`, `ts`, `cpp`, `c`) dan menambah kolom `Bahasa` di tabel progres `README.md`.
- `.gitignore` (`go/**/bin/`, `*.test`).

## [Hari 001d] — 2026-11-19 · Opsi Asia (Rust): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Crate Rust `rust/day001d_asian_options/` (Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel `rayon` + `ChaCha8Rng`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope konvergensi log-log) beserta 5 unit test.
- Output hasil `outputs/day001d/results.json` dan grafik `day001d_asian-convergence-loglog.png`, `day001d_asian-estimators-ci.png` di `outputs/figures/`.
- Dukungan metadata `rust/*/meta.json` pada `scripts/update_readme.py` beserta unit test.

### Diperbarui
- `README.md` (tabel progres baris `1d` dengan `Topik No = bonus`) dan `.gitignore` (`target/`).

## [Hari 001c] — 2026-11-19 · Opsi Asia: Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Modul `src/qb100/options.py` (Black-Scholes call/put, geometric Asian call diskret Kemna–Vorst, simulasi GBM tervektorisasi, estimator Naive/Antithetic/Control Variate, serta Delta *pathwise* & *bump-and-revalue* CRN) beserta unit test di `tests/test_options.py`.
- Dukungan `--bonus` (`topic_no: null` / tampilan `"bonus"` di tabel progres) pada `scripts/new_day.py`, `scripts/update_readme.py`, dan `tests/test_scripts.py`.
- Notebook off-roadmap `days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb` beserta grafik `day001c_asian-convergence-loglog.png` dan `day001c_asian-estimators-ci.png`.

### Diperbarui
- `README.md` (tabel progres baris `1c` dengan `Topik No = bonus`).

## [Hari 001b] — 2026-11-19 · Estimasi π dengan Monte Carlo

### Ditambahkan
- Helper `running_mean` dan `running_se` di `src/qb100/simulation.py` beserta unit test di `tests/test_simulation.py`.
- Notebook `days/day001b_sim_estimasi-dengan-monte-carlo.ipynb` (topik 002) beserta 4 grafik (`day001b_pi-scatter.png`, `day001b_pi-convergence.png`, `day001b_pi-error-scaling.png`, `day001b_pi-histogram.png`) di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1b`) dan `ROADMAP.md` (centang topik 002).

## [Hari 001] — 2026-11-19 · Return harian & kumulatif

### Ditambahkan
- Aturan dan dukungan multi-sesi per hari (`day001`, `day001b`, ...) berbasis tanggal WITA (`Asia/Makassar`) pada `scripts/new_day.py`, `scripts/update_readme.py`, `scripts/run_all.py`, dan `templates/day_template.ipynb`.
- Fungsi `cumulative_returns` di `src/qb100/finance.py` beserta unit test di `tests/test_finance.py`.
- Notebook `days/day001_quant_return-harian-kumulatif.ipynb` (topik 001) beserta grafik `outputs/figures/day001_return-harian-kumulatif.png`.

### Diperbarui
- `README.md` (tabel progres dengan kolom Hari & Topik No, bagian scripts, dan konvensi penamaan).
- `ROADMAP.md` (centang topik 001).

## [Hari 000] — 2026-11-19 · Setup baseline

### Ditambahkan
- Struktur repo: `days/`, `templates/`, `src/qb100/`, `scripts/`, `tests/`, `data/`, `outputs/figures/`.
- Package `qb100`: `plotting` (style seragam, `save_figure`), `data` (yfinance + cache CSV),
  `finance` (returns, log returns, volatilitas, Sharpe, drawdown), `simulation` (`get_rng`,
  random walk, Monte Carlo, ringkasan MC).
- Template notebook harian `templates/day_template.ipynb` dengan metadata `qb100`.
- Scripts `new_day.py`, `update_readme.py`, `run_all.py`.
- Pytest untuk `finance`, `simulation`, `data`/`plotting`, dan scripts.
- `ROADMAP.md` berisi 100 topik (50 quant, 50 sim) dalam 4 fase.
- `README.md`, `LICENSE` (MIT), `requirements.txt`, `pyproject.toml`, `.env.example`.
