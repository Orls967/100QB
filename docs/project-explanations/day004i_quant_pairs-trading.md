# Hari 004i — Cointegration dan Trading OU Optimal: Dickey–Fuller, Engle–Granger, Johansen, dan Ambang Bertram (Rust)

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day004i`
- **Tanggal**: `2026-11-22`
- **Topik**: Cointegration dan trading Ornstein–Uhlenbeck (OU) optimal: uji Dickey–Fuller, Engle–Granger, Johansen, dan ambang Bertram
- **Jenis**: `quant`
- **Bahasa**: Rust (*pure standard library*, `edition = "2021"`, nol *crate* eksternal)
- **Lokasi Proyek**: [projects/rust/day004i_pairs_trading](../../projects/rust/day004i_pairs_trading/README.md)

Proyek ini membangun kerangka simulasi, ekonometrika, dan strategi perdagangan statistik berbasis proses **Ornstein–Uhlenbeck (OU)** pada data sintetis sepenuhnya (tanpa klaim tentang pasar nyata). Proyek mencakup empat bagian: (1) transisi eksak vs Euler, MLE bentuk tertutup AR(1) vs numerik, dan studi replikasi $R = 20.000$ bias sampel terbatas (Kendall/Marriott–Pope); (2) waktu lintasan pertama (*first-passage time*) tertutup, BVP beda hingga, Monte Carlo ($R = 200.000$), dan ambang perdagangan optimal Bertram $a^*$; (3) simulasi distribusi kritis Dickey–Fuller ($T = 1000, R = 10^6$), Engle–Granger ($R = 20.000$), regresi palsu (*spurious*), dan uji *trace* Johansen bivariat; serta (4) *backtest* bergulir ($W \in \{1, 2, 5\}$ tahun, $R = 2.000 \times 10$ tahun), filter uji EG, *structural break* + *stop-loss*, dan sensitivitas biaya.

## 2. Masalah dan Tujuan

1. **Bias estimasi kecepatan *mean-reversion* $\theta$**: Mengukur seberapa besar bias positif sampel terbatas $\hat{\theta}$ pada $T \in \{1, 2, 5, 10\}$ tahun dibanding aproksimasi Kendall/Marriott–Pope, serta kemiringan distribusi *half-life* $\ln 2 / \hat{\theta}$.
2. **Verifikasi rumus waktu lintasan pertama dan ambang Bertram**: Menurunkan dan menguji rumus tertutup $\mathbb{E}[\tau(a \to m)]$ (termasuk suku koreksi integral `erf` untuk batas asimetris) terhadap solver BVP beda hingga dan simulasi MC, serta membuktikan bahwa kondisi stasioneritas analitik $R'(a^*) = 0$ cocok dengan optimisasi langsung dan *backtest* MC 200 tahun.
3. **Kalibrasi uji akar unit dan cointegration**: Membuktikan mengapa nilai kritis Dickey–Fuller standar (`-2,864`) tidak boleh dipakai pada residual regresi Engle–Granger (`-3,353`), mendemonstrasikan bahaya regresi palsu (`93,26%` $|t| > 1{,}96$), dan mengukur kuantil kritis serta *power* uji *trace* Johansen bivariat.
4. **Dampak galat estimasi dan *structural break* pada perdagangan**: Mengkuantifikasi "penalti galat estimasi" (laba *oracle* minus estimasi) per panjang jendela $W$, efektivitas filter EG pada pasangan *spurious*, serta kerugian dan waktu deteksi saat $\theta$ runtuh dari `16` ke `2`.

## 3. Apa yang Dikerjakan Proyek

Runner biner `run` mengeksekusi keempat bagian secara deterministik (`main_seed = 20261122`, `8` *worker*) memakai PRNG `xoshiro256++` per indeks lintasan dan reduksi blok tetap (`BLOCK = 64`), lalu mengekspor:
- Berkas hasil berpresisi 17 digit [outputs/day004i/results.json](../../outputs/day004i/results.json).
- Tujuh grafik SVG diagnostik dan perbandingan ke [outputs/figures/](../../outputs/figures/).

## 4. Uraian Kode

Seluruh kode sumber berada di [projects/rust/day004i_pairs_trading/src/](../../projects/rust/day004i_pairs_trading/src/lib.rs):
- [rng.rs](../../projects/rust/day004i_pairs_trading/src/rng.rs) & [parallel.rs](../../projects/rust/day004i_pairs_trading/src/parallel.rs): PRNG `xoshiro256++` + `SplitMix64` dan eksekusi paralel deterministik `1` vs banyak *worker*.
- [stats.rs](../../projects/rust/day004i_pairs_trading/src/stats.rs) & [linalg.rs](../../projects/rust/day004i_pairs_trading/src/linalg.rs): Statistik sampel, kuantil+SE, OLS sederhana ($t$-stat, $R^2$, Durbin–Watson), solver tridiagonal Thomas, Cholesky/invers/eigen $2\times 2$, LU, dan Jacobi simetris.
- [special.rs](../../projects/rust/day004i_pairs_trading/src/special.rs): Fungsi $\operatorname{erfi}(x)$ (deret Taylor + kuadratur Simpson) dan $\operatorname{erfi}'(x) = \frac{2}{\sqrt{\pi}}e^{x^2}$.
- [ou.rs](../../projects/rust/day004i_pairs_trading/src/ou.rs): `simulate_ou_path`, `exact_loglik`, `euler_loglik`, `mle_closed_form`, `mle_numerical`, dan `run_ou_replication_study`.
- [fpt.rs](../../projects/rust/day004i_pairs_trading/src/fpt.rs): `expected_fpt_exact`, `solve_fpt_bvp` ($\mathbb{E}[\tau]$ dan $\operatorname{Var}(\tau)$), dan `simulate_fpt_mc` (dengan koreksi jembatan Brownian).
- [strategy.rs](../../projects/rust/day004i_pairs_trading/src/strategy.rs): `expected_return_rate`, `optimal_threshold_direct`, `optimal_threshold_foc`, `theoretical_sharpe_at_a`, dan `verify_strategy_mc_200yr`.
- [unitroot.rs](../../projects/rust/day004i_pairs_trading/src/unitroot.rs) & [johansen.rs](../../projects/rust/day004i_pairs_trading/src/johansen.rs): `simulate_dickey_fuller_dist` ($R=10^6$ + limit Brownian), `simulate_engle_granger_dist`, `simulate_spurious_regression`, `run_eg_power_grid`, dan `run_johansen_study`.
- [backtest.rs](../../projects/rust/day004i_pairs_trading/src/backtest.rs): `trade_one_year` (akuntansi *mark-to-market* harian) dan `run_backtest_study`.
- [jsonout.rs](../../projects/rust/day004i_pairs_trading/src/jsonout.rs), [svg.rs](../../projects/rust/day004i_pairs_trading/src/svg.rs), & [bin/run.rs](../../projects/rust/day004i_pairs_trading/src/bin/run.rs): Serialisasi JSON, grafik SVG, dan runner utama.

## 5. Metode dan Alur Kerja

1. **Penurunan & Verifikasi Rumus FPT**: Dari BVP Dynkin $\frac{1}{2}\sigma^2 u'' - \theta x u' = -1$, $u(m)=0, u'(-\infty)=0$, integrasi langsung menghasilkan:
   $$\mathbb{E}[\tau(a \to m)] = \frac{\pi}{2\theta}\bigl[\operatorname{erfi}(z_m) - \operatorname{erfi}(z_a)\bigr] + \frac{\sqrt{\pi}}{\theta}\int_{z_a}^{z_m} e^{z^2}\operatorname{erf}(z)\,dz, \quad z = \frac{x\sqrt{\theta}}{\sigma}.$$
   Untuk batas simetris $a = -m$, integral ganjil $\operatorname{erf}(z)$ lenyap sehingga $\mathbb{E}[\tau(-m \to m)] = \frac{\pi}{\theta}\operatorname{erfi}(z_m)$ (faktor $\pi/\theta$ pada $\operatorname{erfi}(z_m)$, bukan $\frac{\pi}{\theta}[\operatorname{erfi}(z_m)-\operatorname{erfi}(-z_m)]$ yang akan dua kali lipat terlalu besar). Untuk batas asimetris $a = 0 \to m = \mathrm{SD}$, suku integral $\operatorname{erf}$ esensial dan terverifikasi terhadap BVP serta MC (`z = 0,43`).
2. **Alur Komputasi**: Parameter acuan ($\theta=16, \mu=0, \sigma=0{,}20, \Delta=1/252, c=0{,}002$) $\to$ replikasi OU $R=20.000$ $\to$ verifikasi FPT & $a^*$ $\to$ simulasi distribusi DF ($10^6$), EG ($20.000$), *spurious*, & Johansen ($20.000$) $\to$ *backtest* bergulir ($2.000 \times 10$ thn) $\to$ JSON & 7 SVG.

## 6. Output yang Ditargetkan

1. `c1_ou_mle_and_stationary_moments`: MLE tertutup vs numerik cocok (`< 1e-8`), gradien relatif `< 1e-8`, Hessian definit negatif, dan ACF(1) empiris cocok terhadap $e^{-\theta\Delta}$ (`< 0,01`).
2. `c2_fpt_closed_bvp_mc_agree`: $\mathbb{E}[\tau]$ tertutup cocok dengan BVP (`rel < 1e-3`) dan MC (`<= 3,5 SE`) pada seluruh pasangan simetris maupun asimetris.
3. `c3_optimal_threshold_and_200yr_mc`: $a^*$ langsung vs kondisi stasioner cocok (`< 1e-6`), serta MC 200 tahun cocok dengan $1/(2\mathbb{E}[\tau])$ dan $R(a^*)$ (`<= 3 SE`).
4. `c4_df_eg_size_and_spurious_rejection`: Ukuran uji EG pada kuantil kritisnya sendiri $\approx 5\%$ (`± 0,015`) dan proporsi $|t| > 1{,}96$ pada regresi *spurious* $> 50\%$.
5. `c5_johansen_size_and_power`: Ukuran uji Johansen pada `q95` $\approx 5\%$ (`± 0,02`) dan *power* pada pasangan berkointegrasi ($\theta=16, T=2$ thn) $> 90\%$.
6. `c6_rolling_penalty_shrinks_with_w`: Penalti galat estimasi (`oracle - est`) menyusut dari $W=1$ ke $W=5$ tahun.
7. `c7_spurious_filtered_profit_nonpos`: Laba pasangan *spurious* setelah filter EG tidak positif bermakna (`<= 3 SE`).
8. `c8_runtime_under_300s`: Total waktu runner di bawah anggaran `300s`.

## 7. Status Pencapaian

Seluruh runner penuh selesai dalam **`3,79s`** pada `8` *worker* ([results.json](../../outputs/day004i/results.json)):

| Kriteria | Status | Bukti Terukur |
| --- | --- | --- |
| `c1_ou_mle_and_stationary_moments` | `tercapai` | `mle_max_diff = 8,17e-9 < 1e-8`; `mle_rel_score = 3,02e-9 < 1e-8`; `hess_neg_def = true`; `acf1_emp = 0,9370` vs `0,9385` (`diff = 0,00149`) |
| `c2_fpt_closed_bvp_mc_agree` | `tercapai` | `(-SD,+SD)`: `0,18721` vs BVP `0,18721` vs MC `0,18672 ± 0,00036` (`z=1,36`); `(0,+SD)`: `0,13084` vs `0,13084` vs `0,13098` (`z=0,43`) |
| `c3_optimal_threshold_and_200yr_mc` | `tercapai` | `a* = -0,0157411` (`diff = 4,9e-10`, `0,4452 SD`); `R(a*) = 0,2044`; MC 200 thn: `z_trades = 0,637`, `z_profit = 0,637` |
| `c4_df_eg_size_and_spurious_rejection` | `tercapai` | DF (`R=10^6`): `q05 = -2,8641 ± 0,0017`; EG (`R=20k`): `q05 = -3,3531 ± 0,0123`, *out-of-sample size* `4,765%` (dalam `3 SE`); Spurious `|t|>1,96`: `93,26%` |
| `c5_johansen_size_and_power` | `belum tercapai` | `q95 = 17,997`, `q95_h0_drift = 18,000 ± 0,097`, *out-of-sample size* `4,870%` (`10,985%` pada `15,41`), `coint_power_2yr = 81,555% < 90%` (`93,460%` pada `15,41`) |
| `c6_rolling_penalty_shrinks_with_w` | `tercapai` | Penalti estimasi menyusut: `0,02427 ± 0,00093` ($W=1$) $\to$ `0,01343` ($W=2$) $\to$ `0,00598 ± 0,00062` ($W=5$) |
| `c7_spurious_filtered_profit_nonpos` | `tercapai` | Laba *spurious* berfilter EG: `-0,00061 ± 0,00038` (`trade_frac = 4,74%`); tanpa filter `-0,01365` |
| `c8_runtime_under_300s` | `tercapai` | Total waktu `3,79s < 300s` (`8` *worker*) |

Ringkasan terukur tambahan (termasuk audit Johansen & ukuran independen):
- **Audit kuantil Johansen (`17,997` vs `15,41`)**: Simulasi fungsional Brownian asimtotik ($T=5.000, R=20.000$) menghasilkan `case1_no_const = 12,334` (literatur `12,53`), `case2_rest_const = 20,266` (`19,96`), `case3_unrest_const` ber-drift nol `17,987`, dan `case3_nonzero_drift` ($[B_1-\bar{B}_1, u-0{,}5]$) $= 15{,}407$ (literatur Osterwald–Lenum `15,41`). Sapuan $T \in \{250, 500, 1000, 5000\}$ menghasilkan `q95 = 18,207; 18,161; 18,162; 18,136`. Artinya `17,997` bukan bug maupun distorsi sampel kecil, melainkan batas asimtotik eksak Case 3 dengan drift nol pada DGP; menggunakan `15,41` pada DGP drift nol menaikkan ukuran *spurious* menjadi `10,985%`.

Ringkasan terukur tambahan:
- **Replikasi OU ($R=20.000$)**: Bias $\hat{\theta}$ pada $T \in \{1, 2, 5, 10\}$ tahun adalah `4,436`, `2,133`, `0,844`, `0,444` (rasio terhadap Kendall `(1+3\phi)/(T\phi)` $= 1{,}091; 1{,}049; 1{,}038; 1{,}092$); *skewness* half-life turun dari `1,257` ($T=1$) menjadi `0,377` ($T=10$), dan `prob_theta_nonpos = 0`.
- **Kuantil DF vs literatur**: `q01 = -3,4301 ± 0,0031`, `q05 = -2,8641 ± 0,0017`, `q10 = -2,5698 ± 0,0014` (`brownian_q05 = -2,8613`), konsisten dengan tabel MacKinnon (`-3,43; -2,86; -2,57`). Bila EG memakai kuantil DF, ukuran uji membengkak dari `4,87%` menjadi `14,56%`.
- **Structural break ($\theta: 16 \to 2$)**: Filter EG mendeteksi putus hubungan pada rata-rata `103,72` hari (`100%` terdeteksi dalam 5 tahun). Tanpa *stop-loss*, laba pasca-break `0,1471/yr` (`max_dd = 0,2469`); dengan *stop-loss* $4\,\mathrm{SD}$, terlalu sering *stopped-out* pada rezim baru yang melebar sehingga laba turun ke `0,0717/yr` (`max_dd = 0,4489`).
- **Biaya**: Pada $c \in \{0; 0{,}002; 0{,}005\}$, laba *oracle* turun `0,4062` $\to$ `0,3940` $\to$ `0,3573` (estimasi $W=2$: `0,3913` $\to$ `0,3822` $\to$ `0,3440`).

## 8. Diagnosis Bagian yang Belum Tercapai (`c5`)

- **FAKTA**:
  1. Pada [results.json](../../outputs/day004i/results.json), kuantil asimtotik fungsional Brownian ($T=5000, R=20.000$) menghasilkan Case 1 `12,334` (literatur `12,53`), Case 2 `20,266` (`19,96`), Case 3 drift-nol `17,987`, dan Case 3 drift non-nol `15,407` (Osterwald–Lenum `15,41`).
  2. Ketika H0 disimulasikan dengan *drift* yang sama persis seperti DGP aplikasi (`dx = 0,05*dt + 0,20*sqrt(dt)*z`), `q95_h0_drift = 18,000 ± 0,097` (hampir identik dengan H0 tanpa drift `17,997 ± 0,089`, karena rasio sinyal/derau drift selama 2 tahun hanya `0,10 / 0,283 ≈ 0,35`, sehingga belum memasuki rezim dominasi tren linear).
  3. Pada sampel *out-of-sample* independen ($R=20.000$), memakai `q95_h0_drift` memberi `spurious_size = 0,04870` (tepat dalam `3 SE` binomial) dan `coint_power_2yr = 0,81555 < 0,90`; memakai kuantil literatur `15,41` memberi `coint_power_at_lit_1541 = 0,93460`, tetapi menaikkan `spurious_size_at_lit_1541 = 0,10985` (> `2x` ukuran nominal `5%`).
- **HIPOTESIS**: Akar masalah rendahnya *power* (`81,555%` vs `90%`) pada $T=2$ tahun bukan bug implementasi maupun "hilangnya derajat kebebasan", melainkan karena pada DGP ber-drift kecil (`0,05/tahun`) di horizon pendek $T=2$ tahun, distribusi null Johansen berada di rezim Case 3 drift-nol (`q95 ≈ 18,00`) alih-alih rezim tren dominan (`15,41`). Menggunakan `15,41` menaikkan *power* ke `93,46%`, tetapi menyebabkan salah tolak *spurious* sebesar `10,99%`. Ambang `90%` sengaja tidak diturunkan sesuai aturan integritas.

## 9. Landasan Teori (11 Poin)

1. **Proses Ornstein–Uhlenbeck & Transisi Eksak**: Solusi $dX_t = \theta(\mu-X_t)dt + \sigma dW_t$ memiliki distribusi kondisional Gaussian eksak dengan peluruhan eksponensial $\phi = e^{-\theta\Delta}$ dan variansi $\frac{\sigma^2(1-\phi^2)}{2\theta}$; diimplementasikan di [ou.rs](../../projects/rust/day004i_pairs_trading/src/ou.rs).
2. **Maximum Likelihood Eksak & Korespondensi AR(1)**: Log-likelihood Gaussian eksak OU terpetakan secara analitik ke regresi AR(1) $X_{i+1} = c + \phi X_i + \varepsilon_i$ dengan pemulihan tertutup $(\hat{\theta}, \hat{\mu}, \hat{\sigma})$; di [ou.rs](../../projects/rust/day004i_pairs_trading/src/ou.rs).
3. **Bias Sampel Terbatas Kendall / Marriott–Pope**: Estimasi $\hat{\phi}$ dengan intersep tak diketahui bias ke bawah sebesar $-\frac{1+3\phi}{N}$, menyebabkan bias ke atas pada $\hat{\theta}$ sebesar $\approx \frac{4}{T}$; diuji di [ou.rs](../../projects/rust/day004i_pairs_trading/src/ou.rs).
4. **Persamaan Kolmogorov Mundur (Dynkin BVP) untuk First-Passage Time**: Waktu tunggu rata-rata $\mathbb{E}[\tau(a \to m)]$ memenuhi $\frac{1}{2}\sigma^2 u'' - \theta x u' = -1$ dan momen keduanya memenuhi $\frac{1}{2}\sigma^2 v'' - \theta x v' = -2u$; diselesaikan dengan metode Thomas tridiagonal di [fpt.rs](../../projects/rust/day004i_pairs_trading/src/fpt.rs) dan [linalg.rs](../../projects/rust/day004i_pairs_trading/src/linalg.rs).
5. **Fungsi Error Imajiner erfi(x) & Rumus Bertram**: Solusi tertutup waktu lintasan pertama OU simetris $\frac{\pi}{\theta}\operatorname{erfi}(m\sqrt{\theta}/\sigma)$ dan optimisasi laju laba $R(a) = (-2a-c)/(2\mathbb{E}[\tau])$; di [special.rs](../../projects/rust/day004i_pairs_trading/src/special.rs) dan [strategy.rs](../../projects/rust/day004i_pairs_trading/src/strategy.rs).
6. **Koreksi Jembatan Brownian (*Brownian Bridge*) pada Exit Diskret**: Peluang lintasan melampaui batas di antara dua langkah $\Delta t$ berbentuk $\exp(-2(b-x_t)(b-x_{t+\Delta t})/(\sigma^2\Delta t))$; dipakai di [fpt.rs](../../projects/rust/day004i_pairs_trading/src/fpt.rs) dan [strategy.rs](../../projects/rust/day004i_pairs_trading/src/strategy.rs).
7. **Uji Akar Unit Dickey–Fuller & Limit Fungsional Brownian**: Statistik-$t$ DF dengan intersep konvergen ke $\frac{\int_0^1 \widetilde{W}\,dW}{\sqrt{\int_0^1 \widetilde{W}^2 dt}}$ yang tak berdistribusi normal; disimulasikan $R=10^6$ di [unitroot.rs](../../projects/rust/day004i_pairs_trading/src/unitroot.rs).
8. **Regresi Palsu (Spurious Regression) Granger–Newbold**: Regresi OLS antar *random walk* independen menghasilkan $t$-stat divergen dan statistik Durbin–Watson mendekati nol; diukur di [unitroot.rs](../../projects/rust/day004i_pairs_trading/src/unitroot.rs) dan [stats.rs](../../projects/rust/day004i_pairs_trading/src/stats.rs).
9. **Uji Cointegration Dua-Tahap Engle–Granger**: Uji ADF pada residual $\hat{e}_t = y_t - \hat{\alpha} - \hat{\beta}x_t$ memerlukan nilai kritis lebih ketat (`-3,353` vs `-2,864`) karena estimasi $\hat{\beta}$ menyerap variansi; di [unitroot.rs](../../projects/rust/day004i_pairs_trading/src/unitroot.rs).
10. **Uji Trace VECM Johansen**: Masalah nilai eigen tergeneralisasi $S_{01}S_{11}^{-1}S_{10}v = \lambda S_{00}v$ direduksi lewat Cholesky $S_{00} = LL^\top$ menjadi matriks simetris untuk menguji rank kointegrasi; di [johansen.rs](../../projects/rust/day004i_pairs_trading/src/johansen.rs).
11. **Akuntansi Mark-to-Market & Risiko Structural Break**: Evaluasi strategi bergulir menghitung PnL MTM harian agar posisi terbuka pada pasangan *spurious* atau pasca-break tidak menyembunyikan kerugian belum terealisasi; di [backtest.rs](../../projects/rust/day004i_pairs_trading/src/backtest.rs).

## 10. Kesimpulan, Batasan, dan Deviasi dari Prompt

1. **Kesimpulan**: Seluruh penurunan analitik (MLE AR(1) vs numerik hingga `< 1e-8`, koreksi bias Kendall $\approx 4/T$, rumus FPT simetris & asimetris vs BVP & MC, optimalitas Bertram $a^*$, kuantil DF/EG & Johansen drift nol vs non-nol, serta pengurangan penalti estimasi seiring $W$) terkonfirmasi secara kuantitatif dalam `3,79s`. Akuntansi MTM harian juga membuktikan bahwa filter EG berhasil menekan laba *spurious* menjadi non-positif (`-0,00061 ± 0,00038`).
2. **Batasan**: Seluruh data adalah sintetis; `c5_johansen_size_and_power` berstatus `belum tercapai` (`81,57% < 90%` pada horizon 2 tahun); `run_completed_without_errors = true` menandakan eksekusi bebas galat numerik, bukan bahwa semua 8 kriteria terpenuhi.
3. **Deviasi dari Prompt**: Tidak ada pelonggaran toleransi atau seed setelah perbaikan `day004j` Bagian A: toleransi `c1` dikembalikan ketat ke `< 1e-8` (dilengkapi cek skor gradien dan Hessian definit negatif), ukuran uji `c4`/`c5` dievaluasi pada sampel independen *out-of-sample* (`3 SE` binomial), dan `c5` tetap dilaporkan `belum tercapai` apa adanya (`81,565% < 90%`).

## 11. Referensi dan Jejak Bukti

- Bertram, W. K. (2010). *Analytic solutions for optimal statistical arbitrage trading*. Physica A, 389(11), 2234–2243.
- Dickey, D. A., & Fuller, W. A. (1979). *Distribution of the estimators for autoregressive time series with a unit root*. JASA, 74(366), 427–431.
- Engle, R. F., & Granger, C. W. J. (1987). *Co-integration and error correction: representation, estimation, and testing*. Econometrica, 55(2), 251–276.
- Johansen, S. (1995). *Likelihood-based inference in cointegrated vector autoregressive models*. Oxford University Press.
- Kendall, M. G. (1954). *Note on bias in the estimation of autocorrelation*. Biometrika, 41(3/4), 403–404.
- Artefak hasil & tes: [outputs/day004i/results.json](../../outputs/day004i/results.json), [README.md](../../projects/rust/day004i_pairs_trading/README.md), [integration_tests.rs](../../projects/rust/day004i_pairs_trading/tests/integration_tests.rs).
- Grafik SVG (`outputs/figures/`): [spread-and-thresholds](../../outputs/figures/day004i_spread-and-thresholds.svg), [ou-theta-bias-vs-t](../../outputs/figures/day004i_ou-theta-bias-vs-t.svg), [bertram-return-rate](../../outputs/figures/day004i_bertram-return-rate.svg), [eg-power-vs-theta](../../outputs/figures/day004i_eg-power-vs-theta.svg), [rolling-profit-vs-window](../../outputs/figures/day004i_rolling-profit-vs-window.svg), [cost-sensitivity](../../outputs/figures/day004i_cost-sensitivity.svg), [fpt-comparison](../../outputs/figures/day004i_fpt-comparison.svg).
