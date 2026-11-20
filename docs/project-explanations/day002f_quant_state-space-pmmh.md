# Penjelasan Proyek — day002f: Inferensi Bayesian State-Space — Kalman, Particle Filter & PMMH (Rust)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002f` (Hari 2, Sesi f)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Rust (`rust`)
- **Folder Proyek**: [rust/day002f_state_space_pmmh](../../rust/day002f_state_space_pmmh) · **Metadata**: [meta.json](../../rust/day002f_state_space_pmmh/meta.json)
- **Ringkasan**: Inferensi Bayesian model ruang-keadaan linear-Gaussian (Model A) dan volatilitas stokastik non-linear (Model B) menggunakan Kalman Filter eksak, Bootstrap Particle Filter (systematic resampling), kuadratur grid, Metropolis–Hastings, serta Particle Marginal Metropolis–Hastings (PMMH) dengan diagnostik split-$\hat{R}$, ESS Geyer, dan MCSE.

## 2. Masalah dan Tujuan
Pada model ruang-keadaan non-linear seperti volatilitas stokastik ($y_t = \exp(x_t/2)\varepsilon_t$), likelihood marginal $p(y_{1:T}|\theta)$ tidak memiliki bentuk tertutup. Proyek ini membuktikan ketakbiasan estimator likelihood *bootstrap particle filter* $\hat{Z}$ terhadap Kalman filter eksak pada Model A, lalu memasukkannya ke dalam skema pseudo-marginal PMMH (Andrieu–Doucet–Holenstein 2010) untuk mengambil sampel posterior eksak pada Model A dan Model B.

## 3. Apa yang Dikerjakan Proyek
Crate `rust/day002f_state_space_pmmh` menjalankan verifikasi Kalman vs matriks kovarians penuh $T\times T$, menguji ketakbiasan $\hat{Z}$ dan penskalaan $\mathrm{Var}(\log\hat{Z}) \propto 1/N_p$, menghitung posterior acuan via integrasi grid deterministik, menjalankan rantai jamak MH dan PMMH secara paralel dengan `rayon`, serta menulis [outputs/day002f/results.json](../../outputs/day002f/results.json) dan 4 grafik SVG.

## 4. Uraian Kode
- [rust/day002f_state_space_pmmh/src/kalman.rs](../../rust/day002f_state_space_pmmh/src/kalman.rs): Kalman filter 1D eksak dan pembanding log-likelihood Gaussian multivariat $T\times T$ (`Cholesky`).
- [rust/day002f_state_space_pmmh/src/pf.rs](../../rust/day002f_state_space_pmmh/src/pf.rs): bootstrap particle filter dengan *systematic resampling* dan pemicu ESS $< N_p/2$.
- [rust/day002f_state_space_pmmh/src/grid.rs](../../rust/day002f_state_space_pmmh/src/grid.rs): kuadratur grid posterior eksak (Model A) dan filter HMM terdiskretisasi 2D (Model B).
- [rust/day002f_state_space_pmmh/src/mcmc.rs](../../rust/day002f_state_space_pmmh/src/mcmc.rs): sampler random-walk Metropolis–Hastings dan PMMH multi-rantai.
- [rust/day002f_state_space_pmmh/src/diagnostics.rs](../../rust/day002f_state_space_pmmh/src/diagnostics.rs): split-$\hat{R}$ Gelman–Rubin, ESS *initial positive sequence* Geyer, dan MCSE.

## 5. Metode dan Alur Kerja
Verifikasi Kalman vs brute-force $\to$ simulasi $R$ replikasi particle filter untuk menguji $E[\hat{Z}/Z] = 1$ dan kemiringan log-log $\mathrm{Var}(\log \hat{Z})$ vs $N_p$ $\to$ kuadratur grid posterior $\to$ eksekusi 4 rantai paralel MH & PMMH pada Model A dan PMMH pada Model B $\to$ perbandingan momen posterior dan diagnostik MCMC.

## 6. Output yang Ditargetkan
Kalman vs brute-force `< 1e-10`; rasio ketakbiasan PF $\hat{Z}/Z$ dalam 3 SE dari 1.0; slope $\mathrm{Var}(\log \hat{Z})$ vs $N_p$ sekitar `-1.0`; split-$\hat{R} < 1.05$ dan ESS memadai; momen posterior MH, PMMH, dan Grid konsisten.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Kalman filter vs brute-force Gaussian $T\times T$ | `tercapai` | `kalman_vs_bruteforce_max_diff = 7.105e-15 < 1e-10` |
| Ketakbiasan $\hat{Z}$ & penskalaan $\mathrm{Var}(\log\hat{Z})$ | `tercapai` | `mean_ratio = 0.98239 ± 0.02043` (`z = 0.862 < 3`), `loglog_slope = -0.9632` |
| Konsistensi posterior Model A (Grid vs MH vs PMMH) | `tercapai` | Massa tepi grid `< 1e-12`; MH & PMMH cocok dengan momen grid dengan split-$\hat{R} < 1.02$ |
| Inferensi PMMH Model B (Stochastic Volatility) | `tercapai` | Semua rantai konvergen dengan split-$\hat{R} < 1.03$ dan selang kredibel mencakup/dekat acuan grid |
| Unit test `cargo test --release` | `tercapai` | `cargo test --release -q` dijalankan pada sesi `day003c`: 5/5 tes lulus dalam 0.09s |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `cargo test --release` serta `outputs/day002f/results.json`.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Ruang-Keadaan Linear-Gaussian & Dekomposisi Prediksi-Koreksi Kalman**: Rekursi mean/variansi prediktif dan inovasi $v_t = y_t - x_{t|t-1}$, $F_t = P_{t|t-1} + \tau^2$ menghitung log-likelihood marginal $\sum_t -\frac{1}{2}(\ln(2\pi F_t) + v_t^2/F_t)$ dalam waktu $O(T)$ (`kalman.rs`).
2. **Model Volatilitas Stokastik Log-AR(1) (Taylor 1986)**: $x_t = \mu + \phi(x_{t-1}-\mu) + \sigma\eta_t$, $y_t = \exp(x_t/2)\varepsilon_t$ menangkap *volatility clustering* non-Gaussian.
3. **Sequential Importance Sampling & Bootstrap Particle Filter (Gordon et al. 1993)**: Menyebarkan partikel dari transisi prior $p(x_t|x_{t-1})$ dan membobotinya dengan likelihood observasi $g(y_t|x_t)$ (`pf.rs`).
4. **Systematic Resampling $O(N_p)$ (Kitagawa 1996)**: Menggunakan satu angka acak $U \sim U(0, 1/N_p)$ dengan kisi berjarak $1/N_p$ untuk meminimalkan variansi Monte Carlo saat resampling.
5. **Ukuran Sampel Efektif Partikel $N_{\text{eff}} = 1/\sum_i (w_t^{(i)})^2$**: Memicu resampling adaptif hanya ketika degenerasi bobot menyebabkan $N_{\text{eff}} < N_p / 2$.
6. **Ketakbiasan Estimator Likelihood Particle Filter (Del Moral 2004)**: $E[\hat{Z}] = Z$ secara eksak untuk setiap $N_p \ge 1$, sementara ketaksamaan Jensen menyebabkan $E[\log \hat{Z}] < \log Z$ dengan bias $\approx -\frac{1}{2}\mathrm{Var}(\log \hat{Z}) = O(1/N_p)$.
7. **Teorema Pseudo-Marginal / PMMH (Andrieu, Doucet & Holenstein 2010)**: Mengganti likelihood eksak $p(y|\theta)$ di dalam rasio penerimaan Metropolis–Hastings dengan estimator positif tak-bias $\hat{Z}(\theta, u)$ tetap memiliki distribusi stasioner marginal $p(\theta|y)$ yang **eksak** (`mcmc.rs`).
8. **Trade-off Variansi Log-Likelihood $\mathrm{Var}(\log\hat{Z}) \approx 1\text{--}2$ (Doucet et al. 2015)**: Menjaga variansi log-likelihood di sekitar $1\text{--}2$ menyeimbangkan biaya per langkah partikel dengan efisiensi pencampuran rantai MCMC (*sticky chains*).
9. **Diagnostik Konvergensi Split-$\hat{R}$ (Gelman et al. 2013)**: Memecah tiap rantai menjadi dua paruh untuk mendeteksi baik perbedaan antar-rantai maupun ketidakstasioneran di dalam rantai (`diagnostics.rs`).
10. **Effective Sample Size (ESS) via Initial Positive Sequence (Geyer 1992)**: Menjumlahkan pasangan autokorelasi berurutan $\hat{\rho}_{2k} + \hat{\rho}_{2k+1} > 0$ untuk mengestimasi efisiensi sampel independen dan $\text{MCSE} = \text{SD}/\sqrt{\text{ESS}}$.

## 10. Kesimpulan dan Batasan
Proyek `day002f` mendemonstrasikan inferensi Bayesian eksak pada model volatilitas stokastik non-linear menggunakan PMMH di Rust.

## 11. Referensi dan Jejak Bukti
- Andrieu, C., Doucet, A., & Holenstein, R. (2010). *Particle Markov chain Monte Carlo methods*. Journal of the Royal Statistical Society: Series B, 72(3), 269–342.
- Artefak: [rust/day002f_state_space_pmmh/README.md](../../rust/day002f_state_space_pmmh/README.md), [outputs/day002f/results.json](../../outputs/day002f/results.json).
