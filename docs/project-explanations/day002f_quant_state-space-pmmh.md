# Penjelasan Proyek — day002f: Model State-Space, Particle Filter & PMMH (Rust)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002f` (Hari 2, Sesi f)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Rust
- **Proyek**: [projects/rust/day002f_state_space_pmmh/](../../projects/rust/day002f_state_space_pmmh/README.md)
- **Output JSON**: [outputs/day002f/results.json](../../outputs/day002f/results.json)
- **Ringkasan**: Inferensi Bayesian pada model *state-space* linear dan non-linear (Model A: Gaussian AR(1) vs Model B: *Stochastic Volatility* Taylor dengan $\sigma_t = \beta \exp(x_t/2)$) menggunakan *Kalman Filter* eksak, *Bootstrap Particle Filter* (*systematic resampling*), integrasi deterministik *grid*, *Metropolis–Hastings*, serta *Particle Marginal Metropolis–Hastings* (PMMH) paralel di Rust tanpa dependensi eksternal.

## 2. Masalah dan Tujuan
Model deret waktu keuangan sering melibatkan variabel laten $x_t$ yang berevolusi secara Markovian, seperti volatilitas logaritmik stokastik pada Model B ($y_t = \beta \exp(x_t / 2)\varepsilon_t$). Karena distribusi observasi bersifat non-linear dan non-Gaussian terhadap $x_t$, likelihood marginal $p(y_{1:T} \mid \theta) = \int p(y_{1:T} \mid x_{1:T}, \theta) p(x_{1:T} \mid \theta)\,dx_{1:T}$ tidak memiliki bentuk tertutup seperti pada *Kalman Filter*. Pendekatan *Particle Filter* menghasilkan estimator likelihood stokastik $\hat{Z}(\theta)$ yang tak-bias ($E[\hat{Z}(\theta)] = p(y_{1:T} \mid \theta)$), sehingga dapat disisipkan ke dalam algoritma MCMC *pseudo-marginal* (Andrieu, Doucet, & Holenstein 2010) untuk menarik sampel dari posterior eksak $p(\theta \mid y_{1:T})$. Tujuan proyek `day002f` adalah: (1) memvalidasi `kalman.rs` hingga presisi mesin terhadap solusi Gaussian $T \times T$ *brute-force*, (2) membuktikan ketakbiasan $\hat{Z}$ serta hukum penskalaan variansi $\operatorname{Var}(\log \hat{Z}) \propto N_p^{-1}$ lewat $R=2.000$ replikasi, (3) memverifikasi kesesuaian posterior Model A antara *grid*, MH eksak, dan PMMH, serta (4) menjalankan inferensi PMMH multi-rantai pada Model B *Stochastic Volatility* dengan diagnostik $\hat{R}$, ESS, dan MCSE.

## 3. Apa yang Dikerjakan Proyek
Crate Rust [projects/rust/day002f_state_space_pmmh/](../../projects/rust/day002f_state_space_pmmh/Cargo.toml) dibangun dari nol (hanya `std`, 10 modul di `src/`) dan menyelesaikan empat eksperimen utama:
1. **Validasi Kalman Filter**: membandingkan log-likelihood rekursif $O(T)$ di `kalman.rs` terhadap evaluasi matriks kovarians penuh $T \times T$ secara *brute-force*.
2. **Studi Ketakbiasan & Variansi Particle Filter (`pf.rs`)**: menjalankan $R=2.000$ replikasi pada $N_p = 500$ untuk menguji $E[\hat{Z}/Z] = 1$, serta menyapu $N_p \in \{50, 100, 200, 500, 1000, 2000\}$ guna mengukur slope *log-log* $\operatorname{Var}(\log \hat{Z})$ dan *throughput* langkah partikel per detik.
3. **Kalibrasi Silang Model A (`grid.rs`, `mcmc.rs`)**: membandingkan posterior 2-parameter $(\theta_1, \theta_2)$ antara kuadratur *grid* deterministik, MH dengan likelihood Kalman eksak, dan PMMH ($N_p = 500$).
4. **Inferensi PMMH Model B (`models.rs`, `mcmc.rs`, `diagnostics.rs`)**: mengestimasi tiga parameter $(\mu, \phi, \sigma_\eta)$ pada model *Stochastic Volatility* dengan $N_p = 100$ terkalibrasi, menghitung *split*-$\hat{R}$, *Effective Sample Size* (ESS), *Monte Carlo Standard Error* (MCSE), dan mengekspor 3 grafik SVG.

## 4. Uraian Kode
- [src/lib.rs](../../projects/rust/day002f_state_space_pmmh/src/lib.rs) & [src/main.rs](../../projects/rust/day002f_state_space_pmmh/src/main.rs): titik masuk *pipeline* eksperimen, penulis JSON ke [outputs/day002f/results.json](../../outputs/day002f/results.json), dan generator SVG.
- [src/models.rs](../../projects/rust/day002f_state_space_pmmh/src/models.rs): definisi Model A (Gaussian linear state-space) dan Model B (*Stochastic Volatility* non-linear) beserta generator data sintetik.
- [src/kalman.rs](../../projects/rust/day002f_state_space_pmmh/src/kalman.rs): *Kalman filter* skalar/vektor dan evaluator *brute-force* matriks kovarians penuh $T \times T$ untuk uji silang.
- [src/pf.rs](../../projects/rust/day002f_state_space_pmmh/src/pf.rs): *Bootstrap Particle Filter* dengan *log-sum-exp* bobot ternormalisasi, *systematic resampling* ber-variansi rendah, dan penghitung *Effective Sample Size* partikel.
- [src/grid.rs](../../projects/rust/day002f_state_space_pmmh/src/grid.rs): evaluator posterior *grid* deterministik 2D untuk Model A beserta pemeriksaan massa probabilitas di tepi *grid*.
- [src/mcmc.rs](../../projects/rust/day002f_state_space_pmmh/src/mcmc.rs): *sampler* Metropolis–Hastings dan PMMH multi-rantai deterministik berbasis *thread pool* dan `rng.rs` per-rantai.
- [src/diagnostics.rs](../../projects/rust/day002f_state_space_pmmh/src/diagnostics.rs): penghitung *split*-$\hat{R}$ (Gelman–Rubin), fungsi autokorelasi, ESS, dan MCSE.
- [src/rng.rs](../../projects/rust/day002f_state_space_pmmh/src/rng.rs) & [src/plots.rs](../../projects/rust/day002f_state_space_pmmh/src/plots.rs): PRNG deterministik (`SplitMix64` / `Xoroshiro`) dan penulis grafik SVG mandiri.

## 5. Metode dan Alur Kerja
Arsitektur PMMH memisahkan pembangkitan sampel parameter di ruang $\theta$ dari estimasi likelihood laten di ruang partikel $x_{1:T}$:

```mermaid
flowchart LR
    A["Proposal theta* ~ q(.|theta)"] --> B["Bootstrap Particle Filter (Np partikel)"]
    B --> C["Systematic Resampling & Log-Sum-Exp"]
    C --> D["Estimator Tak-Bias log Z_hat(theta*)"]
    D --> E["Rasio Penerimaan MH: p(theta*) Z_hat(theta*) / p(theta) Z_hat(theta)"]
    E --> F["Rantai Paralel Multi-Seed -> split-R_hat, ESS, MCSE"]
```

Pada setiap iterasi rantai PMMH, kandidat $\theta^*$ dievaluasi oleh `pf.rs` dengan $N_p$ partikel independen menggunakan benih RNG deterministik. Karena $E[\hat{Z}(\theta^*)] = p(y_{1:T} \mid \theta^*)$, rantai Markov yang menerima $\theta^*$ dengan peluang $\min(1, \frac{\pi(\theta^*)\hat{Z}(\theta^*)q(\theta\mid\theta^*)}{\pi(\theta)\hat{Z}(\theta)q(\theta^*\mid\theta)})$ memiliki distribusi stasioner persis sama dengan posterior eksak $\pi(\theta \mid y_{1:T})$, tanpa aproksimasi Laplace ataupun linierisasi EKF.

## 6. Output yang Ditargetkan
1. Selisih log-likelihood antara *Kalman Filter* dan evaluasi Gaussian $T \times T$ *brute-force* `< 1e-10`.
2. Rasio rata-rata $\hat{Z}/Z$ pada $R=2.000$ replikasi berada dalam `3 SE` dari `1.0`, dan slope log-log $\operatorname{Var}(\log \hat{Z})$ terhadap $N_p$ mendekati `-1.0`.
3. Momen posterior (mean, SD, kuantil 2.5%, 50%, 97.5%) Model A dari MH dan PMMH berimpit dengan posterior *grid* deterministik, dengan massa tepi *grid* `< 1e-10` dan `split-R_hat < 1.02`.
4. Rantai PMMH Model B (*Stochastic Volatility*) konvergen dengan `split-R_hat < 1.03` dan ESS `> 1.000`.
5. Seluruh 5 unit test Rust (`cargo test --release -q`) lulus dan 3 grafik SVG terekspor ke [outputs/figures/](../../outputs/figures/).

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Kalman filter vs *brute-force* Gaussian $T \times T$ | `tercapai` | `kalman_vs_bruteforce_max_diff = 7.105e-15 < 1e-10` di [outputs/day002f/results.json](../../outputs/day002f/results.json) |
| Ketakbiasan $\hat{Z}$ & penskalaan $\operatorname{Var}(\log \hat{Z})$ | `tercapai` | `mean_ratio = 0.982391 ± 0.020425` (`z = 0.862 < 3`), `var_ll` menyusut `6.3991` ($N_p=50$) ke `0.1798` ($N_p=2000$), `loglog_slope = -0.9632`, `68.35M particle-steps/s` |
| Konsistensi posterior Model A (*Grid* vs MH vs PMMH) | `tercapai` | `grid_edge_mass = 4.958e-13`; Param 1 mean: Grid `-0.5843`, MH `-0.5849` (`r_hat=1.0003, ess=20555.6`), PMMH `-0.5837` (`r_hat=1.0002, ess=11404.6`, acc `22.83%`) |
| Inferensi PMMH Model B (*Stochastic Volatility*) | `tercapai` | `np_calibrated = 100`, `acceptance = 15.43%`, Param 0 (`mean=-0.8313, r_hat=1.0016, ess=1813.7`), Param 1 (`mean=0.9558, r_hat=1.0017, ess=1536.7`) |
| Unit test `cargo test --release -q` | `tercapai` | Dijalankan pada sesi ini (`2026-11-21`): `cargo test --release -q` lulus seluruh tes |

Visualisasi hasil tersimpan pada:
- [day002f_pf-unbiased-and-var.png](../../outputs/figures/day002f_pf-unbiased-and-var.png)
- [day002f_model_a_posteriors.png](../../outputs/figures/day002f_model_a_posteriors.png)
- [day002f_model_b_sv.png](../../outputs/figures/day002f_model_b_sv.png)

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
**FAKTA**: Seluruh kriteria kuantitatif di [outputs/day002f/results.json](../../outputs/day002f/results.json) dan suite `cargo test --release -q` lulus. Pada `pf_unbiased`, tampak bahwa estimator $\log \hat{Z}$ memiliki bias negatif (`bias_ll = -0.3668` di $N_p=500$ dan `-0.1081` di $N_p=2000$); ini merupakan konsekuensi matematis ketaksamaan Jensen $E[\log \hat{Z}] < \log E[\hat{Z}] = \log Z$, bukan cacat pada ketakbiasan $\hat{Z}$ itu sendiri (`mean_ratio = 0.9824`, `z = 0.862`).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model State-Space Markovian & Filtering Marginal**: Sistem *state-space* didefinisikan oleh distribusi awal $p(x_1 \mid \theta)$, transisi laten $p(x_t \mid x_{t-1}, \theta)$, dan emisi observasi $p(y_t \mid x_t, \theta)$. Tujuannya menghitung likelihood marginal $Z(\theta) = p(y_{1:T} \mid \theta)$ dengan mengintegralkan seluruh lintasan laten $x_{1:T}$, yang hanya memiliki solusi analitik tertutup pada kasus linear-Gaussian (`kalman.rs`).
2. **Rekursi Kalman Filter & Verifikasi $T \times T$**: Untuk Model A linear-Gaussian, `kalman.rs` memperbarui mean dan variansi prediktif dalam waktu $O(T)$ sambil mengakumulasi $\log p(y_t \mid y_{1:t-1})$. Kebenarannya dibuktikan oleh `kalman_vs_bruteforce_max_diff = 7.105e-15` terhadap faktorisasi Cholesky langsung matriks kovarians penuh $T \times T$.
3. **Sequential Importance Resampling (Bootstrap Particle Filter)**: Pada Model B non-linear (`y_t = beta * exp(x_t/2) * eps_t`), `pf.rs` mempropagasi $N_p$ partikel $x_t^{(i)} \sim p(x_t \mid x_{t-1}^{(i)}, \theta)$ dan menghitung bobot inkremental $w_t^{(i)} = p(y_t \mid x_t^{(i)}, \theta)$. Estimator likelihood $\hat{Z} = \prod_{t=1}^T \left(\frac{1}{N_p}\sum_{i=1}^{N_p} w_t^{(i)}\right)$ terbukti secara matematis tak-bias terhadap $Z$ untuk sembarang $N_p \ge 1$ (Del Moral 2004).
4. **Systematic Resampling & Degenerasi Partikel**: Tanpa *resampling*, variansi bobot meledak secara eksponensial terhadap $T$ (*weight degeneracy*). Metode *systematic resampling* di `pf.rs` menggunakan satu variabel acak seragam $U \sim \mathcal{U}(0, 1/N_p)$ untuk menarik $N_p$ keturunan secara terstratifikasi dalam waktu $O(N_p)$, meminimalkan variansi Monte Carlo tambahan dibanding *multinomial resampling*.
5. **Ketaksamaan Jensen pada $\log \hat{Z}$ & Penskalaan Variansi $O(N_p^{-1})$**: Meskipun $E[\hat{Z}] = Z$, kecekungan fungsi logaritma menyebabkan $E[\log \hat{Z}] \approx \log Z - \frac{1}{2}\operatorname{Var}(\log \hat{Z})$. Pada `pf_unbiased`, variansi $\operatorname{Var}(\log \hat{Z})$ turun hampir persis berbanding terbalik dengan $N_p$ (`loglog_slope = -0.9632`, dari `6.399` ke `0.180`), sehingga `bias_ll` menyusut sebanding $N_p^{-1}$ dari `-2.662` menjadi `-0.108`.
6. **Teorema Pseudo-Marginal MCMC (Andrieu & Roberts 2009)**: Mengganti likelihood eksak $Z(\theta)$ dengan estimator acak non-negatif $\hat{Z}(\theta; U)$ yang tak-bias di dalam algoritma Metropolis–Hastings mendefinisikan rantai Markov pada ruang diperluas $(\theta, U)$ yang distribusi marjinal stasionernya tetap persis sama dengan posterior eksak $\pi(\theta \mid y_{1:T})$.
7. **Kalibrasi Jumlah Partikel Optimal $N_p$ (Pitt et al. 2012)**: Jika $\operatorname{Var}(\log \hat{Z}) \gg 1$, rantai PMMH mudah "terjebak" (*sticky chain*) saat $\hat{Z}$ beruntung bernilai terlalu tinggi; sebaliknya $N_p$ terlalu besar memboroskan komputasi. Teori menunjukkan target $\operatorname{Var}(\log \hat{Z}) \approx 1.0\text{–}1.5$ (tercapai di sekitar $N_p = 100\text{–}500$ pada Model A & B) memaksimalkan ESS per detik.
8. **Diagnostik Konvergensi Split-$\hat{R}$ (Gelman–Rubin)**: Untuk memastikan rantai MCMC paralel mencapai stasioneritas tanpa bergantung pada titik awal, `diagnostics.rs` membagi tiap rantai menjadi dua paruhan dan membandingkan variansi antar-rantai terhadap variansi dalam-rantai. Nilai `r_hat <= 1.0017` pada Model A dan B jauh di bawah ambang batas `1.05`.
9. **Autokorelasi, Effective Sample Size (ESS) & MCSE**: Sampel MCMC berkorelasi serial sehingga jumlah sampel efektif dihitung sebagai $\text{ESS} = N / (1 + 2\sum_{k=1}^\infty \hat{\rho}_k)$ menggunakan *initial positive sequence* Geyer, dan ketidakpastian estimator rata-rata posterior diukur melalui $\text{MCSE} = \hat{\sigma} / \sqrt{\text{ESS}}$ (`0.00072` pada PMMH Model A).
10. **Paralelisme Deterministik Tanpa Dependensi Eksternal**: Untuk menjaga reproduksibilitas penuh pada eksekusi multi-rantai, `rng.rs` dan `mcmc.rs` mengalokasikan aliran PRNG terpisah per-rantai menggunakan `std::thread::scope`, mencapai *throughput* `68.35 juta particle-steps/detik` pada CPU lokal.

## 10. Kesimpulan dan Batasan
Proyek `day002f` membuktikan implementasi *Particle Marginal Metropolis–Hastings* di Rust murni mampu mereproduksi posterior eksak *Kalman/Grid* pada Model A dan mengkalibrasi model *Stochastic Volatility* non-linear Model B dengan konvergensi $\hat{R} \approx 1.001$. Batasan saat ini adalah penggunaan *Bootstrap Particle Filter* yang mengusulkan partikel dari transisi prior $p(x_t \mid x_{t-1})$ alih-alih *guided/auxiliary particle filter* yang memanfaatkan observasi $y_t$, serta skala $T$ moderat di mana $\operatorname{Var}(\log \hat{Z})$ masih terkendali dengan $N_p = 100$.

## 11. Referensi dan Jejak Bukti
- Andrieu, C., Doucet, A., & Holenstein, R. (2010). *Particle Markov chain Monte Carlo methods*. Journal of the Royal Statistical Society: Series B, 72(3), 269–342.
- Pitt, M. K., Silva, R. S., Giordani, P., & Kohn, R. (2012). *On some properties of Markov chain Monte Carlo simulation methods based on the particle filter*. Journal of Econometrics, 171(2), 134–151.
- Artefak: [projects/rust/day002f_state_space_pmmh/README.md](../../projects/rust/day002f_state_space_pmmh/README.md), [outputs/day002f/results.json](../../outputs/day002f/results.json).
