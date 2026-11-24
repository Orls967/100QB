# Penjelasan Proyek `day006e`: HMC/NUTS dari Nol untuk Volatilitas Stokastik Bayesian dan Kalibrasi Berbasis Simulasi (SBC)

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day006e` (Hari 6, sesi `e`, topik `bonus` di luar `ROADMAP.md`)
- **Tanggal**: `2026-11-24` (WITA)
- **Topik**: Implementasi mandiri *Hamiltonian Monte Carlo* (HMC) dan *No-U-Turn Sampler* (NUTS) multinomial dengan adaptasi *dual averaging* serta matriks massa diagonal berjendela untuk model **Volatilitas Stokastik (SV)** Bayesian, diverifikasi melalui **Simulation-Based Calibration (SBC)**
- **Jenis**: `sim` (simulasi komputasi statistik MCMC ruang berdimensi tinggi & ekonometrika keuangan Bayesian)
- **Bahasa**: `python` (NumPy + SciPy murni di dalam Jupyter Notebook [`days/day006e_sim_bayesian-sv-nuts.ipynb`](../../days/day006e_sim_bayesian-sv-nuts.ipynb), tanpa pustaka pemrograman probabilistik eksternal seperti Stan, PyMC, `emcee`, atau TensorFlow Probability)
- **Artefak Utama**: Notebook tereksekusi penuh [`days/day006e_sim_bayesian-sv-nuts.ipynb`](../../days/day006e_sim_bayesian-sv-nuts.ipynb), ringkasan terstruktur [`outputs/day006e/results.json`](../../outputs/day006e/results.json), dan 6 pasang grafik `PNG` + `SVG` di `outputs/figures/day006e_*`.

Proyek ini membangun mesin inferensi MCMC berdimensi tinggi ($D = T + 3 \ge 1003$) dari prinsip pertama mekanika Hamiltonian, membuktikan secara kuantitatif mengapa parameterisasi non-terpusat (*non-centered parameterization* / NCP) dengan *backpropagation* manual melalui operator adjoin resolven AR(1) memecahkan patologi *Neal's funnel* pada volatilitas stokastik, serta menguji ketepatan matematis dan keseimbangan terperinci (*detailed balance*) sampler melalui 200 replikasi kalibrasi prior-prediktif SBC.

## 2. Masalah dan Tujuan

Inferensi Bayesian penuh pada model Volatilitas Stokastik (SV) waktu diskret mengharuskan penarikan sampel bersama atas parameter struktural $(\mu, \phi, \sigma_\eta)$ dan seluruh lintasan log-variansi laten tak-teramati $h_{1:T} = (h_1, \dots, h_T)$ dari distribusi posterior berdimensi $T + 3$:

$$p(\mu, \phi, \sigma_\eta, h_{1:T} \mid y_{1:T}) \propto p(\mu)\,p(\phi)\,p(\sigma_\eta)\,p(h_1 \mid \mu, \phi, \sigma_\eta)\prod_{t=2}^{T} p(h_t \mid h_{t-1}, \mu, \phi, \sigma_\eta)\prod_{t=1}^{T} p(y_t \mid h_t).$$

Masalah ini menjadi batu ujian klasik bagi algoritma MCMC modern karena tiga alasan mendasar:
1. **Kutukan Dimensi & Korelasi Temporal**: Proposal *random-walk* Metropolis–Hastings gagal berkembang di dimensi $D > 1000$, sementara *particle MCMC* (seperti PMMH yang telah dipelajari terpisah pada sesi Rust [`day002f`](day002f_quant_state-space-pmmh.md)) menanggung variansi penaksir *likelihood* partikel yang tumbuh seiring $T$. HMC memanfaatkan informasi gradien $\nabla \log p$ untuk bergerak koheren di sepanjang himpunan tipikal (*typical set*).
2. **Geometri *Neal's Funnel***: Dalam parameterisasi alami/terpusat $(h_{1:T}, \sigma_\eta)$, skala fluktuasi bersyarat $h_t - \phi h_{t-1}$ dikendalikan langsung oleh $\sigma_\eta$. Saat $\sigma_\eta$ kecil, kerapatan prior menyempit tajam dengan kelengkungan Hessian $\propto \sigma_\eta^{-2}$, sedangkan saat $\sigma_\eta$ besar ruang keadaan melebar — menghentikan integrator simpatik berlangkah tunggal $\varepsilon$.
3. **Kerentanan Bug Halus pada Sampler Dinamis**: Kesalahan tanda Jacobian transformasi, pembatalan numerik dekat $|\phi| \to 1$, ketidaktepatan peluang penerimaan multinomial pada pohon biner NUTS, atau kebocoran adaptasi pasca-*warmup* sering kali tetap menghasilkan kurva posterior yang tampak masuk akal secara kasatmata, sehingga wajib diuji lewat invarian prior-prediktif **SBC**.

## 3. Apa yang Dikerjakan Proyek

Notebook [`days/day006e_sim_bayesian-sv-nuts.ipynb`](../../days/day006e_sim_bayesian-sv-nuts.ipynb) disusun dalam 10 bagian berurutan (25 sel: penjelasan Markdown berbahasa Indonesia, fungsi bernama modular, dan sel uji `assert`):
1. **Pemeriksaan Lingkungan, *Seed* Hierarkis & Kriteria A Priori**: Penurunan 9 aliran independen dari `MASTER_SEED = 20261124` via `np.random.SeedSequence.spawn` dan penetapan ambang `c1`–`c10` sebelum komputasi.
2. **Model SV Non-Terpusat vs Terpusat & Gradien Analitik**: Transformasi ruang bebas batas $(\mu, u_\phi = \mathrm{logit}\,\tilde{\phi}, u_\sigma = \ln\sigma_\eta)$ beserta koreksi log-Jacobian eksak, penurunan *backpropagation* manual lewat operator adjoin resolven AR(1), dan verifikasi selisih hingga terpusat (`c1`).
3. **Integrator Leapfrog & Skema Warmup Berjendela**: Implementasi *leapfrog* tervektorisasi, heuristik pencarian $\varepsilon_0$, *dual averaging* Nesterov–Hoffman–Gelman, dan jendela adaptasi variansi massa diagonal ala Stan, disertai uji reversibilitas, pelestarian determinan Jacobian simplektik $\det J = 1$, dan orde galat energi $O(\varepsilon^2)$ (`c2`).
4. **Multinomial NUTS & Verifikasi Gaussian 50-D**: Pohon biner rekursif berkedalaman maksimum `10` dengan kriteria *generalized U-turn* dalam metrik $M$ dan pencabangan progresif terbias (*biased progressive sampling*), diuji pada Gaussian berkorelasi 50 dimensi dengan kovarians heterogen (`c3`).
5. **Pemulihan Parameter & Jalur $h_{1:1000}$ pada SV Sintetis**: Simulasi $T = 1000$ dengan $(\mu^\star, \phi^\star, \sigma_\eta^\star) = (-9{,}0, 0{,}97, 0{,}25)$, penarikan `4 × 1000` sampel pasca-*warmup* (`n_warmup = 800`), dan evaluasi liputan pita kredibel 90% jalur $h_t$ (`c5`).
6. **Diagnostik Konvergensi & Kontras Terpusat vs Non-Terpusat**: Perhitungan *rank-normalized split-$\hat{R}$*, `ESS bulk`, `ESS tail`, divergensi, kedalaman pohon, dan `E-BFMI` (`c4`), serta perbandingan langsung pada data yang sama terhadap parameterisasi terpusat (`c7`).
7. **Simulation-Based Calibration (SBC)**: `200` replikasi prior-prediktif ($T = 250$, `4 × 300` sampel di-*thin* ke $L = 99$), uji $\chi^2_{19}$ 20-bin, dan pita kepercayaan simultan 99% selisih ECDF (`c6`).
8. **Aplikasi Seri Sintetis Berlompat 4-Rezim & Skor Prediktif Bergulir**: Estimasi volatilitas pada seri lonjakan krisis ($T = 500$, tanpa jaringan) dan evaluasi kausal *log predictive score* (LPS) 1-hari ke depan di `30` titik asal bergulir terhadap **EWMA ($\lambda = 0{,}94$)**.
9. **Rekayasa Performa & Ablasi Efisiensi**: Pengukuran mikrodetik gradien loop Python vs resolven C `scipy.signal.lfilter` serta tabel median $\pm$ MAD (5 pengulangan) `ESS/detik` terhadap variasi adaptasi $(\varepsilon, M)$.
10. **Validasi Artefak & Bukti Reproduksibilitas Dua Eksekusi**: Pemeriksaan validitas XML SVG (`c10`), sidik jari SHA-256 lintas dua eksekusi beruntun dari nol (`c8`), dan ekspor [`outputs/day006e/results.json`](../../outputs/day006e/results.json).

## 4. Uraian Kode

Seluruh implementasi berada di dalam [`days/day006e_sim_bayesian-sv-nuts.ipynb`](../../days/day006e_sim_bayesian-sv-nuts.ipynb) dengan pemisahan tanggung jawab per fungsi:

- **Transformasi & Model SV (Sel 3–4)**:
  - `unpack_sv_theta(theta3)`: Memetakan $(\mu, u_\phi, u_\sigma) \in \mathbb{R}^3$ secara bebas pembatalan numerik memakai identitas $\ln\tilde{\phi} = -\ln(1 + e^{-|u_\phi|}) - \max(0, -u_\phi)$, $\phi = \tilde{\phi} - (1 - \tilde{\phi})$, dan $\sqrt{1 - \phi^2} = 2\exp\bigl(\tfrac{1}{2}(\ln\tilde{\phi} + \ln(1 - \tilde{\phi}))\bigr)$.
  - `reconstruct_h_ncp(q)`: Membangun lintasan $h_{1:T} = \mu + \sigma_\eta R_\phi \eta$ menggunakan `scipy.signal.lfilter`.
  - `sv_ncp_neg_log_post_and_grad(q, log_y2)` & `sv_ncp_grad_python_loop(q, log_y2)`: Menghitung energi potensial $U(q) = -\ln p(q \mid y)$ dan gradien analitik penuh lewat resolven maju $z = R_\phi\eta$ dan resolven adjoin mundur $\lambda = R_\phi^\top(\sigma_\eta g_h)$.
  - `sv_cp_neg_log_post_and_grad(q, log_y2)`: Formulasi parameterisasi terpusat langsung pada koordinat $(h_1, \dots, h_T)$ beserta gradien tridiagonal analitik.
  - `central_fd_grad(fn, q, eps_fd)`: Pemeriksa selisih hingga orde dua untuk sel uji `assert` kriteria `c1`.
- **Integrator & Adaptasi Warmup (Sel 6–7)**:
  - `leapfrog_step` & `hmc_integrate`: Integrator Velocity Verlet bervektor dalam metrik diagonal `inv_m` dengan pengaman luapan numerik.
  - `find_reasonable_epsilon` & `build_stan_warmup_windows`: Penentuan $\varepsilon_0$ awal dan pembagian jadwal *warmup* menjadi *init fast buffer* (15%), jendela kovarians lambat mengganda (`25, 50, 100, ...`), serta *terminal fast buffer* (20%).
- **Mesin Multinomial NUTS & Diagnostik (Sel 9–10)**:
  - `check_merged_no_uturn`: Menguji ketiadaan putar-balik baik pada bentang pohon penuh maupun pada dua sub-lintasan sambungan antarsubpohon secara simetris waktu.
  - `build_nuts_tree` & `nuts_step`: Membangun pohon biner rekursif (`max_depth = 10`) dengan penarikan multinomial dalam sub-pohon dan transisi terbias progresif $\min(1, w_{\text{ext}}/w_{\text{old}})$ antar-penggandaan kedalaman.
  - `run_nuts_chain`: Pengendali utama *warmup* (*dual averaging* + pembaruan `inv_m` ter-regularisasi) dan pengumpulan sampel pasca-*warmup* dengan $(\varepsilon, M)$ yang dibekukan ketat.
  - `split_rhat_and_ess` & `compute_ebfmi`: Implementasi *rank-normalized split-$\hat{R}$* (maksimum komponen *bulk* dan *folded*), `ESS bulk` via FFT autokovarians dan barisan positif monoton awal Geyer, `ESS tail` kuantil 5%/95%, serta diagnostik energi E-BFMI.
- **Eksperimen, SBC, Prediksi Rezim, Profil & Ekspor (Sel 12–23)**:
  - `save_and_show`, `empirical_acf_mean`, `_run_sbc_replication`, `_rolling_forecast_worker`, `median_mad`, dan `verify_all_figures`.

## 5. Metode dan Alur Kerja

### 5.1 Alur Inferensi dan Penurunan Gradien Resolven Adjoin

Untuk inovasi baku $\eta_{1:T} \sim \mathcal{N}(0, I_T)$, barisan AR(1) $z_1 = (1-\phi^2)^{-1/2}\eta_1$, $z_t = \phi z_{t-1} + \eta_t$ ($t \ge 2$), dan $h_t = \mu + \sigma_\eta z_t$, log-likelihood pengamatan adalah $\ell(h) = -\frac{T}{2}\ln(2\pi) - \frac{1}{2}\sum_{t=1}^T \bigl(h_t + y_t^2 e^{-h_t}\bigr)$. Dengan mendefinisikan $g_t = \frac{\partial \ell}{\partial h_t} = \frac{1}{2}(y_t^2 e^{-h_t} - 1)$ dan variabel adjoin mundur $\lambda_T = \sigma_\eta g_T$, $\lambda_t = \sigma_\eta g_t + \phi\lambda_{t+1}$ ($t = T-1, \dots, 1$), aturan rantai menghasilkan gradien eksak berbiaya $O(T)$:

$$\frac{\partial \ell}{\partial \mu} = \sum_{t=1}^T g_t, \qquad \frac{\partial \ell}{\partial u_\sigma} = \sigma_\eta \sum_{t=1}^T g_t z_t, \qquad \frac{\partial \ell}{\partial u_\phi} = \frac{1-\phi^2}{2}\sum_{t=2}^T \lambda_t z_{t-1} + \frac{1}{2}\phi\lambda_1 z_1, \qquad \frac{\partial \ell}{\partial \eta_1} = \frac{\lambda_1}{\sqrt{1-\phi^2}}, \qquad \frac{\partial \ell}{\partial \eta_t} = \lambda_t \;\;(t \ge 2).$$

```mermaid
flowchart LR
    S["SeedSequence(20261124) .spawn(9)"] --> G["Verifikasi Gradien FD (c1) & Leapfrog Simplektik (c2)"]
    G --> N50["Multinomial NUTS Gaussian 50-D (c3)"]
    N50 --> SV["SV Sintetis T=1000: NCP vs CP Funnel & Liputan h_t (c4, c5, c7)"]
    SV --> SBC["SBC 200 Replikasi Prior-Prediktif T=250, L=99 (c6)"]
    SBC --> APP["Aplikasi 4-Rezim & LPS Bergulir vs EWMA + Ablasi Efisiensi"]
    APP --> V["Validasi 2x Eksekusi SHA-256 (c8) + 6 SVG/PNG (c10) -> results.json"]
```

### 5.2 Ringkasan Kuantitatif Diagnostik, Kontras NCP vs CP, dan Efisiensi

Perbandingan langsung parameterisasi Non-Terpusat (`NCP`, waktu `26,16 s`) terhadap Terpusat (`CP`, waktu `57,62 s`) pada data sintetis $T = 1000$ (`4` rantai $\times$ `1000` sampel pasca-*warmup*):

| Parameter | NCP Split-$\hat{R}$ | NCP ESS Bulk | NCP ESS Tail | NCP ESS/s | CP Split-$\hat{R}$ | CP ESS Bulk | CP ESS/s | Rasio ESS/s (NCP / CP) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `mu` | `1,0033` | `3209,3` | `1502,9` | `122,7` | `1,0015` | `3862,9` | `67,0` | `1,83×` |
| `phi` | `1,0040` | `1047,8` | `436,6` | `40,1` | **`1,0573`** | `66,2` | `1,15` | **`34,9×`** |
| `sigma_eta` | `1,0033` | `1441,0` | `2217,3` | `55,1` | **`1,1436`** | `23,2` | `0,40` | **`136,8×`** |

Selain kegagalan $\hat{R} > 1{,}05$ pada $(\phi, \sigma_\eta)$, rantai terpusat memiliki `E-BFMI` sangat buruk (`0,0396`–`0,0935` $\ll 0{,}30$, menandakan difusi energi yang macet di leher *funnel*), sedangkan rantai non-terpusat memiliki `E-BFMI` sehat di keempat rantai (`0,9030`–`0,9830`).

Hasil ablasi efisiensi konfigurasi HMC/NUTS ($T = 300$, median $\pm$ MAD dari 5 pengulangan independen):

| Konfigurasi Sampler | Median $\pm$ MAD `min ESS_bulk / detik` | Median $\varepsilon$ Terpakai |
|---|---:|---:|
| `Fixed eps=0.025, M=I` | `57,84 ± 5,68` | `0,02500` |
| `Fixed eps=0.11, M=I` | `16,84 ± 12,87` | `0,11000` |
| `DualAvg eps, M=I` | `100,06 ± 10,53` | `0,05000` |
| **`DualAvg eps + Diag M` (Stan Warmup)** | **`119,34 ± 8,80`** | `0,08728` |

## 6. Output yang Ditargetkan

1. **Kebenaran Numerik Inti (`c1`, `c2`, `c3`)**: Galat relatif gradien analitik terhadap selisih hingga terpusat $\le 10^{-6}$; reversibilitas *leapfrog* $\le 10^{-9}$ dengan kemiringan galat energi $\in [1{,}95, 2{,}05]$; pemulihan rata-rata Gaussian 50-D dalam $\le 4\,\text{SE}$ dan variansi diagonal relatif $\le 0{,}15$.
2. **Kualitas Inferensi SV & Kontras Parameterisasi (`c4`, `c5`, `c7`)**: Pada SV sintetis $T = 1000$, `split-Rhat < 1,01`, `ESS bulk >= 400`, `divergensi == 0`; liputan pita kredibel 90% jalur $h_t$ di dalam $[0{,}85, 0{,}95]$; serta `ESS/detik` NCP melampaui CP.
3. **Kalibrasi SBC (`c6`)**: Pada 200 replikasi prior-prediktif ($T = 250, L = 99$), uji keseragaman histogram $\chi^2$ memiliki $p > 0{,}01$ untuk $(\mu, \phi, \sigma_\eta)$ dan kurva ECDF-difference berada di dalam pita simultan 99%.
4. **Keandalan Rekayasa Notebook & Grafik (`c8`, `c9`, `c10`)**: Eksekusi penuh dari nol selesai $\le 10$ menit dengan hasil numerik identik pada dua eksekusi beruntun, seluruh `assert` lulus, dan 6 pasang berkas grafik `PNG` + `SVG` valid XML.

## 7. Status Pencapaian

Tabel berikut memakai nama kriteria yang persis sama dengan [`outputs/day006e/results.json`](../../outputs/day006e/results.json) (tanpa menyetel ambang maupun parameter pasca-eksekusi):

| Kriteria | Status | Nilai Terukur | Ambang A Priori & Bukti |
|---|---|---|---|
| `c1_analytic_grad_vs_fd` | `tercapai` | `9,577e-9` (NCP: `9,577e-9`, CP: `8,799e-10`; cocok eksak loop vs `lfilter` `< 1e-12`) | Maks galat relatif vs FD terpusat $\le 10^{-6}$ (Sel 4) |
| `c2_leapfrog_reversible_and_order2` | `tercapai` | Reversibilitas `4,441e-16`; kemiringan energi `1,9985`; `\|det J - 1\| = 1,594e-10` | Reversibilitas $\le 10^{-9}$ dan kemiringan $\in [1{,}95, 2{,}05]$ (Sel 7) |
| `c3_nuts_gaussian_50d` | `tercapai` | Maks `\|Z_mean\| = 2,8759`; galat rel var diag `0,0642`; median `1819,6 ESS/s` (`div = 0`) | `\|mean_hat - mean\| <= 4 SE` & rel diag $\le 0{,}15$ (Sel 10) |
| `c4_sv_synthetic_convergence` | `sebagian` | Maks `split-Rhat = 1,0040` (`< 1,01`), min `ESS bulk = 1047,8` (`>= 400`), **`divergences = 15`** (`0,375%`) | `Rhat < 1,01` & `ESS >= 400` lulus; syarat `divergensi == 0` meleset `15/4000` (Sel 15) |
| `c5_sv_latent_h_coverage_90` | `tercapai` | `0,9050` (`905 / 1000` hari tercakup pita kredibel `[q05, q95]`) | Liputan berada tepat di dalam `[0,85, 0,95]` (Sel 13) |
| `c6_sbc_uniformity` | `sebagian` | ECDF simultan 99% **lulus ketiganya** (`<= 0,055 < 0,110`); $\chi^2_{19}$: `phi` `p=0,6169`, `sigma_eta` `p=0,7599`, **`mu` `p=0,0011`** | Pita ECDF & 2 parameter lulus; $\chi^2$ 20-bin `mu` `< 0,01` (Sel 18) |
| `c7_ncp_beats_centered_ess_per_sec` | `tercapai` | Min `ESS/s` NCP `40,06` vs CP `0,40`; pada $\sigma_\eta$ unggul **`136,8×`** (`55,1` vs `0,40 ESS/s`) | `ESS/detik` NCP `>` CP untuk seluruh parameter (Sel 15) |
| `c8_clean_reproducible_under_10min` | `tercapai` | Run 1 = `295,46 s`, Run 2 = `318,23 s` (`<= 600 s`); SHA-256 identik `aaff5c7979c0fe0c...` (`true`) | Dua eksekusi bersih dari nol $\le 600\text{ s}$ & identik persis (Sel 23) |
| `c9_unit_asserts_pass` | `tercapai` | Seluruh sel uji `assert` (`c1, c2, c3, c5, c7, c10`) lulus tanpa galat | `scripts/run_all.py --day 6` keluar dengan kode `0` |
| `c10_figures_valid` | `tercapai` | `6` pasang `PNG` + `SVG` ada, tak kosong, dan valid `ET.parse` | `day006e_{traces_acf,pair_posterior,centered_vs_ncp,sbc_ecdf,sv_latent_path,efficiency}` |

Daftar tautan grafik terverifikasi:
- Jejak rantai & fungsi autokorelasi (NCP vs CP): [day006e_traces_acf.svg](../../outputs/figures/day006e_traces_acf.svg) ([PNG](../../outputs/figures/day006e_traces_acf.png))
- *Pair plot* posterior $(\mu, \phi, \sigma_\eta)$: [day006e_pair_posterior.svg](../../outputs/figures/day006e_pair_posterior.svg) ([PNG](../../outputs/figures/day006e_pair_posterior.png))
- Geometri *Neal's funnel* terpusat vs non-terpusat & `ESS/detik`: [day006e_centered_vs_ncp.svg](../../outputs/figures/day006e_centered_vs_ncp.svg) ([PNG](../../outputs/figures/day006e_centered_vs_ncp.png))
- Histogram peringkat SBC & kurva selisih ECDF berpita simultan 99%: [day006e_sbc_ecdf.svg](../../outputs/figures/day006e_sbc_ecdf.svg) ([PNG](../../outputs/figures/day006e_sbc_ecdf.png))
- Pemulihan jalur laten $h_t$ & aplikasi 4-rezim vs EWMA: [day006e_sv_latent_path.svg](../../outputs/figures/day006e_sv_latent_path.svg) ([PNG](../../outputs/figures/day006e_sv_latent_path.png))
- Profil latensi gradien, adaptasi $\varepsilon$, & ablasi $(\varepsilon, M)$: [day006e_efficiency.svg](../../outputs/figures/day006e_efficiency.svg) ([PNG](../../outputs/figures/day006e_efficiency.png))

## 8. Jika Belum Tercapai atau Belum Diketahui

Sesuai aturan integritas repositori, parameter, ambang, bin, maupun *seed* **tidak disetel** setelah hasil keluar. Dua kriteria berstatus `sebagian` (`c4` dan `c6`) serta satu temuan menarik pada aplikasi rezim diuraikan secara tajam dengan memisahkan **FAKTA** dari **HIPOTESIS**:

### 8.1 Analisis `c4_sv_synthetic_convergence` (`sebagian`: `15` Divergensi dari `4000` Sampel)
- **FAKTA**: Ketiga parameter struktural pada SV sintetis $T = 1000$ mencapai konvergensi sangat kuat (`split-Rhat` berturut-turut `1,0033, 1,0040, 1,0033 < 1,01`; `ESS bulk` `3209,3, 1047,8, 1441,0 >= 400`; `E-BFMI` `0,903`–`0,983`), dan 98,5% transisi berhenti normal di kedalaman pohon `5`–`6` (`117` di `d=5`, `3874` di `d=6`, tidak pernah menyentuh batas `max_depth = 10`). Namun tercatat **`15` transisi divergen** (`0,375%` dari 4000 sampel), sehingga syarat keras `divergensi == 0` berstatus `sebagian`.
- **HIPOTESIS**: Dua mekanisme saling berinteraksi memicu 15 divergensi langka tersebut: (a) ketika rantai menjelajahi ekor kanan $\phi \to 1$ ($u_\phi > 4{,}5$), faktor stasioner $z_1 = \eta_1 / \sqrt{1 - \phi^2}$ dan sensitivitas operator resolven $(I - \phi L)^{-1}$ terhadap $\phi$ menaikkan kelengkungan lokal secara tajam sementara metrik massa hanya bersifat diagonal (mengabaikan korelasi silang antara $u_\phi$, $u_\sigma$, dan $\eta$), dan (b) penyetelan ulang pusat penyusutan *dual averaging* `mu_da = log(5 * eps)` seusai jendela massa terakhir pada *terminal buffer* 160 langkah menyisakan ukuran langkah $\varepsilon \approx 0{,}087$ yang sedikit agresif di ekor $\phi \to 1$. Bukti untuk memastikannya adalah membandingkan matriks massa padat blok $3\times 3$ pada $(\mu, u_\phi, u_\sigma)$ atau menaikkan *terminal buffer* ke $\ge 300$ iterasi dengan $\mu_{\text{da}} = \ln(2\varepsilon)$.

### 8.2 Analisis `c6_sbc_uniformity` (`sebagian`: ECDF Simultan 99% Lulus Penuh, $\chi^2_{19}$ pada `mu` $p = 0{,}0011$)
- **FAKTA**: Dari `200` replikasi SBC ($T = 250, L = 99$, median `ESS` minimum `408,75`), simpangan maksimum kurva ECDF empiris terhadap diagonal seragam adalah **`0,0500` untuk `mu`**, **`0,0550` untuk `phi`**, dan **`0,0550` untuk `sigma_eta`** — ketiganya jauh di bawah batas kritis pita simultan 99% (`0,1100`), tanpa bentuk huruf-U (dispersi terlalu sempit), huruf-$\cap$ (dispersi terlalu lebar), maupun kemiringan monoton (bias lokasi). Uji $\chi^2_{19}$ 20-bin lulus dengan mudah pada `phi` ($\chi^2 = 16{,}60, p = 0{,}6169$) dan `sigma_eta` ($\chi^2 = 14{,}40, p = 0{,}7599$), tetapi menghasilkan $\chi^2 = 43{,}60$ ($p = 0{,}0011 < 0{,}01$) pada `mu` akibat osilasi gergaji antar-bin tetangga (`[..., 5, 19, 3, 8, 15, 17, 6, 15, 3, 18, ...]`, di mana rata-rata tiap pasang bin bersebelahan tetap $\approx 10$–$11$). Tercatat pula `139/200` replikasi prior mengalami minimal satu divergensi ketika prior $\tilde{\phi} \sim \mathrm{Beta}(20, 1{,}5)$ menarik nilai ekstrem dekat $\phi \to 1$ dengan *warmup* pendek `200` langkah.
- **HIPOTESIS**: Kegagalan $\chi^2_{19}$ khusus pada `mu` ketika deviasi sup-norma ECDF-nya justru paling kecil (`0,0500`) merupakan manifestasi langsung kelemahan klasik uji histogram ber-bin terhadap fluktuasi batas bin diskret dan autokorelasi sisa rantai pendek, sebagaimana didokumentasikan oleh Säilynoja, Bürkner & Vehtari (2022) sebagai alasan mengganti $\chi^2$ dengan uji ECDF simultan. Bukti untuk memastikannya adalah menguji ulang peringkat yang sama dengan statistik Cramér–von Mises / Anderson–Darling bebas-bin serta menaikkan *warmup* SBC pada replikasi ber-$\phi > 0{,}99$.

### 8.3 Temuan Aplikasi Lonjakan Rezim: SV Stasioner vs EWMA ($\lambda = 0{,}94$)
- **FAKTA**: Pada seri sintetis 4-rezim ($T = 500$) yang melompat secara tangga deterministik dari $\mu_1 = -9{,}6$ ke $\mu_2 = -6{,}7$, pita kredibel 90% SV menangkap `88,6%` volatilitas sejati, namun pada evaluasi kausal 1-langkah ke depan di 30 titik asal bergulir, rata-rata *log predictive score* EWMA (`2,53479`) sedikit mengungguli SV stasioner (`2,49554`, selisih kumulatif `-1,1773` nats; `RMSE` volatilitas `0,00938` vs `0,01107`).
- **HIPOTESIS**: Model SV yang dipasang mengasumsikan $\mu$ tunggal konstan sepanjang sampel sehingga sesaat sesudah patahan struktural $t = 180$, gaya pemulih (*mean reversion*) menarik prediksi $h_{t+1\mid t} = \mu + \phi(h_t - \mu)$ turun menuju rata-rata pra-krisis, sedangkan EWMA ($\lambda = 0{,}94$, setara IGARCH akar-unit tanpa *mean reversion*) langsung mengikuti level variansi baru tanpa penalti tarik-balik.

### 8.4 Pemeriksaan Bagian A (`day006d`) & Deviasi dari Prompt
- **Pemeriksaan Bagian A (`day006d`)**: Pemeriksaan pada [day006d_sim_heston-adi-cpp.md](day006d_sim_heston-adi-cpp.md) (`HEAD 9e5344d`) mengonfirmasi bahwa ketiga butir audit — (1) kontribusi utama 10 *thread* pada `c5` (`12,49×` vs `V4` 1-*thread* `2,82×`), (2) pelambatan `V2` terhadap `V1` di `801×401` (`1,56×` vs `1,83×`), dan (3) datarnya pemindaian `tile` (`0,9829`–`0,9928 s`) — telah dinyatakan sebagai **FAKTA** terukur beserta penjelasan mekanismenya sebagai **HIPOTESIS** di §5.3 dan §8.1–§8.3, sehingga sesuai instruksi tidak dibuat commit tambahan untuk Bagian A.
- **Deviasi Kecil & Pilihan Teknis**:
  - Untuk mempercepat `200` replikasi SBC dan `30` jendela prediksi bergulir agar seluruh notebook berjalan dalam $\approx 5$ menit ($\le 10$ menit), replikasi independen dieksekusi paralel menggunakan `multiprocessing` (`fork`) pustaka standar dengan *seed* per-replikasi yang diturunkan di muka dari `SeedSequence.spawn` dan diurutkan kembali menurut indeks replikasi, sehingga hasilnya **100% identik bit-demi-bit** antara eksekusi sekuensial maupun paralel serta antar-eksekusi berulang.
  - Pada perbandingan parameterisasi terpusat (CP) di Bagian 6, kedalaman pohon maksimum CP dibatasi pada `max_depth = 8` (`256` langkah *leapfrog*/iterasi) agar kemacetan leher *funnel* tidak melampaui anggaran waktu 10 menit.
  - Item `090` di `docs/ROADMAP.md` (*MCMC Metropolis–Hastings*) sengaja **tidak dicentang** karena sesuai `RULES.md` sesi berhuruf `day006e` berstatus `topic_no: null` (*bonus/off-roadmap*) dan urutan pencentangan topik reguler masih berada di `003`.

## 9. Landasan Teori

1. **Inferensi Posterior Bayesian pada Model Ruang Keadaan Volatilitas Stokastik**: Model Taylor (1986) memperlakukan log-variansi $h_t$ sebagai proses AR(1) Gaussian laten yang memodulasi skala pengamatan non-linier $y_t = e^{h_t/2}\varepsilon_t$, menghasilkan kerapatan pengamatan log-konkaf dalam $h_t$ karena $\frac{\partial^2}{\partial h_t^2}\bigl(-\frac{1}{2}h_t - \frac{1}{2}y_t^2 e^{-h_t}\bigr) = -\frac{1}{2}y_t^2 e^{-h_t} \le 0$ — diimplementasikan pada `sv_ncp_neg_log_post_and_grad` di [`day006e_sim_bayesian-sv-nuts.ipynb`](../../days/day006e_sim_bayesian-sv-nuts.ipynb).
2. **Geometri Corong Neal (Neal's Funnel) dan Parameterisasi Non-Terpusat (NCP)**: Pada model hierarkis $h \mid \theta \sim \mathcal{N}(\mu\mathbf{1}, \sigma_\eta^2 \Sigma_\phi)$, volume himpunan tipikal menyusut secara eksponensial ketika $\sigma_\eta \to 0$ (Neal, 2003; Papaspiliopoulos, Roberts & Sköld, 2007); memetakan ulang variabel laten ke inovasi putih $\eta \sim \mathcal{N}(0, I_T)$ melalui $h = \mu\mathbf{1} + \sigma_\eta R_\phi \eta$ memutus ketergantungan kelengkungan prior terhadap $\sigma_\eta$ maupun $\phi$ — dibuktikan oleh kenaikan `ESS/detik` $\sigma_\eta$ sebesar `136,8×` pada [`day006e_centered_vs_ncp.svg`](../../outputs/figures/day006e_centered_vs_ncp.svg).
3. **Transformasi Ruang Tanpa Batas dan Koreksi Determinan Jacobian**: Agar dinamika Hamiltonian bergerak bebas tanpa pantulan dinding pada domain terbatas $\tilde{\phi} \in (0, 1)$ dan $\sigma_\eta \in (0, \infty)$, diterapkan difeomorfisme $u_\phi = \ln\frac{\tilde{\phi}}{1-\tilde{\phi}}$ dan $u_\sigma = \ln\sigma_\eta$ dengan penambahan suku $\ln\left|\det \frac{\partial(\tilde{\phi}, \sigma_\eta)}{\partial(u_\phi, u_\sigma)}\right| = \ln\tilde{\phi} + \ln(1-\tilde{\phi}) + u_\sigma$ ke dalam log-target — diimplementasikan di `unpack_sv_theta` dan diverifikasi tak-bias oleh SBC di [`outputs/day006e/results.json`](../../outputs/day006e/results.json).
4. **Diferensiasi Mode-Balik Manual lewat Operator Adjoin Resolven AR(1)**: Pemetaan kausal $z = (I - \phi L)^{-1}\tilde{\eta}$ merupakan filter rekursif orde-1 bawah-segitiga; transpos adjoinnya $(I - \phi L)^{-\top}$ adalah filter anti-kausal mundur $\lambda_t = \bar{z}_t + \phi\lambda_{t+1}$ yang menghitung seluruh vektor $\nabla_\eta \ell$ dan skalar $\nabla_\phi \ell$ dalam satu sapuan balik $O(T)$ — diverifikasi hingga `9,58e-9` terhadap selisih hingga terpusat di Sel 4 dan dipercepat `8,88×` dengan `scipy.signal.lfilter` pada [`day006e_efficiency.svg`](../../outputs/figures/day006e_efficiency.svg).
5. **Dinamika Hamiltonian dan Integrator Simplektik Leapfrog**: Mengangkat ruang target $\mathbb{R}^D$ ke ruang fasa $\mathbb{R}^{2D}$ dengan Hamiltonian $H(q,p) = U(q) + \frac{1}{2}p^\top M^{-1}p$ mendefinisikan aliran konservatif $\dot{q} = \nabla_p H$, $\dot{p} = -\nabla_q H$; skema belah-tiga *leapfrog* merupakan integrator simplektik eksak untuk *shadow Hamiltonian* terdekat $\tilde{H} = H + O(\varepsilon^2)$ sehingga melestarikan ukuran Liouville $\mathrm{d}q\,\mathrm{d}p$ dan menjaga galat energi tetap terbatas tanpa hanyutan sekuler — dibuktikan oleh `|det J - 1| = 1,59e-10` dan kemiringan `1,9985` di Sel 7.
6. **Keseimbangan Terperinci (Detailed Balance) dan Invariansi Distribusi Kanonik**: Karena pemetaan $\Phi_\varepsilon^L$ bersifat involutif bersama pembalikan momentum $\mathcal{R}(q,p) = (q,-p)$ ($\mathcal{R}\circ\Phi_\varepsilon^L\circ\mathcal{R}\circ\Phi_\varepsilon^L = \mathrm{id}$, terukur `4,44e-16`) dan melestarikan volume, pembobotan keadaan dengan faktor Boltzmann $e^{-H(q,p)}$ menjamin distribusi stasioner marginal atas $q$ persis sama dengan posterior target $\pi(q) \propto e^{-U(q)}$.
7. **Multinomial NUTS, Kriteria U-Turn Generalisasi, dan Pencabangan Progresif Terbias**: Perluasan biner simetris maju/mundur secara acak membangun himpunan kandidat lintasan hingga kondisi $\rho^\top M^{-1} p_\pm \le 0$ mendeteksi bahwa lintasan mulai berbalik arah dalam metrik $M$ (Hoffman & Gelman, 2014; Betancourt, 2017); penarikan multinomial internal dikombinasikan dengan peluang terima luar $\min(1, w_{\text{ext}}/w_{\text{old}})$ meningkatkan efisiensi sampel efektif di atas jumlah sampel mentah sebagaimana terlihat pada Gaussian 50-D (`ESS > 3000` dari `4000` sampel, `1819,6 ESS/s`) di `build_nuts_tree` dan `nuts_step`.
8. **Adaptasi Ukuran Langkah Dual Averaging Stokastik**: Skema optimisasi stokastik primal-dual Nesterov yang diadaptasi Hoffman & Gelman (2014) memperbarui rata-rata galat penerimaan $\bar{H}_m = (1 - \frac{1}{m+t_0})\bar{H}_{m-1} + \frac{1}{m+t_0}(\delta - \alpha_m)$ dan $\ln\varepsilon_m = \mu_{\text{da}} - \frac{\sqrt{m}}{\gamma}\bar{H}_m$ dengan penghalusan Polyak–Ruppert $\overline{\ln\varepsilon}_m$, menjamin konvergensi stabil $\alpha_m \to \delta$ tanpa penyetelan tangan — diimplementasikan di `run_nuts_chain`.
9. **Prekondisi Metrik Euclidean Diagonal Berjendela**: Memilih $M^{-1} = \mathrm{diag}(\hat{\sigma}_1^2, \dots, \hat{\sigma}_D^2)$ setara dengan merotasi skala koordinat sehingga variansi posterior marginal tiap komponen mendekati satuan; estimasi variansi sampel ber-penyusutan (*shrinkage* ke $10^{-3}$ dengan bobot $5/(n_w+5)$) pada jendela yang melebar diikuti pembekuan ketat pasca-*warmup* meningkatkan `min ESS/s` dari `100,06` ($M=I$) menjadi `119,34` pada tabel ablasi Sel 21.
10. **Diagnostik Split-Rhat Ternormalisasi Peringkat, ESS, dan E-BFMI**: Pemecahan tiap rantai menjadi dua paruh mendeteksi ketidakstasioneran intra-rantai, sedangkan transformasi peringkat kuantil normal Blom dan pelipatan sekitar median (*folded split-$\hat{R}$*) mendeteksi perbedaan skala ekor bahkan ketika distribusi ber-ekor berat (Vehtari et al., 2021); di sisi lain $\text{E-BFMI} = \frac{\mathbb{E}[(E_n - E_{n-1})^2]}{\mathrm{Var}(E)}$ mendiagnosis kegagalan transisi momentum menembus penghalang energi *funnel* (`0,040` pada CP vs `0,983` pada NCP) — diimplementasikan pada `split_rhat_and_ess` dan `compute_ebfmi`.
11. **Teorema Kalibrasi Berbasis Simulasi (SBC) dan Uji Selisih ECDF Simultan**: Identitas $\int p(\theta \mid y)\,p(y, \theta^\star)\,\mathrm{d}y\,\mathrm{d}\theta^\star = \int p(\theta \mid y)p(\theta^\star \mid y)p(y)\,\mathrm{d}y\,\mathrm{d}\theta^\star$ membuktikan bahwa $\theta^\star$ dan sampel posterior $\theta^{(1:L)}$ bersifat dapat dipertukarkan (*exchangeable*) bersyarat pada $y$, sehingga peringkat $r(\theta^{(1:L)}, \theta^\star)$ wajib seragam diskret (Talts et al., 2018); amplop simultan kuantil ke-99 dari $\sup_z |\hat{F}(z) - z|$ mengendalikan laju positif palsu uji ganda di seluruh titik evaluasi (Säilynoja et al., 2022) — diimplementasikan pada Bagian 7 dan [`day006e_sbc_ecdf.svg`](../../outputs/figures/day006e_sbc_ecdf.svg).
12. **Aturan Skor Wajar Ketat (Strictly Proper Scoring Rule) Log Predictive Score**: Skor logaritmik prediktif satu-langkah $\ln \int p(y_{t+1} \mid h_{t+1})p(h_{t+1}, \theta \mid y_{1:t})\,\mathrm{d}h_{t+1}\mathrm{d}\theta$ mengevaluasi ketajaman dan kalibrasi distribusi prediktif penuh secara serentak (Gneiting & Raftery, 2007), mengungkap kompromi nyata antara model struktural ber-*mean-reversion* tetap (SV) dan penyaring eksponensial akar-unit (EWMA) ketika seri mengalami lompatan rezim mendadak pada Bagian 8 dan [`day006e_sv_latent_path.svg`](../../outputs/figures/day006e_sv_latent_path.svg).

## 10. Kesimpulan dan Batasan

1. **Kesimpulan Utama**: Implementasi mandiri HMC/NUTS multinomial dengan gradien adjoin resolven AR(1), *dual averaging*, dan matriks massa diagonal berjendela terbukti akurat secara matematis (`c1` `9,58e-9`, `c2` reversibilitas `4,44e-16` & orde `1,9985`), sangat cepat pada dimensi tinggi (`1819,6 ESS/s` di Gaussian 50-D; eksekusi penuh notebook 25 sel termasuk 200 replikasi SBC selesai dalam $\approx 5$ menit dan 100% deterministik lintas eksekusi), memulihkan lintasan volatilitas laten dengan liputan `90,5%` (target `90%`), serta memberikan lonjakan efisiensi **`136,8×`** atas parameterisasi terpusat.
2. **Temuan Diagnostik Jujur**: Kriteria `c4` dan `c6` tercatat `sebagian` secara transparan tanpa rekayasa pasca-uji: terdapat `15/4000` divergensi di ekor stasioner $\phi \to 1$ pada $T=1000$ meskipun $\hat{R} \le 1{,}004$ dan $\text{ESS} \ge 1048$, serta osilasi bin tetangga $\chi^2_{19}$ pada `mu` ($p = 0{,}0011$) meskipun seluruh tiga kurva ECDF SBC berada jauh di dalam pita simultan 99% (`<= 0,055 < 0,110`).
3. **Pekerjaan Lanjutan**: Menambahkan blok metrik massa padat pada subruang hiperparameter $(\mu, u_\phi, u_\sigma)$, mengganti inovasi pengamatan Gaussian dengan Student-$t$ ber-efek *leverage* $\rho\,\mathrm{d}W^{(1)}\mathrm{d}W^{(2)}$, dan memasukkan dinamika perpindahan rezim Markov pada $\mu_t$ untuk mengatasi bias tarik-balik saat terjadi krisis struktural.

## 11. Referensi dan Jejak Bukti

### Literatur Eksternal
- Neal, R. M. (2003). *Slice sampling* (termasuk contoh *Neal's funnel*). Annals of Statistics, 31(3), 705–767.
- Neal, R. M. (2011). *MCMC using Hamiltonian dynamics*. Dalam *Handbook of Markov Chain Monte Carlo* (hlm. 113–162). Chapman & Hall/CRC.
- Hoffman, M. D., & Gelman, A. (2014). *The No-U-Turn Sampler: Adaptively setting path lengths in Hamiltonian Monte Carlo*. Journal of Machine Learning Research, 15, 1593–1623.
- Betancourt, M. (2017). *A conceptual introduction to Hamiltonian Monte Carlo*. arXiv:1701.02434.
- Papaspiliopoulos, O., Roberts, G. O., & Sköld, M. (2007). *A general framework for the parametrization of hierarchical models*. Statistical Science, 22(1), 59–73.
- Kim, S., Shephard, N., & Chib, S. (1998). *Stochastic volatility: Likelihood inference and comparison with ARCH models*. Review of Economic Studies, 65(3), 361–393.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). *Validating Bayesian inference algorithms with simulation-based calibration*. arXiv:1804.06788.
- Säilynoja, T., Bürkner, P.-C., & Vehtari, A. (2022). *Graphical test for discrete uniformity and its applications in goodness-of-fit evaluation and multiple sample comparison*. Statistics and Computing, 32(2), 32.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). *Rank-normalization, folding, and localization: An improved $\hat{R}$ for assessing convergence of MCMC*. Bayesian Analysis, 16(2), 667–718.
- Gneiting, T., & Raftery, A. E. (2007). *Strictly proper scoring rules, prediction, and estimation*. Journal of the American Statistical Association, 102(477), 359–378.

### Jejak Artefak Lokal Terverifikasi
- Notebook tereksekusi 2$\times$ identik (`SHA-256: aaff5c7979c0fe0c...`, `0` sel galat): [`days/day006e_sim_bayesian-sv-nuts.ipynb`](../../days/day006e_sim_bayesian-sv-nuts.ipynb).
- Ringkasan metrik & evaluasi kriteria: [`outputs/day006e/results.json`](../../outputs/day006e/results.json).
- Enam pasang grafik `PNG` + `SVG`: [`day006e_traces_acf.svg`](../../outputs/figures/day006e_traces_acf.svg), [`day006e_pair_posterior.svg`](../../outputs/figures/day006e_pair_posterior.svg), [`day006e_centered_vs_ncp.svg`](../../outputs/figures/day006e_centered_vs_ncp.svg), [`day006e_sbc_ecdf.svg`](../../outputs/figures/day006e_sbc_ecdf.svg), [`day006e_sv_latent_path.svg`](../../outputs/figures/day006e_sv_latent_path.svg), dan [`day006e_efficiency.svg`](../../outputs/figures/day006e_efficiency.svg).
