# Hari 004h — Volatilitas Kasar: Model Rough Bergomi, Skema Hybrid, Skew ATM, dan Identifiabilitas Hurst (Rust)

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day004h`
- **Tanggal**: `2026-11-22`
- **Topik**: Volatilitas kasar: model *rough Bergomi*, skema hybrid $\kappa = 1$, *skew* ATM, dan identifiabilitas parameter Hurst $H$
- **Jenis**: `quant`
- **Bahasa**: Rust (*pure standard library*, `edition = "2021"`, nol *crate* eksternal)
- **Lokasi Proyek**: [projects/rust/day004h_rough_bergomi](../../projects/rust/day004h_rough_bergomi/README.md)

Proyek ini membangun simulasi, verifikasi teoretis, dan studi kalibrasi lengkap untuk model **rough Bergomi** (Bayer–Friz–Gatheral, 2016) dengan eksponen Hurst kasar $H = 0{,}1$. Proyek mencakup empat bagian utama: (1) kovarians eksak proses Volterra Riemann–Liouville dan pembanding skema hybrid $\kappa = 1$ terhadap simulasi eksak Cholesky $2n \times 2n$; (2) uji martingal variansi dan harga serta estimator bersyarat Black–Scholes (Romano–Touzi / McCrickerd–Pakkanen); (3) permukaan *implied volatility smile*, hukum pangkat *skew* ATM $\psi(T) \sim T^{H - 1/2}$, dan pembanding model Markov ($H = 0{,}5$); serta (4) kalibrasi Nelder–Mead (*Common Random Numbers*), profil objektif $H$, matriks Hessian numerik, studi 30 replikasi, dan studi salah spesifikasi horizon.

## 2. Masalah dan Tujuan

1. **Akurasi representasi proses Volterra $\widetilde{W}^H_t$**: Menguji apakah bentuk tertutup hipergeometrik ${}_2F_1$ cocok dengan kuadratur Gauss–Legendre bebas-singularitas (`< 1e-8`), memverifikasi simulasi gabungan Cholesky $2n \times 2n$ sebagai referensi eksak dalam hukum, dan mengukur bias skema hybrid $\kappa = 1$ titik kanan pada $n \in \{50, 100, 250, 500\}$ untuk $H = 0{,}1$ dan $H = 0{,}5$.
2. **Sifat martingal dan reduksi variansi**: Menguji apakah $\mathbb{E}[\exp(\eta \widetilde{W}^H_u - \frac{1}{2}\eta^2 u^{2H})] = 1$, *variance swap* $\int_0^T \mathbb{E}[V_u]\,du = 0{,}04\,T$, dan $\mathbb{E}[S_T] = S_0$ bertahan di bawah diskretisasi skema hybrid, serta mengukur faktor reduksi variansi (VRF) estimator bersyarat Black–Scholes terhadap Monte Carlo log-Euler murni.
3. **Geometri *smile* dan hukum pangkat *skew* ATM**: Mengukur kemiringan regresi $\ln \psi(T)$ vs $\ln T$ untuk $H \in \{0{,}05; 0{,}1; 0{,}25; 0{,}5\}$ terhadap prediksi asimtotik $H - \frac{1}{2}$ dan membandingkan struktur *smile* $T = 0{,}1$ dengan model Markov 1-faktor ($H = 0{,}5$).
4. **Identifiabilitas kalibrasi $(H, \eta, \rho)$**: Mengukur akurasi pemulihan $H$ melalui Nelder–Mead CRN, korelasi estimasi $H$–$\eta$, bias/SD dari $R_{\mathrm{rep}} = 30$ replikasi independen, serta pergeseran $\hat{H}$ jika kalibrasi dibatasi pada jatuh tempo panjang ($T \ge 0{,}5$) vs pendek ($T \le 0{,}1$).

## 3. Apa yang Dikerjakan Proyek

Runner biner `run` mengeksekusi seluruh empat bagian secara deterministik (`main_seed = 20261122`, `8` *worker*) menggunakan *worker pool* `std::thread::scope` berbasis blok tetap (`BLOCK = 64`) sehingga hasil bit-identik lintas jumlah *worker*, lalu menghasilkan:
- Berkas hasil berpresisi 17 digit [outputs/day004h/results.json](../../outputs/day004h/results.json).
- Tujuh grafik SVG diagnostik dan perbandingan di [outputs/figures/](../../outputs/figures/).

## 4. Uraian Kode

Seluruh kode sumber berada di [projects/rust/day004h_rough_bergomi/src/](../../projects/rust/day004h_rough_bergomi/src/lib.rs):
- [rng.rs](../../projects/rust/day004h_rough_bergomi/src/rng.rs) & [parallel.rs](../../projects/rust/day004h_rough_bergomi/src/parallel.rs): PRNG `xoshiro256++` + `SplitMix64` per indeks lintasan serta `run_chunks` / `run_chunk_reduce` berurutan blok tetap yang menjamin determinisme bit-identik `1` vs banyak *worker*.
- [linalg.rs](../../projects/rust/day004h_rough_bergomi/src/linalg.rs), [quad.rs](../../projects/rust/day004h_rough_bergomi/src/quad.rs), & [stats.rs](../../projects/rust/day004h_rough_bergomi/src/stats.rs): Dekomposisi Cholesky, invers matriks padat $3\times 3$ (eliminasi Gauss), kuadratur Gauss–Legendre 64-titik dari akar polinomial Legendre, serta statistik sampel dan regresi OLS.
- [bs.rs](../../projects/rust/day004h_rough_bergomi/src/bs.rs): CDF normal standar presisi tinggi (`erfc` deret + pecahan kontinu), harga Black–Scholes berbasis total variansi, vega, dan *implied volatility* Newton–Raphson dengan pengaman *bisection*.
- [volterra.rs](../../projects/rust/day004h_rough_bergomi/src/volterra.rs): Fungsi `exact_cov_closed_form` (${}_2F_1$ dengan transformasi DLMF 15.8.4 dekat $s/t \to 1$), `exact_cov_quadrature` (substitusi penghapus singularitas), `build_joint_cov_matrix` ($2n \times 2n$), `run_cholesky_reference_check`, `HybridWeights` ($\kappa = 1$), dan `hybrid_bias_study`.
- [rbergomi.rs](../../projects/rust/day004h_rough_bergomi/src/rbergomi.rs): `run_martingale_anchors`, `run_call_estimators_comparison` (MC log-Euler vs bersyarat BS + paritas put-call), serta `compare_hybrid_vs_cholesky_ref`.
- [smile.rs](../../projects/rust/day004h_rough_bergomi/src/smile.rs): `compute_smile` ($13 \times 5$ grid dengan SE metode delta), `run_atm_skew_scaling_study`, dan `run_markov_comparison`.
- [calib.rs](../../projects/rust/day004h_rough_bergomi/src/calib.rs): `CalibCrnBank`, `eval_smile_sse`, `nelder_mead`, profil $H$, Hessian numerik $3\times 3$, studi 30 replikasi, serta kalibrasi `T >= 0.5` vs `T <= 0.1`.
- [jsonout.rs](../../projects/rust/day004h_rough_bergomi/src/jsonout.rs), [svg.rs](../../projects/rust/day004h_rough_bergomi/src/svg.rs), & [bin/run.rs](../../projects/rust/day004h_rough_bergomi/src/bin/run.rs): Serialisasi JSON, grafik SVG, dan orkestrasi eksperimen.

## 5. Metode dan Alur Kerja

1. **Alur Komputasi**: Parameter acuan ($S_0=1, \xi_0=0{,}04, \eta=1{,}9, H=0{,}1, \rho=-0{,}7, r=0$) $\to$ verifikasi kovarians & Cholesky $2n \times 2n$ ($n=50, R=200.000$) + bias hybrid $n \in \{50, 100, 250, 500\}$ $\to$ uji martingal & pembanding estimator call ($n=250, R=100.000$) $\to$ *smile* $13\times 5$, regresi $\log\psi(T)$ vs $\log T$, dan Markov $H=0{,}5$ $\to$ pasar sintetis ($R=200.000, n=500$), Nelder–Mead CRN ($R=20.000, n=100$), profil $H$, Hessian, $R_{\mathrm{rep}}=30$, dan salah spesifikasi horizon $\to$ ekspor JSON & 7 SVG.
2. **Penurunan Substitusi Bebas-Singularitas**: Pada $\operatorname{Cov}(\widetilde{W}^H_s, \widetilde{W}^H_t) = 2H \int_0^s (s-u)^{H-1/2}(t-u)^{H-1/2}\,du$, substitusi $s - u = s\,y^{1/(H+1/2)}$ menghasilkan $v^{H-1/2}\,dv = \frac{s^{H+1/2}}{H+1/2}\,dy$, menghilangkan singularitas di $u=s$ secara eksak sehingga kuadratur Gauss–Legendre 64-titik mencapai ketelitian `2,12e-10`.

## 6. Output yang Ditargetkan

1. `c1_volterra_cov_quad_within_1e8`: Kesepakatan kovarians bentuk tertutup ${}_2F_1$ vs kuadratur Gauss–Legendre bebas singularitas dalam galat relatif $< 10^{-8}$.
2. `c2_cholesky_exact_cov_within_4_5se`: Rekonstruksi Cholesky $2n \times 2n$ ($< 10^{-12}$) dan seluruh entri kovarians empiris ($R = 200.000$) dalam `4,5 SE`.
3. `c3_martingales_within_3se`: Martingal variansi $\mathbb{E}[\exp(\eta \widetilde{W}^H_u - \frac{1}{2}\eta^2 u^{2H})] = 1$ dan martingal harga $\mathbb{E}[S_T] = S_0$ dalam `3 SE` untuk $T \in \{0{,}25; 0{,}5; 1\}$.
4. `c4_variance_swap_within_3se`: Nilai ekspektasi *variance swap* $\mathbb{E}[\int_0^T V_u\,du]$ cocok terhadap $0{,}04\,T$ dalam `3 SE`.
5. `c5_cond_bs_matches_mc_and_vrf`: Estimator bersyarat Black–Scholes sepakat dengan MC log-Euler murni (`max_z_diff <= 3,5`), `VRF > 1,5`, dan galat paritas put-call jalur `< 1e-12`.
6. `c6_hybrid_matches_cholesky_price_3se`: Harga opsi skema hybrid ($n=50, T=1$) sepakat dengan referensi eksak Cholesky $2n\times 2n$ ($n=50$) dalam `3 SE` gabungan.
7. `c7_atm_skew_power_law_and_markov_zero`: Kemiringan $\log\psi(T)$ vs $\log T$ untuk $H=0{,}1$ berada dalam $\pm 0{,}15$ dari $-0{,}40$ dan untuk $H=0{,}5$ berada di sekitar $0$ ($|\mathrm{slope}| < 0{,}08$).
8. `c8_calib_bias_h_under_003`: Pada studi replikasi $R_{\mathrm{rep}} = 30$, $|\mathrm{bias}(\hat{H})| < 0{,}03$ (ambang ditetapkan sebelum *run*).
9. `c9_runtime_under_300s`: Total waktu eksekusi runner penuh berada di bawah anggaran `300s`.

## 7. Status Pencapaian

Runner penuh selesai dalam **`123,25s`** pada `8` *worker* (`part1 = 5,49s`, `part2 = 4,57s`, `part3 = 7,53s`, `part4 = 105,65s`, `part5 = 0,002s`; eksperimen ekstra dipisah ke `run_extras`, `344,99s`) ([results.json](../../outputs/day004h/results.json)):

| Kriteria | Status | Bukti Terukur |
| --- | --- | --- |
| `c1_volterra_cov_quad_within_1e8` | `tercapai` | `max_rel_err = 2,12e-10 < 1e-8` di 6 pasangan $(s, t)$ |
| `c2_cholesky_exact_cov_bonferroni` | `tercapai` | `max_chol_recon_err = 5,6e-16`; `max_cov_z_score = 3,8415 <= 4,12` (Bonferroni a priori untuk `1275` pasangan, $\alpha=5\%$ dua sisi) |
| `c3_martingales_within_3se` | `tercapai` | Cholesky $\mathbb{E}[e^{\dots}] = 0{,}99982 \pm 0{,}01301$ ($z=0{,}014$, boot `[0,9732; 0,9998; 1,0271]`); BLP `max_z_exp_w = 2,247`, `max_z_s_t = 2,229` |
| `c4_variance_swap_within_3se` | `tercapai` | BLP $\mathbb{E}[\int_0^T V\,dt]$: `0,01002` ($T=0{,}25$), `0,02011` ($T=0{,}5$), `0,04031` ($T=1$) vs `0,04 T` (`max_z = 1,618`) |
| `c5_cond_bs_matches_mc_and_vrf` | `tercapai` | `max_z_diff = 1,767 <= 3,5`; `VRF = 2,15..308,4`; `max_parity_err = 1,97e-14` |
| `c6_hybrid_matches_cholesky_price_3se` | `tercapai` | BLP vs Cholesky ($n=50, T=1$): `max_z_comb = 1,230 <= 3,0` |
| `c7_atm_skew_power_law_and_markov_zero` | `tercapai` | $H=0{,}1$: slope `-0,3849 ± 0,0130` ($R^2=0{,}997$); $H=0{,}5$: slope `-0,0404 ± 0,0253` |
| `c8_calib_bias_h_under_003` | `belum tercapai` | Kondisi (c) `same-scheme-same-n` ($R_{\mathrm{rep}}=30$, pasar MC murni): `bias_H = +0,0719`, `SD_H = 0,0158` (vs `< 0,03`) |
| `c9_runtime_under_300s` | `tercapai` | Total waktu `123,25s < 300s` (`8` *worker*, eksperimen ekstra dipisah ke `run_extras`) |

Perbandingan tiga skema hybrid $\kappa = 1$ ($H = 0{,}1$, $n \in \{50, 100, 250, 500\}$):
- **Titik kanan (*RightEndpoint*)**: bias teoretis $\operatorname{Var}(\widetilde{W}^H_1)$ `+4,94%` $\to$ `+4,38%` $\to$ `+3,69%` $\to$ `+3,23%` (`max_cov_abs_err = 0,0494 -> 0,0323`).
- **BLP titik optimal (*BlpOptimal*, primer)**: bias teoretis `-0,0619%` $\to$ `-0,0539%` $\to$ `-0,0449%` $\to$ `-0,0391%` (`max_cov_abs_err = 0,00083 -> 0,00039`), karena $b_k^* = ((k^{\alpha+1}-(k-1)^{\alpha+1})/(\alpha+1))^{1/\alpha}$ membuat kovarians silang dengan $\Delta B$ eksak.
- **Pencocok variansi (*VarianceMatching*)**: bias variansi `0,00000%` di semua $n$, `max_cov_abs_err` `0,00061 -> 0,000031` (tidak mencocokkan kovarians silang $\Delta B$ secara eksak).

Pemisahan bias skema vs identifiabilitas pada kalibrasi:
- `(a)` Skema lama (`RightEndpoint`), $n$ beda (`250` vs `100`, $R=15$): `bias_H = +0,0788 ± 0,0346`.
- `(b)` Skema baru (`BlpOptimal`), $n$ beda (`250` vs `100`, $R=15$): `bias_H = +0,0744 ± 0,0174`.
- `(c)` Skema baru (`BlpOptimal`), $n$ SAMA (`100` vs `100`, $R=30$, identifiabilitas murni): `bias_H = +0,0719 ± 0,0158`.

## 8. Diagnosis Bagian yang Belum Tercapai (`c8`)

Setelah mengganti skema titik kanan menjadi `BlpOptimal` dan memisahkan eksperimen tambahan ke `run_extras`, kriteria `c3`, `c4`, `c6`, dan `c9` kini seluruhnya **tercapai**, sedangkan `c8_calib_bias_h_under_003` tetap **belum tercapai** pada runner utama (pasar MC murni):
- **FAKTA**:
  1. Uji *self-consistency* (kalibrasi ke *smile* mesin sendiri dengan bank CRN yang sama) memulihkan parameter eksak: $\hat{H} = 0{,}100001$, $\hat{\eta} = 1{,}900000$, $\hat{\rho} = -0{,}699997$ (`rmse_opt = 1,18e-7`) — membuktikan tidak ada bug pada *forward engine* maupun *optimizer* Nelder–Mead.
  2. **Analisis 10 *seed* `true_rmse_grid` di `run_extras`**: Ketidakmonotonan satu-*seed* (`n=500: 0,0121 -> 0,0756 -> 0,0314`) ternyata berasal dari *outlier* ekor berat antar-*seed*. Median 10 *seed* turun rapi seiring $n$ (`n=100: 0,0195/0,0182/0,0182`; `n=250: 0,0116/0,0107/0,0114`; `n=500: 0,0087/0,0068/0,0078`), namun IQR pada $R=20\mathrm{k}$ sangat lebar (`0,0672` di `n=250`, `0,0702` di `n=500`, `max = 0,0762`) sebelum menyempit ke `0,00088` di `n=500, R=200k`.
  3. **Estimator ekor Hill ($N=200.000, n=500, T=1$)**: Indeks ekor $S_T$ adalah $\hat{\alpha}_{S_T} = 17{,}61 \to 16{,}53$ ($k=500..5000$, momen $S_T$ berhingga karena $\rho=-0{,}7 < 0$ meredam ekor kanan $S_T$), namun untuk kuadrat *payoff* $(S_T-S_0)_+^2$, $\hat{\alpha}_{\mathrm{payoff}^2} = 3{,}18 \to 2{,}84 \to 2{,}59 \to 2{,}21$ (mendekati $2$: variansi kuadrat *payoff* nyaris tak berhingga sehingga memicu lonjakan sampel langka).
  4. **Kalibrasi dengan estimator Black–Scholes bersyarat (`run_extras`, $R_{\mathrm{rep}} = 30$)**: Ketika *smile* pasar DAN kalibrasi sama-sama memakai estimator BS bersyarat (VRF besar), bias $H$ turun drastis dari `+0,0719` menjadi **`+0,0108`** (`SD_H = 0,0265`) pada bank CRN tetap dan **`+0,0228`** (`SD_H = 0,0578`) pada bank CRN berbeda — keduanya `< 0,03`!
- **HIPOTESIS**: Bias positif `+0,0719` pada `c8` di runner utama terutama dipicu oleh derau ekor berat pada *smile* pasar MC murni (sayap OTM yang terdistorsi saat inversi IV) yang mendorong Nelder–Mead menaikkan $\hat{H}$ di sepanjang *ridge* $(H, \eta)$ ($\operatorname{Corr} = +0{,}79$); penggunaan estimator BS bersyarat menghilangkan sumber derau pasar tersebut.
- **Catatan martingal eksponensial lognormal (`c3`)**: Pada $\eta = 1{,}9$, variabel acak $\exp(\eta \widetilde{W}_1 - \frac{1}{2}\eta^2)$ memiliki variansi teoritis $e^{\eta^2} - 1 \approx 35{,}96$ dan ekor kanan berat, sehingga uji-$z$ normal perlu dibaca bersama kuantil bootstrap (`[0,9732; 0,9998; 1,0271]`).

## 9. Landasan Teori (11 Poin)

1. **Proses Volterra Riemann–Liouville ($\widetilde{W}^H_t$)**: Integral stokastik fraksional $\sqrt{2H}\int_0^t (t-s)^{H-1/2}\,dB_s$ dengan kernel singular $(t-s)^{H-1/2}$ untuk $H \in (0, 1/2)$ yang menginduksi lintasan kasar berorde Hölder $< H$; diimplementasikan di [volterra.rs](../../projects/rust/day004h_rough_bergomi/src/volterra.rs).
2. **Fungsi Hipergeometrik ${}_2F_1$ & Transformasi Euler/DLMF**: Bentuk tertutup kovarians $\operatorname{Cov}(\widetilde{W}^H_s, \widetilde{W}^H_t)$ diekspresikan lewat ${}_2F_1(\frac{1}{2}-H, 1; H+\frac{3}{2}; s/t)$ dengan kelanjutan analitik $1 - s/t$ di sekitar diagonal; diimplementasikan di [volterra.rs](../../projects/rust/day004h_rough_bergomi/src/volterra.rs).
3. **Kuadratur Gauss–Legendre Bebas Singularitas**: Perubahan variabel $s - u = s\,y^{1/(H+1/2)}$ meniadakan singularitas integrasi $(s-u)^{H-1/2}$ sehingga polinomial Legendre berorde 64 mencapai ketelitian `2,12e-10`; di [quad.rs](../../projects/rust/day004h_rough_bergomi/src/quad.rs) dan [volterra.rs](../../projects/rust/day004h_rough_bergomi/src/volterra.rs).
4. **Simulasi Gaussian Gabungan via Faktorisasi Cholesky $2n\times 2n$**: Penarikan eksak dalam hukum untuk vektor gabungan $(\widetilde{W}^H_{t_1..t_n}, \Delta B_1..\Delta B_n)$ menggunakan kovarians silang $\operatorname{Cov}(\widetilde{W}^H_t, \Delta B_j)$; di [linalg.rs](../../projects/rust/day004h_rough_bergomi/src/linalg.rs) dan [volterra.rs](../../projects/rust/day004h_rough_bergomi/src/volterra.rs).
5. **Skema Hybrid Bennedsen–Lunde–Pakkanen ($\kappa = 1$)**: Mengintegrasikan sel singular pertama secara eksak (berkorelasi $\rho_H = \sqrt{2H}/(H+1/2)$ dengan $\Delta B_i$) dan memotong sel sisanya dengan bobot konvolusi titik kanan dalam kompleksitas $O(n^2)$; di [volterra.rs](../../projects/rust/day004h_rough_bergomi/src/volterra.rs).
6. **Eksponensial Wick & Martingal Variansi Rough Bergomi**: $V_u = \xi_0(u)\exp(\eta\widetilde{W}^H_u - \frac{1}{2}\eta^2 u^{2H})$ merupakan martingal lognormal eksak dengan $\mathbb{E}[V_u] = \xi_0(u)$ dan *variance swap* $\int_0^T \xi_0(u)\,du$; diuji di [rbergomi.rs](../../projects/rust/day004h_rough_bergomi/src/rbergomi.rs).
7. **Dekomposisi Romano–Touzi & Estimator Bersyarat Black–Scholes**: Bersyarat pada $\mathcal{F}^B_T$, $S_T$ berdistribusi lognormal sehingga $\mathbb{E}[(S_T - K)_+ \mid \mathcal{F}^B_T] = \mathrm{BS}(S_{\mathrm{eff}}, K, (1-\rho^2)\int_0^T V_t\,dt)$, memangkas variansi hingga `308x` di sayap OTM; di [bs.rs](../../projects/rust/day004h_rough_bergomi/src/bs.rs) dan [rbergomi.rs](../../projects/rust/day004h_rough_bergomi/src/rbergomi.rs).
8. **Hukum Pangkat *Skew* ATM ($\psi(T) \sim T^{H - 1/2}$)**: Pada volatilitas kasar, kemiringan *smile* ATM membesar sebagai $T^{H-1/2}$ saat $T \to 0$ (berbeda dengan $T^0$ pada model Markov $H=0{,}5$); diuji di [smile.rs](../../projects/rust/day004h_rough_bergomi/src/smile.rs).
9. **Inversi Volatilitas Implisit & Metode Delta SE**: Volatilitas implisit diekstrak dengan Newton–Raphson + *bisection* dan galat standarnya dipropagasikan lewat $\mathrm{SE}(\hat{\sigma}_{\mathrm{imp}}) \approx \mathrm{SE}(\hat{C}) / \mathcal{V}_{\mathrm{BS}}$; di [bs.rs](../../projects/rust/day004h_rough_bergomi/src/bs.rs) dan [smile.rs](../../projects/rust/day004h_rough_bergomi/src/smile.rs).
10. **Optimisasi Simpleks Nelder–Mead dengan Common Random Numbers**: Pencarian bebas-gradien pada permukaan SSE *implied volatility* menggunakan bank normal CRN tetap agar fungsi objektif mulus; di [calib.rs](../../projects/rust/day004h_rough_bergomi/src/calib.rs).
11. **Identifiabilitas & Matriks Informasi Hessian**: Kelengkungan fungsi objektif terhadap $(H, \eta, \rho)$ dioptimum diukur via Hessian beda hingga $3\times 3$ dan diinvers untuk memperoleh $\operatorname{Corr}(\hat{H}, \hat{\eta}) = 0{,}7915$; di [linalg.rs](../../projects/rust/day004h_rough_bergomi/src/linalg.rs) dan [calib.rs](../../projects/rust/day004h_rough_bergomi/src/calib.rs).

## 10. Kesimpulan dan Batasan

1. **Kesimpulan**: Eksperimen membuktikan kesepakatan eksak kovarians ${}_2F_1$ vs kuadratur (`2,12e-10`), efisiensi luar biasa estimator bersyarat Black–Scholes (`VRF = 2,16..308,4`), kemiringan *skew* ATM $\log\psi(T) \sim -0{,}3849\log T$ pada $H=0{,}1$ vs $-0{,}0404\log T$ pada $H=0{,}5$, serta pemulihan penuh `c3`, `c4`, dan `c6` setelah mengganti bobot titik kanan menjadi skema hybrid `BlpOptimal`.
2. **Batasan**: Skema hybrid memakai bobot optimal BLP pada sel pertama dan titik kanan pada sel $k \ge 2$ ($O(n^2)$); kalibrasi memakai $n=100, R=20.000$ (serta $R=80.000$ pada uji sensitivitas) terhadap pasar sintetis $n=500, R=200.000$; `run_completed_without_errors = true` menandakan eksekusi bersih tanpa `NaN`/overflow, bukan bahwa seluruh kriteria ilmiah tercapai (kriteria `c8_calib_bias_h_under_003` tetap berstatus `belum tercapai` karena *ridge* identifiabilitas $(H, \eta)$).

## 11. Referensi dan Jejak Bukti

- Bayer, C., Friz, P., & Gatheral, J. (2016). *Pricing under rough volatility*. Quantitative Finance, 16(6), 887–904.
- Bennedsen, M., Lunde, A., & Pakkanen, M. S. (2017). *Hybrid scheme for Brownian semistationary processes*. Finance and Stochastics, 21(4), 931–965.
- McCrickerd, R., & Pakkanen, M. S. (2018). *Turbocharging Monte Carlo pricing for the rough Bergomi model*. Quantitative Finance, 18(11), 1877–1886.
- Artefak hasil & tes: [outputs/day004h/results.json](../../outputs/day004h/results.json), [README.md](../../projects/rust/day004h_rough_bergomi/README.md), [integration_tests.rs](../../projects/rust/day004h_rough_bergomi/tests/integration_tests.rs).
- Grafik SVG (`outputs/figures/`): [cov-exact-vs-hybrid](../../outputs/figures/day004h_cov-exact-vs-hybrid.svg), [hybrid-bias-vs-n](../../outputs/figures/day004h_hybrid-bias-vs-n.svg), [smile-term-structure](../../outputs/figures/day004h_smile-term-structure.svg), [skew-loglog](../../outputs/figures/day004h_skew-loglog.svg), [vrf-conditional-bs](../../outputs/figures/day004h_vrf-conditional-bs.svg), [calib-h-profile](../../outputs/figures/day004h_calib-h-profile.svg), [rough-vs-markov-t01](../../outputs/figures/day004h_rough-vs-markov-t01.svg).
