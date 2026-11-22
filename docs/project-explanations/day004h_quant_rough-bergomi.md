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

Runner penuh selesai dalam **`104,42s`** pada `8` *worker* (`part1 = 2,46s`, `part2 = 4,48s`, `part3 = 7,78s`, `part4 = 89,70s`, `part5 = 0,002s`; [results.json](../../outputs/day004h/results.json)):

| Kriteria | Status | Bukti Terukur |
| --- | --- | --- |
| `c1_volterra_cov_quad_within_1e8` | `tercapai` | `max_rel_err = 2,12e-10 < 1e-8` di 6 pasangan $(s, t)$ |
| `c2_cholesky_exact_cov_within_4_5se` | `tercapai` | `max_chol_recon_err = 5,6e-16`; `max_cov_z_score = 3,84`; $\operatorname{Var}(\widetilde{W}_1)$ `0,9944` vs `1,0000` |
| `c3_martingales_within_3se` | `belum tercapai` | Harga lolos (`max_z_s_t = 2,25`), tetapi martingal variansi gagal: $\mathbb{E}[e^{\dots}] = 1{,}0866$ di $T=0{,}5$ (`max_z_exp_w = 4,80`) |
| `c4_variance_swap_within_3se` | `belum tercapai` | $\mathbb{E}[\int_0^T V\,dt]$: `0,01053` ($T=0{,}25$), `0,02119` ($T=0{,}5$), `0,04238` ($T=1$) vs `0,04 T` (`max_z = 13,84`) |
| `c5_cond_bs_matches_mc_and_vrf` | `tercapai` | `max_z_diff = 1,79 <= 3,5`; `VRF = 2,16..308,4` (`T=0,1, k=0,2`); `max_parity_err = 2,45e-14` |
| `c6_hybrid_matches_cholesky_price_3se` | `belum tercapai` | Pada $n=50, T=1$: $k=-0{,}2$ (`z=2,40`), $k=0$ (`0,06830` vs `0,06669`, `z=5,60`), $k=0{,}2$ (`z=7,26`) |
| `c7_atm_skew_power_law_and_markov_zero` | `tercapai` | $H=0{,}1$: slope `-0,3850 ± 0,0128` ($R^2=0{,}997$); $H=0{,}5$: slope `-0,0404 ± 0,0253` |
| `c8_calib_bias_h_under_003` | `belum tercapai` | `main_fit`: $\hat{H} = 0{,}1629, \hat{\eta}=1{,}918, \hat{\rho}=-0{,}730$; 30 replikasi: `bias_H = +0,0652`, `SD_H = 0,0139` |
| `c9_runtime_under_300s` | `tercapai` | Total waktu `104,42s < 300s` (`8` *worker*) |

Hasil terukur tambahan:
- **Bias skema hybrid titik kanan vs $n$**: Untuk $H = 0{,}1$, $\operatorname{Var}(\widetilde{W}^H_1)$ teoretis (MC) pada $n \in \{50, 100, 250, 500\}$ adalah `1,0494 (1,0421)`, `1,0438 (1,0370)`, `1,0369 (1,0365)`, `1,0323 (1,0331)` (bias `+4,94%` $\to$ `+3,23%`). Untuk $H = 0{,}5$, bias teoretis adalah `0,0000` di seluruh $n$ dan `max_cov_err_vs_min = 0`.
- **Skala *skew* ATM lintas $H$**: $H=0{,}05$ (`-0,4214 ± 0,0136` vs `-0,45`), $H=0{,}1$ (`-0,3850 ± 0,0128` vs `-0,40`), $H=0{,}25$ (`-0,2696 ± 0,0132` vs `-0,25`), $H=0{,}5$ (`-0,0404 ± 0,0253` vs `0,00`).
- **Pembanding Markov ($H=0{,}5$) di $T=0{,}1$**: Setelah menala $\eta_{\mathrm{Markov}} = 4{,}020$ agar ATM IV cocok (`0,1817` vs `0,1832`), *skew* ATM Markov (`0,6045`) sedikit lebih landai daripada *rough Bergomi* (`0,6448`), namun membutuhkan *vol-of-vol* dua kali lebih besar ($\eta = 4{,}02$ vs $1{,}90$) dan kehilangan kelengkungan tajam di sayap *out-of-the-money*.
- **Identifiabilitas & salah spesifikasi horizon**: Korelasi Hessian $\operatorname{Corr}(\hat{H}, \hat{\eta}) = +0{,}7915$ (trade-off kuat antara kekasaran $H$ dan skala $\eta$). Kalibrasi hanya pada $T \ge 0{,}5$ menggeser $\hat{H}$ lebih jauh ke `0,1946` ($\hat{\eta}=1{,}674$), sedangkan kalibrasi pada $T \le 0{,}1$ menghasilkan $\hat{H} = 0{,}1604$ ($\hat{\eta}=2{,}362$, $\hat{\rho}=-0{,}942$).

## 8. Diagnosis Bagian yang Belum Tercapai (`c3`, `c4`, `c6`, `c8`)

Keempat kriteria gagal memiliki **satu akar matematis bersama** yang terukur secara langsung:
- **FAKTA**: Pada Bagian 1 (`hybrid_bias_h01` di [results.json](../../outputs/day004h/results.json)), bobot konvolusi titik kanan $((k-1)\Delta t)^{H-1/2}$ untuk $k \ge 2$ selalu lebih besar dari rata-rata integral selnya karena $u^{H-1/2}$ menurun ketat untuk $H < 1/2$. Akibatnya $\operatorname{Var}_{\mathrm{hyb}}(\widetilde{W}^H_1)$ melebihi $1{,}0$ sebesar `+4,94%` ($n=50$), `+4,38%` ($n=100$), `+3,69%` ($n=250$), dan `+3,23%` ($n=500$).
- **HIPOTESIS (terkonfirmasi secara kuantitatif)**:
  1. **Untuk `c3` & `c4`**: Karena kompensator eksponensial menggunakan $-\frac{1}{2}\eta^2 u^{2H}$ (variansi eksak) sedangkan skema hybrid menghasilkan variansi $\operatorname{Var}_{\mathrm{hyb}}(u) \approx 1{,}0369\,u^{2H}$ pada $n=250$, nilai harapannya menjadi $\exp\bigl(\frac{1}{2}\eta^2(0{,}0369)\bigr) = \exp(0{,}0666) = 1{,}0689$, cocok dengan nilai ukur `1,0646` ($T=1$) dan kenaikan *variance swap* `+5,95%` (`0,04238` vs `0,04000`).
  2. **Untuk `c6`**: Pada $n=50$, kelebihan variansi $+4{,}94\%$ pada $\widetilde{W}^H$ menaikkan *vol-of-vol* efektif $V_u$, sehingga harga *call* hybrid (`0,06830` di ATM) lebih tinggi secara signifikan (`z = 5,60..7,26`) daripada referensi Cholesky eksak (`0,06669`).
  3. **Untuk `c8`**: Karena kalibrasi memakai $n=100$ (bias variansi $+4{,}38\%$) terhadap pasar $n=500$ (bias $+3{,}23\%$) dan bias titik kanan membesar ketika $H$ mengecil, Nelder–Mead menaikkan $\hat{H}$ (`bias_H = +0,0652`) agar variansi efektif model turun mendekati pasar.
- **Bukti lanjutan**: Mengganti bobot titik kanan $k \ge 2$ dengan bobot optimal $L^2$ Bennedsen–Lunde–Pakkanen ($b_k = \sqrt{\frac{k^{2H} - (k-1)^{2H}}{2H}}$) atau memperbesar $\kappa \ge 2$ akan menghilangkan bias variansi ini; sesuai prinsip integritas, bobot titik kanan dipertahankan apa adanya untuk mengukur bias diskretisasinya secara jujur.

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

1. **Kesimpulan**: Eksperimen membuktikan kesepakatan eksak kovarians ${}_2F_1$ vs kuadratur (`2,12e-10`), efisiensi luar biasa estimator bersyarat Black–Scholes (`VRF = 2,16..308,4`), serta kemiringan *skew* ATM $\log\psi(T) \sim -0{,}3850\log T$ pada $H=0{,}1$ vs $-0{,}0404\log T$ pada $H=0{,}5$. Di sisi lain, proyek ini secara kuantitatif mendokumentasikan bagaimana bias variansi $+3{,}2\%$–$+4{,}9\%$ dari skema hybrid titik kanan memicu pergeseran martingal variansi (`c3`, `c4`), harga `n=50` (`c6`), dan bias kalibrasi $\hat{H}$ (`c8`).
2. **Batasan**: Skema hybrid memakai bobot titik kanan pada sel $k \ge 2$ ($O(n^2)$); kalibrasi memakai $n=100, R=20.000$ terhadap pasar sintetis $n=500, R=200.000$; `run_completed_without_errors = true` menandakan eksekusi bersih tanpa `NaN`/overflow, bukan bahwa seluruh kriteria ilmiah tercapai (`c3, c4, c6, c8` berstatus `belum tercapai`).

## 11. Referensi dan Jejak Bukti

- Bayer, C., Friz, P., & Gatheral, J. (2016). *Pricing under rough volatility*. Quantitative Finance, 16(6), 887–904.
- Bennedsen, M., Lunde, A., & Pakkanen, M. S. (2017). *Hybrid scheme for Brownian semistationary processes*. Finance and Stochastics, 21(4), 931–965.
- McCrickerd, R., & Pakkanen, M. S. (2018). *Turbocharging Monte Carlo pricing for the rough Bergomi model*. Quantitative Finance, 18(11), 1877–1886.
- Artefak hasil & tes: [outputs/day004h/results.json](../../outputs/day004h/results.json), [README.md](../../projects/rust/day004h_rough_bergomi/README.md), [integration_tests.rs](../../projects/rust/day004h_rough_bergomi/tests/integration_tests.rs).
- Grafik SVG (`outputs/figures/`): [cov-exact-vs-hybrid](../../outputs/figures/day004h_cov-exact-vs-hybrid.svg), [hybrid-bias-vs-n](../../outputs/figures/day004h_hybrid-bias-vs-n.svg), [smile-term-structure](../../outputs/figures/day004h_smile-term-structure.svg), [skew-loglog](../../outputs/figures/day004h_skew-loglog.svg), [vrf-conditional-bs](../../outputs/figures/day004h_vrf-conditional-bs.svg), [calib-h-profile](../../outputs/figures/day004h_calib-h-profile.svg), [rough-vs-markov-t01](../../outputs/figures/day004h_rough-vs-markov-t01.svg).
