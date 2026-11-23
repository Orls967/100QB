# Day 005d — Jump-Diffusion Merton: Seri Eksak, PIDE Crank–Nicolson, Kalibrasi Smile & Risiko Ekor (Rust)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day005d` (Hari 5, sesi `d` pada `2026-11-23`)
- **Topik**: Jump-diffusion Merton (`ROADMAP 075`, topik cadangan): *pricing* PIDE Crank–Nicolson vs seri eksak Merton, kalibrasi *smile*, dan analisis risiko ekor
- **Jenis**: `quant` (topik cadangan sesuai Langkah 0.6 karena topik utama *Marchenko–Pastur & Ledoit–Wolf* telah dikerjakan penuh pada `day002b`)
- **Bahasa**: `rust` (Rust stable, edition 2021, tanpa crate eksternal)
- **Lokasi Proyek**: [`projects/rust/day005d_merton_jump_diffusion/`](../../projects/rust/day005d_merton_jump_diffusion/)

Proyek ini membangun implementasi lengkap model penetapan harga opsi **Merton Jump-Diffusion (1976)** dalam Rust murni tanpa dependensi eksternal. Empat komponen utama dikerjakan dan divalidasi secara independen: (1) rumus seri eksak Merton dengan verifikasi *put-call parity* serta simulasi Monte Carlo paralel deterministik ($R = 400{,}000$ lintasan), (2) pemecah *Partial Integro-Differential Equation* (PIDE) Crank–Nicolson dengan langkah awal Rannacher dan evaluasi integral lompatan *predictor-corrector*, (3) ekstraksi *implied volatility smile* jangka pendek/menengah/panjang beserta kalibrasi parameter $(\sigma, \lambda, \mu_J, \delta_J)$ melalui optimisasi Nelder–Mead 4D, dan (4) kuantifikasi risiko ekor pada opsi *put* OTM dalam serta *Expected Shortfall* 99% posisi *short put* 1 bulan dibandingkan model Black–Scholes dengan volatilitas ATM yang sama.

---

## 2. Masalah dan Tujuan
Model Black–Scholes mengasumsikan lintasan harga aset bergerak kontinu mengikuti gerak Brownian geometris (GBM) dengan volatilitas konstan. Dalam pasar keuangan nyata, harga aset mengalami loncatan diskret (*jumps*) akibat pengumuman makroekonomi, kejutan laba, atau krisis likuiditas. Akibat ketiadaan komponen lompatan, Black–Scholes memiliki dua kelemahan fundamental:
1. **Permukaan Volatilitas Datar**: Black–Scholes gagal menghasilkan *implied volatility skew* dan *smile* yang sangat curam pada jatuh tempo pendek ($T = 0{,}25$ tahun) sebagaimana diamati di pasar opsi ekuitas.
2. **Peremehan Risiko Ekor (*Tail Risk*)**: Pada horizon pendek, difusi kontinu hampir tidak memiliki massa peluang untuk pergerakan ekstrem $-20\%$ hingga $-30\%$, sehingga opsi *put* *out-of-the-money* (OTM) dalam dihargai jauh terlalu murah dan risiko *Expected Shortfall* (ES 99%) posisi jual *put* diremehkan secara drastis.

Tujuan proyek `day005d` adalah menjawab secara kuantitatif:
- Apakah seri eksak Merton (1976) konsisten dengan simulasi Monte Carlo independen ($\mathbb{Q}$-martingal dan harga opsi dalam $3\,\mathrm{SE}$)?
- Apakah pemecah PIDE Crank–Nicolson + Rannacher mampu mereproduksi seri eksak pada tingkat kesalahan $< 0{,}01$ dengan konvergensi spasial-temporal orde kedua ($\sim 16\times$ saat grid diperhalus $4\times$)?
- Seberapa akurat algoritma Nelder–Mead 4D memulihkan parameter $(\sigma, \lambda, \mu_J, \delta_J)$ dari data harga sintetis?
- Seberapa besar selisih harga *put* OTM dalam dan ES 99% antara Merton dan Black–Scholes yang dicocokkan pada ATM IV?

---

## 3. Apa yang Dikerjakan Proyek
Seluruh eksperimen dijalankan secara deterministik menggunakan `MAIN_SEED = 20261123` pada `10` *worker* paralel (`std::thread::scope` dengan reduksi blok tetap ukuran `64`), menghasilkan:
- Berkas hasil numerik berpresisi 17 digit signifikan di [`outputs/day005d/results.json`](../../outputs/day005d/results.json).
- Empat grafik SVG buatan sendiri di `outputs/figures/day005d_*.svg`: perbandingan PIDE vs seri eksak, kurva konvergensi grid PIDE, *implied volatility smile* lintas jatuh tempo, dan perbandingan harga *put* OTM Merton vs Black–Scholes.
- Suite tes lengkap terdiri dari `9` tes unit di [`src/lib.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/lib.rs) (dan `src/rng.rs`) serta `2` tes integrasi di [`tests/integration_tests.rs`](../../projects/rust/day005d_merton_jump_diffusion/tests/integration_tests.rs).
- Biner eksperimen tambahan [`src/bin/run_extras.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/bin/run_extras.rs) untuk uji grid halus `801x400` dan kalibrasi multi-start.

---

## 4. Uraian Kode
Seluruh kode sumber ditempatkan secara modular di [`projects/rust/day005d_merton_jump_diffusion/src/`](../../projects/rust/day005d_merton_jump_diffusion/src/):

1. **[`src/rng.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/rng.rs)**: Implementasi PRNG `xoshiro256++` dengan *seeder* `splitmix64` dan `derive_seed` untuk aliran acak deterministik per indeks lintasan. Menyediakan sampler `uniform`, `normal` (Box–Muller), `poisson` (metode produk Knuth dengan pemecahan sub-interval untuk $\lambda$ besar), dan `gamma` (Marsaglia–Tsang).
2. **[`src/special.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/special.rs) & [`src/stats.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/stats.rs)**: Evaluasi `norm_pdf`, `norm_cdf` (aproksimasi rasional Chebyshev Cody untuk $\operatorname{erfc}$ dengan akurasi $\sim 10^{-15}$), `inv_norm_cdf`, serta statistik sampel `mean`, `sd`, `se`, dan `quantile`.
3. **[`src/parallel.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/parallel.rs)**: Utilitas paralel deterministik `default_workers`, `reduce_vec_fixed_blocks` (penjumlahan berurutan per blok ukuran `64` yang menjamin keluaran bit-identik untuk `1` maupun `N` *worker*), dan `run_chunks`.
4. **[`src/linalg.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/linalg.rs)**: `solve_tridiagonal` yang mengimplementasikan algoritma Thomas $O(N)$ untuk menyelesaikan sistem linear tridiagonal $a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i$ pada setiap langkah waktu PIDE.
5. **[`src/bs.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/bs.rs)**: `bs_call`, `bs_put`, `bs_vega`, dan `implied_vol_call` (inversi Newton–Raphson dengan *fallback* biseksi yang dijamin konvergen dalam batas arbitrase).
6. **[`src/merton.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/merton.rs)**:
   - `MertonParams`: menyimpan $S_0=100, r=0{,}05, q=0{,}01, \sigma=0{,}15, \lambda=0{,}80, \mu_J=-0{,}10, \delta_J=0{,}25$ beserta `kappa()` $= e^{\mu_J + \frac{1}{2}\delta_J^2} - 1$.
   - `merton_call_series` & `merton_put_series`: mengevaluasi seri eksak Merton hingga $n = 80$ (berhenti lebih awal saat bobot Poisson $< 10^{-17}$) dan *put* via *put-call parity*.
   - `simulate_merton_mc`: menyimulasikan nilai terminal $S_T$ secara eksak menggunakan $N_T \sim \mathrm{Poisson}(\lambda T)$ dan jumlah lompatan log $\sum_{i=1}^{N_T} J_i \sim \mathcal{N}(N_T \mu_J, N_T \delta_J^2)$, memverifikasi sifat martingal $e^{-(r-q)T}S_T$ serta harga *call* dan *put* lintas *strike*.
7. **[`src/pide.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/pide.rs)**:
   - `price_call_pide_cn`: menyelesaikan PIDE Merton pada grid log-spot seragam $x \in [\ln S_0 - 2{,}5,\; \ln S_0 + 2{,}5]$. Bobot probabilitas lompatan diskret $w_m = \Phi((y_{m+1/2} - \mu_J)/\delta_J) - \Phi((y_{m-1/2} - \mu_J)/\delta_J)$ dihitung sekali di awal; nilai di luar batas kanan grid ($j \ge N_x$) dievaluasi dengan kelanjutan asimtotik eksak $(e^{x_j - q\tau} - K e^{-r\tau})^+$. Empat *half-step* Rannacher ($\theta = 1$) dipakai di awal untuk meredam diskontinuitas turunan *payoff* di $S = K$, diikuti langkah Crank–Nicolson ($\theta = 0{,}5$) dengan skema *predictor-corrector* pada integral lompatan.
   - `compare_pide_vs_exact` & `pide_grid_convergence`: mengukur galat absolut/relatif PIDE terhadap seri eksak dan rasio konvergensi grid $(101\times 50) \to (201\times 100) \to (401\times 200)$.
8. **[`src/calib.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/calib.rs)**: `build_smile_table` untuk mengekstrak *implied volatility* pada $T \in \{0{,}25; 1{,}0; 2{,}0\}$ dan $K \in [70, 130]$, serta `calibrate_merton_nelder_mead` untuk meminimalkan jumlah kuadrat galat harga (*SSE*) di ruang 4D $(\sigma, \lambda, \mu_J, \delta_J)$.
9. **[`src/tail_risk.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/tail_risk.rs)**: `compare_otm_puts` untuk membandingkan harga *put* OTM ($K \in \{70, 80, 90\}, T = 0{,}25$) terhadap Black–Scholes pada volatilitas ATM, serta `simulate_short_put_tail_risk` untuk menghitung 99% VaR dan *Expected Shortfall* posisi jual *put* 1 bulan ($K=90$, $R = 200{,}000$ lintasan).
10. **[`src/bin/run.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/bin/run.rs) & [`src/bin/run_extras.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/bin/run_extras.rs)**: Orkestrator utama dan eksperimen tambahan.

---

## 5. Metode dan Alur Kerja
Alur kerja eksekusi pada [`src/bin/run.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/bin/run.rs) berjalan berurutan melalui empat tahap:
1. **Bagian 1 (Seri Eksak vs Monte Carlo)**: Menentukan kompensator martingal $\kappa = e^{-0{,}10 + 0{,}5(0{,}25)^2} - 1 = -0{,}066444$, menghitung harga *call* dan *put* seri eksak untuk $K \in \{80, 90, 100, 110, 120\}$ pada $T=1{,}0$, lalu menjalankan `400.000` lintasan MC paralel untuk memeriksa apakah $\mathbb{E}[e^{-(r-q)T}S_T] = S_0$ dan seluruh harga opsi berada dalam $3\,\mathrm{SE}$.
2. **Bagian 2 (PIDE Crank–Nicolson)**: Menyelesaikan PIDE pada grid `401 x 200` untuk setiap *strike* $K \in \{80, 90, 100, 110, 120\}$ serta mengukur konvergensi grid pada $(101\times 50)$, $(201\times 100)$, dan $(401\times 200)$.
3. **Bagian 3 (Smile & Kalibrasi)**: Membangun tabel IV pada $T \in \{0{,}25; 1{,}0; 2{,}0\}$ dan $K \in \{70..130\}$, kemudian menjalankan Nelder–Mead dari tebakan awal terganggu $(\sigma=0{,}20, \lambda=0{,}45, \mu_J=-0{,}04, \delta_J=0{,}18)$ untuk memulihkan parameter sejati.
4. **Bagian 4 (Risiko Ekor)**: Menghitung rasio harga *put* OTM Merton terhadap Black–Scholes ATM IV pada $T=0{,}25$, lalu menyimulasikan distribusi kerugian posisi *short put* 1 bulan ($K=90$, $R=200{,}000$) menggunakan *Common Random Numbers* (CRN) untuk komponen difusi agar perbandingan VaR 99% dan ES 99% tajam dan adil.

---

## 6. Output yang Ditargetkan
1. Pemeriksaan kode: `cargo fmt --check` dan `cargo clippy -- -D warnings` lolos tanpa galat maupun peringatan.
2. Pengujian otomatis: seluruh `9` tes unit dan `2` tes integrasi pada `cargo test --release` lolos.
3. Anggaran waktu: eksekusi penuh `cargo run --release --bin run` selesai di bawah `240s`.
4. Kriteria kuantitatif (`c1`–`c7`) dievaluasi secara objektif dengan ambang yang ditetapkan sebelum eksekusi (*a priori*), menghasilkan `outputs/day005d/results.json` dan 4 grafik SVG di `outputs/figures/`.

---

## 7. Status Pencapaian

Seluruh runner `day005d` selesai dalam **`0,49s`** (`0,488s`) pada `10` *worker* (`part1 = 0,01s`, `part2 = 0,46s`, `part3 = 0,00s`, `part4 = 0,01s`), dengan `run_completed_without_errors = true`:

| Kriteria | Status | Nilai Terukur | Ambang A Priori |
|---|---|---|---|
| `c1_merton_series_vs_mc_3se` | `tercapai` | `martingale_z = 0,49`, `max_call_z = 1,10`, `max_put_z = 0,75` | `martingale_z, max_call_z, max_put_z <= 3,0` (`R=400k`) |
| `c2_pide_cn_matches_exact_and_converges` | `tercapai` | `max_pide_abs_err = 0,00303` (`rel = 0,00066`), `conv_ratio_101_to_401 = 15,43x` | `max_pide_abs_err < 0,01` & `conv_ratio > 8,0` |
| `c3_short_maturity_jump_skew` | `tercapai` | `IV(T=0,25; K=80) = 35,40%` vs `ATM = 22,55%` | `IV(K=80) > IV(K=100) + 2%` |
| `c4_nelder_mead_calibration_recovery` | `tercapai` | `price_rmse = 1,15e-8`, `iv_rmse_bps = 9,17e-6`, `iters = 214` | `price_rmse < 1e-2` & `iv_rmse_bps < 5,0` |
| `c5_tail_risk_otm_put_and_es99` | `tercapai` | `put(K=70)` rasio `191,14x` BS; `ES99 = 28,156` (Merton) vs `4,017` (BS) | `ratio(K=70) > 2,0` & `ES99_merton > ES99_bs` |
| `c6_worker_invariance_1_vs_n` | `tercapai` | `bit_diffs = 0` (`1` vs `10` *worker*) | `bit_diffs == 0` |
| `c7_runtime_under_240s` | `tercapai` | `0,49s` (`10` *worker*) | `< 240,0s` |

### Rincian Hasil Terukur
- **Bagian 1 (Seri vs MC)**: Rata-rata martingal terdiskonto adalah `100,0198 ± 0,0404` (`z = 0,49`). Untuk $K=100, T=1$: harga *call* eksak `10,8036` vs MC `10,8296 ± 0,0235` (`z = 1,10`); harga *put* eksak `6,9211` vs MC `6,9283 ± 0,0234` (`z = 0,31`).
- **Bagian 2 (PIDE Crank–Nicolson)**: Galat absolut maksimum turun dari `0,04676` (`101x50`) menjadi `0,01190` (`201x100`) dan `0,00303` (`401x200`), menunjukkan rasio konvergensi `15,43x` (sangat dekat dengan $4^2 = 16\times$ untuk orde kedua).
- **Bagian 3 (Smile & Kalibrasi)**: Nelder–Mead berkonvergensi dalam `214` iterasi ke $\hat{\sigma} = 0{,}150000$, $\hat{\lambda} = 0{,}800000$, $\hat{\mu}_J = -0{,}100000$, $\hat{\delta}_J = 0{,}250000$ dengan `price_rmse = 1,15e-8`.
- **Bagian 4 (Risiko Ekor)**: Pada $T=0{,}25$ (`BS ATM IV = 22,55%`), harga *put* $K=70$ Merton adalah `0,2732` dibanding BS `0,00143` (`191,14x`); pada $K=80$, `0,6879` vs `0,0707` (`9,72x`). Untuk posisi *short put* 1 bulan ($K=90$, `BS ATM IV = 19,97%`), VaR 99% Merton adalah `19,426` vs BS `2,330`, dan ES 99% Merton adalah `28,156` vs BS `4,017`.

---

## 8. Pemisahan Fakta vs Hipotesis
- **FAKTA**:
  1. Seluruh `7` kriteria (`c1`–`c7`), `9` tes unit, dan `2` tes integrasi lolos setelah suntingan kode terakhir.
  2. Pada percobaan awal `test_pide_cn_matches_bs_at_zero_lambda`, domain spasial lebar `[-4,5; +4,5]` dengan `n_x = 301` (`dx = 0,03`) menghasilkan galat `0,0279` di titik *strike* $K = S_0$ akibat diskontinuitas gradien *payoff*; penyempitan domain ke `[-2,5; +2,5]` (yang tetap mencakup $> 8\sigma$) dan grid `401 x 200` (`dx = 0,0125`) menurunkan galat maksimum lintas seluruh *strike* menjadi `0,00303` dengan konvergensi orde kedua `15,43x`.
  3. Kalibrasi Nelder–Mead 4D pada permukaan harga sintetis bebas-derau memulihkan seluruh empat parameter hingga 7–8 digit signifikan.
- **HIPOTESIS**:
  1. Pada data pasar riil yang mengandung *bid-ask noise* dan keterbatasan likuiditas di sayap ekstrem, pemulihan simultan $(\lambda, \mu_J, \delta_J)$ dapat mengalami korelasi parameter yang lebih tinggi daripada kasus sintetis bebas-derau ini; verifikasi lebih lanjut dapat dilakukan dengan menambahkan *noise* pengamatan pada `run_extras`.

---

## 9. Landasan Teori (10 Poin Terkait Kode)
1. **Proses Jump-Diffusion Merton (1976)**: Dinamika harga aset menggabungkan difusi Brownian kontinu dan proses Poisson majemuk dengan ukuran lompatan log-normal $J \sim \mathcal{N}(\mu_J, \delta_J^2)$; dimodelkan dalam `MertonParams` pada [`src/merton.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/merton.rs).
2. **Kompensator Martingal Netral-Risiko ($\kappa$)**: Agar proses harga terdiskonto $e^{-(r-q)t}S_t$ menjadi $\mathbb{Q}$-martingal, laju *drift* harus dikompensasi sebesar $-\lambda\kappa$ dengan $\kappa = \mathbb{E}[e^J - 1] = e^{\mu_J + \frac{1}{2}\delta_J^2} - 1$; di `MertonParams::kappa` dan `simulate_merton_mc` (`martingale_z = 0,49`).
3. **Ekspansi Seri Eksak Campuran Poisson**: Bersyarat pada $N_T = n$, distribusi $\ln S_T$ tetap Gaussian dengan variansi $\sigma^2 T + n\delta_J^2$ dan rata-rata tergeser $n(\mu_J + \frac{1}{2}\delta_J^2) - \lambda\kappa T$, sehingga harga opsi adalah penjumlahan tertimbang Poisson dari rumus Black–Scholes; di `merton_call_series`.
4. **Identitas Put-Call Parity pada Jump-Diffusion**: Karena *put-call parity* $C - P = S_0 e^{-qT} - K e^{-rT}$ hanya bergantung pada sifat martingal dan bebas arbitrase tanpa bergantung pada distribusi, ia berlaku eksak pada model Merton; di `merton_put_series` dan `test_merton_put_call_parity`.
5. **Partial Integro-Differential Equation (PIDE)**: Penerapan lema Itô untuk proses lompatan-difusi menghasilkan persamaan diferensial parsial orde dua dengan suku integral konvolusi non-lokal $\lambda\int V(\tau, x+y)\phi(y)\,dy$; di [`src/pide.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/pide.rs).
6. **Skema Crank–Nicolson & Langkah Awal Rannacher**: Metode Crank–Nicolson ($\theta = 0{,}5$) memberikan akurasi orde kedua dalam waktu dan ruang, sementara 4 *half-step* Euler implisit ($\theta = 1$) di $\tau = 0$ meredam osilasi akibat patahan *payoff* di $S = K$; di `price_call_pide_cn` (`conv_ratio = 15,43x`).
7. **Algoritma Thomas ($O(N)$) & Skema Predictor-Corrector**: Matriks tridiagonal dari operator difusi-adveksi diselesaikan dalam waktu linear via [`solve_tridiagonal`](../../projects/rust/day005d_merton_jump_diffusion/src/linalg.rs), sedangkan suku integral non-lokal dievaluasi secara eksplisit-implisit dua tahap agar struktur tridiagonal tetap terjaga.
8. **Implied Volatility Skew & Teorema Limit Pusat (CLT)**: Lompatan negatif ($\mu_J = -0{,}10$) menghasilkan distribusi return miring ke kiri (*negative skewness*) dan leptokurtik pada $T$ pendek ($T=0{,}25$: `35,40%` di $K=80$ vs `22,55%` ATM), yang perlahan mendatar pada $T=2{,}0$ sesuai CLT; di [`src/calib.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/calib.rs).
9. **Optimisasi Simpleks Nelder–Mead 4D**: Pencarian bebas-turunan melalui refleksi, ekspansi, kontraksi, dan penyusutan simpleks meminimumkan galat kuadrat harga opsi terhadap parameter Merton; di `calibrate_merton_nelder_mead`.
10. **Value-at-Risk (VaR 99%) & Expected Shortfall (ES 99%)**: Pengukuran risiko ekor pada posisi *short put* 1 bulan menunjukkan bahwa pencocokan volatilitas ATM Black–Scholes gagal menangkap kerugian ekstrem akibat lompatan mendadak (`ES99 = 28,156` vs `4,017`); di [`src/tail_risk.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/tail_risk.rs).
11. **Common Random Numbers (CRN) & Determinisme Paralel**: Pemakaian aliran acak `xoshiro256++` terindeks per lintasan dan reduksi blok tetap ukuran `64` menjamin hasil bit-identik lintas jumlah *worker* (`bit_diffs = 0`); di [`src/parallel.rs`](../../projects/rust/day005d_merton_jump_diffusion/src/parallel.rs).

---

## 10. Kesimpulan, Batasan, dan Deviasi dari Prompt
1. **Kesimpulan**: Proyek `day005d` berhasil mengimplementasikan dan memverifikasi model Merton Jump-Diffusion dari empat sudut pandang (seri eksak, Monte Carlo, PIDE Crank–Nicolson, kalibrasi *smile*, dan risiko ekor) dalam `0,49s` dengan seluruh `7` kriteria berstatus `tercapai`, `9` tes unit, dan `2` tes integrasi lolos.
2. **Batasan**: Seluruh data harga dan *smile* bersifat sintetis; kalibrasi diuji pada permukaan bebas-derau; suku bunga, *dividend yield*, dan intensitas lompatan diasumsikan konstan sepanjang horizon.
3. **Deviasi dari Prompt**:
   - **Pemakaian Topik Cadangan (`ROADMAP 075`)**: Topik utama (*kovarians sampel, Marchenko–Pastur, dan penyusutan Ledoit–Wolf*) telah dikerjakan penuh pada sesi [`day002b`](day002b_quant_estimasi-kovarians-portofolio-varians.md) (Python), sehingga sesuai instruksi Langkah 0.6 proyek ini beralih ke topik cadangan **Jump-Diffusion Merton (`ROADMAP 075`)** dalam Rust murni tanpa mencentang `docs/ROADMAP.md`.

---

## 11. Referensi dan Jejak Bukti
- Merton, R. C. (1976). *Option Pricing When Underlying Stock Returns Are Discontinuous*. Journal of Financial Economics, 3(1–2), 125–144.
- Cont, R., & Tankov, P. (2004). *Financial Modelling with Jump Processes*. Chapman & Hall/CRC.
- Tavella, D., & Randall, C. (2000). *Pricing Financial Instruments: The Finite Difference Method*. Wiley.
- Artefak lokal: [`projects/rust/day005d_merton_jump_diffusion/`](../../projects/rust/day005d_merton_jump_diffusion/), [`outputs/day005d/results.json`](../../outputs/day005d/results.json), dan 4 grafik `outputs/figures/day005d_*.svg`.
