# Penjelasan Proyek 004k — AMM, Impermanent Loss, LVR, Biaya, dan Likuiditas Terkonsentrasi (Rust)

## 1. Identitas Proyek

- **Kode Sesi**: `4k` (`day004k`) — Hari 4, sesi k (22 November 2026, WITA).
- **Bahasa & Lingkungan**: Rust stabil (`rustc 1.93.1`, *edition* `2021`, *zero external crates*), eksekusi paralel deterministik menggunakan `std::thread::scope` (`8` *worker*) dengan reduksi blok tetap (`BLOCK = 64`) dan `MAIN_SEED = 20261122`.
- **Lokasi Proyek**: [projects/rust/day004k_amm_lvr/](../../projects/rust/day004k_amm_lvr/)
- **Keluaran Utama**: [results.json](../../outputs/day004k/results.json) dan 8 grafik SVG di [outputs/figures/](../../outputs/figures/).
- **Catatan Ruang Lingkup**: Seluruh lintasan harga dan aliran pesanan pada proyek ini bersifat **sintetis** (PRNG deterministik `xoshiro256++`) dan tidak memuat klaim empiris mengenai pasar atau protokol nyata, serta bukan nasihat investasi.

---

## 2. Permasalahan yang Dipecahkan

Penyedia likuiditas (*Liquidity Provider*, LP) pada *Automated Market Maker* produk-konstan (*Constant Product Market Maker*, CPMM $x \cdot y = L^2$) menanggung kerugian struktural terhadap pasar acuan yang bergerak kontinu karena harga pool hanya diperbarui ketika pedagang atau arbitraseur bertransaksi:

1. **Impermanent Loss vs Loss-Versus-Rebalancing (LVR)**: *Impermanent loss* ($\mathrm{IL}(r) = 2\sqrt{r}/(1+r) - 1$) hanya membandingkan nilai terminal pool terhadap portofolio *hodl* statis; sedangkan **LVR** (Milionis et al., 2022) memisahkan komponen arah pasar dari biaya konveksitas negatif murni dengan membandingkan pool terhadap strategi replikasi delta kontinu, menghasilkan laju biaya eksak $\frac{\sigma^2}{8}V$.
2. **Blok diskret, biaya pool $\varphi$, dan limit refleksi**: Ketika blok berjarak $\Delta > 0$ dan biaya $\varphi$ dibayar pada input, arbitrase hanya terjadi jika $|\ln(P'/p)| > -\ln(1-\varphi)$. Diperlukan integrasi eksak per blok dan analisis limit gerak Brownian terefleksi ($\Delta \to 0$) untuk memisahkan LVR bersih arbitraseur dari pendapatan biaya LP.
3. **Likuiditas terkonsentrasi & salah spesifikasi**: Rentang $[P_a, P_b]$ meningkatkan efisiensi modal $\mathrm{CE}$ tetapi memperbesar LVR/TVL secara proporsional serta memicu risiko keluar rentang dan biaya gas *rebalancing*. Selain itu, asumsi difusi murni dapat meleset di bawah lompatan harga, volatilitas stokastik, *drift* tak nol, dan ukuran blok besar.

---

## 3. Daftar Fungsi & Modul Utama

1. `rng.rs` — `splitmix64`, `derive_seed`, `Rng::normal`, `Rng::poisson`: PRNG deterministik `xoshiro256++` + Box–Muller + Knuth Poisson.
2. `parallel.rs` — `run_chunks`, `reduce_tuple4_fixed_blocks`: paralelisme `std::thread::scope` berurutan blok tetap (`BLOCK = 64`) agar hasil `1` vs `8` *worker* identik bit-demi-bit.
3. `stats.rs` — `mean`, `sd`, `se`, `sharpe`, `corr`: metrik sampel dan korelasi lintasan.
4. `special.rs` — `norm_pdf`, `norm_cdf` ($\approx 10^{-15}$), `gauss_legendre(n)`: fungsi khusus dan kuadratur Legendre buatan sendiri.
5. `cpmm.rs` — `reserves`, `pool_value`, `pool_delta`, `pool_gamma`, `impermanent_loss`, `expected_pool_value`, `variance_pool_value`, `mc_pool_moments`, `pathwise_lvr_ito`.
6. `blocks.rs` — `block_arb_outcome`, `exact_block_expectations_gl`, `mc_block_annualized`.
7. `reflect.rs` — `reflected_limit_rates`, `mc_reflected_block_sim`, `mc_noise_trader_breakeven`.
8. `conc.rs` — `ConcRange` (`reserves`, `value`, `gamma`, `lvr_rate`, `capital_efficiency`), `expected_exit_time_gbm`, `mc_conc_metrics`.
9. `strategies.rs` — `compare_lp_strategies`: pembanding CRN 5 strategi 30 hari (*Hodl*, *Full-range LP*, *Passive Conc ±10%*, *Rebalanced Conc ±10%*, *Delta-hedged Conc ±10%*).
10. `misspec.rs` — `misspec_merton_jumps`, `misspec_stoch_vol`, `misspec_drift`.
11. `jsonout.rs`, `svg.rs`, `bin/run.rs` — penulis JSON 17 digit, generator 8 grafik SVG, dan runner penuh.

---

## 4. Uraian Kode & Penurunan Matematis

### 4.1 CPMM, Momen GBM, dan LVR Pathwise (`cpmm.rs`)

Dari $x \cdot y = L^2$ dan $p = y/x$, diperoleh $x(p) = L/\sqrt{p}$ dan $y(p) = L\sqrt{p}$, sehingga nilai pool $V(P) = P x(P) + y(P) = 2L\sqrt{P}$ dengan $V_P = L/\sqrt{P}$ dan $V_{PP} = -\frac{1}{2}LP^{-3/2}$. Untuk $P_0 = 1, V_0 = 1 \implies L = 0{,}5$.

Di bawah GBM $\mu = 0$ ($P_T = P_0 \exp(-\frac{1}{2}\sigma^2 T + \sigma W_T)$), $\sqrt{P_T/P_0} = \exp(-\frac{1}{4}\sigma^2 T + \frac{1}{2}\sigma W_T)$, sehingga:
$$\mathbb{E}[V(P_T)] = V_0 e^{-\sigma^2 T / 8}, \qquad \operatorname{Var}[V(P_T)] = \mathbb{E}[V_0^2 P_T/P_0] - \mathbb{E}[V(P_T)]^2 = V_0^2\bigl(1 - e^{-\sigma^2 T/4}\bigr).$$
Dengan lema Itô pada portofolio replikasi delta ($dV_{\mathrm{rebal}} = V_P\,dP$), selisihnya memenuhi $d(V_{\mathrm{rebal}} - V_{\mathrm{pool}}) = -\frac{1}{2}V_{PP}\,(dP)^2 = \frac{\sigma^2}{8}V_t\,dt$.

### 4.2 Mekanik Blok Diskret & Limit Refleksi (`blocks.rs`, `reflect.rs`)

Ketika biaya $\varphi$ dipotong dari input, pembelian $x$ (saat $P' > p$) memasukkan $(1-\varphi)\Delta y$ ke invarian $x' y' = L^2$ dan menggeser harga pool menjadi $p' = y'/x' = p(1 + (1-\varphi)\Delta y / y)^2$. Memaksimalkan laba arbitrase $\Pi(\Delta y) = P'(x - x') - \Delta y$ memberikan syarat orde pertama $(1-\varphi)P'/p' = 1 \iff p' = (1-\varphi)P'$. Karena arbitrase hanya menguntungkan jika $p' > p$, ambang perdagangan adalah $z = \ln(P'/p) > -\ln(1-\varphi) =: \delta(\varphi) \approx \varphi$. Dengan substitusi $u = \sqrt{(1-\varphi)e^z} > 1$, diperoleh masukan $\Delta y_{\mathrm{in}} = \frac{L\sqrt{p}}{1-\varphi}(u-1)$, keluaran $\Delta x_{\mathrm{out}} = \frac{L}{\sqrt{p}}(1 - \frac{1}{u})$, biaya LP $\varphi \Delta y_{\mathrm{in}} = \frac{\varphi}{1-\varphi}L\sqrt{p}(u-1)$, serta laba arbitrase tertutup $\Pi(z) = \frac{L\sqrt{p}}{1-\varphi}(u-1)^2$. Secara simetris, saat $P' < p(1-\varphi)$ ($z < -\delta$), arbitraseur menjual $\Delta x_{\mathrm{in}} = \frac{L}{(1-\varphi)\sqrt{p}}(u-1)$ ($u = \sqrt{(1-\varphi)e^{-z}}$) untuk menggeser harga pool ke $p' = P'/(1-\varphi)$ dengan laba $L\sqrt{p}\frac{(u-1)^2}{u^2}$ dan biaya $\varphi P' \Delta x_{\mathrm{in}}$.

Pada limit $\Delta \to 0$ dengan $\varphi > 0$, proses $X_t = \ln(P_t/p_t)$ menjadi gerak Brownian terefleksi di $[-\delta, +\delta]$. Dengan menerapkan rumus Itô pada $X_t^2$ di keadaan stasioner ($\mathbb{E}[d(X_t^2)] = 0$), laju waktu lokal total pada kedua batas memenuhi $2\delta\,\mathbb{E}[dL_t^+ + dL_t^-]/dt = \sigma^2 \implies \mathbb{E}[dL_t]/dt = \frac{\sigma^2}{2\delta}$. Karena $dy = \frac{1}{4}V\frac{dp}{p}$, laju volume arbitrase adalah $\frac{\sigma^2}{8\delta}V$ dan laju pendapatan biaya LP adalah $\varphi\frac{\sigma^2}{8\delta}V \approx \frac{\sigma^2}{8}V$.

### 4.3 Likuiditas Terkonsentrasi & Salah Spesifikasi (`conc.rs`, `strategies.rs`, `misspec.rs`)

Untuk posisi terkonsentrasi pada $[P_a, P_b]$ dengan likuiditas $L$, cadangan di dalam rentang berbentuk $x(P) = L(P^{-1/2} - P_b^{-1/2})$ dan $y(P) = L(P^{1/2} - P_a^{1/2})$, sehingga nilai pool bernilai piecewise:
$$V(P) = \begin{cases} P L(P_a^{-1/2} - P_b^{-1/2}), & P \le P_a, \\ L(2\sqrt{P} - \sqrt{P_a} - P/\sqrt{P_b}), & P_a < P < P_b, \\ L(\sqrt{P_b} - \sqrt{P_a}), & P \ge P_b. \end{cases}$$
Di dalam rentang, $V_{PP} = -\frac{1}{2}LP^{-3/2}$ tidak bergantung pada batas $[P_a, P_b]$, sehingga $\mathrm{LVR}_{\mathrm{rate}} = \frac{\sigma^2}{4}L\sqrt{P}$. Ketika dinormalisasi terhadap modal awal $V_0 = 1$ pada $P_0 = 1$, rasio $\mathrm{LVR}_0 / V_0 = \frac{\sigma^2}{8} \cdot \mathrm{CE}$ meningkat secara linear terhadap efisiensi modal $\mathrm{CE} = L / 0{,}5$; di luar rentang, posisi berubah menjadi $100\%$ satu aset dengan gamma dan LVR nol. Waktu keluar pertama dari $x_0 = 0$ ke $\{a = \ln P_a, b = \ln P_b\}$ untuk $\ln P_t$ ber-*drift* $\nu = -\sigma^2/2$ memenuhi ODE $\frac{1}{2}\sigma^2 u'' + \nu u' = -1$ dengan syarat batas $u(a)=u(b)=0$, yang diselesaikan secara analitik lewat $u(x) = C_1 + C_2 e^{x} - x/\nu$.

Pada skenario salah spesifikasi, `misspec.rs` menguji empat pelanggaran asumsi difusi standar: (1) proses lompatan-difusi Merton ($\lambda_J = 20/\text{thn}, \mu_J = -0{,}02, \sigma_J = 0{,}05$) dengan kompensator martingal eksak; (2) varians stokastik CIR/Heston ($\theta = \sigma^2 = 0{,}64, \kappa = 3, \xi = 1$) yang mengukur kesesuaian $\mathrm{LVR} \approx \mathbb{E}[v_t V_t/8]$ serta korelasi lintasan $\operatorname{Corr}(\mathrm{LVR}, \bar{v})$; (3) *drift* tak nol $\mu \in \{-0{,}5; +0{,}5\}$ untuk memverifikasi invariansi rasio $\mathrm{LVR}/\int V_t\,dt = \sigma^2/8$ terhadap $\mu$; serta (4) keterlambatan blok $\Delta = 600\text{s}$ terhadap prediksi kontinu.

---

## 5. Metode Eksperimen & Parameter

1. **Bagian 1**: $R = 10^6$ lintasan GBM ($\sigma = 0{,}8, V_0 = 1$) pada $T \in \{1/365, 7/365, 30/365, 1\}$; LVR Itô 30 hari (`1000` langkah, `50.000` lintasan); sapuan $\sigma \in \{0{,}3; 0{,}5; 0{,}8; 1{,}2\}$.
2. **Bagian 2**: Grid $\Delta \in \{12\text{s}, 60\text{s}, 600\text{s}\} \times \varphi \in \{0; 0{,}0005; 0{,}003; 0{,}01\}$ membandingkan kuadratur Gauss–Legendre 64 titik vs MC (`500.000` blok); simulasi multi-langkah refleksi di $\Delta \in \{1\text{s}, 12\text{s}, 60\text{s}\}$ pada $\varphi = 0{,}003$; dan aliran *noise trader* Poisson untuk volume harian/TVL $v \in \{0; 0{,}05; 0{,}1; 0{,}3; 1{,}0\}$.
3. **Bagian 3**: Rentang terkonsentrasi $\pm 5\%$ (`[0,95; 1,05]`), $\pm 10\%$ (`[0,90; 1,10]`), $\pm 30\%$ (`[0,70; 1,30]`), dan rentang penuh (`40.000` lintasan, 14 hari); perbandingan CRN 5 strategi 30 hari (`50.000` lintasan, $\varphi = 0{,}003$, volume harian/TVL `0,15`, gas `0,0005`, *rebalance* tiap 3 hari jika keluar rentang).
4. **Bagian 4**: Lompatan Merton ($\lambda_J = 20/\text{thn}, \mu_J = -0{,}02, \sigma_J = 0{,}05$), volatilitas stokastik CIR ($\theta = 0{,}64, \kappa = 3, \xi = 1$), *drift* $\mu \in \{-0{,}5; +0{,}5\}$, dan ukuran blok $\Delta = 600\text{s}$ (`40.000` lintasan).

---

## 6. Target & Kriteria A Priori

1. `c1_gbm_pool_moments_3se`: $\mathbb{E}[V(P_T)]$ MC ($R = 10^6$) cocok dengan $V_0 e^{-\sigma^2 T/8}$ dalam $3\,\mathrm{SE}$ pada seluruh 4 horizon.
2. `c2_pathwise_lvr_ito_matches`: Selisih LVR replikasi diskret vs integral Itô 30 hari `< 2e-4`.
3. `c3_block_gl_matches_mc_3se`: Integral Gauss–Legendre per blok cocok dengan MC pada seluruh 12 kombinasi $(\Delta, \varphi)$ dalam $3\,\mathrm{SE}$.
4. `c4_reflected_limit_matches_1s_block`: Laju biaya LP simulasi blok $\Delta = 1\text{s}$ ($\varphi = 0{,}003$) mendekati limit Brownian terefleksi dalam `< 25%`.
5. `c5_concentrated_lvr_scales_with_ce`: LVR MC tiap rentang cocok dengan integral Itô dalam $3\,\mathrm{SE}$ dan $\mathrm{LVR}_{5\%} > \mathrm{LVR}_{10\%} > \mathrm{LVR}_{\mathrm{full}}$.
6. `c6_merton_jumps_increase_lvr`: LVR realisasi di bawah lompatan Merton melampaui komponen difusi murni $\sigma^2 V/8$.
7. `c7_stoch_vol_lvr_matches_ev_over_8`: LVR realisasi CIR cocok dengan $\mathbb{E}[v_t V_t/8]$ dalam $3\,\mathrm{SE}$.
8. `c8_drift_invariance_lvr_per_v`: Rasio tahunan $\mathrm{LVR}/\int V_t\,dt$ tetap `0,0800 ± 0,005` untuk $\mu = -0{,}5$ dan $\mu = +0{,}5$.
9. `c9_runtime_under_480s`: Total waktu eksekusi `< 480s`.

---

## 7. Status Pencapaian Berbasis Bukti

Seluruh runner penuh selesai dalam **`86,58s`** pada `8` *worker* (`part1 = 0,98s`, `part2 = 79,93s`, `part3 = 1,85s`, `part4 = 3,81s`), dengan `run_completed_without_errors = true` dan **seluruh 9 kriteria `tercapai`** ([results.json](../../outputs/day004k/results.json)):

| Kriteria | Status | Bukti Terukur | Ambang A Priori |
| --- | --- | --- | --- |
| `c1_gbm_pool_moments_3se` | `tercapai` | `max_z = 2,31` (1d `0,25`, 7d `2,31`, 30d `0,95`, 1y `0,65`) | `<= 3 SE` |
| `c2_pathwise_lvr_ito_matches` | `tercapai` | Rebal `0,0065501` vs Itô `0,0065516` (`diff = 1,50e-6`) | `< 2e-4` |
| `c3_block_gl_matches_mc_3se` | `tercapai` | `max_z = 1,65` pada seluruh 12 pasangan `(Delta, phi)` | `<= 3 SE` |
| `c4_reflected_limit_matches_1s_block` | `tercapai` | Limit `0,07988` vs `1s` MC `0,07761` (`diff = 2,84%`) | `< 25%` |
| `c5_concentrated_lvr_scales_with_ce` | `tercapai` | `+-5%: 3,238`, `+-10%: 1,635`, `+-30%: 0,559`, `Full: 0,080` (`max_z = 1,13`) | `<= 3 SE` & monoton |
| `c6_merton_jumps_increase_lvr` | `tercapai` | Lompatan `0,08678 ± 0,00005` > difusi `0,07973` (`+0,00705`) | `LVR_jump > LVR_diff` |
| `c7_stoch_vol_lvr_matches_ev_over_8` | `tercapai` | `0,07968 ± 0,00008` vs `0,07968` (`z = 0,07`, `corr = 0,923`) | `<= 3 SE` |
| `c8_drift_invariance_lvr_per_v` | `tercapai` | `0,08001` ($\mu=-0{,}5$) & `0,07999` ($\mu=+0{,}5$) | `0,0800 ± 0,005` |
| `c9_runtime_under_480s` | `tercapai` | `86,58s` (`8` *worker*) | `< 480s` |

### Rincian Temuan Kuantitatif

- **Bagian 1 (Bentuk tertutup & $\sigma$-sweep)**: Untuk $r \in \{0{,}25; 0{,}5; 1; 2; 4\}$, nilai pool memenuhi $x\cdot y = L^2$ hingga `1e-14` dan $\mathrm{IL}(r)$ cocok dengan rasio terhadap *hodl*. Pada $\sigma \in \{0{,}3; 0{,}5; 0{,}8; 1{,}2\}$, laju LVR teoritis $\sigma^2/8 = \{0{,}01125; 0{,}03125; 0{,}08000; 0{,}18000\}$ tercapai oleh replikasi diskret hingga galat absolut `< 4e-5` (`0,011245; 0,031248; 0,079961; 0,179994`).
- **Bagian 2 (Blok diskret & limit refleksi)**: Pada grid 12 kombinasi $(\Delta \in \{12\text{s}, 60\text{s}, 600\text{s}\}, \varphi \in \{0; 0{,}0005; 0{,}003; 0{,}01\})$, kuadratur Gauss–Legendre 64 titik cocok sempurna dengan MC (`500.000` blok, `max_z = 1,65`). Pada $\varphi = 0$, LVR bersih sama persis dengan `0,0800` untuk semua $\Delta$ (mis. `12s: 0,08000` GL vs `0,07991 ± 0,00016` MC; `600s: 0,08000` GL vs `0,07994 ± 0,00016` MC), sementara volume arbitrase tahunan turun seiring $\Delta^{-1/2}$ (`258,7 -> 115,7 -> 36,5`). Pada simulasi berkesinambungan $\varphi = 0{,}003$, saat $\Delta$ mengecil `60s -> 12s -> 1s`, LVR bersih arbitraseur menyusut `0,01411 -> 0,00699 -> 0,00215`, sementara pendapatan biaya LP dari arbitraseur naik `0,06582 -> 0,07290 -> 0,07761` mendekati limit teoritis `0,07988`. Jadi **pendapatan biaya dari arbitrase jauh lebih besar daripada LVR bersih** pada $\Delta$ kecil (kebalikan dari blok tunggal i.i.d. tanpa memori harga).
- **Bagian 2.4 & 3 (Titik impas, rentang terkonsentrasi & strategi 30 hari)**: Syarat impas volume/TVL harian pada $\varphi = 0{,}003$ adalah `0,0731/hari` (full range), naik menjadi `0,510` ($\pm 30\%$), `1,493` ($\pm 10\%$), dan `2,957` ($\pm 5\%$). Pada aliran *noise trader* Poisson ($\Delta = 12\text{s}, \varphi = 0{,}003$), laba bersih tahunan LP naik linear dari `0,06583` ($v=0$) menjadi `0,12067` ($v=0{,}05$), `0,17549` ($v=0{,}1$), `0,39430` ($v=0{,}3$), dan `1,15885` ($v=1{,}0$). Pada perbandingan CRN 30 hari (`v = 0,15/hari`, gas `0,0005`, *rebalance* tiap 3 hari), *delta-hedged ±10%* memberi `PnL = 0,07541, SD = 0,03753, Sharpe = 2,009` (selisih berpasangan vs *passive ±10%* `-0,00010 ± 0,00060`, memangkas SD `4,16x`); *passive ±10%* memberi `PnL = 0,07551, SD = 0,15608, Sharpe = 0,484` (unggul `+0,07600 ± 0,00056` atas *hodl* `-0,00049`); *rebalanced ±10%* memberi `PnL = 0,12724, SD = 0,11955, Sharpe = 1,064` dengan rata-rata `2,11` *rebalance* (`gas = 0,00106`, waktu dalam rentang `86,71%` vs `53,13%` pasif); dan *full-range LP* memberi `PnL = 0,00646, SD = 0,11419, Sharpe = 0,057`.
- **Bagian 4 (Salah spesifikasi)**: Lompatan Merton menaikkan LVR tahunan dari `0,07973` menjadi `0,08678 ± 0,00005` (`kontribusi lompatan = +0,00705`); volatilitas stokastik CIR memberi `0,07968 ± 0,00008` vs prediksi $\mathbb{E}[v_t V_t/8] = 0,07968$ (`corr = 0,9228`); *drift* $\mu = -0{,}5$ (`LP ret = -0,02718, IL = -0,00663`) dan $\mu = +0{,}5$ (`LP ret = +0,01332, IL = -0,00681`) mempertahankan rasio `LVR / int V dt = 0,08001` dan `0,07999` (invarian terhadap $\mu$); dan blok `600s` pada $\varphi = 0{,}003$ hanya menghasilkan LVR `0,01631` (jauh lebih kecil dari `0,0800` bila LP salah mengira $\varphi=0$ atau $\Delta \to 0$).

---

## 8. Pemisahan Fakta vs Hipotesis

- **FAKTA**: Seluruh 9 kriteria a priori berstatus `tercapai` pada simulasi deterministik `MAIN_SEED = 20261122` (`86,58s`), dan seluruh uji silang (momen $R=10^6$, Itô *pathwise*, kuadratur GL vs MC, limit refleksi, serta 4 skenario salah spesifikasi) terkonfirmasi secara numerik.
- **HIPOTESIS**: Pada model riil dengan *gas fee* diskret dan *adverse selection* terinformasi (*toxic flow*), selisih antara biaya yang dikutip dari arbitraseur dan LVR bersih pada blok $\Delta = 12\text{s}$ juga dipengaruhi oleh *priority gas auction* / MEV yang tidak dimodelkan pada simulasi sintetis tanpa antrean mempool ini.

---

## 9. Landasan Teori (10 Poin Terkait Kode)

1. **Invarian CPMM ($x\cdot y = L^2$) & Nilai Pool Cekung**: $V(P) = 2L\sqrt{P}$ dengan gamma negatif; di [cpmm.rs](../../projects/rust/day004k_amm_lvr/src/cpmm.rs).
2. **Impermanent Loss & Replikasi Hodl**: $\mathrm{IL}(r) = 2\sqrt{r}/(1+r) - 1$; di [cpmm.rs](../../projects/rust/day004k_amm_lvr/src/cpmm.rs).
3. **Momen Tertutup $\sqrt{P_T}$ di bawah GBM**: $\mathbb{E}[V_T] = V_0 e^{-\sigma^2 T/8}$ dan $\operatorname{Var}[V_T] = V_0^2(1 - e^{-\sigma^2 T/4})$; di [cpmm.rs](../../projects/rust/day004k_amm_lvr/src/cpmm.rs).
4. **Loss-Versus-Rebalancing (Milionis et al., 2022)**: $d(V_{\mathrm{rebal}} - V_{\mathrm{pool}}) = -\frac{1}{2}\sigma^2 P^2 V_{PP}\,dt = \frac{\sigma^2}{8}V_t\,dt$; di [cpmm.rs](../../projects/rust/day004k_amm_lvr/src/cpmm.rs).
5. **Mekanik Biaya pada Input & Zona Tanpa-Arbitrase**: Ambang $-\ln(1-\varphi)$ dan harga marginal pasca-dagang $p' = (1-\varphi)P'$; di [blocks.rs](../../projects/rust/day004k_amm_lvr/src/blocks.rs).
6. **Kuadratur Gauss–Legendre 64-Titik**: Integrasi deterministik ekspektasi laba, biaya, dan volume per blok; di [special.rs](../../projects/rust/day004k_amm_lvr/src/special.rs) dan [blocks.rs](../../projects/rust/day004k_amm_lvr/src/blocks.rs).
7. **Gerak Brownian Terefleksi & Waktu Lokal**: Laju refleksi $\sigma^2/(2\delta)$ pada pita $[-\delta, +\delta]$; di [reflect.rs](../../projects/rust/day004k_amm_lvr/src/reflect.rs).
8. **Likuiditas Terkonsentrasi Uniswap v3**: Kurva rezerve *piecewise* $[P_a, P_b]$, efisiensi modal $\mathrm{CE}$, dan waktu keluar ODE GBM; di [conc.rs](../../projects/rust/day004k_amm_lvr/src/conc.rs).
9. **Proses Lompatan-Difusi Merton & Volatilitas Stokastik CIR**: Pengaruh lompatan diskret dan varians stokastik pada LVR; di [misspec.rs](../../projects/rust/day004k_amm_lvr/src/misspec.rs).
10. **Invariansi LVR terhadap Drift $\mu$**: Melalui replikasi delta dinamis, komponen *drift* terhapus oleh posisi lindung nilai sehingga laju LVR/V tetap $\sigma^2/8$; di [misspec.rs](../../projects/rust/day004k_amm_lvr/src/misspec.rs).

---

## 10. Kesimpulan dan Batasan

1. **Kesimpulan**: Kerangka analitik dan simulasi AMM/LVR dari nol di Rust memverifikasi kesesuaian penuh antara bentuk tertutup, kuadratur Gauss–Legendre, dan Monte Carlo paralel deterministik dalam `86,58s`.
2. **Batasan & Deviasi dari Prompt**: Seluruh data sintetis; tidak ada pelonggaran ambang atau deviasi parameter dari spesifikasi prompt.

---

## 11. Referensi dan Jejak Bukti

- Milionis, J., Moallemi, C. C., Roughgarden, T., & Zhang, A. L. (2022). *Automated market making and loss-versus-rebalancing*. arXiv:2208.06046.
- Milionis, J., Moallemi, C. C., & Roughgarden, T. (2023). *Automated market making and arbitrage profits in the presence of fees*. arXiv:2305.14604.
- Adams, H., Zinsmeister, N., Salem, M., Keefer, R., & Robinson, D. (2021). *Uniswap v3 Core*.
- Artefak hasil & tes: [outputs/day004k/results.json](../../outputs/day004k/results.json), [README.md](../../projects/rust/day004k_amm_lvr/README.md), [integration_tests.rs](../../projects/rust/day004k_amm_lvr/tests/integration_tests.rs).
- Grafik SVG (`outputs/figures/`): [vp-and-il](../../outputs/figures/day004k_vp-and-il.svg), [ev-vs-mc](../../outputs/figures/day004k_ev-vs-mc.svg), [lvr-vs-delta](../../outputs/figures/day004k_lvr-vs-delta.svg), [block-profit-vs-z](../../outputs/figures/day004k_block-profit-vs-z.svg), [reflect-vs-block](../../outputs/figures/day004k_reflect-vs-block.svg), [conc-range-lvr](../../outputs/figures/day004k_conc-range-lvr.svg), [lp-strategies](../../outputs/figures/day004k_lp-strategies.svg), [misspec-lvr](../../outputs/figures/day004k_misspec-lvr.svg).
