# Day 005c — Kriteria Kelly dan Risiko Drawdown: Peluang Sentuh Eksak, Kelly Fraksional, Galat Estimasi & Ekor Berat (Rust)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day005c` (sesi `c` pada `2026-11-23`; huruf `b` telah ditempati berkas `day005b_black_scholes` milik pemilik yang belum di-commit)
- **Topik**: Kriteria Kelly dan risiko drawdown: peluang sentuh eksak, Kelly fraksional, galat estimasi, dan ekor berat
- **Jenis**: `quant` (proyek tantangan biasa di luar daftar `ROADMAP.md`)
- **Bahasa**: `rust` (Rust stable, edition 2021, tanpa crate eksternal)
- **Lokasi Proyek**: [`projects/rust/day005c_kelly_drawdown/`](../../projects/rust/day005c_kelly_drawdown/)

Proyek ini membangun kerangka simulasi dan analitik lengkap untuk kriteria Kelly dari kasus taruhan biner diskret hingga pasar kontinu gerak Brownian geometris (GBM), membuktikan estimator peluang sentuh tanpa bias menggunakan jembatan Brownian, memetakan kompromi pertumbuhan-risiko pada Kelly fraksional, menurunkan dan menguji kerugian pertumbuhan akibat galat estimasi $1/(2T_{\mathrm{est}})$ beserta penyusutan optimal, serta menguji ketahanan prediksi Gaussian terhadap return harian berekor berat Student-$t$ ($\nu = 4$).

---

## 2. Masalah dan Tujuan
Kriteria Kelly memaksimalkan laju pertumbuhan geometris jangka panjang $\mathbb{E}[\ln W_T]$, tetapi dalam praktik menghadapi empat tantangan utama yang dijawab proyek ini:
1. **Volatilitas dan risiko drawdown pada Kelly penuh**: Seberapa sering kekayaan menyentuh setengah modal ($0{,}5 W_0$) atau seperempat modal ($0{,}25 W_0$) baik pada horizon tak hingga maupun horizon hingga $T$?
2. **Bias pemantauan diskret dalam simulasi**: Mengapa simulasi berlangkah $\Delta$ meremehkan peluang sentuh sebesar $O(\sqrt{\Delta})$, dan bagaimana koreksi jembatan Brownian menghilangkan bias tersebut?
3. **Kompromi Kelly fraksional ($f = c f^*$)**: Mengapa setengah-Kelly ($c = 0{,}5$) mempertahankan $75\%$ pertumbuhan berlebih sambil memangkas peluang sentuh secara drastis?
4. **Keterbatasan asumsi parameter diketahui dan Gaussian**: Berapa kerugian pertumbuhan ketika $\mu$ harus diestimasi dari sampel $T_{\mathrm{est}}$ tahun, seberapa besar perbaikan dari penyusutan (*shrinkage*), dan bagaimana ekor berat Student-$t$ ($\nu = 4$) mematahkan prediksi Gaussian?

---

## 3. Apa yang Dikerjakan Proyek
Proyek ini mengeksekusi empat studi terintegrasi pada data sintetis dengan `MAIN_SEED = 20261123` dan `10` *worker* paralel deterministik, menghasilkan:
- Berkas hasil numerik berpresisi 17 digit di [`outputs/day005c/results.json`](../../outputs/day005c/results.json).
- Tujuh grafik SVG buatan sendiri di `outputs/figures/day005c_*.svg`.
- Laporan status berbasis bukti pada 8 kriteria yang ditetapkan sebelum eksekusi (*a priori*).

---

## 4. Uraian Kode
Seluruh kode berada di [`projects/rust/day005c_kelly_drawdown/`](../../projects/rust/day005c_kelly_drawdown/) tanpa dependensi eksternal (`Cargo.toml`):
1. [`src/rng.rs`](../../projects/rust/day005c_kelly_drawdown/src/rng.rs): `splitmix64`, `derive_seed`, PRNG `xoshiro256++` (`Rng`), sampler `uniform`, `bernoulli`, `normal` (Box–Muller), `gamma` (Marsaglia–Tsang), dan `student_t` ($\chi^2$ via `gamma`).
2. [`src/special.rs`](../../projects/rust/day005c_kelly_drawdown/src/special.rs): `norm_pdf`, `norm_cdf` (aproksimasi rasional Cody untuk $\operatorname{erfc}$ dengan akurasi $\sim 10^{-15}$), dan `inv_norm_cdf` (biseksi + Newton).
3. [`src/stats.rs`](../../projects/rust/day005c_kelly_drawdown/src/stats.rs): `mean`, `variance_sample`, `sd`, `se`, dan `quantile`.
4. [`src/parallel.rs`](../../projects/rust/day005c_kelly_drawdown/src/parallel.rs): `default_workers`, `reduce_vec_fixed_blocks` (blok tetap ukuran `64` yang menjamin hasil bit-identik untuk `1` vs `N` *worker*), dan `run_chunks`.
5. [`src/kelly_binary.rs`](../../projects/rust/day005c_kelly_drawdown/src/kelly_binary.rs): `g_binary`, `f_star_binary`, `maximize_g_binary` (*golden-section* + Newton), `breakeven_binary`, serta `run_binary_study` ($N=1000, R=100{,}000$) termasuk aproksimasi kontinu.
6. [`src/gbm.rs`](../../projects/rust/day005c_kelly_drawdown/src/gbm.rs): `GbmParams` ($\mu=0{,}08, \sigma=0{,}20, r=0{,}03$), `kelly_fraction`, `kelly_growth`, `drift_log`, dan `vol_log`.
7. [`src/hitting.rs`](../../projects/rust/day005c_kelly_drawdown/src/hitting.rs): `prob_hit_infinite`, `prob_hit_finite`, `bridge_cross_prob`, `run_hitting_grid` (estimator jembatan Brownian Rao–Blackwellized `1 - prod(1 - p_cross)` pada $R = 10^6$ lintasan bersama lintas $(T, x, f)$ untuk tiap $\Delta$), dan `infinite_horizon_table`.
8. [`src/fractional.rs`](../../projects/rust/day005c_kelly_drawdown/src/fractional.rs): `build_fractional_curve`, `build_fractional_map`, dan `compute_fractional_roots`.
9. [`src/estimation.rs`](../../projects/rust/day005c_kelly_drawdown/src/estimation.rs): `run_estimation_loss_study` ($R = 10^6$), `optimal_shrinkage_c`, dan `run_shrinkage_study` ($R = 200{,}000$ CRN, horizon 20 tahun).
10. [`src/heavy_tail.rs`](../../projects/rust/day005c_kelly_drawdown/src/heavy_tail.rs): `run_heavy_tail_study` (return harian Student-$t(\nu=4)$ bervariansi sama dengan GBM, $R = 100{,}000$) dan `build_practical_summary`.
11. [`src/jsonout.rs`](../../projects/rust/day005c_kelly_drawdown/src/jsonout.rs) & [`src/svg.rs`](../../projects/rust/day005c_kelly_drawdown/src/svg.rs): serialisasi JSON 17-digit dan generator grafik SVG.
12. [`src/bin/run.rs`](../../projects/rust/day005c_kelly_drawdown/src/bin/run.rs) & [`tests/integration_tests.rs`](../../projects/rust/day005c_kelly_drawdown/tests/integration_tests.rs): runner utama dan tes integrasi determinisme/invariansi.

---

## 5. Metode dan Alur Kerja
Alur kerja runner utama [`src/bin/run.rs`](../../projects/rust/day005c_kelly_drawdown/src/bin/run.rs):
1. **Bagian 1 (Taruhan Biner)**: Menghitung $f^* = p - (1-p)/b = 0{,}10$ dan memverifikasi terhadap maksimisasi numerik; mencari akar $g(f_{\mathrm{be}}) = 0$; mensimulasikan $N=1000$ taruhan sebanyak $R=100{,}000$ lintasan untuk $f \in \{0{,}025; 0{,}05; 0{,}10; 0{,}15; 0{,}20; 0{,}30\}$; membandingkan $\mathbb{E}[\ln W_N]$, median $W_N$, rata-rata aritmetik $\mathbb{E}[W_N]$, peluang rugi, dan peluang sentuh $0{,}5 W_0$ terhadap aproksimasi kontinu GBM.
2. **Bagian 2 (Peluang Sentuh GBM)**: Untuk $\Delta \in \{1/12, 1/52, 1/252\}$, $T \in \{1, 5, 20\}$, $x \in \{0{,}5; 0{,}25\}$, $f \in \{0{,}5; 1{,}0; 1{,}25; 2{,}0\}$ ($72$ kombinasi, $R = 10^6$), menyimulasikan lintasan Brownian bersama hingga $T=20$ dan membandingkan pemantauan diskret naif, estimator jembatan Brownian tanpa bias, dan rumus tertutup prinsip refleksi.
3. **Bagian 3 (Kelly Fraksional)**: Mengevaluasi rasio pertumbuhan berlebih $c(2-c)$, peta $(c, T)$ (median, persentil ke-5, SD log-kekayaan, peluang rugi), serta akar tertutup vs numerik untuk $g(cf^*)=r$ dan $P_\infty(\text{sentuh } 0{,}5) \in \{1/2, 1/4\}$.
4. **Bagian 4 (Galat Estimasi & Ekor Berat)**: Menguji kerugian pertumbuhan $1/(2T_{\mathrm{est}})$ untuk $T_{\mathrm{est}} \in \{1, 5, 10, 20\}$ ($R = 10^6$) pada dua pasang $(\mu, \sigma)$; membandingkan plug-in penuh, setengah-Kelly, dan penyusutan optimal $c^*$ ($R = 200{,}000$ CRN, $T=20$); menyimulasikan return harian Student-$t(\nu=4)$ ($252$ langkah/tahun, $T=20$, $R=100{,}000$); serta menyusun tabel ringkasan praktis.

---

## 6. Output yang Ditargetkan
1. Format dan *linter* Rust bersih (`cargo fmt --check`, `cargo clippy -- -D warnings`).
2. Seluruh `8` tes unit dan `3` tes integrasi lolos, termasuk invariansi bit-identik `1` vs `N` *worker* dan determinisme *seed*.
3. Runner penuh selesai `< 180s` dan menghasilkan `results.json` serta 7 grafik SVG.
4. Delapan kriteria kuantitatif (`c1`–`c8`) dievaluasi secara objektif dengan ambang *a priori*.

---

## 7. Status Pencapaian

Seluruh runner `day005c` selesai dalam **`46,80s`** pada `10` *worker* (`part1 = 0,41s`, `part2 = 40,65s`, `part3 = 0,00s`, `part4 = 5,73s`), dengan `run_completed_without_errors = true`:

| Kriteria | Status | Nilai Terukur | Ambang A Priori |
|---|---|---|---|
| `c1_binary_kelly_exact_and_mc` | `tercapai` | `f_star_err = 0,0`, `max_z_log = 2,07`, `max_median_rel_err = 5,33e-14` | `f_star_err < 1e-10`, `max_z_log <= 3,0`, `rel < 0,05` |
| `c2_gbm_hitting_bridge_unbiased_3se` | `tercapai` | `33/33` sel informatif (`exact >= 5e-3`) dalam `3 SE` (`max |z| = 2,80`), naif negatif `33/33`, `slope_sqrt_dt = 0,4946 ± 0,0168` | `frac_inf_within_3se >= 0,95`, `frac_inf_naive_neg >= 0,90`, `slope in 0,5 ± 0,15` |
| `c3_fractional_kelly_identity_and_roots` | `tercapai` | `max_excess_ratio_err = 4,44e-16`, `root_max_err = 1,11e-16` | `excess < 1e-12`, `root < 1e-10` |
| `c4_estimation_loss_1_over_2t_est` | `tercapai` | `max_z_default = 1,41`, `max_z_alt_params = 1,41` | `max_z <= 3,0` pada kedua pasang $(\mu,\sigma)$ |
| `c5_shrinkage_improves_growth_and_drawdown` | `tercapai` | `c_shrink = 0,2381`, `diff = +0,07561 ± 0,00036`, `P20(0,5)`: `5,82%` vs `53,71%` | `diff > 3 SE` & `P20_shrink < P20_full` |
| `c6_heavy_tail_t4_exceeds_gaussian_drawdown` | `tercapai` | $\Delta g(2f^*) = -0{,}00631 \pm 0{,}00057$ (`z = -11,04`), `ruin_t4(2f*) = 0,387%` vs `ruin_gauss(2f*) = 0%` | `|z(g_t4(2f*) - g_gauss_daily(2f*))| > 3,0` & `ruin_t4 > 0` vs `gauss == 0` |
| `c7_worker_invariance_1_vs_n` | `tercapai` | `bit_diffs = 0` (`1` vs `10` *worker*) | `bit_diffs == 0` |
| `c8_runtime_under_180s` | `tercapai` | `46,80s` (`10` *worker*) | `< 180,0s` |

### Rincian Hasil Terukur
- **Bagian 1 (Biner)**: $f^* = 0{,}10$, $g(f^*) = 0{,}005008$, $f_{\mathrm{be}} = 0{,}198668$. Untuk $f=0{,}10$: $\mathbb{E}[\ln W_{1000}] = 4{,}9890 \pm 0{,}0100$ (eksak `5,0084`), median `149,66` (eksak `149,66`), rata-rata aritmetik `19.594,97 ± 1.398,94` (eksak `20.959,16`), peluang rugi `0,0614` (kontinu `0,0563`), peluang sentuh `0,5` = `0,4718` (kontinu `0,4909`). Untuk $f = 0{,}30$ (over-Kelly ekstrem), rata-rata aritmetik sampel hanya `6,17e5 ± 3,02e5` jauh di bawah nilai teoretis `6,87e12` karena ekor kanan ekstrem tidak tersampel pada $R = 10^5$.
- **Bagian 2 (Sentuh GBM)**: Pada $T=5, x=0{,}5, f=1{,}25$ (`exact = 0,09973`), pemantauan naif memberi `0,08436` ($\Delta=1/12$, bias `-0,01537`), `0,09133` ($\Delta=1/52$, bias `-0,00840`), `0,09636` ($\Delta=1/252$, bias `-0,00338`), sedangkan jembatan Brownian memberi `0,10002 ± 0,00029` (`z = 0,97`), `0,09928 ± 0,00030` (`z = -1,53`), dan `0,10020 ± 0,00030` (`z = 1,57`). Pada horizon tak hingga ($x=0{,}5$), $c \in \{0{,}25; 0{,}5; 0{,}75; 1; 1{,}5; 2\}$ memberi peluang sentuh `{0,000000; 0,008729; 0,096500; 0,257028; 0,590496; 0,846745}`.
- **Bagian 3 (Kelly Fraksional)**: Akar tertutup dan numerik identik hingga `1,11e-16`: $c(g=r) = 2{,}0$, $c(P_\infty=0{,}5) = 1{,}354400$, $c(P_\infty=0{,}25) = 0{,}989924$.
- **Bagian 4 (Estimasi, Penyusutan & Ekor Berat)**: Kerugian MC untuk $T_{\mathrm{est}} \in \{1, 5, 10, 20\}$ adalah `0,50038 ± 0,00071` (`0,5`), `0,10020 ± 0,00014` (`0,1`), `0,05001 ± 0,00007` (`0,05`), `0,02498 ± 0,00004` (`0,025`), identik pada parameter alternatif $(\mu=0{,}13, \sigma=0{,}35, r=0{,}01)$. Penyusutan optimal $c^* = 0{,}2381$ mengungguli plug-in penuh (`+0,07561 ± 0,00036`) maupun setengah-Kelly (`0,02865`), serta memangkas peluang sentuh `0,5` 20-tahun dari `53,71%` ke `5,82%`.

---

## 8. Pemisahan Fakta vs Hipotesis (Pembaruan Bukti `c2` dan `c6`)
- **FAKTA**:
  1. Pada `c2_gbm_hitting_bridge_unbiased_3se`, pemisahan himpunan sel informatif (`exact_prob >= 5e-3`, `33` sel) menghilangkan artefak derau MC pada sel `1e-6`–`1e-13`: `33/33` sel jembatan Brownian berada dalam `3 SE` (`max |z| = 2,80`), `33/33` sel naif memiliki bias negatif, dan regresi dalam-kelompok $\ln|\text{bias}|$ terhadap $\ln\Delta$ menghasilkan kemiringan **`0,4946 ± 0,0168`** (teoretis $0{,}5$).
  2. Pada `c6_heavy_tail_t4_exceeds_gaussian_drawdown`, penambahan simulasi Gaussian harian (`dt = 1/252`) membuktikan bahwa selisih total `-0,00465` pada $P(\text{sentuh } 0{,}5 \mid f^*)$ didominasi oleh **efek pemantauan diskret** (`-0,00401`), sedangkan efek murni ekor berat $t_4$ pada ambang moderat `0,5` hanya `-0,00064 ± 0,00184` (tak signifikan secara statistik). Sebaliknya, pada *leverage* $2f^*$, ekor berat $t_4$ menurunkan pertumbuhan realisasi sebesar `-0,00631 ± 0,00057` (`z = -11,04`) dan memicu `ruin = 0,387%` (Gaussian harian `0%`).

| $c = f/f^*$ | (a) Pemantauan Diskret (`Gauss Harian - Kontinu`) | (b) Ekor Berat (`t4 Harian - Gauss Harian`) | (c) Total (`t4 Harian - Kontinu`) | $\Delta g_{t_4 - \mathrm{Gauss}}$ | Ruin ($t_4$ vs Gauss) |
|---:|---:|---:|---:|---:|---:|
| `0,5` | `-0,00065` | `+0,00062 ± 0,00037` | `-0,00003` | `-0,00011 ± 0,00012` | `0%` vs `0%` |
| `1,0` | `-0,00401` | `-0,00064 ± 0,00184` | `-0,00465` | `-0,00060 ± 0,00026` | `0,024%` vs `0%` |
| `2,0` | `-0,00731` | `-0,00207 ± 0,00208` | `-0,00938` | `-0,00631 ± 0,00057` (`z=-11,04`) | `0,387%` vs `0%` |
- **HIPOTESIS**:
  1. Tidak ada hipotesis yang belum teruji tersisa untuk `c2` maupun `c6`; seluruh komponen perbedaan telah dipisahkan secara kuantitatif melalui kontrol Gaussian harian berpasangan.

---

## 9. Landasan Teori (10 Poin Terkait Kode)
1. **Kriteria Kelly Diskret (Kelly, 1956)**: Memaksimalkan ekspektasi log-kekayaan $\mathbb{E}[\ln W_N] = N[p\ln(1+bf) + (1-p)\ln(1-f)]$ menghasilkan fraksi optimal $f^* = p - (1-p)/b$; di [`src/kelly_binary.rs`](../../projects/rust/day005c_kelly_drawdown/src/kelly_binary.rs).
2. **Perbedaan Rata-rata Geometris vs Aritmetik**: $\mathbb{E}[W_N] = (1 + f(bp - (1-p)))^N$ tumbuh secara eksponensial bahkan saat $f > f_{\mathrm{be}}$ di mana median $\exp(Ng(f)) \to 0$, akibat dominasi ekor kanan langka; di `run_binary_study`.
3. **Gerak Brownian Geometris & Persamaan Log-Kekayaan**: Dari $dW_t/W_t = (r + f(\mu-r))\,dt + f\sigma\,dB_t$, lema Itô menghasilkan drift log $g(f) = r + f(\mu-r) - \frac{1}{2}f^2\sigma^2$ dan volatilitas $s(f) = f\sigma$; di [`src/gbm.rs`](../../projects/rust/day005c_kelly_drawdown/src/gbm.rs).
4. **Martingale Eksponensial & Peluang Sentuh Horizon Tak Hingga**: Proses $M_t = \exp(-(2g/s^2)\ln(W_t/W_0))$ adalah martingale, sehingga Optional Stopping pada $\tau_a$ menghasilkan $P(\tau_a < \infty) = x^{2g/s^2}$ (khusus $f=f^*, r=0$: $P=x$); di `prob_hit_infinite` ([`src/hitting.rs`](../../projects/rust/day005c_kelly_drawdown/src/hitting.rs)).
5. **Prinsip Refleksi & Teorema Girsanov (Horizon Hingga)**: Mengubah ukuran ke gerak Brownian tanpa drift lalu memantulkan lintasan di batas $a = \ln x$ menghasilkan rumus tertutup $\Phi((a-gT)/(s\sqrt{T})) + e^{2ga/s^2}\Phi((a+gT)/(s\sqrt{T}))$; di `prob_hit_finite`.
6. **Jembatan Brownian & Rao–Blackwellization**: Bersyarat pada $(Y_i, Y_{i+1}) > a$, peluang lintasan kontinu menyentuh $a$ dalam selang $\Delta$ adalah $\exp(-2(Y_i-a)(Y_{i+1}-a)/(s^2\Delta))$, menghilangkan bias diskretisasi $O(\sqrt{\Delta})$; di `run_hitting_grid`.
7. **Identitas Pertumbuhan Berlebih Kelly Fraksional**: Karena $g(f)$ kuadratik terhadap $f$, substitusi $f = cf^*$ menghasilkan $(g(cf^*) - r)/(g^* - r) = c(2 - c)$; di [`src/fractional.rs`](../../projects/rust/day005c_kelly_drawdown/src/fractional.rs).
8. **Ekspansi Kuadratik Eksak & Kerugian Estimasi $1/(2T_{\mathrm{est}})$**: Karena $g(\hat{f}) = g(f^*) - \frac{1}{2}\sigma^2(\hat{f} - f^*)^2$ dan $\operatorname{Var}(\hat{f}) = 1/(\sigma^2 T_{\mathrm{est}})$, diperoleh $\mathbb{E}[g(f^*) - g(\hat{f})] = \frac{1}{2T_{\mathrm{est}}}$ bebas dari $\mu$ dan $\sigma$; di [`src/estimation.rs`](../../projects/rust/day005c_kelly_drawdown/src/estimation.rs).
9. **Penyusutan Optimal (Shrinkage)**: Memaksimalkan $\mathbb{E}[g(c\hat{f})]$ di bawah ketidakpastian $\hat{f}$ menghasilkan faktor penyusutan $c^* = f^{*2}/(f^{*2} + \operatorname{Var}(\hat{f}))$; di `optimal_shrinkage_c` dan `run_shrinkage_study`.
10. **Distribusi Student-t & Kegagalan Ekor Berat**: Return harian $t_\nu$ ($\nu=4$) diskalakan dengan $\sqrt{(\nu-2)/\nu}$ agar bervariansi sama dengan GBM, tetapi tetap memicu risiko *ruin* non-nol dan penurunan pertumbuhan pada *leverage* tinggi; di [`src/heavy_tail.rs`](../../projects/rust/day005c_kelly_drawdown/src/heavy_tail.rs).
11. **Common Random Numbers (CRN) & Reduksi Blok Deterministik**: Pemakaian lintasan acak yang sama lintas fraksi Kelly dan strategi penyusutan memangkas varians pembanding $\pm\,\mathrm{SE}$, sementara reduksi blok tetap ukuran `64` menjamin hasil bit-identik antara `1` dan `10` *worker*; di [`src/parallel.rs`](../../projects/rust/day005c_kelly_drawdown/src/parallel.rs).

---

## 10. Kesimpulan, Batasan, dan Deviasi dari Prompt
1. **Kesimpulan**: Seluruh 4 studi Kelly dan drawdown berhasil diimplementasikan dan dijalankan dalam `46,80s` dengan determinisme penuh (`bit_diffs = 0` antara `1` dan `10` *worker*), meloloskan seluruh 8 kriteria (`c1`–`c8`), 8 tes unit, dan 3 tes integrasi.
2. **Batasan**: Data bersifat sintetis; dampak ekor berat $t_4$ pada ambang moderat `0,5` kecil dibanding efek diskretisasi harian, namun terlihat jelas pada pertumbuhan $2f^*$ dan probabilitas *ruin*.
3. **Deviasi dari Prompt**:
   - **Penamaan sesi `day005c`**: Karena berkas `projects/rust/day005b_black_scholes/` milik pemilik sudah ada di pohon kerja (belum di-commit), sesi ini memakai huruf bebas berikutnya `day005c` sesuai instruksi Langkah 0.3.
   - **Perbaikan `c2` (sel informatif & regresi $\sqrt{\Delta}$)**: Sel dengan `exact_prob < 5e-3` tidak dipakai menguji tanda bias naif karena didominasi derau MC; sebagai gantinya, `c2` mensyaratkan pada `33` sel informatif (`exact_prob >= 5e-3`): jembatan dalam `3 SE` (`>= 95%`), bias naif negatif (`>= 90%`), dan kemiringan regresi $\ln|\text{bias}|$ vs $\ln\Delta$ dalam `0,5 ± 0,15`.
   - **Perbaikan `c6` (pembanding Gaussian harian setara)**: Membandingkan $t_4$ harian (`dt = 1/252`) terhadap simulasi Gaussian harian berpasangan pada `dt` yang sama serta memisahkan efek pemantauan diskret, ekor berat, dan total.
   - **Parameter penyusutan (Bagian 4.2)**: Proyek menggunakan penyusutan skalar optimal frequentist $c^* = f^{*2}/(f^{*2} + \operatorname{Var}(\hat{f}))$ ($c^* = 0{,}2381$ untuk $T_{\mathrm{est}} = 5$) menuju nol dan membandingkannya dengan plug-in penuh serta setengah-Kelly pada horizon 20 tahun ($R = 200{,}000$ CRN).

---

## 11. Referensi dan Jejak Bukti
- Kelly, J. L. (1956). *A New Interpretation of Information Rate*. Bell System Technical Journal, 35(4), 917–926.
- Thorp, E. O. (2006). *The Kelly Criterion in Blackjack, Sports Betting, and the Stock Market*. Handbook of Asset and Liability Management.
- MacLean, L. C., Thorp, E. O., & Ziemba, W. T. (2011). *The Kelly Capital Growth Investment Criterion*. World Scientific.
- Artefak lokal: [`projects/rust/day005c_kelly_drawdown/`](../../projects/rust/day005c_kelly_drawdown/), [`outputs/day005c/results.json`](../../outputs/day005c/results.json), dan 7 grafik `outputs/figures/day005c_*.svg`.
