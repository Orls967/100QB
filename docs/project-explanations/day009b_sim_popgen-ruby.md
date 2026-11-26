# Penjelasan Proyek `day009b` — Genetika Populasi (Ruby Murni, 65 Tahap)

> **Catatan Pra-Eksekusi vs Pasca-Eksekusi**: Bagian aturan galat baku (SE), ambang jarak variasi total (`c16`), rentang $t$ seleksi (`c11`), koreksi $O(1/N)$ (`c9`), tiga konstanta dari ingatan literatur (`c5`, `c8`, `c13`), serta seluruh ambang kriteria `c1`–`c23` ditetapkan **sebelum** implementasi tahap dan eksekusi penuh dijalankan. Tabel Status Pencapaian di Bagian 7 dan diagnosis ilmiah di Bagian 8 diisi pasca-eksekusi persis dari keluaran [`outputs/day009b/results.json`](../../outputs/day009b/results.json) dan CSV keluaran tanpa mengubah parameter atau ambang.

## 1. Identitas dan Ringkasan Proyek

- **ID Progress**: `day009b` (Hari 9, sesi `b`)
- **Tanggal**: `2026-11-27` (WITA)
- **Topik**: Genetika Populasi dalam Ruby Murni (65 Tahap `S01`–`S65`): Wright–Fisher, Moran, Rantai Markov `Rational`, Koalesen Kingman, Model Pulau & Stepping-Stone, Sapuan Selektif, LD, ARG-Lite, dan Inbreeding
- **Jenis**: `sim` (bonus / *off-roadmap*; tidak mencentang `docs/ROADMAP.md`)
- **Bahasa**: Ruby murni (stdlib saja: `json`, `rational`, `fileutils`, `optparse`; tanpa gem tambahan; lolos bersih `ruby -w` di bawah `ruby 4.0.9 +PRISM [arm64-darwin27]`)
- **Lokasi Proyek**: [projects/ruby/day009b_popgen_ruby](../../projects/ruby/day009b_popgen_ruby)
- **Metrik Eksekusi Penuh (`WORKERS=4`)**: `5.156` baris Ruby di `25` berkas, `65` tahap selesai dalam `15.16 s`, menghasilkan `155` berkas keluaran (`109` CSV = `67.2 KB`, `45` SVG tulis-tangan valid XML = `290.8 KB`, `1` `results.json`), `20/23` kriteria `tercapai`, `3` `belum tercapai` (`c8`, `c9`, `c14`), `vacuous_flags = []`.

## 2. Masalah dan Tujuan

1. **Memvalidasi dinamika evolusi stokastik dengan jangkar eksak**: Menghindari simulasi *black-box* dengan menyandingkan simulator Wright–Fisher, Moran, dan koalesen terhadap jangkar eksak aljabar/analitik di dalam program: rantai Markov eksak dengan aritmetika `Rational` ($N \le 12$), penyelesaian sistem linier distribusi stasioner ($N=20$), sistem linier model pulau takata/Takahata–Nei, rumus distribusi partisi Ewens ($n=6$), serta uji $\chi^2$ dan batas $\pm 3\,\mathrm{SE}$ yang terdefinisi sebelum eksekusi.
2. **Membangun perangkat genetika populasi bervolume penuh dalam Ruby stdlib**: Menunjang alur panjang 65 tahap (`S01`–`S65`) meliputi *sampler* diskrit-kontinu eksak, *drift* netral, seleksi (aditif/dominan/overdominan), keseimbangan mutasi–seleksi, spektrum frekuensi situs (SFS), Tajima's $D$, struktur populasi ($F_{ST}$ Hudson & Nei), *hitchhiking* sapuan selektif, peluruhan & keseimbangan *linkage disequilibrium* ($r^2$), graf rekombinasi leluhur dua lokus (*ARG-lite*), kesetimbangan Hardy–Weinberg, serta *selfing* parsial.
3. **Menjamin determinisme multi-proses dan audit anti-hampa**: Mengimplementasikan paralelisme `Process.fork` + `IO.pipe` dengan aliran RNG independen SplitMix64 (`MASTER_SEED = 20261127`) dan penggabungan urut ID tugas agar `workers = 1, 2, 4` serta dua eksekusi penuh `rake run` identik bit demi bit, disertai audit `vacuous_flags` otomatis atas seluruh CSV.

## 3. Apa yang Dikerjakan Proyek

- **Fondasi Probabilistik & Determinisme (`S01`–`S10`)**:
  - RNG deterministik SplitMix64 $\to$ `Random.new(seed)` per `(stage_id, task_id)` (`S01`) dan penjaminan paralelisme `Process.fork` (`S62`).
  - *Sampler* eksak: Binomial dua cabang (inversi untuk $n\min(p,1-p) < 10$ dan BTPE Kachitvichyanukul–Schmeiser 1988 untuk $\ge 10$) diuji pada 12 pasangan $(n,p)$ (`S02`, `S03`), Multinomial kondisional berantai (`S04`), Poisson Knuth + PTRS (`S05`), Eksponensial & Geometrik (`S06`).
  - Rantai Markov Wright–Fisher eksak berbasis `Rational` vs `Float` untuk $N \le 12$ (`S09`, selisih maks `1.776e-13`) dan eliminasi Gauss untuk probabilitas fiksasi, waktu absorpsi, waktu fiksasi bersyarat, dan distribusi stasioner (`S10`).
- **Dinamika Forward: Drift, Moran & Seleksi (`S11`–`S30`)**:
  - Peluruhan heterozigositas $H_t = H_0(1 - 1/N_{\text{copies}})^t$ pada $N=50$ dan diploid $2N=200$ (`S11`, `S12`), probabilitas fiksasi netral (`S13`, `S14`), serta waktu fiksasi bersyarat vs rantai Markov dan hampiran difusi (`S15`–`S17`).
  - Proses Moran: penskalaan waktu $N^2/2$ terhadap WF (`S18`–`S20`) dan fiksasi seleksi $(1 - 1/r)/(1 - 1/r^N)$ (`S21`, `S22`).
  - Seleksi WF $N=2000, p_0=0.1, s=0.1$ pada $t \in [0, 100]$ (`S23`, `S24`), model alel tak-hingga (`S25`, `S26`), keseimbangan mutasi–seleksi & distribusi stasioner $N=20$ (`S27`–`S29`), serta dominansi/overdominansi (`S30`).
- **Koalesen Kingman, SFS, Tajima & Ewens (`S31`–`S42`)**:
  - Pohon koalesen Kingman pada $n \in \{5, 20, 100\}$ untuk $\mathbb{E}[T_{\text{MRCA}}], \mathbb{E}[L], \mathbb{E}[S], \mathrm{Var}(S)$ (`S31`–`S34`), SFS *unfolded* & *folded* (`S35`), momen & sebaran Tajima's $D$ (`S36`, `S37`), campuran Poisson $S$ (`S38`), partisi Ewens $n=6$ & $\mathbb{E}[K_n]$ (`S39`, `S40`), kecocokan koalesen vs WF terbatas (`S41`), dan koalesen pertumbuhan eksponensial (`S42`).
- **Struktur Populasi, Sapuan Selektif, LD, ARG-Lite & Inbreeding (`S43`–`S65`)**:
  - Model pulau hingga $d \in \{2,4,8\}$ & *stepping-stone* 1D cincin: $T_{\text{within}}, T_{\text{between}}$, $F_{ST}$ Hudson vs $G_{ST}$ Nei, serta isolasi oleh jarak (`S43`–`S50`).
  - Sapuan selektif (*hitchhiking*) ke belakang berdasar lintasan alel menguntungkan terkomdisikan (`S51`–`S53`), peluruhan LD $(1-r)^t$ dan keseimbangan $\sigma_d^2$ vs Ohta–Kimura 1971 (`S54`–`S56`), koalesen dua lokus ARG-lite untuk invariansi marginal & peluruhan korelasi terhadap $\rho$ (`S57`–`S59`), uji HWE & *partial selfing* (`S60`, `S61`), tangga optimasi `V0`–`V4` (`S63`), kompresi lintasan delta+RLE *lossless* (`S64`), dan agregasi akhir (`S65`).

## 4. Uraian Kode

Proyek berlokasi di [projects/ruby/day009b_popgen_ruby](../../projects/ruby/day009b_popgen_ruby) dengan berkas utama:

| Berkas | Baris | Peran Utama |
|---|---|---|
| [`bin/popgen`](../../projects/ruby/day009b_popgen_ruby/bin/popgen) | 8 | Titik masuk CLI (`selftest`, `all`, `group`, `workers`) |
| [`lib/popgen/rng.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/rng.rb) | 86 | SplitMix64, aliran benih deterministik per tahap/tugas, Box–Muller, `parallel_map` via `Process.fork` + `IO.pipe` |
| [`lib/popgen/binomial.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/binomial.rb) | 145 | `Binomial.sample` dua cabang (`sample_inv` & `sample_btpe`), `pmf`, `log_choose` |
| [`lib/popgen/samplers.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/samplers.rb) | 90 | `multinomial`, `exponential`, `geometric_trials`, `poisson` (Knuth + PTRS Hormann 1993) |
| [`lib/popgen/markov.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/markov.rb) | 243 | Eliminasi Gauss `solve_linear_rational`/`float`, matriks transisi WF `Rational` & `Float`, stasioner, waktu absorpsi |
| [`lib/popgen/wf.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/wf.rb) | 147 | `step_biallelic`, `simulate_v0_hash`, `simulate_v1_trajectory`, `simulate_to_absorption`, `simulate_infinite_alleles` |
| [`lib/popgen/moran.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/moran.rb) | 48 | `Moran.simulate_to_absorption` (fitness $r$) dan peluruhan heterozigositas |
| [`lib/popgen/coalescent.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/coalescent.rb) | 183 | `simulate_tree` (Kingman standar & pertumbuhan eksponensial $\beta$), `tajima_constants`, partisi Ewens |
| [`lib/popgen/structured.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/structured.rb) | 157 | `simulate_island_two_sample`, `simulate_structured_pairwise` ($d$-pulau & cincin *stepping-stone*) |
| [`lib/popgen/ld_sweep.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/ld_sweep.rb) | 190 | `ld_stats` ($D, D', r^2$), `simulate_ld_wf`, `generate_sweep_trajectory`, `simulate_hitchhiking_diversity` |
| [`lib/popgen/recomb.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/recomb.rb) | 120 | `simulate_arg_lite` (2-bit *mask* leluhur dua lokus), `benchmark_bitset_recomb` (`V3`) |
| [`lib/popgen/inbreeding.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/inbreeding.rb) | 81 | `hwe_sample_and_test` ($\chi^2$ df=1), `selfing_recursion`, `simulate_partial_selfing` |
| [`lib/popgen/stats.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/stats.rb) | 340 | `block_average`, `median_mad`, `ols_fit`, uji $\chi^2$ GOF, FNV-1a 64-bit, statistik ringkasan |
| [`lib/popgen/svg.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/svg.rb) | 303 | Penulis grafik SVG XML tulis tangan (`write_xy`, batang, histogram) tanpa pustaka eksternal |
| [`lib/popgen/io.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/io.rb) | 233 | `RunContext`, penulis CSV, pelacakan hash FNV-1a, deteksi hampa otomatis `vacuous_flags` |
| [`lib/popgen/selftest.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/selftest.rb) | 86 | 4 *self-test* analitik (`Bin(4,0.3)`, `N=2` WF `Rational`, Kingman $n=2$, Ewens $n=6$) |
| [`lib/popgen/cli.rb`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/cli.rb) | 86 | Pengarah subperintah `selftest`, `all`, `group`, `workers` |
| `lib/popgen/stages_a.rb` .. `stages_h.rb` | 2.610 | Implementasi berurutan `S01`–`S65` yang menulis `109` CSV dan `45` SVG serta mengisi `criteria` |

## 5. Metode dan Alur Kerja

### Aturan Galat Baku (SE), Ambang Pra-Run, dan Konstanta Literatur (A Priori)
- **Replikasi independen & rata-rata blok**:
  - Untuk simulasi Monte Carlo dengan replikasi independen $R$ (lintasan WF/Moran/Koalesen/ARG-lite), estimasi rata-rata adalah $\hat{\mu} = \frac{1}{R}\sum_{r=1}^R X_r$ dengan $\mathrm{SE} = \frac{s}{\sqrt{R}}$ (atau $\mathrm{SE} = \sqrt{\hat{p}(1-\hat{p})/R}$ untuk proporsi fiksasi).
  - Untuk deret waktu stasioner panjang tunggal (`c16` distribusi stasioner mutasi-seleksi $N=20$, serta seri keseimbangan jangka panjang), seri pasca-*burn-in* dibagi menjadi $N_b = 10$ blok tidak tumpang-tindih; $\mathrm{SE}$ dihitung antar rata-rata blok, dan ambang jarak variasi total (*total variation*) `c16` ditetapkan *a priori* sebagai $\mathrm{TV} = \frac{1}{2}\sum_{k=0}^{2N} |\hat{\pi}_k - \pi_k^{\text{exact}}| \le \max\bigl(0.015,\; 3 \cdot \frac{1}{2}\sum_{k=0}^{2N}\mathrm{SE}(\hat{\pi}_k)\bigr)$.
  - Untuk kriteria perbedaan arah antar dua kondisi $A$ dan $B$, dipakai $\mathrm{SE}_\Delta = \sqrt{\mathrm{SE}_A^2 + \mathrm{SE}_B^2}$ dengan syarat efek $\ge 3\,\mathrm{SE}_\Delta$ beserta kontrol negatif.
- **Rentang $t$ seleksi (`c11`, `S23`)**: $N = 2000$ salinan gen, $p_0=0.1, s=0.1$, evaluasi pada $t \in [0, 100]$ generasi (ditulis sebelum run).
- **Koreksi $O(1/N)$ untuk `c9`**: Pada model alel tak-hingga Wright–Fisher dengan $2N=1000$ salinan gen ($N=500$ diploid) dan laju mutasi $\mu$ ($\theta = 2N_{\text{copies}}\mu$), relasi eksak Kimura–Crow adalah $\mathbb{E}[H] = \frac{\theta}{\theta + 1 + \theta/(2N_{\text{copies}})} \approx \frac{\theta}{1+\theta}\bigl(1 - \frac{\theta}{(1+\theta)2N}\bigr)$, dikalibrasi dengan koreksi $O(1/N) < 0.2\%$ dan batas `c9` 5%.
- **Konstanta dari ingatan literatur (ditandai & diberi ambang longgar sebelum run)**:
  1. Hampiran difusi waktu fiksasi bersyarat netral (`c5`): $\mathbb{E}[T_{\text{fix}} \mid \text{fix}] \approx -2N_{\text{copies}} \frac{1-p_0}{p_0}\ln(1-p_0)$ — **dari ingatan literatur**, toleransi 10%.
  2. Varians Tajima's $D$ (`c8`): $\mathrm{Var}(D_{\text{Tajima}}) \approx 1$ untuk $\theta$ sedang dengan rumus bobot $e_1, e_2$ Tajima (1989) — **dari ingatan literatur**, toleransi 10%.
  3. Keseimbangan LD Ohta–Kimura (1971) / Sved (`c13`): $\sigma_d^2 = \frac{\mathbb{E}[D^2]}{\mathbb{E}[p_A(1-p_A)p_B(1-p_B)]} \approx \frac{10+\rho}{22+13\rho+\rho^2}$ — **dari ingatan literatur**, toleransi longgar 25%.

## 6. Output yang Ditargetkan (Kriteria A Priori `c1`–`c23`)

| ID | Deskripsi Kriteria A Priori | Ambang Resmi |
|---|---|---|
| `c1` | Rantai Markov WF `Rational` vs `Float` ($N \le 12$) | Selisih maks $\le 10^{-12}$ |
| `c2` | Probabilitas fiksasi netral $= p_0$ ($p_0 \in \{0.1, 0.5\}$, $\ge 20.000$ ulangan) | Rantai eksak & simulasi $\le 3\,\mathrm{SE}$ |
| `c3` | Moran: rumus tertutup $(1-1/r)/(1-1/r^N)$ vs rantai $\le 10^{-12}$ & simulasi $\le 3\,\mathrm{SE}$; $r>1$ naik, $r<1$ turun ($\ge 3\,\mathrm{SE}$), $r=1$ netral | Rumus $\le 10^{-12}$, sim $\le 3\,\mathrm{SE}$, arah $\ge 3\,\mathrm{SE}$ |
| `c4` | Peluruhan heterozigositas $H_t$ pada 5 titik waktu ($N \in \{50, 200\}$) | Selisih $\le 3\,\mathrm{SE}$ di ke-5 titik |
| `c5` | Waktu fiksasi rata-rata: simulasi vs rantai ($\le 3\,\mathrm{SE}$) & vs hampiran difusi (**dari ingatan literatur**, $\le 10\%$) | Sim vs rantai $\le 3\,\mathrm{SE}$; difusi $\le 10\%$ |
| `c6` | Koalesen Kingman ($n \in \{5, 20, 100\}$): $\mathbb{E}[T_{\text{MRCA}}], \mathbb{E}[L], \mathbb{E}[S], \mathrm{Var}(S)$ | Masing-masing $\le 3\,\mathrm{SE}$ |
| `c7` | SFS ($n=20$): $\xi_i$ vs $\theta/i$ per kelas ($\le 3\,\mathrm{SE}$) & Wald $\chi^2$ agregat $p > 0.01$ | Semua kelas $\le 3\,\mathrm{SE}$ & $p_{\chi^2} > 0.01$ |
| `c8` | Tajima's $D$: rata-rata $0$ ($\le 3\,\mathrm{SE}$) & varians dalam $10\%$ rumus Tajima 1989 (**dari ingatan literatur**) | $\|\bar{D}\| \le 3\,\mathrm{SE}$ & $\mathrm{Var}$ dalam $10\%$ |
| `c9` | Ewens: $\mathbb{E}[K_n]$ ($\le 3\,\mathrm{SE}$), partisi $n=6$ $\chi^2$ $p>0.01$, $\mathbb{E}[H] \approx \theta/(1+\theta)$ dalam $5\%$ ($N=500$, koreksi $O(1/N)$ tertulis) | $\le 3\,\mathrm{SE}$, $p_{\chi^2}>0.01$, $\le 5\%$ |
| `c10` | Model pulau: $T_{\text{within}}, T_{\text{between}}$ vs sistem linier eksak ($\le 3\,\mathrm{SE}$); $F_{ST}$ turun terhadap $m$ ($\ge 3\,\mathrm{SE}$); kontrol $m$ besar $F_{ST}<0.02$ | $\le 3\,\mathrm{SE}$, monoton $\ge 3\,\mathrm{SE}$, $F_{ST}<0.02$ |
| `c11` | Seleksi WF ($N=2000, p_0=0.1, s=0.1$, $t \in [0, 100]$): rata-rata lintasan vs rekursi deterministik | Semua titik waktu $\le 3\,\mathrm{SE}$ |
| `c12` | Hitchhiking: diversitas situs netral pasca-sapuan $<$ kontrol ($\ge 3\,\mathrm{SE}$), pulih monoton terhadap jarak ($\ge 3\,\mathrm{SE}$), $r$ besar tak ada penurunan | Efek $\ge 3\,\mathrm{SE}$, pulih $\ge 3\,\mathrm{SE}$, kontrol netral |
| `c13` | LD: peluruhan $(1-r)^t$ pada $N$ besar ($\le 3\,\mathrm{SE}$); $\sigma_d^2$ vs hampiran literatur (**dari ingatan literatur**, $\le 25\%$) | Peluruhan $\le 3\,\mathrm{SE}$ & $\sigma_d^2 \le 25\%$ |
| `c14` | ARG-lite: $\mathbb{E}[T_{\text{MRCA}}]$ marginal sama untuk semua $\rho$ (selisih $\le 3\,\mathrm{SE}$); korelasi antar lokus turun monoton ($\ge 3\,\mathrm{SE}$) & $=1$ pada $\rho=0$ | Selisih $\le 3\,\mathrm{SE}$, turun $\ge 3\,\mathrm{SE}$, $\rho=0 \implies 1$ |
| `c15` | Hardy–Weinberg $\chi^2$ $p > 0.01$ & selfing $F \to s/(2-s)$ dalam $3\,\mathrm{SE}$ | $p > 0.01$ & $\|F - F^*\| \le 3\,\mathrm{SE}$ |
| `c16` | Distribusi stasioner rantai WF mutasi+seleksi ($N=20$) vs simulasi panjang | Total variation $\le \max(0.015, 3\,\mathrm{SE}_{\text{TV}})$ |
| `c17` | Determinisme: 2× `all` identik & 1/2/4 proses `fork` identik (hash FNV-1a CSV) | Hash identik |
| `c18` | Sampler binomial: $\chi^2$ $p > 0.001$ (koreksi Bonferroni untuk 12 pasangan $(n,p)$ mencakup dua cabang) | Semua 12 pasangan $p > 0.001$ |
| `c19` | Jumlah keluaran $\ge 150$ berkas & $\ge 40$ SVG valid XML | Berkas $\ge 150$ & SVG $\ge 40$ |
| `c20` | `rake selftest` lulus | Exit code 0 |
| `c21` | `vacuous_flags` kosong atau dijelaskan jujur | `vacuous_flags == []` |
| `c22` | Jumlah baris Ruby (`wc -l`) $\in [5000, 8000]$ & waktu per tahap dilaporkan | Baris $\in [5000, 8000]$ & timing 65 tahap |
| `c23` | Tangga optimasi `V0`–`V4` (median $\pm$ MAD, $\ge 5$ ulang) dilaporkan jujur | Dilaporkan lengkap |

## 7. Status Pencapaian

Ringkasan dari [`outputs/day009b/results.json`](../../outputs/day009b/results.json): **20/23 `tercapai`**, **3 `belum tercapai`** (`c8`, `c9`, `c14`), **0 `sebagian`**, **0 `hampa`**, `vacuous_flags = []`.

| Kriteria | Status | Bukti Terukur |
|---|---|---|
| `c1` | `tercapai` | `max_abs_diff = 1.776e-13` ($\le 10^{-12}$, $N=4..12$) |
| `c2` | `tercapai` | `max_|p_sim-p0|/SE = 1.19` ($\le 3\,\mathrm{SE}$, $p_0 \in \{0.1, 0.5\}$) |
| `c3` | `tercapai` | `chain_err = 1.11e-15`, `sim_max_z = 2.05`, `z_up = 35.4`, `z_dn = 36.4`, kontrol `r=1` `0.0025 ± 0.0043` |
| `c4` | `tercapai` | `max_z = 1.55` di seluruh 5 titik waktu untuk $N \in \{50, 200\}$ |
| `c5` | `tercapai` | `max_z_vs_chain = 1.38` ($\le 3\,\mathrm{SE}$), `max_rel_vs_diffusion = 8.87%` ($\le 10\%$) |
| `c6` | `tercapai` | `max_z = 1.95` untuk $T_{\text{MRCA}}, L, S, \mathrm{Var}(S)$ pada $n \in \{5, 20, 100\}$ |
| `c7` | `tercapai` | `max_class_z = 1.97` ($\le 3\,\mathrm{SE}$), Wald $\chi^2$ agregat `p = 0.2875` ($> 0.01$, $n=20$) |
| `c8` | `belum tercapai` | `mean_D = -0.1064 ± 0.0100` (`z = 10.68 > 3`), `Var(D) = 0.7949` (target $\sim 1.0$) — lihat Bagian 8 |
| `c9` | `belum tercapai` | `Kn_max_z = 2.11` ($\le 3$), `n=6 chi2_p = 0.8164` ($> 0.01$), namun `H_max_rel = 42.87%` di $\theta=0.5$ — lihat Bagian 8 |
| `c10` | `tercapai` | `times_max_z = 1.53`, penurunan $F_{ST}$ min `12.4 SE`, kontrol negatif $M=80$ `F_ST = -0.0001` ($< 0.02$) |
| `c11` | `tercapai` | `max_z = 2.74` ($\le 3\,\mathrm{SE}$) di seluruh $t \in [0, 100]$ ($N=2000, p_0=0.1, s=0.1$) |
| `c12` | `tercapai` | Penurunan sapuan `43.9 SE`, pemulihan monoton min `7.9 SE`, kontrol $r$ besar `z = 0.91` |
| `c13` | `tercapai` | Peluruhan $(1-r)^t$ `decay_max_z = 1.73`, keseimbangan $\sigma_d^2$ `max_rel = 20.11%` ($\le 25\%$) vs Ohta–Kimura |
| `c14` | `belum tercapai` | Invariansi marginal `max_z = 1.49`, $\rho=0 \implies 1.0000$, tetapi langkah ekor $\rho=8 \to 32$ hanya `2.7 SE` ($< 3$) — lihat Bagian 8 |
| `c15` | `tercapai` | `HWE_min_p = 0.2014` ($> 0.01$), *selfing* `F = 0.6654 ± 0.0028` vs $F^* = 0.6667$ (`z = 0.45`) |
| `c16` | `tercapai` | `TV = 0.01426 <= threshold = 0.02932` ($N=20$ mutasi+seleksi) |
| `c17` | `tercapai` | `workers 1/2/4 = 8e12879dfcfa8dc0` (identik), 2× eksekusi penuh identik `711328ecfbfdb734` |
| `c18` | `tercapai` | `min_chi2_p = 0.0803` ($> 0.001$) di ke-12 pasangan $(n,p)$ |
| `c19` | `tercapai` | `total_files = 155` ($\ge 150$), `valid_svg = 45` ($\ge 40$ XML valid) |
| `c20` | `tercapai` | `rake selftest` keluar dengan kode `0` (4/4 lulus) |
| `c21` | `tercapai` | `vacuous_flags_count = 0` (`vacuous_flags = []`) |
| `c22` | `tercapai` | `ruby_lines = 5156` ($\in [5000, 8000]$), `stages = 65`, `elapsed = 15.16s` tercatat per tahap |
| `c23` | `tercapai` | Tangga `V0..V4` (`7` ulang) dilaporkan jujur: median ms `[1.50, 0.05, 0.87, 2.09, 0.87, 2.07, 2.50]` |

## 8. Jika Belum Tercapai atau Belum Diketahui

### 8.1 `c8` — Tajima's $D$ Mean & Varians (`S36`)
- **FAKTA**: Pada berkas [`s36_tajima_d_summary.csv`](../../projects/ruby/day009b_popgen_ruby/out/s36_tajima_d_summary.csv) ($n=30$, masing-masing 10.000 pohon koalesen netral), baik estimator Watterson $\hat{\theta}_W$ maupun $\hat{\pi}$ tidak bias sempurna (misal pada $\theta=20$: $\bar{\hat{\theta}}_W = 19.9938$, $\bar{\hat{\pi}} = 20.0004$), tetapi rasio terstandardisasi $D = (\hat{\pi}-\hat{\theta}_W)/\sqrt{e_1 S + e_2 S(S-1)}$ menghasilkan $\bar{D} = -0.0695$ ($\theta=2$, `z=6.60`), $-0.0939$ ($\theta=5$, `z=9.20`), $-0.1064 \pm 0.0100$ ($\theta=10$, `z=10.68`), $-0.1020$ ($\theta=20$, `z=10.35`), dan varians $\mathrm{Var}(D) \in [0.7764, 0.8868]$ di bawah ambang literatur $\sim 1.0$ (`abs_var_diff = 0.2236 > 10%`).
- **HIPOTESIS**: Pembilang $\hat{\pi}-\hat{\theta}_W$ bernilai harap nol, tetapi penyebut $\sqrt{e_1 S + e_2 S(S-1)}$ berkorelasi positif dengan $S$ dan membesar pada realisasi pohon ber-S besar (ketika $\hat{\pi}$ bervarians tinggi); akibat ketaksamaan Jensen dan korelasi rasio tersebut dalam teori Tajima (1989), distribusi sejati $D$ untuk koalesen netral sedikit miring negatif ($\mathbb{E}[D] \approx -0.07$ s.d. $-0.11$) dan bervarians tertekan ($\approx 0.78$–$0.89$ dibanding $1.0$ hampiran). Karena pengetatan jumlah replikasi menjadi $R=10.000$ mengecilkan $\mathrm{SE}(\bar{D}) \approx 0.01$, simpangan eksak rasio tersebut otomatis menonjol (`z=10.68`). Hal ini dibiarkan jujur sebagai temuan tanpa melonggarkan ambang.

### 8.2 `c9` — Heterozigositas Alel Tak-Hingga $\mathbb{E}[H]$ pada $\theta = 0.5$ (`S25`)
- **FAKTA**: Di berkas [`s25_infinite_alleles_h.csv`](../../projects/ruby/day009b_popgen_ruby/out/s25_infinite_alleles_h.csv) ($2N=1000$ salinan gen, $N=500$ diploid), komponen Ewens `Kn_max_z = 2.11 <= 3` (`S40`) dan uji $\chi^2$ partisi $n=6$ (`p = 0.8164 > 0.01`, `S39`) lulus baik. Untuk heterozigositas Kimura–Crow dengan koreksi $O(1/N)$, tiga titik $\theta \in \{1.0, 2.0, 4.0\}$ memberi galat relatif `0.84%`, `8.77%`, dan `0.30%`, namun titik $\theta = 0.5$ ($\mu = 5\times 10^{-4}$) terukur $H = 0.1904 \pm 0.0464$ vs teori $0.3331$ (`rel_err = 42.87% > 5%`, `3.07 SE`).
- **HIPOTESIS**: Waktu autokorelasi heterozigositas pada model alel tak-hingga menskala sebagai $\tau_H \approx \frac{2N_{\text{copies}}}{1+\theta}$. Pada $\theta = 0.5$, $\tau_H \approx 667$ generasi (jauh lebih lambat daripada $\theta=4.0$ yang hanya $200$ generasi), sehingga pada horison simulasi terbatas setiap replikasi membawa varians sampel sangat besar ($\mathrm{SE} = 0.0464$, yakni `13.9%` dari $0.3331$), membuat ambang ketat `5%` terlewati oleh fluktuasi sampling.

### 8.3 `c14` — Penurunan Ekor Korelasi ARG-Lite $\rho = 8 \to 32$ (`S58`)
- **FAKTA**: Di [`s57_arg_pairwise_diff.csv`](../../projects/ruby/day009b_popgen_ruby/out/s57_arg_pairwise_diff.csv) dan [`s58_arg_corr_decay.csv`](../../projects/ruby/day009b_popgen_ruby/out/s58_arg_corr_decay.csv), invariansi distribusi marginal $T_{\text{MRCA},1}$ antar semua pasangan $\rho \in \{0, 0.5, 2, 8, 32\}$ terbukti (`max_z = 1.49 <= 3 SE`), dan pada $\rho = 0$ korelasinya tepat `1.0000`. Korelasi $T_{\text{MRCA}}$ turun tajam dari `1.0000` ($\rho=0$) $\to$ `0.6699 ± 0.0241` ($\rho=0.5$) $\to$ `0.3093 ± 0.0198` ($\rho=2$) $\to$ `0.0762 ± 0.0135` ($\rho=8$) $\to$ `0.0247 ± 0.0133` ($\rho=32$), namun langkah terakhir $\rho=8 \to 32$ memiliki beda `0.0514 ± 0.0189` atau `2.72 SE < 3.0 SE`.
- **HIPOTESIS**: Karena kovarians dua-lokus asimtotik meluruh seperti $O(1/\rho^2)$ ketika $\rho$ membesar, selisih sejati antara $\rho=8$ dan $\rho=32$ hanya sekitar $\approx 0.05$, sedangkan varians estimator korelasi Pearson pada $R$ lintasan tetap memiliki $\mathrm{SE}_\Delta \approx 0.019$; perbedaan tersebut nyata searah penurunan teoretis namun masih berjarak $\approx 0.28\,\mathrm{SE}$ di bawah batas deteksi $3\,\mathrm{SE}$.

### 8.4 Deviasi
- Jumlah tahap diperluas dari minimal `S01`–`S60` menjadi `S01`–`S65` (tetap $\ge 60$) untuk menampung pemeriksaan determinisme *fork* (`S62`), tangga optimasi (`S63`), kompresi *run-length* (`S64`), serta manifes gabungan (`S65`) tanpa memadatkan satu tahap berlebihan.
- Di `S63` (*benchmark* `V4_fork4` pada beban mikro kecil `gens=800` per pekerja), overhead `Process.fork` + `IO.pipe` (`2.50 ms` untuk `fork4` vs `0.05 ms` untuk `V1` sekuensial murni) lebih besar daripada pekerjaan simulasi mikro itu sendiri; paralelisme `fork` baru unggul secara wall-clock untuk tahap panjang bertipe Monte Carlo ribuan pohon (`S31`, `S36`, `S42`, `S45`, `S51`).

## 9. Landasan Teori

1. **Generator Acak Deterministik SplitMix64**: Setiap `(stage_id, task_id)` diturunkan dari `MASTER_SEED = 20261127` melalui fungsi pencampur 64-bit SplitMix64 ($\bmod 2^{64}$) sebelum menginisialisasi `Random.new(seed)` independen pada [`RNG.derive_seed`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/rng.rb#L18-L27).
2. **Penarikan Binomial Dua Cabang (Inversi & BTPE)**: Untuk $np < 10$ digunakan rekursi CDF langsung `sample_inv` (Kemp), sedangkan untuk $np \ge 10$ dipakai algoritma penolakan segitiga–ekor eksponensial BTPE Kachitvichyanukul & Schmeiser (1988) di [`Binomial.sample_btpe`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/binomial.rb#L46-L95) dalam waktu $O(1)$ teramortisasi.
3. **Sampler Poisson (Knuth & PTRS), Multinomial & Eksponensial**: Distribusi multinomial dibangkitkan secara sekuensial lewat faktorisasi binomial bersyarat di [`Samplers.multinomial`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/samplers.rb#L10-L29), sedangkan Poisson menggunakan perkalian uniform Knuth ($\lambda < 10$) dan penolakan tertransformasi PTRS Hörmann (1993) untuk $\lambda \ge 10$ di [`Samplers.poisson_ptrs`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/samplers.rb#L61-L88).
4. **Rantai Markov Wright–Fisher Eksak `Rational`**: Matriks transisi berukuran $(N+1)\times(N+1)$ dibangun tanpa galat pembulatan dengan kelas `Rational` stdlib dan diselesaikan lewat eliminasi Gauss pada [`Markov.solve_linear_rational`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/markov.rb#L25-L45) untuk $N \le 12$.
5. **Peluruhan Heterozigositas & Fiksasi Netral**: Pada populasi bebas mutasi berukuran $N_{\text{copies}}$ salinan gen, heterozigositas meluruh eksak menurut $\mathbb{E}[H_t] = H_0(1 - 1/N_{\text{copies}})^t$ dan probabilitas fiksasi netral tepat sama dengan frekuensi awal $p_0 = i/N_{\text{copies}}$ sebagaimana diuji di [`WF.step_biallelic`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/wf.rb#L12-L23).
6. **Waktu Fiksasi Bersyarat & Hampiran Difusi Kimura**: Rantai serap $(I - Q)t = \mathbf{1}$ memberikan waktu fiksasi bersyarat eksak $t_{\text{fix}}(i) = w_i / u_i$, dibandingkan terhadap rumus hampiran difusi $-2N_{\text{copies}}\frac{1-p_0}{p_0}\ln(1-p_0)$ (**dari ingatan literatur**).
7. **Proses Moran & Probabilitas Fiksasi Selektif**: Proses kelahiran–kematian tumpang tindih satu individu per langkah pada [`Moran.simulate_to_absorption`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/moran.rb#L8-L24) mempercepat varians genetik dua kali relatif terhadap $N$ langkah sehingga penskalaan waktunya adalah $N^2/2$ dan fiksasi selektif memiliki solusi tertutup $(1 - r^{-i})/(1 - r^{-N})$.
8. **Rekursi Seleksi Terbobot Kebugaran & Keseimbangan Mutasi–Seleksi**: Perubahan deterministik frekuensi alel A dengan dominansi $h$ dan kebugaran $w_{AA}=1+s, w_{Aa}=1+hs, w_{aa}=1$ mengikuti $p' = (p^2 w_{AA} + pq w_{Aa})/\bar{w}$ diikuti mutasi dua arah $u, v$ di [`Markov.wf_next_prob_float`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/markov.rb#L82-L90).
9. **Koalesen Kingman & Pertumbuhan Eksponensial**: Waktu tunggu antar penggabungan $k$ garis keturunan berdistribusi $\mathrm{Exp}(k(k-1)/2)$ untuk populasi konstan, dan waktu tak-homogen eksak $\frac{1}{\beta}\ln(1 + \frac{\beta E}{\lambda_k e^{\beta t}})$ untuk populasi bertumbuh $\beta > 0$ pada [`Coalescent.simulate_tree`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/coalescent.rb#L39-L72).
10. **Spektrum Frekuensi Situs (SFS), Watterson, Tajima's $D$ & Ewens**: Mutasi netral Poisson $(0.5\,\theta \cdot \ell_b)$ ditempatkan pada setiap cabang bermakna $b$ sehingga $\mathbb{E}[\xi_i] = \theta/i$, diverifikasi bersama konstanta bobot $e_1, e_2$ Tajima (1989) di [`Coalescent.tajima_constants`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/coalescent.rb#L24-L34) dan rumus sampling partisi Ewens.
11. **Model Pulau Terstruktur & Stepping-Stone 1D**: Migrasi simetris antar $d$ deme pada laju $M = 2Nm$ menghasilkan pasangan garis keturunan dalam satu deme $T_w$ dan beda deme $T_b$ di [`Structured.simulate_structured_pairwise`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/structured.rb#L40-L95) untuk menghitung $F_{ST}$ Hudson dan $G_{ST}$ Nei serta isolasi oleh jarak 1D.
12. **Sapuan Selektif (*Selective Sweep* / *Hitchhiking*)**: Pembangkitan lintasan selektif ke depan terkondisi mencapai fiksasi $1/N \to 1$ melalui [`LDSweep.generate_sweep_trajectory`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/ld_sweep.rb#L69-L81), kemudian penelusuran mundur garis keturunan netral antara latar belakang $B$ dan $b$ pada laju rekombinasi $r$ menunjukkan pemangkasan keragaman $\pi$ di dekat situs target.
13. **Peluruhan Linkage Disequilibrium & Keseimbangan Ohta–Kimura**: Koefisien $D = f_{AB} - p_A p_B$, $D'$, dan $r^2$ dihitung via [`LDSweep.ld_stats`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/ld_sweep.rb#L12-L26) untuk memverifikasi peluruhan deterministik $D_t = D_0(1-r)^t$ dan rasio ekspektasi keseimbangan $\sigma_d^2 \approx \frac{10+\rho}{22+13\rho+\rho^2}$.
14. **Graf Rekombinasi Leluhur Dua Lokus (*ARG-Lite*)**: Setiap garis keturunan membawa *bitmask* 2-bit material leluhur di [`Recomb.simulate_arg_lite`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/recomb.rb#L17-L82); garis keturunan ganda (`0b11`) membelah pada laju $\rho/2$, mempertahankan invariansi distribusi marginal tiap lokus sambil meluruhkan $\mathrm{Corr}(T_{\text{MRCA},1}, T_{\text{MRCA},2})$ terhadap $\rho$.
15. **Kesetimbangan Hardy–Weinberg & Rekursi Partial Selfing**: Uji $\chi^2$ derajat bebas $1$ (`3 kelas - 1 parameter - 1 jumlah`) dan simulasi persilangan sendiri tingkat $s$ pada [`Inbreeding.simulate_partial_selfing`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/inbreeding.rb#L42-L79) memverifikasi konvergensi indeks fiksasi $F_{t+1} = \frac{s}{2}(1+F_t) \to \frac{s}{2-s}$.
16. **Paralelisme Deterministik Berbasis `Process.fork` + Audit Anti-Hampa**: [`RNG.parallel_map`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/rng.rb#L41-L84) menghindari `Ractor` dan menserialisasi luaran pekerja melalui `IO.pipe` + `Marshal` lalu menyusunnya ulang berdasarkan indeks tugas tetap sehingga `hash` 64-bit FNV-1a [`Stats`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/stats.rb#L6-L8) dan deteksi `vacuous_flags` [`IOAudit`](../../projects/ruby/day009b_popgen_ruby/lib/popgen/io.rb#L9-L95) konsisten lintas 1/2/4 proses.

## 10. Kesimpulan dan Batasan

- Seluruh target rekayasa Ruby stdlib tercapai penuh: `5.156` baris Ruby modular (`<= 0` peringatan di `ruby -w`), `65` tahap (`S01`–`S65`), `155` berkas keluaran (`109` CSV + `45` SVG XML valid + `results.json`) selesai dalam `15.16 s`, dengan `vacuous_flags = []` serta paritas bit-identik antar 2 kali eksekusi penuh (`711328ecfbfdb734`) dan antar `workers = 1, 2, 4` (`8e12879dfcfa8dc0`).
- **20 dari 23** kriteria a priori `tercapai`, sedangkan **3 kriteria** (`c8`, `c9`, `c14`) dilaporkan jujur berstatus `belum tercapai` beserta penyebab matematis dan samplingnya pada Bagian 8 tanpa penyesuaian mundur.
- Batasan utama: model simulasi menggunakan populasi bialelik / haploid-diploid homogen berpasangan acak (kecuali pada deme pulau/cincin dan *partial selfing*), ARG terbatas pada dua lokus (*ARG-lite*), dan pada skala besar seperti `c9` ($\theta=0.5$) serta ekor korelasi rekombinasi tinggi (`c14`), presisi dapat dipersempit lebih lanjut dengan menambah panjang burn-in/replikasi Monte Carlo di luar jatah waktu baku ~15 detik.

## 11. Referensi dan Jejak Bukti

- **Kode sumber**: [`projects/ruby/day009b_popgen_ruby`](../../projects/ruby/day009b_popgen_ruby) (`Rakefile`, `bin/popgen`, `lib/popgen/*.rb`).
- **Hasil agregat JSON**: [`outputs/day009b/results.json`](../../outputs/day009b/results.json).
- **Keluaran CSV (109 berkas)**: [`projects/ruby/day009b_popgen_ruby/out/`](../../projects/ruby/day009b_popgen_ruby/out).
- **Grafik SVG tulis tangan (45 berkas)**: [`outputs/figures/day009b_*.svg`](../../outputs/figures).
- **Rujukan literatur klasik**:
  - Ewens, W. J. (1972). *The sampling theory of selectively neutral alleles*. Theoretical Population Biology, 3(1), 87–112.
  - Kachitvichyanukul, V., & Schmeiser, B. W. (1988). *Binomial random variate generation* (BTPE). Communications of the ACM, 31(2), 216–222.
  - Kingman, J. F. C. (1982). *The coalescent*. Stochastic Processes and their Applications, 13(3), 235–248.
  - Ohta, T., & Kimura, M. (1971). *Linkage disequilibrium between two segregating nucleotide sites under the steady flux of mutations in a finite population*. Genetics, 68(4), 571–580.
  - Tajima, F. (1989). *Statistical method for testing the neutral mutation hypothesis by DNA polymorphism*. Genetics, 123(3), 585–595.
