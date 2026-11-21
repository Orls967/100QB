# day004b — AFNS: Bias Sampel Terbatas, Rentang Kalender (*Long-Span*), dan Frekuensi Observasi (*Infill*)

Modul Go 1.22+ mandiri (`zero external dependencies`) untuk menguji hipotesis bias Kendall / *finite-sample* pada estimasi *Maximum Likelihood* (MLE) ruang-keadaan filter Kalman model *Arbitrage-Free Nelson–Siegel* (AFNS) 11 parameter, serta memisahkan efek **panjang rentang kalender ($H$)** dari **kerapatan frekuensi observasi ($\Delta$)** pada rentang kalender tetap.

## 1. Intuisi Masalah

Pada `day004a`, studi replikasi $R=50$ dengan $T = 1{,}000$ observasi mingguan ($H \approx 19{,}23$ tahun) menemukan bahwa:
- Parameter peluruhan lintasan Q ($\lambda$), volatilitas faktor ($\sigma_1, \sigma_2, \sigma_3$), dan derau observasi ($\sigma_\varepsilon$) diidentifikasi sangat presisi dari penampang *cross-sectional* kurva imbal hasil pada 10 jatuh tempo.
- Sebaliknya, kecepatan *mean reversion* di bawah ukuran fisik $\mathbb{P}$ ($K^{\mathbb{P}}_{11}, K^{\mathbb{P}}_{22}, K^{\mathbb{P}}_{33}$) mengalami bias ke atas sebesar $+0{,}21\text{ hingga }+0{,}27\text{ thn}^{-1}$, dan cakupan interval Wald 95% bagi $\theta^{\mathbb{P}}_1$ ($0{,}045$, dengan *level* paling persisten $K^{\mathbb{P}}_{11} = 0{,}05\text{ thn}^{-1}$, waktu paruh $\approx 13{,}9$ tahun) merosot ke $68\%$.

Pertanyaan utama `day004b`: apakah bias ke atas $K^{\mathbb{P}}_{ii}$ dan cakupan rendah $\theta^{\mathbb{P}}_1$ dikendalikan oleh **rentang kalender $H$** (sebagaimana diprediksi teori asimtotik *long-span* difusi Ornstein–Uhlenbeck / bias Kendall–Tang–Chen) atau dapat dihilangkan oleh **pengambilan sampel lebih rapat** (*infill* bulanan $\to$ mingguan $\to$ harian pada $H=20$ tahun tetap)?

## 2. Spesifikasi Model dan Satuan

Seluruh eksperimen mempertahankan spesifikasi AFNS `day004a`:
- **Tingkat suku bunga & rata-rata faktor $\theta^{\mathbb{P}}$**: desimal per tahun ($\theta^{\mathbb{P}} = (0{,}045, -0{,}020, 0{,}010)$ sesuai $4{,}5\%$, $-2{,}0\%$, $1{,}0\%$).
- **Parameter peluruhan Q & mean-reversion P**: $\lambda = 0{,}6\text{ thn}^{-1}$, $K^{\mathbb{P}} = \operatorname{diag}(0{,}05,\; 0{,}5,\; 0{,}4)\text{ thn}^{-1}$.
- **Matriks difusi faktor**: $\Sigma = \operatorname{diag}(0{,}0060,\; 0{,}0100,\; 0{,}0200)\text{ thn}^{-1/2}$.
- **Derau observasi**: $\sigma_\varepsilon = 0{,}0005$ ($5\text{ bp}$), berdistribusi $\mathcal{N}(0, \sigma_\varepsilon^2)$ independen antar-tenor dan antar-waktu ($\Delta t$). Pada skema *infill*, setiap observasi tambahan membawa satu realisasi baru $\varepsilon_{t,j}$, bukan subsampling tanpa derau.
- **Jatuh tempo (10 tenor)**: $\tau \in \{0{,}25,\; 0{,}5,\; 1,\; 2,\; 3,\; 5,\; 7,\; 10,\; 20,\; 30\}\text{ tahun}$.
- **Kondisi awal**: ditarik eksak dari distribusi stasioner $X_0 \sim \mathcal{N}(\theta^{\mathbb{P}}, Q_\infty)$, lalu dimajukan eksak menggunakan $\Phi(\Delta) = e^{-K^{\mathbb{P}}\Delta}$ dan $Q_\Delta$.
- **Seed & determinisme**: `MainSeed = 2026112202`, seed turunan deterministik per `(salt, repIdx)` dengan `xoshiro256++`, reduksi berurutan independen terhadap jumlah goroutine (`workers = 8` pada hasil penuh; bit-identik dengan `workers = 1`).

## 3. Desain Eksperimen (`R = 50` per Konfigurasi)

1. **Baseline terpisah (`baseline_1000w`)**: $T = 1{,}000$ minggu ($\Delta = 1/52$, $H \approx 19{,}2308$ tahun) untuk menguji reproduktibilitas terhadap `day004a` dalam batas MCSE gabungan.
2. **Long-span mingguan (`longspan_{10,20,40,80}y_weekly`)**: $H \in \{10, 20, 40, 80\}$ tahun ($T \in \{520, 1040, 2080, 4160\}$ minggu). Keempat horizon berbagi prefiks lintasan 80 tahun dan tebakan awal yang sama untuk mengisolasi efek $H$.
3. **Infill 20 tahun (`infill_20y_{monthly,weekly,daily}`)**: horizon tetap $H = 20$ tahun dengan $\Delta \in \{1/12, 1/52, 1/252\}$ ($T \in \{240, 1040, 5040\}$). Konfigurasi mingguan (`infill_20y_weekly`) memakai ulang hasil `longspan_20y_weekly` secara eksak dan terpisah dari `baseline_1000w`.

Seluruh ringkasan tetap menyertakan ke-$50$ replikasi apa adanya (termasuk kasus gagal konvergen atau solusi menyentuh batas pengaman).

## 4. Ringkasan Temuan Utama

### Validasi Metode & Jangkar Eksak (`c1..c4 = tercapai`, waktu total `41,93 detik`)
- **Jangkar Gaussian gabungan ($T=5$, 3 tenor)**: selisih absolut log-likelihood Kalman penuh `3,41e-13`, `FastLogLikelihood` tunak Riccati `3,55e-13`, rata-rata smoother RTS `2,83e-16`, kovarians smoother `6,39e-17` (seluruhnya `< 1e-9`).
- **Gradien beda hingga vs Richardson (koordinat terskala)**: galat relatif L2 maksimum sepanjang studi di titik awal `2,71e-9` dan di nilai benar `4,34e-7` (`< 1e-5`).
- **Determinisme worker**: identik bit-demi-bit antara `1` dan `8` worker pada fixture deterministik (`T=100, R=4`).
- **Baseline `1000w` vs `day004a`**: `11/11` parameter selisih biasnya berada dalam $3\text{ MCSE}$ gabungan (`max |z| = 2,82` pada `kp33`; `10/11` dalam $2\text{ MCSE}$).

### Long-Span (`10y -> 80y`): Peluruhan Bias $K^{\mathbb{P}}_{ii} \propto 1/H$
| Parameter | Bias $10\text{y}$ | Bias $20\text{y}$ | Bias $40\text{y}$ | Bias $80\text{y}$ | Tren OLS vs $1/H$ ($\pm$ prefix MCSE) | Tren WLS vs $1/H$ ($\pm$ prefix MCSE) |
|---|---:|---:|---:|---:|---:|---:|
| `kp11` (`0.05`) | $+0{,}4657 \pm 0{,}0802$ | $+0{,}2241 \pm 0{,}0290$ | $+0{,}0983 \pm 0{,}0130$ | $+0{,}0451 \pm 0{,}0081$ | $+4{,}839 \pm 0{,}945$ | $+4{,}701 \pm 0{,}719$ |
| `kp22` (`0.50`) | $+0{,}4254 \pm 0{,}0730$ | $+0{,}1998 \pm 0{,}0424$ | $+0{,}1007 \pm 0{,}0241$ | $+0{,}0731 \pm 0{,}0235$ | $+4{,}120 \pm 0{,}839$ | $+3{,}864 \pm 0{,}784$ |
| `kp33` (`0.40`) | $+0{,}3590 \pm 0{,}0725$ | $+0{,}2026 \pm 0{,}0425$ | $+0{,}0981 \pm 0{,}0278$ | $+0{,}0452 \pm 0{,}0184$ | $+3{,}562 \pm 0{,}878$ | $+3{,}778 \pm 0{,}791$ |

Produk `bias × H` berada pada rentang $+3{,}59\text{ hingga }+5{,}85$ di semua horizon (rata-rata dekat heuristik skalar $4/H$), namun korespondensinya dengan kasus skalar univariat berskala perkiraan karena faktor AFNS diekstrak melalui kurva lintas-tenor berderau oleh filter Kalman multivariat.

### Infill pada Horizon Tetap `H = 20 tahun` (`1/12 -> 1/52 -> 1/252`)
| Kelas Parameter | Rasio RMSE Harian / Bulanan ($T=5040$ vs $240$) | Pola Empiris |
|---|---|---|
| Drift $K^{\mathbb{P}}_{ii}$ (`kp11, kp22, kp33`) | `0,980`, `1,181`, `0,936` | Stagnan ($\sim +0{,}20\text{–}0{,}23$ bias tetap persisten) |
| Mean jangka panjang $\theta^{\mathbb{P}}_1$ | `0,503` (mingguan ke harian `0,819`) | Terbatas panjang rentang kalender ($\text{RMSE}$ $0{,}0181 \to 0{,}0148$) |
| Penampang Q (`lambda`) & Derau (`sigma_eps`) | `0,012`, `0,003` | Sangat tajam oleh informasi penampang & banyak observasi |
| Difusi $\sigma_1, \sigma_2, \sigma_3$ | `0,217`, `0,100`, `0,128` | Konsisten dengan laju variasi kuadratik $\sim \sqrt{240/5040} \approx 0{,}218$ |

### Cakupan Wald 95% & Rasio Hessian SE
- Parameter persisten $\theta_1$ tetap mengalami *undercoverage* berat di $H = 10\text{–}40\text{y}$ (`34%–48%`), membaik menjadi `72%` [Wilson 95%: `58,3%–82,5%`] di $H=80\text{y}$ dan `60%` [Wilson: `46,2%–72,4%`] pada *infill* harian 20y.
- Mekanisme teramati: varians asimtotik rata-rata sampel OU berbanding terbalik dengan $(K^{\mathbb{P}}_{11})^2 H$; saat $\hat{K}^{\mathbb{P}}_{11}$ terbias ke atas $3\text{–}5\times$, matriks informasi Hessian teramati *overconfident* (rasio `Median Hessian SE / Empirical SD` untuk $\theta_1$ hanya `0,06–0,39`), sehingga interval Wald menjadi terlalu sempit. Sebaliknya, pada 1 replikasi batas aktif (`Rep 31` 10y, `Rep 13` 40y/80y), invers Hessian membengkak dan menarik `Mean Hessian SE` jauh melampaui `Empirical SD`.

## 5. Cara Build, Run, dan Test

Jalankan dari direktori proyek:

```bash
cd projects/go/day004b_afns_inference
gofmt -l .
go vet ./...
go test -race -count=1 ./...
go run ./cmd/run           # Runner penuh R=50 (~42s di 8 core) -> outputs/day004b/results.json & outputs/figures/
go run ./cmd/run -quick    # Mode cepat pengembangan (R=4, menulis hanya ke direktori sementara sistem)
```

## 6. Struktur Kode

- `internal/rng/` — Pembangkit deterministik `xoshiro256++` + `SplitMix64` (`MainSeed = 2026112202`).
- `internal/afns/` — Kurva Nelson–Siegel bebas pembatalan, koreksi konveksitas CDR 2011, transisi diskret eksak $\Phi(\Delta)$, $Q_\Delta$, $Q_\infty$, dan simulator stasioner $\mathbb{P}$.
- `internal/optim/` — Gradien beda hingga pusat, verifikasi Richardson, Hessian pusat, dan BFGS terkendala kotak dengan pelacakan batas per parameter.
- `internal/kalman/` — `FastLogLikelihood` dengan akselerasi steady-state Riccati, `FilterAndSmooth` (Joseph + RTS), `VerifyJointGaussianAnchor`, `EstimateMLE`, serta penggerak studi *baseline*, *long-span*, dan *infill*.
- `internal/stats/`, `internal/parallel/`, `internal/jsonout/`, `internal/svg/` — Statistik (Wilson CI, tren OLS/WLS + prefix MCSE), paralelisme deterministik per chunk, serialisasi JSON presisi 17 digit, dan generator SVG.
- `cmd/run/` — Runner CLI utama dan ekspor hasil JSON serta 5 grafik SVG.

## 7. Grafik & Artefak Terkait

- JSON keluaran: [`../../outputs/day004b/results.json`](../../../outputs/day004b/results.json)
- Grafik SVG:
  - [`day004b_kp-bias-vs-inv-h.svg`](../../../outputs/figures/day004b_kp-bias-vs-inv-h.svg)
  - [`day004b_infill-vs-longspan-rmse.svg`](../../../outputs/figures/day004b_infill-vs-longspan-rmse.svg)
  - [`day004b_wald-coverage-intervals.svg`](../../../outputs/figures/day004b_wald-coverage-intervals.svg)
  - [`day004b_hessian-se-to-emp-sd-ratio.svg`](../../../outputs/figures/day004b_hessian-se-to-emp-sd-ratio.svg)
  - [`day004b_convergence-and-bound-hits.svg`](../../../outputs/figures/day004b_convergence-and-bound-hits.svg)
- Penjelasan rinci 11 bagian: [`day004b_quant_afns-inference.md`](../../../docs/project-explanations/day004b_quant_afns-inference.md)
