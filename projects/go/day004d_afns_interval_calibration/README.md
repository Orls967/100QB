# day004d — AFNS: Kalibrasi Interval Bootstrap yang Peka Bias dan Batas Parameter

Modul Go 1.22+ mandiri (`zero external dependencies`) untuk membandingkan lima metode konstruksi selang kepercayaan 95% (**Wald Hessian**, **persentil bootstrap**, **basic / reversed percentile**, **bootstrap-t / studentized**, dan **log-scale basic** untuk parameter *mean-reversion* positif) secara berpasangan pada dataset luar dan *draw* bootstrap yang sama untuk model ruang-keadaan *Arbitrage-Free Nelson–Siegel* (AFNS) 11 parameter.

## 1. Intuisi Masalah

Pada `day004c`, bootstrap parametrik persentil naif gagal memulihkan cakupan 95% pada dua parameter kunci horizon $H = 20$ tahun mingguan ($T = 1{,}040$):
- Cakupan `kp11` jatuh ke `2%` karena estimator $\hat{\kappa}_1$ mengalami bias ke atas ganda (pertama pada data luar, lalu tergeser ke atas lagi dalam simulasi bootstrap).
- Cakupan `theta1` hanya mencapai `42%` (mirip Wald `46%`) karena simulasi bootstrap di bawah $\hat{\kappa}_1$ yang terlalu tinggi mempersempit sebaran $\hat{\theta}_1^*$.

Proyek `day004d` memperluas studi ke $R = 100$ dataset luar $\times B = 100$ refit bootstrap dalam ($10{,}000$ refit lengkap dengan evaluasi matriks Hessian dan *standard error* $\widehat{\text{SE}}^*$ di setiap refit) untuk menjawab lima pertanyaan kalibrasi secara empiris tanpa memangkas batas (*unclamped*).

## 2. Spesifikasi Model, Satuan, dan Desain Eksperimen

Spesifikasi model identik dengan `day004b` dan `day004c`:
- **Yield & $\theta^{\mathbb{P}}$**: desimal per tahun ($\theta^{\mathbb{P}} = (0{,}045,\; -0{,}020,\; 0{,}010)$ sesuai $4{,}5\%$, $-2{,}0\%$, $1{,}0\%$).
- **Peluruhan $\mathbb{Q}$ & *mean-reversion* $\mathbb{P}$**: $\lambda = 0{,}6\text{ thn}^{-1}$, $K^{\mathbb{P}} = \operatorname{diag}(0{,}05,\; 0{,}5,\; 0{,}4)\text{ thn}^{-1}$.
- **Volatilitas difusi faktor**: $\Sigma = \operatorname{diag}(0{,}0060,\; 0{,}0100,\; 0{,}0200)\text{ thn}^{-1/2}$.
- **Derau observasi**: $\sigma_\varepsilon = 0{,}0005$ ($5\text{ bp}$), i.i.d. antar-tenor dan antar-waktu ($\Delta = 1/52\text{ thn}$).
- **Tenor (10 tenor)**: $\tau \in \{0{,}25,\; 0{,}5,\; 1,\; 2,\; 3,\; 5,\; 7,\; 10,\; 20,\; 30\}\text{ tahun}$.
- **Horizon & Kondisi Awal**: $H = 20$ tahun ($T = 1{,}040$ minggu), $X_0 \sim \mathcal{N}(\theta^{\mathbb{P}}, Q_\infty)$ stasioner $\mathbb{P}$.
- **Seed & Determinisme**: `MainSeed = 2026112204` (`xoshiro256++` + `SplitMix64`), seed turunan tetap per `(salt, rep_idx)` dan `(salt, outer_idx, boot_idx)`, reduksi terurut; bit-identik `1` vs `8` worker.

## 3. Metode Interval yang Dibandingkan

Seluruh metode dievaluasi pada 100 fit luar dan $10{,}000$ *draw* bootstrap yang sama, dengan ambang validitas pra-run `MinValidFraction = 0.95` ($\ge 95/100$ *draw* valid):
1. **Wald (`wald`)**: $\hat{\psi} \pm z_{0{,}975}\,\widehat{\text{SE}}_{\text{Hess}}$.
2. **Persentil (`percentile`)**: $[q_{0{,}025}(\psi^*),\; q_{0{,}975}(\psi^*)]$ pada skala asli.
3. **Basic / Reversed Percentile (`basic`)**: $[2\hat{\psi} - q_{0{,}975}(\psi^*),\; 2\hat{\psi} - q_{0{,}025}(\psi^*)]$ tanpa pemangkasan ke batas fisik.
4. **Bootstrap-t / Studentized (`studentized`)**: menggunakan pivot $t^*_b = (\psi^*_b - \hat{\psi}) / \widehat{\text{SE}}^*_b$ dengan selang $[\hat{\psi} - q_{0{,}975}(t^*)\,\widehat{\text{SE}},\; \hat{\psi} - q_{0{,}025}(t^*)\,\widehat{\text{SE}}]$. Draw dengan $\psi^*$ non-finite atau $\widehat{\text{SE}}^* \le 0$ dikeluarkan; bila draw valid $< 95$ atau $\widehat{\text{SE}} \le 0$, interval ditandai *invalid* namun tetap masuk denominator $R = 100$.
5. **Log-Scale Basic (`log_basic`)**: diterapkan pada parameter positif ($\lambda, K^{\mathbb{P}}_{ii}, \sigma_i, \sigma_\varepsilon$) melalui $\eta = \ln\psi$, membangun interval dasar di ruang log, lalu mengeksponensialkan kembali $[\exp(2\hat{\eta} - q_{0{,}975}(\eta^*)),\; \exp(2\hat{\eta} - q_{0{,}025}(\eta^*))]$.

## 4. Ringkasan Hasil Utama (`outputs/day004d/results.json`, `369,51s` pada `8` worker)

### Verifikasi Teknis & Diagnostik Optimizer (`c1..c5 = tercapai`)
- **Jangkar Gaussian gabungan ($T=5$, 3 tenor)**: selisih log-likelihood Kalman `5,68e-13`, `FastLogLikelihood` `5,54e-13`, smoother RTS `3,16e-16` (mean) & `6,39e-17` (kovarians).
- **Skor FD vs Richardson**: galat relatif maksimum di titik awal `1,54e-9` dan di nilai benar `3,38e-7` (`< 1e-5`).
- **Determinisme worker**: identik bit-demi-bit antara `1` dan `8` worker pada fixture cepat (`T=80, R=2, B=4`).
- **Diagnostik Optimizer & Validitas**:
  - Fit luar (`100` total): `100/100` *finite*, `97/100` `scaled_grad <= 1e-6`, `99` `line_search_failed` (`86` `step_underflow_plateau`, `13` `armijo_not_met`, `52` berobjektif identik bit-demi-bit), `1` `max_iterations`, `2/100` menyentuh batas & regularisasi Hessian.
  - Refit bootstrap dalam (`10.000` total): `10.000/10.000` *finite*, `9.875/10.000` (`98,75%`) `scaled_grad <= 1e-6`, `9.997` `line_search_failed` (`9.118` `step_underflow_plateau`, `879` `armijo_not_met`, `5.564` berobjektif identik bit-demi-bit), `3` `max_iterations`, `171` (`1,71%`) batas aktif, `169` (`1,69%`) regularisasi Hessian.
  - Seluruh interval (`100/100` per metode & parameter) valid (`0` invalid).

### Cakupan 95%, Lebar Interval, dan Pelanggaran Domain (`kp11`, `kp22`, `kp33`, `theta1`)

| Param (True) | Metode | Cakupan [Wilson 95%] | Lebar (Mean / Med) | Non-Fisik ($\le 0$) | OOB |
|---|---|---|---:|---:|---:|
| `kp11 (0.05)` | Wald | `85,0%` `[76,7%, 90,7%]` | `0,5466 / 0,5094` | `79/100 (79%)` | `79%` |
| `kp11 (0.05)` | Percentile | `14,0%` `[8,5%, 22,1%]` | `0,8896 / 0,8819` | `0/100 (0%)` | `0%` |
| `kp11 (0.05)` | Basic | `95,0%` `[88,8%, 97,9%]` | `0,8896 / 0,8819` | `100/100 (100%)` | `100%` |
| `kp11 (0.05)` | Studentized | `86,0%` `[77,9%, 91,5%]` | `0,4856 / 0,4405` | `87/100 (87%)` | `87%` |
| `kp11 (0.05)` | Log-Basic | `49,0%` `[39,4%, 58,7%]` | `0,3635 / 0,3178` | `0/100 (0%)` | `2%` |
| `kp22 (0.50)` | Wald / Pct / Basic / Stud / LogBasic | `90% / 67% / 92% / 87% / 86%` | med `0,99 / 1,22 / 1,22 / 0,97 / 0,83` | `12% / 0% / 90% / 36% / 0%` | `12%..90%` |
| `kp33 (0.40)` | Wald / Pct / Basic / Stud / LogBasic | `95% / 62% / 91% / 94% / 80%` | med `0,90 / 1,10 / 1,10 / 0,79 / 0,67` | `20% / 0% / 93% / 47% / 0%` | `20%..93%` |
| `theta1 (0.045)` | Wald | `52,0%` `[42,3%, 61,5%]` | `1,5871 / 0,0215` | `0%` | `2%` |
| `theta1 (0.045)` | Percentile | `46,0%` `[36,6%, 55,7%]` | `0,0302 / 0,0202` | `0%` | `0%` |
| `theta1 (0.045)` | Basic | `47,0%` `[37,5%, 56,7%]` | `0,0302 / 0,0202` | `0%` | `2%` |
| `theta1 (0.045)` | Studentized | `72,0%` `[62,5%, 79,9%]` | `19,1171 / 0,0462` | `0%` | `3%` |

**Kesimpulan Utama & Keterbatasan**:
- **Trade-off `kp11` skala asli vs log**: `Basic` mencapai cakupan `95%` namun `100%` intervalnya berbatas bawah negatif (non-fisik). `Log-Basic` menjamin batas positif (`0%` non-fisik), tetapi cakupan `kp11` turun ke `49%` karena transformasi rasio log mengoreksi kurang kuat saat bias aditif besar relatif terhadap nilai benar `0,05`.
- **Peningkatan parsial `theta1` lewat *studentization***: `Studentized` menaikkan cakupan `theta1` secara signifikan dari `46%–52%` menjadi `72,0%` (median lebar naik lebih dari $2\times$ menjadi `0,0462`), namun tetap di bawah `95%` dan sensitif terhadap nilai $\hat{\kappa}_1$ luar (`94,0%` bila $\hat{\kappa}_1 \le 0{,}208$ vs `50,0%` bila $\hat{\kappa}_1 > 0{,}208$) serta dua outlier batas aktif yang meledakkan rata-rata lebar menjadi `19,12`.

## 5. Cara Build, Run, dan Test

```bash
cd projects/go/day004d_afns_interval_calibration
gofmt -l .
go vet ./...
go test -race -count=1 -v ./...  # 11 paket, 14 tes utama
go run ./cmd/run               # Runner penuh R=100, B=100 (~370s pada 8 worker)
go run ./cmd/run -quick        # Mode cepat pengembangan (R=4, B=6 ke temp dir)
```

## 6. Struktur Kode

- `internal/rng/` — Pembangkit deterministik `xoshiro256++` + `SplitMix64` (`MainSeed = 2026112204`).
- `internal/parallel/` — `DefaultWorkers` dan `RunChunks` deterministik.
- `internal/afns/` — Pemuat kurva Nelson–Siegel, koreksi konveksitas CDR 2011, transisi eksak $\Phi(\Delta)$, $Q_\Delta$, $Q_\infty$, dan simulator stasioner.
- `internal/optim/` — Beda hingga pusat, Richardson, Hessian, dan `MinimizeBFGS` terkendala kotak dengan diagnostik pencarian garis lengkap.
- `internal/kalman/` — `FastLogLikelihood` (*steady-state* Riccati), `FilterAndSmooth`, `VerifyJointGaussianAnchor`, aljabar matriks (`linalg.go`), dan `EstimateMLEWithOptions`.
- `internal/interval/` — Implementasi 5 metode interval (`EvalWald`, `EvalPercentile`, `EvalBasic`, `EvalStudentized`, `EvalLogBasic`) tanpa pemangkasan serta validitas pra-run.
- `internal/bootstrap/` — `RunOuterFits`, `RunCalibrationStudy` (`R x B` dengan Hessian $\widehat{\text{SE}}^*$), agregasi diagnostik optimizer, dan stratifikasi sensitivitas `kp11` vs `theta1`.
- `internal/stats/`, `internal/jsonout/`, `internal/svg/`, `cmd/run/` — Statistik (`WilsonCI`, `CoverageWithValidity`), JSON 17-digit, pembuat SVG, dan runner CLI.

## 7. Grafik & Artefak Terkait

- JSON keluaran: [`outputs/day004d/results.json`](../../../outputs/day004d/results.json)
- Grafik SVG:
  - [`day004d_coverage-by-method-kp11-theta1.svg`](../../../outputs/figures/day004d_coverage-by-method-kp11-theta1.svg)
  - [`day004d_kp11-endpoints-raw-vs-log.svg`](../../../outputs/figures/day004d_kp11-endpoints-raw-vs-log.svg)
  - [`day004d_theta1-interval-widths.svg`](../../../outputs/figures/day004d_theta1-interval-widths.svg)
  - [`day004d_kp11-vs-theta1-width-coverage.svg`](../../../outputs/figures/day004d_kp11-vs-theta1-width-coverage.svg)
- Penjelasan proyek 11 bagian: [`day004d_quant_afns-interval-calibration.md`](../../../docs/project-explanations/day004d_quant_afns-interval-calibration.md)
