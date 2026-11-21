# Hari 3 — Sesi e (`day003e`): CVA Swap Suku Bunga di Bawah Hull–White 1-Faktor, Intensitas Kredit CIR, *Wrong-Way Risk*, Netting, dan Kolateral

Proyek Go 1.22+ murni (*zero external dependencies* — hanya *standard library*) yang membangun mesin valuasi dan *Counterparty Credit Risk* (CCR/CVA) *end-to-end* untuk swap suku bunga 10 tahun *payer ATM* dengan model suku bunga jangka pendek Hull–White 1-faktor dan intensitas gagal stokastik Cox–Ingersoll–Ross (CIR).

## Struktur Paket

- [`cmd/run/main.go`](cmd/run/main.go) — *runner* eksperimen penuh ($200{,}000$ lintasan eksposur/swaption, $100{,}000$ lintasan CVA/WWR/netting/kolateral pada 8 goroutine), penulis [`results.json`](../../../outputs/day003e/results.json) (presisi 17 digit), dan generator 7 grafik vektor SVG.
- [`internal/rng/`](internal/rng/rng.go) — generator `xoshiro256++` deterministik berbenih `SplitMix64` (`DeriveChunkSeed`), normal Box–Muller, eksponensial, Gamma Marsaglia–Tsang + Ahrens–Dieter ($\alpha < 1$), Beta, serta Poisson eksak (Knuth + dekomposisi Gamma–Binomial Ahrens–Dieter).
- [`internal/special/`](internal/special/special.go) — fungsi normal baku (`NormPDF`, `NormCDF`, Acklam+Halley `NormQuantile`) serta kuadratur Golub–Welsch untuk Gauss–Hermite dan Gauss–Legendre.
- [`internal/hw/`](internal/hw/hw.go) — model Hull–White 1-faktor ($r_0 = 0.03$, $a = 0.08$, $\sigma = 0.01$): dekomposisi eksak $r_t = x_t + \alpha(t)$, langkah transisi bivariat bersama eksak $(x_{t+\Delta}, \int x_s\,ds)$ bebas pembatalan numerik untuk $\Delta$ kecil, uji martingale 1D/2D/MC, harga swaption 2y×8y Jamshidian vs kuadratur bebas-*kink*, profil eksposur kuartalan (EE, EE terdiskonto, PFE 95%/99%, EPE), serta simulasi dua-waktu eksak untuk kolateral variasi margin + MPOR.
- [`internal/cir/`](internal/cir/cir.go) — model intensitas gagal CIR ($\lambda_0 = 0.02$, $\kappa = 0.5$, $\theta = 0.025$, $\xi = 0.1$, $R = 0.4$): pemeriksaan syarat Feller ($2\kappa\theta = 0.025 \ge \xi^2 = 0.01$), kurva survival afin eksak $Q(T) = A(T)e^{-B(T)\lambda_0}$ vs integrasi numerik ODE Riccati RK4, transisi eksak campuran Poisson–Gamma ($\chi^2$ tak-sentral), Euler *full-truncation*, dan par spread CDS kuartalan.
- [`internal/cva/`](internal/cva/cva.go) — CVA semi-analitik independen & sensitivitas ($\partial\text{CVA}/\partial\lambda_0$ analitik vs beda hingga pusat, $\text{CS01}$ terhadap $\lambda_0$ dan $\theta$, $\sigma_{\text{HW}}$, strike $K$), CVA MC waktu gagal eksplisit ($\tau$ lewat $u \sim \text{Exp}(1)$) + estimator Rao–Blackwellized dengan sub-langkah Brown bersarang $\{20, 40, 80\}$ terhadap referensi $m = 320$, kurva *Wrong-Way Risk* $\text{CVA}(\rho)$ untuk $\rho \in \{-0.8, -0.4, 0, 0.4, 0.8\}$ via *Common Random Numbers* (CRN), netting portofolio (`5y payer + 10y receiver` dan identitas nol), serta kolateral MPOR $\delta$ dengan threshold $H$.
- [`internal/parallel/`](internal/parallel/parallel.go), [`internal/stats/`](internal/stats/stats.go), [`internal/jsonout/`](internal/jsonout/jsonout.go), [`internal/svg/`](internal/svg/svg.go) — pekerja goroutine deterministik berbasis *chunk*, statistik sampel, penulis JSON 17 digit, dan mesin plot SVG.

## Menjalankan Pengujian dan Eksperimen

```bash
cd projects/go/day003e_cva_wrong_way
go vet ./...
go test -race -count=1 ./...
go run ./cmd/run
```

## Ringkasan Hasil Utama (`results.json`)

Semua angka di bawah ini dilaporkan per `1` unit nosional (dan ekuivalen per `$1.000.000` nosional):

| Komponen / Kriteria | Hasil Kunci | Acuan / Verifikasi |
|---|---|---|
| **Martingale Hull–White** ($t \in \{2,5,8\}, T=10$) | Galat kuadratur maks `8.88e-16` | MC `200k` semua $|z| \le 0.58$ |
| **Swaption 2y×8y Jamshidian** (ITM / ATM / OTM) | `0.071464`, `0.026532`, `0.005768` | Selisih vs kuadratur 2D maks `1.11e-16`; MC $|z| \le 0.95$ |
| **Profil Eksposur 10y Payer** | `EPE = 0.017726` (`$17,726/1M`), puncak `EE = 0.031838` di $t=4$ | Galat MC `200k` vs analitik di 41 tanggal maks `z = 2.01` |
| **Kurva Survival CIR & Par CDS** | $2\kappa\theta = 0.025 \ge \xi^2 = 0.01$; `Q(10) = 0.789127` | RK4 ODE selisih maks `1.33e-15`; par CDS `126.7 / 134.0 / 137.8 / 141.7 bps` |
| **CVA Independen ($\rho = 0$)** | `CVA_anal = 0.00191515` (`$1,915.15/1M`) | $\partial\text{CVA}/\partial\lambda_0$ selisih vs FD pusat `8.08e-14`; RB MC `z = 0.24`, $\tau$ MC `z = 0.16` |
| **Wrong-Way Risk** ($\rho: -0.8 \to +0.8$) | `$1,080.34 -> $1,468.94 -> $1,912.75 -> $2,410.41 -> $2,962.22 / 1M` | Monoton naik (`0.565x -> 1.549x` terhadap $\rho = 0$) |
| **Netting (`5y payer + 10y receiver`)** | Analitik `$2,220.83/1M` vs jumlah mandiri `$3,504.96/1M` (hemat `36.9%` pada $\rho=0$) | Pelanggaran subaditivitas `0`; `10y payer + 10y receiver = 0` eksak |
| **Kolateral ($\delta = 0.04$, $H \in \{0, 0.01, 0.05, \infty\}$)** | Pada $\rho = 0$: `$112.77 -> $383.84 -> $1,376.49 -> $1,912.87 / 1M` | Pada $H = 0$ jadwal bersih, $\delta: 0.04 \to 0.01 \to 0.0025$ menyusut `$112.77 -> $61.26 -> $31.99 / 1M` (*slope* log-log `0.454`) |

## Artefak Grafik SVG

1. [Profil EE & PFE Swap 10Y ATM](../../../outputs/figures/day003e_ee-pfe-profile.svg)
2. [Kurva Survival CIR: Analitik vs RK4 vs MC](../../../outputs/figures/day003e_cir-survival.svg)
3. [Kurva Wrong-Way Risk: CVA terhadap Korelasi $\rho$](../../../outputs/figures/day003e_cva-vs-rho.svg)
4. [Konvergensi Bias Diskretisasi Integral Intensitas](../../../outputs/figures/day003e_trapezoid-discretization-bias.svg)
5. [Perbandingan Eksposur Portofolio Netting](../../../outputs/figures/day003e_netting-exposure.svg)
6. [Profil Eksposur Terkolateralisasi Lintas Threshold $H$](../../../outputs/figures/day003e_collateral-thresholds.svg)
7. [Konvergensi MPOR $\delta \to 0$ pada $H = 0$](../../../outputs/figures/day003e_mpor-convergence.svg)

Penjelasan teori lengkap, audit metode, dan dokumentasi jejak verifikasi tersedia di [day003e_quant_cva-wrong-way.md](../../../docs/project-explanations/day003e_quant_cva-wrong-way.md).
