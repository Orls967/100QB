# Hari 3 — Sesi f (`day003f`): Teori Ruin Cramér–Lundberg dan Simulasi Kejadian Langka

Proyek Go 1.22+ murni (*zero external dependencies* — hanya *standard library*) yang mengimplementasikan teori ruin Cramér–Lundberg klasik dan modern: derivasi bentuk tertutup lewat transformasi Laplace rasional, rekursi batas atas/bawah Panjer/Beekman, *importance sampling* (IS) eksponensial berbasis eksponen Lundberg, estimator kondisional Asmussen–Kroese untuk klaim berekor berat (*heavy-tailed*), serta optimisasi retensi reinsurans *excess-of-loss* dan modal minimum $u^*(\varepsilon)$.

## Intuisi Dasar: Ruin, Eksponen Lundberg, dan Kejadian Langka

Perusahaan asuransi memulai dengan cadangan modal awal $u \ge 0$, menerima premi kontinu dengan laju $c = (1 + \theta)\lambda\mu$ ($\theta = 0.2 > 0$ adalah *safety loading*), dan membayar klaim majemuk Poisson $S_t = \sum_{i=1}^{N_t} X_i$ bermuara ke proses surplus $U_t = u + c t - S_t$. Peluang kebangkrutan (*ruin*) sepanjang waktu $\psi(u) = \mathbb{P}(\inf_{t \ge 0} U_t < 0 \mid U_0 = u)$ tetap positif meski drift bersihnya positif ($\theta\lambda\mu > 0$), karena rangkaian klaim ekstrem dapat menggerus surplus di bawah nol.

- **Ekor ringan (*light-tailed*)**: Ketika fungsi pembangkit momen $M_X(r) = \mathbb{E}[e^{r X}]$ berhingga di sekitar nol dan persamaan Lundberg $\lambda(M_X(R) - 1) = c R$ memiliki akar positif unik $R > 0$, peluang ruin meluruh secara eksponensial ($\psi(u) \sim C e^{-R u}$). Simulasi Monte Carlo biasa mustahil mengukur $\psi(100) \approx 10^{-10}$, sehingga *exponential tilting* menuju ukuran $\mathbb{Q}_R$ memaksa drift surplus menjadi negatif dan menstabilkan galat baku relatif di bawah `0.07%` untuk semua $u$.
- **Ekor berat (*heavy-tailed* / subeksponensial)**: Pada Pareto ($\alpha = 2.5$) dan Weibull ($k = 0.5$), $M_X(r) = \infty$ untuk setiap $r > 0$, sehingga eksponen Lundberg maupun IS eksponensial gagal (*variansi tak hingga*). Kebangkrutan cenderung terjadi akibat **satu klaim tunggal raksasa** (*principle of a single big jump*). Estimator bersyarat Asmussen–Kroese memanfaatkan simetri penjumlahan acak Pollaczek–Khinchine untuk menekan galat jauh melampaui simulasi naif tanpa kehilangan tak-bias.

## Cara Build, Test, dan Run

```bash
cd projects/go/day003f_ruin_theory
gofmt -l .
go vet ./...
go test -race -count=1 ./...
go run ./cmd/run
```

## Struktur Kode

- [`cmd/run/main.go`](cmd/run/main.go) — eksekutor seluruh eksperimen Bagian 1–5 ke [`results.json`](../../../outputs/day003f/results.json) (presisi 17 digit) dan 7 grafik vektor SVG.
- [`internal/claims/`](internal/claims/claims.go) — lima distribusi klaim ($\mu = 1$): (A) Eksponensial(`1`), (B) Erlang(`2`, rate `2`), (C) Hipereksponensial-2 (`p = 0.5`, $\beta_1 = 2/3$, $\beta_2 = 2$, $\text{CV} = \sqrt{1.5} > 1$), (D) Pareto Lomax ($\alpha = 2.5$, skala `1.5`), dan (E) Weibull ($k = 0.5$, skala `0.5`); lengkap dengan distribusi tinggi tangga $F_L$, sampler *tilted*, dan MGF.
- [`internal/lundberg/`](internal/lundberg/lundberg.go) & [`internal/closedform/`](internal/closedform/closedform.go) — akar Brent $R$, koefisien Cramér–Lundberg $C$, serta solusi tertutup rasional lewat transformasi Laplace dan sistem linear LU.
- [`internal/panjer/`](internal/panjer/panjer.go) — rekursi batas bawah/atas Panjer/Beekman untuk campuran geometrik terdiskretisasi (pembulatan bawah/atas), orde konvergensi vs $h$, serta perbandingan formulasi ekor positif vs $1 - \sum g_k$.
- [`internal/process/`](internal/process/process.go) — simulasi proses eksak antar-klaim ($U_t$, defisit saat ruin + uji KS dengan nilai kritis simulasi nol sendiri, dan MC biasa $1.000.000$ lintasan).
- [`internal/is/`](internal/is/is.go) — *importance sampling* eksponensial (optimal $\tilde{r} = R$ maupun *misspecified* $0.7R, 1.3R$) dan pembuktian divergensi $M_X(r)$ pada ekor berat.
- [`internal/ak/`](internal/ak/ak.go) — estimator kondisional Asmussen–Kroese berbobot-ukuran (*size-biased*) untuk klaim Pareto & Weibull hingga $u = 200$.
- [`internal/reinsurance/`](internal/reinsurance/reinsurance.go) — retensi *excess-of-loss* $\min(X, M)$ ($\xi = 0.4$), rumus tertutup langkah berhingga untuk (A), batas Panjer, optimisasi $M^*$, dan kebutuhan modal $u^*(\varepsilon)$ via Brent.
- Utilitas bersama murni stdlib: [`rng`](internal/rng/rng.go) (`xoshiro256++`), [`special`](internal/special/special.go) (`BrentRoot`, `GoldenSectionMin`, `SolveLinearLU`), [`stats`](internal/stats/stats.go) (uji KS & simulasi nol $O(n)$), [`parallel`](internal/parallel/parallel.go), [`jsonout`](internal/jsonout/jsonout.go), dan [`svg`](internal/svg/svg.go).

## Ringkasan Hasil Utama

| Bagian | Temuan & Angka Utama (`results.json`) |
|---|---|
| **Bagian 1 — Bentuk Tertutup & Defisit** | (A) $R = 1/6$ (`abs_err = 4.16e-16`), (B) $R_1 = 0.226765$ (`err = 1.61e-15`), (C) $R_1 = 0.130502$ (`err = 5.83e-16`). Momen `200k` $E[N_{10}], E[U_{10}]$ lulus `3 SE`; defisit (A) di $u \in \{2, 5\}$ lolos uji KS ($p = 0.394, 0.985$). |
| **Bagian 2 — Rekursi Panjer** | Seluruh batas mengapit eksak untuk (A)–(C) pada semua $(h, u)$ dengan *slope* lebar vs $h$ = `0.997–0.999` (orde 1). Rasio $\psi/(C e^{-Ru}) \to 1$; $\psi_{\text{Weibull}}(200)$ `~4.70e-5` masih `217x` di atas asimtotik subeksponensial. |
| **Bagian 3 — MC Biasa vs IS** | MC `1.000.000` cocok dalam `1.47 SE`. IS eksponensial mencapai $\psi_A(100) = 4.817\times 10^{-8}$ (`relSE = 0.053%`, VRF waktu `1.34e8x`) dan $\psi_B(100) = 1.209\times 10^{-10}$ (`relSE = 0.045%`, VRF waktu `1.42e11x`). Tilt $0.7R$ & $1.3R$ tetap tak-bias namun bervariansi `3.3x–10.1x` lebih besar. |
| **Bagian 4 — Ekor Berat & AK** | Integral $M_X(0.1)$ meledak di $x_{\max}=500$ (`1.36e14` Pareto, `4.53e7` Weibull). AK (`200k`) cocok ketat dengan batas Panjer $h=0.01$ hingga $u = 200$: Pareto `0.004144` (`relSE = 0.092%`, VRF waktu `3.49e5x`) dan Weibull `4.526e-5` (`relSE = 5.94%`, VRF waktu `1.13e4x`). |
| **Bagian 5 — Reinsurans & Modal** | Di $u = 10$ ($\xi = 0.4$): retensi optimal $M^*_A = 1.529$ memangkas $\psi_A(10)$ `0.1573 -> 0.0941` dan menghemat modal $u^*(10^{-4})$ sebesar `25.65%` (`54.17 -> 40.28`); $M^*_C = 1.760$ memangkas $\psi_C(10)$ `0.2206 -> 0.1255` dan menghemat modal `33.13%` (`69.00 -> 46.15`). |

## Grafik SVG Keluaran

1. [Peluang ruin eksak vs batas Panjer](../../../outputs/figures/day003f_psi-exact-vs-panjer.svg)
2. [Orde diskretisasi Panjer log-log](../../../outputs/figures/day003f_panjer-discretization-order.svg)
3. [Konvergensi rasio Cramér–Lundberg](../../../outputs/figures/day003f_cramer-lundberg-ratio.svg)
4. [SE relatif IS eksponensial vs MC biasa](../../../outputs/figures/day003f_is-vs-mc-rel-se.svg)
5. [Ekor berat: Asmussen–Kroese vs Panjer vs asimtotik](../../../outputs/figures/day003f_heavy-tail-ak-vs-asymptotic.svg)
6. [Kurva $\psi_M(10)$ terhadap batas retensi $M$](../../../outputs/figures/day003f_reinsurance-psi-vs-m.svg)
7. [Modal awal minimum $u^*(\varepsilon)$ dengan vs tanpa reinsurans](../../../outputs/figures/day003f_capital-u-star.svg)

Analisis teori terperinci dan jejak verifikasi lengkap tersedia di [day003f_quant_ruin-theory.md](../../../docs/project-explanations/day003f_quant_ruin-theory.md).
