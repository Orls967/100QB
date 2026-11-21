# Hari 4 — Sesi a (`day004a`): Struktur Suku Bunga: Nelson–Siegel, AFNS, dan Filter Kalman

Proyek Go 1.22+ murni (*zero external dependencies* — hanya *standard library*) yang mengimplementasikan pemodelan struktur suku bunga (*term structure of interest rates*) secara menyeluruh dari kurva statis bebas-arbitrase hingga estimasi ruang-keadaan (*state-space*) dinamis dan dekomposisi komponen utama (PCA). Seluruh aljabar linear (LU, Cholesky, tridiagonal Thomas, eigen simetris Jacobi, eksponensial matriks), kuadratur Gauss–Legendre, integrator ODE RK4, pembangkit acak `xoshiro256++`, serta optimizer BFGS ditulis sendiri dari nol.

## Intuisi Dasar: Dari Kurva Nelson–Siegel ke AFNS dan Filter Kalman

1. **Nelson–Siegel(–Svensson) & Bootstrap**: Kurva imbal hasil (*zero-coupon yield*) $y(\tau)$ diparameterisasi oleh tiga faktor berinterpretasi ekonomi—*level* ($L$), *slope* ($S$), dan *curvature* ($C$) dengan peluruhan $\lambda$—sehingga setiap titik tenor konsisten dengan kurva *forward* seketika $f(\tau)$ melalui identitas rata-rata $y(\tau) = \frac{1}{\tau}\int_0^\tau f(s)\,ds$. Dari harga kupon par obligasi di pasar, *bootstrapping* berurutan memulihkan faktor diskonto $P(\tau)$ eksak, sementara perbandingan interpolasi log-linear $P$ vs spline kubik $y$ mengungkap pertukaran antara ketepatan titik *yield* dan getaran kurva *forward* ketika data mengandung derau.
2. **AFNS (*Arbitrage-Free Nelson–Siegel*)**: Model Nelson–Siegel statis belum menjamin ketiadaan arbitrase di bawah pergerakan stokastik. Dengan menetapkan dinamika Q-netral-risiko $dX_t = K^{\mathbb{Q}}(\theta^{\mathbb{Q}} - X_t)\,dt + \Sigma\,dW_t^{\mathbb{Q}}$ bersama $r_t = X_{1,t} + X_{2,t}$, ODE Riccati atas harga obligasi $P(t, T) = \exp(A(\tau) + B(\tau)^\top X_t)$ memulihkan pemuatan faktor Nelson–Siegel $-B(\tau)/\tau = (1, g_1(\tau), g_2(\tau))$ secara persis ditambah suku penyesuaian konveksitas $-A(\tau)/\tau < 0$ (Christensen–Diebold–Rudebusch 2011). Pemisahan dinamika $\mathbb{Q}$ (pembentukan harga) vs dinamika historis $\mathbb{P}$ ($K^{\mathbb{P}}, \theta^{\mathbb{P}}$) langsung mendefinisikan harga risiko afinitas dan premi jangka waktu (*term premium*).
3. **Filter Kalman, Smoother RTS & MLE**: Karena faktor laten $X_t$ tidak teramati langsung dan seluruh kurva *yield* tercampur derau pengukuran $\varepsilon_t$, filter Kalman mengekstrak faktor terfilter dan menghitung *log-likelihood* Gaussian eksak melalui dekomposisi inovasi, sementara *smoother* Rauch–Tung–Striebel (RTS) memulihkan lintasan faktor terbaik menggunakan seluruh sampel. Pengoptimalan BFGS atas 11 parameter bebas memverifikasi ketepatan estimasi dan menyingkap batas identifiabilitas parameter peluruhan jangka panjang ($K^{\mathbb{P}}_{11}, \theta_1^{\mathbb{P}}$) pada sampel terbatas.
4. **PCA Perubahan Yield**: Dekomposisi eigen atas matriks kovarians teoretis perubahan yield mingguan $\operatorname{Cov}(\Delta y_t)$ membuktikan bagaimana tiga faktor laten AFNS bersama derau pengukuran memunculkan pola klasik Litterman–Scheinkman (*level*, *slope*, *curvature*) dan efek pemisahan nilai eigen derau terdegenerasi.

## Cara Build, Test, dan Run

```bash
cd projects/go/day004a_term_structure
gofmt -l .
go vet ./...
go test -race -count=1 ./...
go run ./cmd/run
```

## Struktur Kode

- [`cmd/run/main.go`](cmd/run/main.go) — runner utama seluruh eksperimen Bagian 1–4 yang mengekspor [`results.json`](../../../outputs/day004a/results.json) (presisi 17 digit) dan 8 grafik SVG ke [`outputs/figures/`](../../../outputs/figures/).
- [`internal/ns/`](internal/ns/ns.go) — kurva Nelson–Siegel dan Svensson, stabilitas deret Taylor untuk $\tau$ kecil vs evaluasi naif, serta verifikasi integral Gauss–Legendre.
- [`internal/bootstrap/`](internal/bootstrap/bootstrap.go) — pembangkit tingkat kupon par semi-tahunan $0.5\text{–}30$ tahun dan *bootstrap* faktor diskonto berurutan.
- [`internal/interp/`](internal/interp/interp.go) — interpolasi log-linear $P(\tau)$, spline kubik natural atas $y(\tau)$ via solver Thomas, dan studi replikasi derau par $R = 200$.
- [`internal/risk/`](internal/risk/risk.go) — durasi Macaulay/modifikasi/efektif & konveksitas bentuk tertutup vs *bump-and-reprice*, *key-rate duration* (KRD) bobot segitiga, PV01 swap 10y, serta lindung nilai dua obligasi.
- [`internal/afns/`](internal/afns/afns.go) — solusi tertutup $B(\tau)$ dan $-A(\tau)/\tau$ (CDR 2011), pembanding ODE Riccati RK4/Richardson, eksponensial matriks transisi $\mathbb{P}$ segitiga bawah, Van Loan vs Gauss–Legendre untuk $Q_\Delta$, harga risiko afinitas, dan premi jangka waktu 10y.
- [`internal/kalman/`](internal/kalman/kalman.go) — filter Kalman (beserta *fast stack path* bebas alokasi), *smoother* RTS, jangkar eksak Gaussian gabungan $T = 5$, MLE BFGS 11 parameter, studi replikasi $R = 50$, serta uji kesalahan spesifikasi DNS vs AFNS.
- [`internal/pca/`](internal/pca/pca.go) — kovarians teoretis $\operatorname{Cov}(\Delta y_t)$, dekomposisi eigen Jacobi, analisis replikasi $T = 5000, R = 100$, uji *noise bulk*, kemiripan kosinus, dan penghitungan perubahan tanda.
- Utilitas matematika murni stdlib: [`lag`](internal/lag/lag.go) (LU, Cholesky, Thomas, Jacobi eigen, `MatExpSeries`, Gauss–Legendre), [`ode`](internal/ode/ode.go) (RK4 & Richardson), [`optim`](internal/optim/optim.go) (BFGS, beda hingga, gradien Richardson, Hessian), [`rng`](internal/rng/rng.go) (`xoshiro256++` dengan `MainSeed = 20261122`), [`stats`](internal/stats/stats.go), [`parallel`](internal/parallel/parallel.go), [`jsonout`](internal/jsonout/jsonout.go), dan [`svg`](internal/svg/svg.go).

## Ringkasan Hasil Utama

| Bagian | Temuan & Angka Utama (`results.json`, `32.89s` dengan `8` worker) |
|---|---|
| **Bagian 1 — NS/NSS, Bootstrap, Interpolasi & Risiko** | Identitas $y(\tau) = \frac{1}{\tau}\int_0^\tau f(s)ds$ galat maks `6.94e-18` ($< 10^{-12}$); bootstrap par $0.5\text{–}30\text{y}$ pulihkan $P(\tau)$ sampai `2.78e-16`. Pada $R=200$ derau par `0.5 bp`, nol *forward* negatif pada kedua metode (`TV(f)`: `0.300` LogLin vs `0.422` Spline; `yRMSE`: `0.823 bp` LogLin vs `0.688 bp` Spline). Obligasi 10y 4%: $D_{\text{mod}} = 8.1351$ (`relErr = 5.84e-11`), $C_{\text{YTM}} = 78.2980$ (`relErr = 1.12e-7`); jumlah KRD cocok persis dengan $D_{\text{eff}} = 8.2921$ (`err = 0`); swap 10y par `4.3334%` (`PV01 = $811.01/1M`); lindung nilai 5y 4% + 10y 2% menetralkan KRD 5y & 10y hingga `< 5e-16`. |
| **Bagian 2 — AFNS, Riccati & Premi Jangka Waktu** | Pemuatan $-B(\tau)/\tau$ cocok eksak (`0`); penyesuaian $-A(\tau)/\tau$ selalu negatif dan membesar monoton (`-10.09 bps` di 10y, `-59.98 bps` di 30y), cocok dengan RK4/Richardson hingga `3.49e-14` ($A$) dan `1.16e-15` ($-A/\tau$). Matriks $\Phi(\Delta)$ segitiga bawah cocok dengan `MatExpSeries` (`1.11e-16`) dan $Q_\Delta$ Van Loan vs GL (`4.74e-20`); MC $200\text{k}$ transisi 1y `max_z = 2.68`. Premi jangka waktu 10y jangka panjang tertutup `+173.04 bps` vs MC $100\text{k}$ `+172.21 ± 0.38 bps` (`z = 2.19`). |
| **Bagian 3 — Kalman, RTS Smoother & MLE** | Jangkar $T=5$ (3 tenor) cocok dengan Gaussian gabungan eksak: `ll_diff = 2.84e-14`, `mean_diff = 7.27e-16`, `cov_diff = 1.24e-16` ($< 10^{-9}$). Single MLE ($T=1000$) dari tebakan terganggu mencapai `est_ll = 58231.66 > true_ll = 58225.07` dengan verifikasi skor vs Richardson `grad_rel_err_init = 2.32e-10`. Replikasi $R=50$: `50/50` $\ell_{\text{MLE}} \ge \ell_{\text{true}}$; seluruh 5 parameter teridentifikasi kuat ($\lambda, \sigma_{1..3}, \sigma_\varepsilon$) memenuhi `|bias| < 2 SE_rep` (`max 1.74`) dan cakupan 95% `90%–96%`, sedangkan $K^{\mathbb{P}}_{11..33}$ memperlihatkan bias naik Kendall $+\frac{4}{T_{\text{thn}}} \approx +0.21\text{/thn}$. Kesalahan spesifikasi DNS vs AFNS menjatuhkan *log-likelihood* sebesar `767.51` dan membiaskan $\hat{\lambda}$ (`+0.0299`). |
| **Bagian 4 — PCA Perubahan Yield ($T=5000, R=100$)** | Tiga faktor pertama menjelaskan `82.90%` varians teoretis (`69.81%`, `8.65%`, `4.43%`) dengan pola perubahan tanda `[0, 1, 2]`; nilai eigen empiris PC1–PC3 cocok dalam `1.33 SE` (`z = 1.32, 0.72, 0.64`) dan kosinus eigenvektor `0.99995, 0.99923, 0.99748` ($> 0.99$). Rata-rata *bulk* derau PC4–PC10 cocok dengan $2\sigma_\varepsilon^2 = 5.00\times 10^{-7}$ (`z = 0.28`), sementara mode terurut individual terpecah oleh repulsi statistik orde Marchenko–Pastur. |

## Grafik SVG Keluaran

1. [Kurva Nelson–Siegel, Svensson & komponen faktor](../../../outputs/figures/day004a_ns-curves-and-loadings.svg)
2. [Kurva *forward* seketika vs *yield* di bawah derau par 0.5 bp](../../../outputs/figures/day004a_yield-vs-forward-interp.svg)
3. [Profil *Dollar Key-Rate Duration* & lindung nilai dua obligasi](../../../outputs/figures/day004a_key-rate-durations.svg)
4. [Suku penyesuaian bebas-arbitrase AFNS $-A(\tau)/\tau$ dan komponennya](../../../outputs/figures/day004a_afns-yield-adjustment.svg)
5. [Faktor laten AFNS simulasi benar vs *RTS smoother* Kalman](../../../outputs/figures/day004a_kalman-smoothed-factors.svg)
6. [Studi replikasi MLE ($R = 50, T = 1000$): rasio estimasi terhadap nilai benar](../../../outputs/figures/day004a_mle-true-vs-estimated.svg)
7. [*Scree plot* nilai eigen PCA teoretis vs empiris](../../../outputs/figures/day004a_pca-scree-eigenvalues.svg)
8. [Tiga eigenvektor pertama PCA perubahan *yield* (Level, Slope, Kurvatur)](../../../outputs/figures/day004a_pca-eigenvectors.svg)

Analisis matematis lengkap, pemisahan fakta vs hipotesis, serta jejak verifikasi terdokumentasi di [day004a_quant_term-structure.md](../../../docs/project-explanations/day004a_quant_term-structure.md).
