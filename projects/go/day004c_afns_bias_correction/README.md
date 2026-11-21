# day004c — AFNS: Koreksi Bias *Mean-Reversion* dan Inferensi Sampel Terbatas dengan Jackknife serta Bootstrap Parametrik

Modul Go 1.22+ mandiri (`zero external dependencies`) untuk menguji secara empiris apakah metode koreksi bias sampel terbatas—**jackknife temporal dua-paruh** dan **bootstrap parametrik**—mampu mengurangi bias ke atas parameter *mean-reversion* $K^{\mathbb{P}}_{ii}$ serta memperbaiki cakupan interval kepercayaan parameter persisten $\theta^{\mathbb{P}}_1$ pada model ruang-keadaan *Arbitrage-Free Nelson–Siegel* (AFNS) 11 parameter.

## 1. Intuisi Masalah

Pada `day004b`, eksperimen *long-span* dan *infill* menunjukkan bahwa:
- Bias ke atas pada kecepatan *mean reversion* fisik $K^{\mathbb{P}}_{ii}$ meluruh seiring $1/H$ saat rentang kalender $H$ diperpanjang ($10 \to 80$ tahun), tetapi tetap stagnan ($\sim +0{,}20\text{ hingga }+0{,}23\text{ thn}^{-1}$) ketika frekuensi observasi dipadatkan pada $H = 20$ tahun.
- Cakupan interval Wald 95% berbasis invers Hessian untuk rata-rata jangka panjang *level* ($\theta^{\mathbb{P}}_1 = 0{,}045$, waktu paruh $\ln 2 / 0{,}05 \approx 13{,}86$ tahun) jatuh jauh di bawah nominal $95\%$ pada $H = 20$ tahun karena $\hat{K}^{\mathbb{P}}_{11}$ yang terbias ke atas membuat informasi Fisher lokal teramati terlalu optimistis.

Proyek `day004c` menguji apakah koreksi bias non-parametrik berbasis pembagian waktu (*temporal jackknife*) maupun koreksi bias serta interval persentil *bootstrap parametrik* benar-benar memperbaiki akurasi titik dan cakupan pada data AFNS yang sama—tanpa pemangkasan (*clamping*) estimasi terkoreksi dan tanpa memaksa hasil koreksi menjadi lebih baik.

## 2. Spesifikasi Model, Satuan, dan Konfigurasi

Seluruh eksperimen mempertahankan spesifikasi AFNS `day004b`:
- **Tingkat *yield* & $\theta^{\mathbb{P}}$**: desimal per tahun ($\theta^{\mathbb{P}} = (0{,}045,\; -0{,}020,\; 0{,}010)$ sesuai $4{,}5\%$, $-2{,}0\%$, $1{,}0\%$).
- **Parameter peluruhan $\mathbb{Q}$ & *mean-reversion* $\mathbb{P}$**: $\lambda = 0{,}6\text{ thn}^{-1}$, $K^{\mathbb{P}} = \operatorname{diag}(0{,}05,\; 0{,}5,\; 0{,}4)\text{ thn}^{-1}$.
- **Volatilitas difusi faktor**: $\Sigma = \operatorname{diag}(0{,}0060,\; 0{,}0100,\; 0{,}0200)\text{ thn}^{-1/2}$.
- **Derau observasi**: $\sigma_\varepsilon = 0{,}0005$ ($5\text{ bp}$), $\varepsilon_{t,j} \stackrel{\text{iid}}{\sim} \mathcal{N}(0, \sigma_\varepsilon^2)$ independen antar-tenor dan antar-waktu ($\Delta = 1/52\text{ thn}$).
- **Tenor observasi (10 tenor)**: $\tau \in \{0{,}25,\; 0{,}5,\; 1,\; 2,\; 3,\; 5,\; 7,\; 10,\; 20,\; 30\}\text{ tahun}$.
- **Keadaan awal & transisi**: $X_0 \sim \mathcal{N}(\theta^{\mathbb{P}}, Q_\infty)$ dari distribusi stasioner $\mathbb{P}$, bertransisi eksak via $e^{-K^{\mathbb{P}}\Delta}$ dan $Q_\Delta$.
- **Seed & paralelisme deterministik**: `MainSeed = 2026112203` (`xoshiro256++` + `SplitMix64`, deterministik per `(salt, rep_idx)` atau `(salt, outer_idx, boot_idx)`, reduksi terurut; bit-identik `1` vs `8` worker).

### Desain Tiga Studi Berpasangan
1. **Jackknife temporal 20 tahun mingguan (`longspan_20y_weekly`)**: $T = 1{,}040$ (setengah sampel $520$), $R = 50$ menggunakan prefiks awal lintasan 80 tahun.
2. **Jackknife temporal 80 tahun mingguan (`longspan_80y_weekly`)**: $T = 4{,}160$ (setengah sampel $2{,}080$), $R = 50$ pada lintasan 80 tahun penuh yang berbagi prefiks dengan studi 20 tahun.
3. **Bootstrap parametrik 20 tahun mingguan (`bootstrap_study_20y_weekly`)**: $R = 50$ dataset luar (fit `Full` dari `20y`) $\times B = 100$ simulasi & refit penuh 11 parameter ($5{,}000$ refit bootstrap dalam). Inisialisasi refit memakai fit luar yang diganggu deterministik kecil di interior kotak batas, dengan ambang validitas pra-run $\ge 95/100$ refit *finite* dan seluruh $R = 50$ replikasi luar tetap masuk denominator cakupan.

## 3. Metodologi Koreksi dan Asumsi

### A. Jackknife Temporal Dua-Paruh
Setiap lintasan waktu penuh $Y_{0:T-1}$ dibagi menjadi dua segmen waktu berurutan yang sama panjang: paruh pertama $Y_{0:T/2-1}$ dan paruh kedua $Y_{T/2:T-1}$ (keduanya diinisialisasi dari distribusi stasioner model yang sama). Estimator jackknife dihitung untuk `kp11`, `kp22`, `kp33`, dan `theta1` (serta dicatat untuk seluruh 11 parameter) tanpa pemangkasan ke batas optimizer:

$$\hat{\psi}_{\text{JK}} = 2\,\hat{\psi}_{\text{full}} - \frac{1}{2}\left(\hat{\psi}_{\text{half1}} + \hat{\psi}_{\text{half2}}\right)$$

**Asumsi & Keterbatasan**: Rumus ini menghapus tepat suku bias orde pertama $c/H$ *hanya jika* $\mathbb{E}[\hat{\psi}_H] = \psi + c/H + o(1/H)$ dan paruh sampel mengikuti struktur bias yang sama. Karena dua paruh berasal dari satu lintasan proses stasioner yang sama, $\hat{\psi}_{\text{half1}}$ dan $\hat{\psi}_{\text{half2}}$ **saling bergantung secara temporal** (terutama untuk faktor *level* persisten $\kappa_1 = 0{,}05\text{ thn}^{-1}$); oleh karena itu proyek ini tidak menggunakan rumus selang kepercayaan maupun simpangan baku dua-sampel i.i.d., dan mengevaluasi ketidakpastian Monte Carlo dari replikasi luar independen serta selisih berpasangan $\hat{\psi}_{\text{JK},r} - \hat{\psi}_{\text{raw},r}$.

### B. Bootstrap Parametrik
Pada tiap dataset luar $r \in \{0,\dots,49\}$ pada $H=20\text{y}$, model AFNS disimulasikan ulang $B=100$ kali di bawah $\hat{\psi}_{\text{outer},r}$ (keadaan awal stasioner $Q_\infty(\hat{\psi}_r)$ + transisi eksak + derau i.i.d. $\hat{\sigma}_{\varepsilon,r}$) lalu di-refit penuh 11 parameter:
- Bias bootstrap kondisional: $\widehat{\text{bias}}_{\text{boot},r} = \frac{1}{B_{\text{fin}}} \sum_{b=1}^{B_{\text{fin}}} (\hat{\psi}^*_{r,b} - \hat{\psi}_{\text{outer},r})$.
- Estimasi terkoreksi bootstrap (tanpa pemangkasan): $\hat{\psi}_{\text{boot-corr},r} = \hat{\psi}_{\text{outer},r} - \widehat{\text{bias}}_{\text{boot},r}$.
- Interval persentil dua sisi 95%: $[q_{0{,}025}(\hat{\psi}^*_r),\; q_{0{,}975}(\hat{\psi}^*_r)]$ dibandingkan dengan interval Wald Hessian $\hat{\psi}_{\text{outer},r} \pm z_{0{,}975}\,\widehat{\text{SE}}_{\text{Hess},r}$ dan dievaluasi dengan interval skor binomial Wilson 95% atas seluruh $R = 50$ replikasi luar.

## 4. Ringkasan Hasil Utama (`outputs/day004c/results.json`, `180,13s` pada `8` worker)

### Verifikasi Teknis Pra-Spesifikasi (`c1..c5 = tercapai`)
- **Jangkar Gaussian gabungan ($T=5$, 3 tenor)**: selisih log-likelihood Kalman penuh `5,68e-14`, `FastLogLikelihood` Riccati `5,68e-14`, rata-rata smoother RTS `5,55e-16`, kovarians smoother `6,39e-17` (target `< 1e-9`).
- **Gradien FD vs Richardson**: galat relatif L2 maksimum sepanjang studi di titik awal `2,42e-9` dan di nilai benar `2,75e-7` (target `< 1e-5`).
- **Uji aljabar jackknife $c/H$**: galat `1,39e-17` di $H=20$ dan $H=80$ (verifikasi aljabar murni di bawah asumsi $c/H$, bukan bukti berlakunya asumsi pada MLE AFNS).
- **Determinisme worker**: identik bit-demi-bit antara `1` dan `8` worker pada fixture deterministik (`T=80, R=2, B=4`).
- **Kelengkapan & waktu**: `100%` replikasi luar (`50/50`) dan `5.000/5.000` refit bootstrap bernilai *finite*; total waktu `180,13 detik` (`< 600s`): anchor/determinisme `0,90s`, pembangkitan jalur `0,01s`, jackknife 20y `6,39s`, jackknife 80y `21,01s`, bootstrap 20y `151,81s`, SVG `0,003s`.

### Perbandingan MLE Mentah vs Jackknife Temporal (`H = 20y` dan `H = 80y`, `R = 50`)

| Horizon & Param (True) | Bias Raw $\pm$ MCSE | Bias JK $\pm$ MCSE | SD Raw $\to$ SD JK | RMSE Raw $\to$ JK (Rasio) | JK OOB / $\le 0$ | Rasio Hessian SE (Mean / Med) |
|---|---:|---:|---:|---:|---:|---:|
| `20y: kp11 (0.05)` | `+0,3170 ± 0,0415` | `+0,0994 ± 0,0486` | `0,2933 -> 0,3439` | `0,4299 -> 0,3546` (`0,825`) | `20/50 (40%)` | `0,581 / 0,560` |
| `20y: kp22 (0.50)` | `+0,2130 ± 0,0491` | `+0,0463 ± 0,0591` | `0,3475 -> 0,4181` | `0,4046 -> 0,4165` (`1,029`) | `3/50 (6%)` | `0,740 / 0,714` |
| `20y: kp33 (0.40)` | `+0,1821 ± 0,0351` | `-0,0715 ± 0,0523` | `0,2483 -> 0,3698` | `0,3059 -> 0,3730` (`1,219`) | `5/50 (10%)` | `0,927 / 0,949` |
| `20y: theta1 (0.045)` | `-0,0013 ± 0,0025` | `-0,0009 ± 0,0025` | `0,0179 -> 0,0174` | `0,0177 -> 0,0172` (`0,971`) | `0/50 (0%)` | `0,279 / 0,217` |
| `80y: kp11 (0.05)` | `+0,0604 ± 0,0092` | `-0,0063 ± 0,0134` | `0,0650 -> 0,0948` | `0,0882 -> 0,0941` (`1,066`) | `16/50 (32%)` | `0,730 / 0,720` |
| `80y: kp22 (0.50)` | `+0,0379 ± 0,0202` | `+0,0043 ± 0,0212` | `0,1431 -> 0,1497` | `0,1467 -> 0,1482` (`1,011`) | `0/50 (0%)` | `0,797 / 0,794` |
| `80y: kp33 (0.40)` | `+0,0563 ± 0,0167` | `+0,0105 ± 0,0184` | `0,1178 -> 0,1302` | `0,1295 -> 0,1293` (`0,998`) | `0/50 (0%)` | `0,894 / 0,900` |
| `80y: theta1 (0.045)` | `+0,0000 ± 0,0018` | `+0,0107 ± 0,0109` | `0,0124 -> 0,0773` | `0,0123 -> 0,0772` (`6,282`) | `1/50 (2% OOB)` | `0,606 / 0,498` |

### Bootstrap Parametrik 20y (`R = 50, B = 100`) dan Cakupan Interval 95%

| Parameter (True) | Bias Raw | Bias Boot-Corr | RMSE Raw $\to$ Boot (Rasio) | Boot OOB / $\le 0$ | Wald Cov 95% [Wilson 95%] | Boot Persentil Cov 95% [Wilson 95%] |
|---|---:|---:|---:|---:|---|---|
| `kp11 (0.05)` | `+0,3170` | `+0,1014 ± 0,0429` | `0,4299 -> 0,3172` (`0,738`) | `19/50 (38%)` | `76,0%` `[62,6%, 85,7%]` | `2,0%` `[0,35%, 10,5%]` |
| `kp22 (0.50)` | `+0,2130` | `+0,0034 ± 0,0492` | `0,4046 -> 0,3446` (`0,852`) | `0/50 (0%)` | `90,0%` `[78,6%, 95,7%]` | `72,0%` `[58,3%, 82,5%]` |
| `kp33 (0.40)` | `+0,1821` | `-0,0282 ± 0,0351` | `0,3059 -> 0,2476` (`0,810`) | `0/50 (0%)` | `94,0%` `[83,8%, 97,9%]` | `62,0%` `[48,2%, 74,1%]` |
| `theta1 (0.045)` | `-0,0013` | `-0,0013 ± 0,0025` | `0,0177 -> 0,0178` (`1,002`) | `0/50 (0%)` | `46,0%` `[33,0%, 59,6%]` | `42,0%` `[29,4%, 55,8%]` |

**Temuan Kunci & Keterbatasan**:
- Koreksi bootstrap parametrik memangkas RMSE ketiga $K^{\mathbb{P}}_{ii}$ pada $H=20\text{y}$ sebesar $15\%\text{–}26\%$ tanpa pembengkakan SD yang terjadi pada jackknife (`0,3439` vs `0,3036` untuk `kp11`). Namun baik jackknife (`40%` di 20y, `32%` di 80y) maupun bootstrap koreksi aditif (`38%` di 20y) sering memproduksi estimasi $\hat{\kappa}_1 \le 0$ yang secara fisik non-stasioner saat nilai luar awal dekat dengan nilai benar $0{,}05$.
- Pada `80y`, satu kegagalan batas di `half1` & `half2` pada `Rep 30` (menyentuh $\lambda=5$, $\theta_1=-0{,}5$, $\sigma_\varepsilon=0{,}05$) merusak nilai jackknife untuk `theta1` (`RMSE ratio = 6,28`), `lambda` (`263,18`), dan `sigma_eps` (`2085,21`), menunjukkan kerentanan jackknife tanpa penyaring kestabilan paruh-sampel.
- **Interval persentil bootstrap naif gagal berat pada parameter berbias kuat**: cakupan `kp11` jatuh ke `2,0%` (`1/50`) karena bootstrap mensimulasikan di bawah $\hat{\kappa}_1$ yang sudah terbias ke atas dan mengalami bias ke atas *kedua kali* (`+0,216` lagi), menggeser seluruh kuantil $[2{,}5\%, 97{,}5\%]$ menjauhi `0,05`. Sementara itu cakupan bootstrap `theta1` tetap rendah (`42,0%`) karena data bootstrap dibangkitkan dengan *mean reversion* $\hat{\kappa}_1$ yang terlalu cepat sehingga lebar selang persentil menyempit (`mean width = 0,0183`) seperti halnya SE Hessian.
- **Diagnostik Optimizer**: Pada segmen penuh (`20y`, `80y`) dan seluruh `5.000/5.000` refit bootstrap, $100\%$ optimisasi mencapai `scaled_grad_norm <= 1e-6` dengan penghentian `line_search_failed` akibat lantai pembulatan `float64` fungsi log-likelihood bernilai besar ($|\ell| \sim 6\times 10^4\text{–}2{,}4\times 10^5$), sehingga bendera konvergensi gradien tak-terskala (`gradient_tolerance`) bernilai `0%`. Pada segmen paruh-sampel, tercatat batas aktif di `20y half1 Rep 28`, `80y half1 Rep 30`, dan `80y half2 Rep 30`, serta `max_iterations` (`150` iterasi tahap-2) pada `80y half2 Rep 23 & 28`.

## 5. Cara Build, Run, dan Test

Jalankan dari direktori proyek:

```bash
cd projects/go/day004c_afns_bias_correction
gofmt -l .
go vet ./...
go test -race -count=1 -v ./...  # 11 paket, 17 tes utama
go run ./cmd/run               # Runner penuh R=50, B=100 (~180s pada 8 worker)
go run ./cmd/run -quick        # Mode cepat pengembangan (R=3, B=6 ke temp dir)
```

## 6. Struktur Kode

- `internal/rng/` — Generator deterministik `xoshiro256++` + `SplitMix64` (`MainSeed = 2026112203`, `DeriveSeed`, `DeriveBootstrapSeed`).
- `internal/parallel/` — `DefaultWorkers` (`min(8, NumCPU)`) dan `RunChunks` dengan reduksi terurut.
- `internal/afns/` — Pemuat kurva `G1`/`G2`, penyesuaian konveksitas CDR 2011 (`YieldAdjustmentComponents`), matriks transisi eksak $\Phi(\Delta)$, $Q_\Delta$, $Q_\infty$, dan simulator stasioner `SimulateStationaryPath`.
- `internal/optim/` — Gradien beda hingga pusat, verifikasi Richardson, Hessian pusat, serta `MinimizeBFGS` terkendala kotak dengan klasifikasi 5 alasan penghentian eksplisit.
- `internal/kalman/` — `FastLogLikelihood` (*steady-state* Riccati), `FilterAndSmooth` (Joseph + *RTS smoother*), `VerifyJointGaussianAnchor`, aljabar matriks/Jacobi mandiri (`linalg.go`), dan `EstimateMLEWithOptions`.
- `internal/jackknife/` — `CombineSplit`, `CombineVec` (tanpa pemangkasan + penanda batas & fisik), pembangkit jalur prefiks-bersama, dan `RunJackknifeHorizonStudy`.
- `internal/bootstrap/` — Ambang validitas pra-run (`MinValidFraction = 0.95`), `PerturbBootstrapGuess`, `EvalParamBootstrap`, dan `RunParametricBootstrapStudy` (`R x B`).
- `internal/stats/`, `internal/jsonout/`, `internal/svg/` — Statistik (kuantil linier, `WilsonCI`, `CoverageWithValidity`), penulis JSON presisi 17 digit, dan generator SVG buatan sendiri.
- `cmd/run/` — Runner CLI utama, verifikasi jangkar & determinisme worker, ekspor JSON serta 4 grafik SVG.

## 7. Grafik & Artefak Terkait

- JSON keluaran: [`outputs/day004c/results.json`](../../../outputs/day004c/results.json)
- Grafik SVG:
  - [`day004c_jackknife-vs-mle-bias-rmse.svg`](../../../outputs/figures/day004c_jackknife-vs-mle-bias-rmse.svg)
  - [`day004c_bootstrap-vs-jackknife-bias-20y.svg`](../../../outputs/figures/day004c_bootstrap-vs-jackknife-bias-20y.svg)
  - [`day004c_wald-vs-bootstrap-coverage.svg`](../../../outputs/figures/day004c_wald-vs-bootstrap-coverage.svg)
  - [`day004c_bootstrap-distributions-rep0.svg`](../../../outputs/figures/day004c_bootstrap-distributions-rep0.svg)
- Penjelasan proyek 11 bagian: [`day004c_quant_afns-bias-correction.md`](../../../docs/project-explanations/day004c_quant_afns-bias-correction.md)
