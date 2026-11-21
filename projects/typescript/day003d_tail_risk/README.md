# day003d — Risiko Ekor: GARCH–EVT, VaR/ES & Backtesting (TypeScript)

Implementasi kuantitatif risiko ekor terintegrasi dalam **TypeScript** (*package* ESM `"strict": true`, nol dependensi *runtime* eksternal) yang memadukan model volatilitas kondisional *heavy-tail* GARCH(1,1)-$t$, Teori Nilai Ekstrem (*Extreme Value Theory*: POT GPD & estimator Hill), pemodelan VaR dan *Expected Shortfall* (ES) satu hari ke depan, serta validasi tiga uji *backtesting* statistik.

## Intuisi Model

1. **GARCH(1,1)-$t$**: Imbal hasil keuangan menampilkan pengelompokan volatilitas ($\sigma_t^2 = \omega + \alpha \varepsilon_{t-1}^2 + \beta \sigma_{t-1}^2$) sekaligus ekor inovasi lebih tebal dari distribusi Gaussian ($z_t \sim t_{\text{std}}(\nu)$ berskala varians $1$). Mengabaikan $\nu$ berhingga ($\nu=5$) membuat model Gaussian meremehkan kerugian ekor ekstrem ($\text{VaR}_{0.999}$ bias `-32.4%`).
2. **Teori Nilai Ekstrem (Pickands–Balkema–de Haan)**: Distribusi kelebihan di atas ambang tinggi $u$ ($Y = X - u \mid X > u$) dari hampir semua distribusi berekor berat konvergen ke *Generalized Pareto Distribution* (GPD) dengan parameter bentuk $\xi \to 1/\nu$ dan skala $\beta > 0$. Estimator semi-parametrik FHS-EVT mengekstrapolasi ekor ekstrem jauh lebih mulus dibanding kuantil historis mentah.
3. **VaR & Expected Shortfall (Konvensi Kerugian Positif)**: Pada kerugian harian $L_t = -r_t$, ambang satu hari ke depan bernilai positif: $\text{VaR}_p(t) = -\mu + \sigma_t q_p$ dan $\text{ES}_p(t) = -\mu + \sigma_t \text{ES}_p(Z)$.
4. **Backtesting (Kupiec, Christoffersen, Acerbi–Székely $Z_2$)**: Uji Kupiec mengevaluasi frekuensi pelanggaran tanpa syarat; uji Christoffersen menambahkan ketergantungan Markov orde-1 untuk menangkap klaster pelanggaran; uji $Z_2 = 1 - \frac{1}{T\alpha}\sum \frac{L_t \mathbf{1}\{L_t > \text{VaR}_t\}}{\text{ES}_t}$ menguji kalibrasi rata-rata ekor kerugian ES terhadap nol.

## Cara Install, Test & Run

```bash
cd projects/typescript/day003d_tail_risk
npm ci
npm run typecheck
npm run lint
npm run format:check
npm run build
npm test
npm start
```

Keluaran tersimpan ke:
- [outputs/day003d/results.json](../../../outputs/day003d/results.json)
- 7 grafik vektor di `outputs/figures/day003d_*.svg` (`t-vs-gpd-tail`, `xi-and-hill-vs-q`, `garch-mle-dist`, `var-bias-rmse`, `residuals-qq`, `regime-change-violations`, `backtest-size-power`).

## Struktur Kode

- `src/rng.ts` — PRNG 64-bit `Xoshiro256PlusPlus` (aritmetika limb 32-bit, ekivalen bit-demi-bit dengan spesifikasi `BigInt`), `SplitMix64`, Box–Muller, Marsaglia–Tsang Gamma, dan Student-$t$ standar.
- `src/special.ts` — fungsi `logGamma`, `digamma`, `trigamma`, `regularizedIncompleteBeta`, `erfc` Cody, `normalQuantile` Acklam–Halley, `chi2Sf`, dan kuadratur Simpson adaptif.
- `src/distributions.ts` — PDF/CDF/SF/kuantil Student-$t$ (baku dan lokasi-skala), rumus tertutup ES Student-$t$ baku, serta PDF/CDF/kuantil/VaR/ES GPD.
- `src/optim.ts` — optimisasi BFGS (pencarian garis *Strong Wolfe*), Nelder–Mead, inversi matriks Gauss–Jordan, dan beda hingga terpusat.
- `src/garch.ts` — simulasi GARCH(1,1)-$t$, rekonstruksi deterministik $\sigma_t^2$, log-likelihood & skor analitik 5-parameter, reparameterisasi halus tanpa kendala aktif, SE *Information* & *Sandwich*, serta QMLE Gaussian.
- `src/evt.ts` — MLE GPD (gradien analitik dengan *guard* kelayakan dan ekspansi $\xi \to 0$), estimator Hill, dan *pipeline* POT.
- `src/risk.ts` — empat estimator VaR/ES (`gaussian`, `t_parametric`, `fhs`, `fhs_evt`) serta uji sensitivitas ambang $q \in \{0.85, 0.90, 0.95\}$.
- `src/backtest.ts` — Kupiec POF (binomial eksak + LR asimtotik), Christoffersen (Ind & CC), Acerbi–Székely $Z_2$, *rolling backtest* ($W=1000$, estimasi ulang tiap 20 hari), dan simulasi perubahan rezim hari ke-3000.
- `src/core_tasks.ts`, `src/worker.ts`, `src/parallel.ts` — pembagian *chunk* `worker_threads` deterministik (hasil identik untuk 1 maupun banyak *worker*).
- `src/stats.ts`, `src/json.ts`, `src/svg.ts`, `src/main.ts` — statistik ringkasan, penulis JSON 17 digit signifikan, pembuat 7 grafik SVG, dan orkestrator penuh (`seed = 20261121`).
- `tests/*.test.ts` — 14 unit test Vitest di 3 berkas tes.

## Ringkasan Hasil (`seed = 20261121`, 10 Worker, `18.44 detik`)

### Bagian 1 — Distribusi Eksak & Pickands–Balkema–de Haan ($N = 2.000.000$)
- **Student-$t$ CDF vs kuadratur adaptif**: galat maks `2.33e-15`; *roundtrip* kuantil `1.11e-16`.
- **ES Student-$t$ baku vs kuadratur** ($\nu \in \{3.5, 5, 10, 30\}$, $p \in \{0.95, 0.975, 0.99, 0.999\}$): galat maks `3.38e-14`.
- **GPD ES numerik** ($\xi \in \{-0.1, 0, 0.2, 0.4\}$): galat maks `1.78e-15`.
- **PBdH pada $t_{\text{std}}(5)$ ($1/\nu = 0.20$)**: $\hat{\xi}_{\text{POT}}$ naik dari `0.1031` ($q=0.90$) $\to$ `0.1310` ($0.95$) $\to$ `0.1613` ($0.99$) $\to$ `0.1785` ($0.999$); sementara Hill turun dari `0.3870` ($0.90$) $\to$ `0.2147` ($0.999$).

### Bagian 2 — GARCH(1,1)-$t$ ($R = 200, n = 2500$)
- **Rekonstruksi $\sigma_t^2$**: selisih maks `2.71e-20`; skor analitik vs beda hingga pusat: `2.02e-6` (`< 1e-5`); selisih NLL BFGS vs Nelder–Mead: `1.29e-14`.
- **Momen & Batas Momen ke-4**: varians tak bersyarat empiris $N=200\text{k}$ `0.9817e-4 ± 0.0376e-4` vs `1.0000e-4`; kurtosis berlebih $z_t$ (`R=200`): `4.96 ± 0.56` vs `6`. Secara aljabar, $D = 1 - 9(0.08)^2 - 2(0.08)(0.90) - 0.90^2 = -0.0116 < 0 \implies E[\varepsilon_t^4] = \infty$, sehingga rumus ACF $\varepsilon_t^2$ standar tidak berlaku untuk $\nu=5$ ini (anchor $t(12), \alpha=0.10, \beta=0.80$ dengan $D=0.1625>0$ terverifikasi cocok dalam `2.14 SE`).
- **Studi Replikasi $R=200$** (`200/200` konvergen, `0` batas aktif): $\hat{\alpha} = 0.0800 \pm 0.0135$ (cakupan Sand `91.5%`, Info `93.0%`), $\hat{\beta} = 0.8983 \pm 0.0160$ (`93.5%` / `94.0%`), $\hat{\nu} = 5.166 \pm 0.559$ (`95.5%` / `95.0%`); rasio efisiensi MSE QMLE/$t$-MLE `1.53x` ($\mu$), `2.19x` ($\omega$), `2.92x` ($\alpha$), dan `2.67x` ($\beta$).

### Bagian 3 & 4 — VaR/ES Bersyarat ($R=200$) & Backtesting ($40 \times 6000$ hari)

| Metode / Uji | Hasil Utama |
|---|---|
| VaR `99.9%` (Gaussian vs t-Param vs FHS vs FHS-EVT) | Bias rel: `-32.41%` vs `-0.78%` vs `-4.40%` vs `-3.05%` \| RMSE rel: `0.3251` vs `0.0635` vs `0.1304` vs `0.1024` |
| Kupiec $H_0$ (`100k`, $T=1000$) | $p=0.95$: `4.933%` (eksak) / `5.092%` (asymp); $p=0.99$: `3.593%` (eksak, ukuran diskrit teoretis `3.646%`) / `5.414%` (asymp) |
| Christoffersen ($20\text{k}$, $T=1000$) | Daya klaster Markov $P(1\|1)=0.30$: `99.7%` CC ($p=0.95$) & `81.2%` CC ($p=0.99$) |
| Acerbi–Székely $Z_2$ ($10\text{k}$, $T=1000$) | Ukuran $H_0$: `4.95%` ($97.5\%$) / `4.98%` ($99\%$); daya vs Gaussian `64.67%` ($99\%$) |
| *Rolling Backtest* ($R=40$, 5000 OOS) | TrueCtrl rate99 `1.009% ± 0.025%`; Gauss-GARCH `1.600%` (ditolak `97.5%` Kupiec, `100%` $Z_2$); $t$-GARCH `1.088%`; FHS `1.146%`; FHS-EVT `1.069%` (`0%` penolakan Kupiec) |
| Perubahan Rezim ($t=3000, \alpha: 0.08\to 0.15$) | Tetap rate99 `0.862%` (ditolak akhir `27.5%`); Bergulir adaptif $\hat{\alpha}_t: 0.072 \to 0.141$, rate99 pulih `1.042%` (`0%` penolakan akhir) |
