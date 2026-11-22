# Hari 004i — Cointegration dan Trading OU Optimal: Dickey–Fuller, Engle–Granger, Johansen, dan Ambang Bertram (Rust)

Proyek Rust mandiri (*pure standard library*, `edition = "2021"`, nol *crate* eksternal) yang mempelajari **cointegration** dan **trading Ornstein–Uhlenbeck (OU) optimal** pada data sintetis: estimasi MLE eksak vs bias Kendall/Marriott–Pope, waktu lintasan pertama (*first-passage time*) dan ambang optimal Bertram ($a^*$), distribusi simulasi Dickey–Fuller ($R = 10^6$), Engle–Granger, regresi palsu, uji trace Johansen bivariat, serta simulasi perdagangan bergulir dengan filter uji, *structural break*, dan sensitivitas biaya.

> **Catatan**: Seluruh data pada proyek ini bersifat **sintetis** (dihasilkan oleh PRNG deterministik `xoshiro256++`) dan tidak memuat klaim empiris tentang pasar nyata.

Dokumentasi penjelasan proyek lengkap tersedia di [day004i_quant_pairs-trading.md](../../../docs/project-explanations/day004i_quant_pairs-trading.md).

## Intuisi Ekonomi & Matematis

Dua aset berkointegrasi $x_t$ dan $y_t = 0{,}1 + 1{,}2\,x_t + X_t$ membentuk *spread* stasioner $X_t$ mengikuti proses OU $dX_t = \theta(\mu - X_t)\,dt + \sigma\,dW_t$ ($\theta = 16$, half-life $\approx 11$ hari, $\sigma = 0{,}20$, $\mathrm{SD}_\infty = \sigma/\sqrt{2\theta} \approx 0{,}0354$).
- **Bias estimasi $\hat{\theta}$**: Pada sampel terbatas berdurasi $T$ tahun, MLE untuk $\theta$ dengan $\mu$ tak diketahui memiliki bias positif orde $4/T$ (Kendall / Marriott–Pope), yang memengaruhi ambang masuk/keluar $\hat{a}^*$.
- **Ambang Bertram optimal**: Memasang batas simetris $a < 0$ dan $m = -a$ menukar laba per siklus $(-2a - c)$ terhadap waktu tunggu rata-rata $2\,\mathbb{E}[\tau(a \to -a)]$ yang tumbuh eksponensial dalam $a^2$, menghasilkan titik maksimum tunggal $a^*$ untuk laju laba tahunan $R(a)$.
- **Bahaya regresi palsu (*spurious*)**: Dua *random walk* independen menghasilkan $|t| > 1{,}96$ pada `93,26%` regresi OLS biasa, sehingga diperlukan uji akar unit residual (Engle–Granger) atau rank VECM (Johansen) dengan nilai kritis yang tepat.

## Rumus Kunci

- **Transisi eksak & MLE OU ($\Delta = 1/252$)**:
  $$X_{t+\Delta} \mid X_t \sim \mathcal{N}\!\left(\mu + (X_t - \mu)e^{-\theta\Delta},\; \frac{\sigma^2(1 - e^{-2\theta\Delta})}{2\theta}\right), \qquad \hat{\theta} = -\frac{\ln\hat{\phi}}{\Delta}.$$
- **Waktu lintasan pertama (diturunkan dari BVP Dynkin $\frac{1}{2}\sigma^2 u'' - \theta x u' = -1$)**:
  $$\mathbb{E}[\tau(a \to m)] = \frac{\pi}{2\theta}\bigl[\operatorname{erfi}(z_m) - \operatorname{erfi}(z_a)\bigr] + \frac{\sqrt{\pi}}{\theta}\int_{z_a}^{z_m} e^{z^2}\operatorname{erf}(z)\,dz, \qquad z = \frac{x\sqrt{\theta}}{\sigma},$$
  yang tereduksi menjadi $\frac{\pi}{\theta}\operatorname{erfi}(m\sqrt{\theta}/\sigma)$ untuk batas simetris $a = -m$.
- **Laju laba Bertram**: $R(a) = \dfrac{-2a - c}{2\,\mathbb{E}[\tau(a \to -a)]}$, dimaksimalkan saat $2\operatorname{erfi}(z) = \left(2z - \frac{c\sqrt{\theta}}{\sigma}\right)\operatorname{erfi}'(z)$.
- **Johansen trace bivariat**: $\mathrm{LR}_{\mathrm{trace}}(0) = -N\sum_{i=1}^2 \ln(1 - \lambda_i)$ dari masalah eigen $S_{01} S_{11}^{-1} S_{10} v = \lambda S_{00} v$.

## Struktur Kode

```
projects/rust/day004i_pairs_trading/
├── Cargo.toml
├── meta.json
├── README.md
├── src/
│   ├── lib.rs          # Re-export modul
│   ├── rng.rs          # xoshiro256++ + Box-Muller + SplitMix64 (deterministik per indeks path)
│   ├── parallel.rs     # Worker pool std::thread::scope deterministik (reduksi blok tetap bit-identik)
│   ├── stats.rs        # Rata-rata, SD, SE, skewness, kuantil+SE, OLS, R^2, Durbin-Watson
│   ├── linalg.rs       # Thomas tridiagonal, Cholesky 2x2, invers 2x2, LU, Jacobi simetris
│   ├── special.rs      # erfi(x) (deret Taylor + Simpson) dan turunannya
│   ├── ou.rs           # Simulasi OU eksak, exact vs Euler loglik, MLE tertutup vs numerik, R=20k replikasi
│   ├── fpt.rs          # E[tau] tertutup (simetris & asimetris), BVP beda hingga E[tau] & Var(tau), MC fine-step
│   ├── strategy.rs     # Ambang optimal Bertram a* (langsung vs FOC), Sharpe, verifikasi MC 200 tahun
│   ├── unitroot.rs     # Dickey-Fuller (R=10^6) & batas Brownian, Engle-Granger (R=20k), spurious, power
│   ├── johansen.rs     # Uji trace Johansen bivariat (R=20k, T=500), size vs power
│   ├── backtest.rs     # Backtest bergulir W in {1,2,5}, filter EG, structural break + stop-loss, biaya
│   ├── jsonout.rs      # Penulis JSON 17 digit presisi tanpa dependensi eksternal
│   ├── svg.rs          # Generator 7 grafik SVG
│   └── bin/run.rs      # Runner utama yang menghasilkan outputs/day004i/results.json dan SVG
└── tests/
    └── integration_tests.rs
```

## Cara Build, Test, dan Run

Jalankan dari direktori proyek `projects/rust/day004i_pairs_trading`:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --release
cargo run --release --bin run
```

Keluaran utama ditulis ke [../../../outputs/day004i/results.json](../../../outputs/day004i/results.json) dan 7 grafik ke [../../../outputs/figures/](../../../outputs/figures/).

## Ringkasan Hasil Eksperimen (`main_seed = 20261122`, `8` worker, `3,79s`)

| Bagian / Besaran | Nilai Terukur | Acuan / Pembanding | Status Kriteria |
| --- | --- | --- | --- |
| MLE tertutup vs numerik & ACF(1) | `max_diff = 8,17e-9`; `acf1 = 0,9370` | Eksak `0,9385`; `rel_score = 3,02e-9` | `c1 = tercapai` |
| Bias $\hat{\theta}$ vs Kendall ($T=1,2,5,10$) | `4,44 / 2,13 / 0,84 / 0,44` | Rasio Kendall `1,04..1,09` | — |
| FPT tertutup vs BVP vs MC (`-SD -> +SD`) | `0,18721` / `0,18721` / `0,18672` | `z_MC = 1,36` (asimetris `0,43`) | `c2 = tercapai` |
| Ambang optimal $a^*$ & MC 200 tahun | `a* = -0,01574` (`0,445 SD`), `R(a*) = 0,2044` | `z_trades = 0,64`, `z_profit = 0,64` | `c3 = tercapai` |
| Dickey–Fuller ($R=10^6$) & EG ($R=20k$) | `DF q05 = -2,864`; `EG q05 = -3,353` | Size OOS EG `4,77%`; Spurious `93,26%` | `c4 = tercapai` |
| Johansen trace ($T=500, R=20k$) | `q95_drift = 18,000`, size OOS `4,87%` | Power `2yr = 81,56%` (`93,46%` di `15,41`, size `10,99%`) | `c5 = belum tercapai` |
| Penalti estimasi bergulir ($W=1,2,5$ yr) | `0,0243` → `0,0134` → `0,0060` | Oracle `0,3920 ± 0,0009/yr` | `c6 = tercapai` |
| Filter EG pada pasangan *spurious* ($W=2$) | Laba `-0,00061 ± 0,00038` | Tanpa filter `-0,01365` | `c7 = tercapai` |
| Anggaran waktu runner | `9,05s` (`8` worker) | `< 300s` | `c8 = tercapai` |

## Daftar Grafik

- [day004i_spread-and-thresholds.svg](../../../outputs/figures/day004i_spread-and-thresholds.svg)
- [day004i_ou-theta-bias-vs-t.svg](../../../outputs/figures/day004i_ou-theta-bias-vs-t.svg)
- [day004i_bertram-return-rate.svg](../../../outputs/figures/day004i_bertram-return-rate.svg)
- [day004i_eg-power-vs-theta.svg](../../../outputs/figures/day004i_eg-power-vs-theta.svg)
- [day004i_rolling-profit-vs-window.svg](../../../outputs/figures/day004i_rolling-profit-vs-window.svg)
- [day004i_cost-sensitivity.svg](../../../outputs/figures/day004i_cost-sensitivity.svg)
- [day004i_fpt-comparison.svg](../../../outputs/figures/day004i_fpt-comparison.svg)
