# Hari 004h — Volatilitas Kasar: Model Rough Bergomi, Skema Hybrid, Skew ATM, dan Identifiabilitas Hurst (Rust)

Proyek Rust mandiri (*pure standard library*, `edition = "2021"`, nol *crate* eksternal) yang mengimplementasikan **model rough Bergomi (Bayer–Friz–Gatheral, 2016)** untuk volatilitas kasar dengan parameter Hurst $H = 0{,}1$: proses Volterra Riemann–Liouville, pembanding eksak Cholesky $2n\times 2n$, skema hybrid $\kappa = 1$ (Bennedsen–Lunde–Pakkanen), estimator bersyarat Black–Scholes (Romano–Touzi / McCrickerd–Pakkanen), hukum pangkat *skew* ATM $\psi(T) \sim T^{H-1/2}$, serta kalibrasi & identifiabilitas $(H, \eta, \rho)$.

Dokumentasi penjelasan proyek lengkap tersedia di [day004h_quant_rough-bergomi.md](../../../docs/project-explanations/day004h_quant_rough-bergomi.md).

## Intuisi Ekonomi & Matematis

Dalam pasar opsi nyata, *skew* volatilitas implisit ATM untuk jatuh tempo pendek $T \to 0$ meledak mengikuti hukum pangkat $\psi(T) \sim T^{H - 1/2}$ dengan $H \approx 0{,}1$ (jauh lebih curam dibanding model Markov 1-faktor $H = 0{,}5$ yang memiliki $\psi(T) \to \mathrm{konstan}$). Model **rough Bergomi** menangkap fenomena ini dengan menggerakkan variansi sesaat memakai proses Volterra fraksional $\widetilde{W}^H_t$:
- **Kasar ($H < 1/2$)**: lintasan volatilitas lebih tidak teratur (*rough*) dibanding gerak Brown standar, menghasilkan *skew* jangka pendek yang tajam tanpa memerlukan lompatan (*jumps*).
- **Reduksi variansi bersyarat**: bersyarat pada lintasan penggerak volatilitas $B$, komponen ortogonal $W^\perp$ berdistribusi Gaussian murni, sehingga harga opsi bersyarat dapat dihitung analitik melalui rumus Black–Scholes dengan spot efektif $S_{\mathrm{eff}}$ dan variansi tersisa $(1-\rho^2)\int_0^T V_t\,dt$.

## Rumus Kunci

- **Proses Volterra & dinamika rough Bergomi ($r = 0$)**:
  $$\widetilde{W}^H_t = \sqrt{2H}\int_0^t (t - s)^{H - \frac{1}{2}}\,dB_s, \qquad V_u = \xi_0(u)\exp\!\left(\eta \widetilde{W}^H_u - \frac{1}{2}\eta^2 u^{2H}\right), \qquad dS_t = S_t\sqrt{V_t}\,dW_t,$$
  dengan $W = \rho B + \sqrt{1 - \rho^2}\,W^\perp$ dan parameter acuan $S_0 = 1$, $\xi_0 = 0{,}04$, $\eta = 1{,}9$, $H = 0{,}1$, $\rho = -0{,}7$.
- **Kovarians eksak**:
  $$\operatorname{Var}(\widetilde{W}^H_t) = t^{2H}, \qquad \operatorname{Cov}(\widetilde{W}^H_s, \widetilde{W}^H_t) = \frac{2H}{H + \frac{1}{2}} s^{H+\frac{1}{2}} t^{H-\frac{1}{2}} {}_2F_1\!\left(\frac{1}{2}-H,\, 1;\, H+\frac{3}{2};\, \frac{s}{t}\right) \quad (0 < s \le t).$$
- **Skema hybrid $\kappa = 1$ (kuadratur eksak sel pertama + titik kanan sel lain)**:
  $$\widetilde{W}^H_{t_i} \approx \Delta t^H\left(\rho_H U_i + \sqrt{1-\rho_H^2}\,V_i\right) + \sqrt{2H}\sum_{k=2}^{i} \bigl((k-1)\Delta t\bigr)^{H-\frac{1}{2}}\Delta B_{i-k+1}, \qquad \rho_H = \frac{\sqrt{2H}}{H + \frac{1}{2}}.$$
- **Estimator bersyarat Black–Scholes**:
  $$S_{\mathrm{eff}} = S_0 \exp\!\left(\rho \int_0^T \sqrt{V_t}\,dB_t - \frac{1}{2}\rho^2 \int_0^T V_t\,dt\right), \qquad V_{\mathrm{cond}} = (1-\rho^2)\int_0^T V_t\,dt, \qquad C = \mathbb{E}\bigl[\mathrm{BS}(S_{\mathrm{eff}}, K, V_{\mathrm{cond}})\bigr].$$

## Struktur Kode

```
projects/rust/day004h_rough_bergomi/
├── Cargo.toml
├── meta.json
├── README.md
├── src/
│   ├── lib.rs          # Re-export modul
│   ├── rng.rs          # xoshiro256++ + Box-Muller + SplitMix64 (deterministik per indeks path)
│   ├── parallel.rs     # Worker pool std::thread::scope deterministik (reduksi blok tetap bit-identik)
│   ├── stats.rs        # Rata-rata, SD, SE, regresi OLS
│   ├── linalg.rs       # Dekomposisi Cholesky, solver SPD, dan invers matriks padat 3x3
│   ├── quad.rs         # Kuadratur Gauss-Legendre dari akar polinomial Legendre
│   ├── bs.rs           # Harga Black-Scholes, vega, dan implied volatility (Newton + bisection)
│   ├── volterra.rs     # Kovarians eksak (2F1 & kuadratur bebas singularitas), Cholesky 2n x 2n, skema hybrid
│   ├── rbergomi.rs     # Anchor martingal, variance swap, MC murni vs bersyarat BS, pembanding Cholesky
│   ├── smile.rs        # Smile 13x5 (R=100k, n=250), hukum pangkat skew ATM, dan pembanding Markov H=0.5
│   ├── calib.rs        # Pasar sintetis (R=200k, n=500), Nelder-Mead CRN, profil H, Hessian, 30 replikasi
│   ├── jsonout.rs      # Penulis JSON 17 digit presisi tanpa dependensi eksternal
│   ├── svg.rs          # Generator 7 grafik SVG
│   └── bin/run.rs      # Runner utama yang menghasilkan outputs/day004h/results.json dan SVG
└── tests/
    └── integration_tests.rs
```

## Cara Build, Test, dan Run

Jalankan dari direktori proyek `projects/rust/day004h_rough_bergomi`:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --release
cargo run --release --bin run
```

Keluaran utama ditulis ke [../../../outputs/day004h/results.json](../../../outputs/day004h/results.json) dan 7 grafik ke [../../../outputs/figures/](../../../outputs/figures/).

## Ringkasan Hasil Eksperimen (`main_seed = 20261122`, `8` worker, `423,41s`)

| Bagian / Besaran | Nilai Terukur | Acuan / Pembanding | Status Kriteria |
| --- | --- | --- | --- |
| Kovarians Volterra (`2F1` vs Gauss-Legendre) | `max_rel_err = 2,12e-10` | `< 1e-8` | `c1 = tercapai` |
| Referensi Cholesky `2n x 2n` (`n=50, R=200k`) | `max_chol_err = 5,6e-16` | `max_cov_z = 3,8415 <= 4,12` (Bonferroni) | `c2 = tercapai` |
| Bias variansi hybrid BLP optimal ($H=0{,}1$) | `-0,062%` ($n=50$) → `-0,039%` ($n=500$) | Titik kanan: `+4,94%..+3,23%`; VM: `0%` | — |
| Martingal $V_u$ & $S_T$ | Cholesky $\mathbb{E}[e^{\dots}] = 0{,}9998$ ($z=0{,}014$); BLP `max_z=2,25` | `1,0` dalam `3 SE` | `c3 = tercapai` |
| Variance swap ($T=1$) | `0,04031` (`max_z = 1,62`) | `0,04000` | `c4 = tercapai` |
| Estimator Bersyarat BS vs MC | `max_z_diff = 1,77`, `VRF = 2,15..308,4` | `parity_err = 1,97e-14` | `c5 = tercapai` |
| Harga Hybrid BLP vs Cholesky (`n=50, T=1`) | `max_z_comb = 1,23` | `<= 3 SE` | `c6 = tercapai` |
| Hukum pangkat *skew* ATM ($\log \psi$ vs $\log T$) | $H=0{,}1$: `-0,3849 ± 0,0130`; $H=0{,}5$: `-0,0404 ± 0,0253` | `-0,40` dan `0,00` | `c7 = tercapai` |
| Kalibrasi & replikasi ($R_{\mathrm{rep}} = 30$, *same-n*) | `bias_H = +0,0719 ± 0,0158` | `|bias_H| < 0,03` | `c8 = belum tercapai` |
| Anggaran waktu runner | `423,41s` (`8` worker) | `< 300s` | `c9 = belum tercapai` |

## Catatan Diagnostik (`c8`)

Dengan skema BLP titik optimal, `c3`, `c4`, dan `c6` kini seluruhnya `tercapai`. Uji *self-consistency* memulihkan `(0,100001; 1,900000; -0,699997)` (`rmse = 1,18e-7`), membuktikan tidak ada bug *optimizer/engine*; namun `c8` tetap `belum tercapai` pada *same-scheme-same-n* (`bias_H = +0,0719 ± 0,0158`), bank CRN berbeda per replikasi (`rep_diff_crn`: `bias_H = +0,1027 ± 0,0514`), maupun `R_calib = 80.000` (`H_hat = 0,1657`), yang membuktikan akar masalahnya adalah *ridge* identifiabilitas datar antara $(H, \eta)$ ($\operatorname{Corr} = 0{,}778$) pada horizon $T \in [0{,}1; 1{,}0]$.

## Daftar Grafik

- [day004h_cov-exact-vs-hybrid.svg](../../../outputs/figures/day004h_cov-exact-vs-hybrid.svg)
- [day004h_hybrid-bias-vs-n.svg](../../../outputs/figures/day004h_hybrid-bias-vs-n.svg)
- [day004h_smile-term-structure.svg](../../../outputs/figures/day004h_smile-term-structure.svg)
- [day004h_skew-loglog.svg](../../../outputs/figures/day004h_skew-loglog.svg)
- [day004h_vrf-conditional-bs.svg](../../../outputs/figures/day004h_vrf-conditional-bs.svg)
- [day004h_calib-h-profile.svg](../../../outputs/figures/day004h_calib-h-profile.svg)
- [day004h_rough-vs-markov-t01.svg](../../../outputs/figures/day004h_rough-vs-markov-t01.svg)
