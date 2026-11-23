# Day 005d — Jump-Diffusion Merton: Seri Eksak, PIDE Crank–Nicolson, Kalibrasi & Risiko Ekor (Rust)

Proyek Rust mandiri (*Rust stable, edition 2021*, tanpa crate eksternal) menggunakan **topik cadangan (`ROADMAP 075`)** setelah topik utama (*Marchenko–Pastur & Ledoit–Wolf*) terdeteksi telah dikerjakan penuh pada `day002b`. Proyek ini mengimplementasikan model **Merton Jump-Diffusion (1976)** secara lengkap pada data sintetis:

$$\frac{dS_t}{S_{t^-}} = (r - q - \lambda \kappa)\,dt + \sigma\,dW_t + (e^J - 1)\,dN_t, \qquad J \sim \mathcal{N}(\mu_J, \delta_J^2), \qquad \kappa = e^{\mu_J + \frac{1}{2}\delta_J^2} - 1,$$

dengan parameter dasar $S_0 = 100$, $r = 0{,}05$, $q = 0{,}01$, $\sigma = 0{,}15$, $\lambda = 0{,}80$, $\mu_J = -0{,}10$, $\delta_J = 0{,}25$.

## Intuisi & Rumus Kunci

1. **Seri Eksak Merton (1976)**: Bersyarat pada jumlah lompatan $N_T = n \sim \mathrm{Poisson}(\lambda T)$, harga opsi *call* Eropa adalah campuran Poisson dari harga Black–Scholes dengan $S_n = S_0 e^{n(\mu_J + \frac{1}{2}\delta_J^2) - \lambda\kappa T}$ dan $\sigma_n^2 = \sigma^2 + n\delta_J^2/T$.
2. **PIDE Crank–Nicolson + Rannacher**: Pada koordinat log-spot $x = \ln S$ ($\tau = T - t$):
   $$\frac{\partial V}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 V}{\partial x^2} + \left(r - q - \frac{1}{2}\sigma^2 - \lambda\kappa\right)\frac{\partial V}{\partial x} - (r + \lambda)V + \lambda \int_{-\infty}^{\infty} V(\tau, x+y)\phi_{\mu_J,\delta_J}(y)\,dy.$$
   Bagian diferensial diselesaikan secara implisit dengan algoritma Thomas $O(N_x)$ dan 4 *half-step* Rannacher, sedangkan integral lompatan non-lokal dievaluasi via skema *predictor-corrector* dengan kelanjutan asimtotik *far-field*.
3. **Implied Volatility Smile & Kalibrasi**: Lompatan negatif ($\mu_J = -0{,}10$) menciptakan *skew* tajam pada jatuh tempo pendek ($T=0{,}25$: `35,40%` di $K=80$ vs `22,55%` ATM) yang melandai seiring $T$ membesar. Optimisasi Nelder–Mead 4D memulihkan $(\sigma, \lambda, \mu_J, \delta_J)$ hingga `price_rmse = 1,15e-8`.
4. **Risiko Ekor**: Model Black–Scholes dengan volatilitas ATM meremehkan harga *put* OTM $K=70$ hingga `191,14x` dan meremehkan *Expected Shortfall* 99% posisi *short put* 1 bulan (`28,156` vs `4,017`).

## Build, Run & Test

```bash
cargo fmt --manifest-path projects/rust/day005d_merton_jump_diffusion/Cargo.toml --all -- --check
cargo clippy --manifest-path projects/rust/day005d_merton_jump_diffusion/Cargo.toml --all-targets --all-features -- -D warnings
cargo test --release --manifest-path projects/rust/day005d_merton_jump_diffusion/Cargo.toml
cargo run --release --manifest-path projects/rust/day005d_merton_jump_diffusion/Cargo.toml --bin run
cargo run --release --manifest-path projects/rust/day005d_merton_jump_diffusion/Cargo.toml --bin run_extras
```

## Struktur Kode

- [src/rng.rs](src/rng.rs) — `SplitMix64`, `xoshiro256++`, normal, Poisson, Gamma
- [src/special.rs](src/special.rs) & [src/stats.rs](src/stats.rs) — PDF/CDF normal (Cody) dan statistik sampel
- [src/parallel.rs](src/parallel.rs) — Reduksi paralel deterministik blok tetap `64`
- [src/linalg.rs](src/linalg.rs) — Pemecah tridiagonal Thomas $O(N)$
- [src/bs.rs](src/bs.rs) — Harga Black–Scholes, vega, dan *implied volatility* (Newton + bisection)
- [src/merton.rs](src/merton.rs) — Seri eksak Merton (1976), *put-call parity*, dan simulasi MC paralel ($R=400\mathrm{k}$)
- [src/pide.rs](src/pide.rs) — Solver PIDE Crank–Nicolson + Rannacher + uji konvergensi grid
- [src/calib.rs](src/calib.rs) — Permukaan *smile* dan kalibrasi Nelder–Mead 4D
- [src/tail_risk.rs](src/tail_risk.rs) — Perbandingan *put* OTM dan 99% VaR/ES *short put* 1 bulan
- [src/bin/run.rs](src/bin/run.rs) & [src/bin/run_extras.rs](src/bin/run_extras.rs) — Runner utama dan eksperimen tambahan
- [tests/integration_tests.rs](tests/integration_tests.rs) — Tes integrasi invariansi worker, determinisme seed, dan 3 SE

## Ringkasan Hasil (`MAIN_SEED = 20261123`, `10` worker, `0,49s`)

| Bagian | Metrik Utama | Hasil Terukur |
|---|---|---|
| 1. Seri vs MC ($R=400\mathrm{k}$) | Martingal & harga call/put | `martingale_z = 0,49`, `max |call_z| = 1,10`, `max |put_z| = 0,75` |
| 2. PIDE CN vs Eksak | Grid `401x200` & konvergensi | `max_abs_err = 0,00303` (`rel = 0,00066`), rasio `101->401` = `15,43x` |
| 3. Smile & Kalibrasi | Skew $T=0{,}25$ & Nelder–Mead | `IV(K=80) = 35,40% > ATM 22,55%`; `price_rmse = 1,15e-8`, `iv_rmse = 9,17e-6 bps` |
| 4. Risiko Ekor | Put OTM $K=70$ & 1m ES 99% | Put `K=70`: `191,14x` BS; `ES99` = `28,156` (Merton) vs `4,017` (BS) |

Artefak keluaran: [results.json](../../../outputs/day005d/results.json), 4 grafik SVG `outputs/figures/day005d_*.svg`, dan dokumen penjelasan [day005d_quant_merton-jump-diffusion.md](../../../docs/project-explanations/day005d_quant_merton-jump-diffusion.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
