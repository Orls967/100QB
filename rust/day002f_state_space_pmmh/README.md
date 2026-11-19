# day002f — Inferensi Bayesian State-Space: Kalman, Particle Filter & PMMH (Rust)

Implementasi inferensi Bayesian untuk model *state-space* dalam **Rust murni** tanpa crate aljabar linear, statistik, atau MCMC siap pakai (Cholesky, Kalman filter, bootstrap particle filter dengan resampling sistematik tanpa alokasi di loop panas, 2D grid posterior, Haario adaptive MH/PMMH, split $\hat{R}$, ESS *initial positive sequence*, dan MCSE ditulis sendiri):

- **Model A (Local Level, linear-Gaussian)**: $x_0 \sim \mathcal{N}(0, 1)$, $x_t = x_{t-1} + \sigma_\eta \eta_t$, $y_t = x_t + \sigma_\varepsilon \varepsilon_t$ ($T = 200$, $\sigma_\eta = 0.5, \sigma_\varepsilon = 1.0$, `seed = 42`).
- **Model B (Stochastic Volatility, nonlinear)**: $h_0 \sim \mathcal{N}(\mu, \sigma^2/(1-\phi^2))$, $h_t = \mu + \phi(h_{t-1} - \mu) + \sigma \eta_t$, $y_t = \exp(h_t/2)\varepsilon_t$ ($T = 500$, $\mu = -1.0, \phi = 0.95, \sigma = 0.25$, `seed = 43`).

## Intuisi Pseudo-Marginal (PMMH)
Ketika likelihood $p(y_{1:T} \mid \theta)$ tidak tersedia secara tertutup (seperti pada model SV nonlinear), **Particle Marginal Metropolis–Hastings (PMMH)** mengganti $p(y_{1:T} \mid \theta)$ dengan estimator *particle filter* $\hat{p}_N(y_{1:T} \mid \theta)$ yang **tak bias** ($\mathbb{E}[\hat{p}_N] = p$). Secara mengejutkan, rantai Markov pada ruang diperluas $(\theta, u)$ tetap memiliki distribusi stasioner **eksak** $p(\theta \mid y_{1:T})$ untuk setiap jumlah partikel $N_p \ge 1$! Pemilihan $N_p$ sehingga $\text{SD}(\ln \hat{p}_N) \lesssim 1.5$ menjaga *acceptance rate* dan ESS tetap efisien.

## Cara Build, Test & Run

```bash
cd rust/day002f_state_space_pmmh
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
cargo run --release
```

Output dihasilkan ke:
- `outputs/day002f/results.json`
- `outputs/figures/day002f_pf-unbiased-and-var.png`
- `outputs/figures/day002f_model_a_posteriors.png`
- `outputs/figures/day002f_model_b_sv.png`

## Struktur Kode

- `src/rng.rs` — `ChaCha8Rng` deterministik (`derive_seed` SplitMix64), Cholesky dense & *forward solve*.
- `src/models.rs` — simulator sintetis dan *log-prior* Model A & Model B.
- `src/kalman.rs` — Kalman filter eksak dan verifikasi *brute force* via Cholesky ($T = 10$).
- `src/pf.rs` — Bootstrap Particle Filter (`PfWorkspace` pra-alokasi) & resampling sistematik.
- `src/grid.rs` — posterior grid 2D eksak Model A ($200 \times 200$, integrasi trapesium).
- `src/mcmc.rs` — MH & PMMH paralel (`rayon`) dengan kovarians proposal adaptif Haario saat *burn-in*.
- `src/diagnostics.rs` — split $\hat{R}$, ESS (*initial positive sequence*), MCSE, dan kuantil.
- `src/plots.rs` & `src/main.rs` — generator grafik PNG `plotters` dan orkestrator eksperimen.

## Ringkasan Hasil (`cargo run --release`, `146.87s`)

- **Kalman vs Brute Force ($T = 10$)**: selisih maksimum `7.11e-15`.
- **Ketidakbiasan PF ($N_p = 500, R = 2000$)**: $\text{mean}(\exp(\hat{\ell} - \ell)) = 0.9824 \pm 0.0204$ ($z = 0.86$), *slope* $\ln \text{Var}(\hat{\ell})$ vs $\ln N_p$ = `-0.963`, throughput `683.5 PF/detik` (`68.35M partikel-langkah/detik`).
- **Model A (`Grid` vs `MH Eksak` vs `PMMH N_p=500`)**:
  - $a = \ln \sigma_\eta$: Grid `-0.5843 ± 0.1371`, MH `-0.5849 ± 0.1373` ($\hat{R}=1.0003$, ESS `20556`), PMMH `-0.5837 ± 0.1364` ($\hat{R}=1.0002$, ESS `11405`).
  - $b = \ln \sigma_\varepsilon$: Grid `-0.1733 ± 0.0778`, MH `-0.1726 ± 0.0778` ($\hat{R}=1.0003$, ESS `21139`), PMMH `-0.1739 ± 0.0777` ($\hat{R}=1.0002$, ESS `11587`).
- **Model B (SV, PMMH $N_p=100$ vs $N_p=200$)**:
  - $\mu = -0.831 \pm 0.330$ (CI 95% `[-1.398, -0.077]`, sejati `-1.0`, $\hat{R}=1.0016$, ESS `1814`).
  - $\phi = 0.9558 \pm 0.0230$ (CI 95% `[0.903, 0.991]`, sejati `0.95`, $\hat{R}=1.0017$, ESS `1537`).
  - $\sigma = 0.2627 \pm 0.0546$ (CI 95% `[0.172, 0.394]`, sejati `0.25`, $\hat{R}=1.0019$, ESS `1753`).
  - Selisih antara $N_p = 100$ dan $N_p = 200$: `[0.55, 1.43, 1.33]` MCSE gabungan ($< 3$).
