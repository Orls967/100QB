# day002k — Opsi Put Amerika (Ruby): CRR Binomial, LSM Monte Carlo & PSOR PDE

Implementasi **Opsi Put Amerika** dalam **Ruby murni (Ruby 3.3+)** tanpa gem numerik eksternal (nol dependensi *runtime*, hanya *standard library* Ruby):

- Parameter standar: `S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `seed = 42`.
- **CRR Binomial Tree**: pohon binomial Cox–Ross–Rubinstein hingga `n = 4.000` langkah, lengkap dengan harga Eropa, *early exercise premium* (EEP), Greeks (`Delta`, `Gamma`, `Theta`), serta *convergence sweep* paralel (`n = 50..4000`) via `Process.fork`.
- **LSM Monte Carlo (Longstaff–Schwartz)**: `100.000` path (`50.000` pasang *antithetic*) dengan `50` titik exercise, regresi polinomial derajat 3 (`[1, x, x^2, x^3]`) memakai *linear solver* Gauss + *partial pivoting* buatan sendiri, simulasi path paralel deterministik (`Process.fork` + `Random.new(seed)` + Box–Muller polar).
- **PSOR PDE (Crank–Nicolson + Rannacher Smoothing)**: grid log-harga `x = ln(S / K)` ukuran `400 x 400`, Projected Successive Over-Relaxation (`omega = 1.25`, `tol = 1e-9`) untuk batas bebas (*free boundary*) Amerika dan algoritma Thomas tridiagonal untuk pembanding Eropa.

## Cara Install, Test & Run

Butuh **Ruby 3.3+** (diuji pada Ruby 4.0.9 + YJIT) dan Bundler.

```bash
cd ruby/day002k_american_put
bundle install
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
```

Output dihasilkan ke:
- `outputs/day002k/results.json`
- `outputs/figures/day002k_american-convergence.svg`
- `outputs/figures/day002k_american-spot-profile.svg`

## Ringkasan Hasil (`seed = 42`)

- **Black-Scholes Put Eropa (tertutup)**: `5.57352602`
- **CRR Binomial (`n = 4.000`)**: `6.09018072` (Eropa: `5.57302611`, EEP: `0.517155`, waktu `407.1 ms`)
- **PSOR PDE (`400 x 400`)**: `6.08755108` (Eropa: `5.57120004`, rata-rata sweep: `15.5`, waktu `206.5 ms`)
- **LSM Monte Carlo (`100.000 x 50`, antithetic)**: `6.051915 ± 0.013566` (CI 95% `[6.025325, 6.078505]`, waktu `1200.6 ms`)
- **Cross-check**: `|CRR - PSOR| = 2.63e-3`, `|LSM - CRR| = 2.82 SE` (Bermudan 50-step vs kontinu).
