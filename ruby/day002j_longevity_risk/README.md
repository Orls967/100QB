# Day 002j — Longevity Risk: Lee–Carter, Anuitas Stokastik & Modal Solvabilitas (Ruby)

Implementasi Ruby 3.3+ murni (tanpa gem runtime eksternal, YJIT aktif, paralelisme deterministik `Process.fork` + `Marshal`) untuk pemodelan risiko umur panjang (*longevity risk*) dan modal solvabilitas portofolio anuitas.

## Intuisi Model

1. **Model Lee–Carter (1992)**: memfaktorkan permukaan log-mortalitas $\ln m_{x,t} = a_x + b_x k_t + \varepsilon_{x,t}$ menggunakan SVD rank-1 buatan sendiri (rotasi Jacobi satu sisi Hestenes) dengan kendala identifikasi $\sum_x b_x = 1$ dan $\sum_t k_t = 0$.
2. **Kurva Gompertz–Makeham**: memodelkan mortalitas dasar $\mu_x = A + B c^x$ lewat optimizer Nelder–Mead buatan sendiri dalam ruang log.
3. **Proyeksi Stokastik & Ketidakpastian Parameter**: indeks waktu $k_t$ diproyeksikan sebagai *random walk with drift* ($H = 40$ tahun, $N = 20.000$ jalur), termasuk ketidakpastian estimasi drift $\hat{d} \sim \mathcal{N}(\hat{d}, \hat{\sigma}^2 / (T-1))$.
4. **Dekomposisi Risiko Portofolio**: berdasarkan Hukum Total Variansi,
   $$\mathrm{Var}(Y) = N\,\mathbb{E}[\sigma^2(\text{path})] + N^2\,\mathrm{Var}(\mu(\text{path})),$$
   risiko idiosinkratik ($O(N)$) terdiversifikasi saat $N$ membesar, sedangkan risiko sistematik umur panjang ($O(N^2)$) tidak dapat didiversifikasi sehingga $\mathrm{CV}(N)$ mendatar ke batas $\sqrt{\mathrm{Var}(\mu)}/\bar{\mu} > 0$.

## Struktur Kode

- `lib/longevity/rng.rb`: `Random.new(seed)` + Box–Muller + derivasi seed deterministik.
- `lib/longevity/parallel.rb`: paralelisme `Process.fork` + `IO.pipe` (`Marshal`) dengan chunk tetap.
- `lib/longevity/linalg.rb`: aljabar linear `Array<Array<Float>>` dan SVD Jacobi satu sisi.
- `lib/longevity/optimizers.rb`: optimizer Nelder–Mead.
- `lib/longevity/lee_carter.rb`: generator data sintetis & estimasi Lee–Carter.
- `lib/longevity/gompertz_makeham.rb`: fit Gompertz–Makeham.
- `lib/longevity/forecast.rb`: proyeksi $k_t$, bootstrap parametrik, fan chart $e_{65}$.
- `lib/longevity/annuity.rb`: rekursi mundur, jumlah maju, variansi eksak, dan simulasi umur (pencarian biner).
- `lib/longevity/portfolio.rb`: Hukum Total Variansi, simulasi penuh, aproksimasi normal bersyarat, VaR/ES, dan stres 0.5%.
- `lib/longevity/svg.rb`: penulis 6 grafik SVG ke `outputs/figures/`.

## Cara Menjalankan & Menguji

```bash
cd ruby/day002j_longevity_risk
bundle install
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
```

## Ringkasan Hasil (`outputs/day002j/results.json`)

- **SVD (6x4, 10x7)**: galat rekonstruksi & ortogonalitas $< 1.1 \times 10^{-14}$.
- **LC tanpa noise**: galat maksimum $a_x, b_x, k_t$ masing-masing $8.88 \times 10^{-16}, 1.39 \times 10^{-17}, 7.99 \times 10^{-15}$.
- **LC bernoise ($R = 200$)**: rasio rata-rata varians empiris/teoretis $\hat{a}_x = 1.0052$.
- **Anuitas (kohor 2020, usia 65, $r = 3\%$)**: eksak `14.718650` (var `28.0112`) vs simulasi 200.000 umur `14.715535 ± 0.013051` (var `28.1620 ± 0.0978`).
- **CV(N)**: `11.39%` ($N=10$) $\to$ `3.65%` ($N=100$) $\to$ `1.30%` ($N=1.000$) $\to$ `0.72%` ($N=10.000$) $\to$ batas `0.625%` ($N=\infty$).
- **Modal VaR 99.5% ($N = 10.000$)**: `1.83%` (tanpa ketidakpastian parameter) dan `2.20%` (dengan ketidakpastian parameter); stres deterministik 0.5% = `149,939.55` vs VaR 99.5% MC = `149,889.57`.
