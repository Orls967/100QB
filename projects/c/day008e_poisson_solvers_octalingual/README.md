# Day 8e — 2D Poisson Solvers (Octalingual 19-Stage Pipeline)

Implementasi alur **19 tahap** (`s01`–`s19`, plus `all`) untuk menyelesaikan persamaan Poisson 2D $-\Delta u = f$ pada $(0,1)^2$ dengan syarat batas Dirichlet $u = 0$ (beda hingga 5-titik orde 2) secara **independen dan penuh** di **8 bahasa pemrograman**:

| Bahasa | Sumber | Toolchain | $N_{\max}$ | $N_{\text{stat},\max}$ |
| :--- | :--- | :--- | :---: | :---: |
| C11 | `c/poisson.c` | `cc -std=c11 -O3 -fno-fast-math -ffp-contract=off` | 511 | 127 |
| C++20 | `cpp/poisson.cpp` | `c++ -std=c++20 -O3 -fno-fast-math -ffp-contract=off` | 511 | 127 |
| Rust | `rust/src/main.rs` | `cargo build --release --offline` | 511 | 127 |
| Go | `go/main.go` | `go build` (`//go:noinline fm()` anti-FMA) | 511 | 127 |
| TypeScript | `ts/poisson.ts` | `node --experimental-strip-types` | 255 | 63 |
| Ruby | `ruby/poisson.rb` | `ruby` | 127 | 63 |
| Python | `py/poisson.py` | `.venv/bin/python` | 127 | 63 |
| Jupyter Notebook | `../../../days/day008e_sim_poisson-solvers-octalingual.ipynb` | `nbclient` | 63 | 63 |

## Metode (7 Solver)
1. **Jacobi** (`s04`)
2. **Gauss–Seidel leksikografis** (`s05`)
3. **Red–Black Gauss–Seidel** (`s06`)
4. **SOR ($\omega_{\text{opt}}$)** (`s07`)
5. **Conjugate Gradient (CG)** (`s08`)
6. **PCG-SSOR ($\omega = 1$)** (`s09`)
7. **Multigrid V(2,1) + FMG** (`s10`, `s11`, `s15`)

## Cara Menjalankan
```bash
make test     # build 4 bahasa kompilasi + jalankan selftest (s01) 7 bahasa CLI
make run      # jalankan seluruh 8 bahasa (19 tahap) + compare.py
```

Keluaran per-bahasa disimpan di `out/<lang>/`, ringkasan verifikasi di `outputs/day008e/results.json`, dan 4 figur gabungan di `outputs/figures/day008e_*.{png,svg}`.
