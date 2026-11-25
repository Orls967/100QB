# `day008a_ising_trilingual` — Model Ising 2D Lintas Tiga Bahasa (C11, C++20, Python)

Proyek ini mengimplementasikan simulasi Model Ising 2D periodik ($J=1, H=0$) dengan algoritma **Metropolis**, **Wolff**, dan **Swendsen–Wang**, mengikuti spesifikasi tunggal [`SPEC.md`](SPEC.md) secara literal di tiga bahasa: **C11** (`c/`), **C++20** (`cpp/`), dan **Python** ([`days/day008a_sim_ising-trilingual.ipynb`](../../../days/day008a_sim_ising-trilingual.ipynb)).

## Struktur Berkas

- `SPEC.md`: Kontrak tunggal RNG (`SplitMix64` + `xoshiro256++`), urutan situs/ikatan, Union-Find, hash FNV-1a, beban uji, dan kriteria a priori `c1`–`c14`.
- `Makefile`: Kompilasi `-std=c11` / `-std=c++20` dengan `-Wall -Wextra -Werror -O3 -fno-fast-math -ffp-contract=off`, target `all`, `test`, `asan`, `parity`, `clean`.
- `c/`: `rng.h`, `ising.h`, `ising.c`, `main.c`, `test_c.c` (paralel `pthread`, tabel penerimaan bebas cabang + tetangga prekomputasi, serta *multi-spin coding* 64-bit opsional).
- `cpp/`: `rng.hpp`, `ising.hpp`, `ising.cpp`, `main.cpp`, `test_cpp.cpp` (paralel `std::thread`, templat `constexpr L` vs runtime `L`).
- `data/`: Keluaran CSV sementara (dapat dihasilkan ulang via `Makefile` dan CLI; tidak di-commit bila > 5 MB).

## Perintah Cepat

```bash
make all
make test
make parity
./build/ising_c --mode enum --out-dir data
./build/ising_cpp --mode enum --out-dir data
./build/ising_c --mode bench --out-dir data
./build/ising_cpp --mode bench --out-dir data
./build/ising_c --mode full --threads 8 --out-dir data
./build/ising_cpp --mode full --threads 8 --out-dir data
```
