# day002g — AAD untuk Greeks Monte Carlo (C++20)

Implementasi mesin **Adjoint Algorithmic Differentiation (AAD)** mode-balik (*tape* + *mark/rewind*), mode maju (`Dual<T>`), *forward-over-reverse* (`Dual<Number>`), dan Cholesky terdiferensiasi dalam **C++20 murni** tanpa library AD atau numerik eksternal:

- **Tape Mode-Balik + Mark/Rewind**: node menyimpan nilai, *adjoint*, paling banyak dua operan induk, dan turunan parsial lokal. `mark()` merekam bagian bersama semua path (parameter model dan faktorisasi Cholesky $10 \times 10$). Pada tiap path Monte Carlo, rekaman di atas *mark* disapu balik hingga *mark* (`backprop_to_mark`) lalu di-*rewind* (`rewind_to_mark()`), sehingga memori tape tetap konstan kecil ($O(\text{node per path})$) dan graf bersama tidak pernah direkam ulang.
- **Mode Maju & Forward-over-Reverse**: `Dual<double>` memverifikasi AAD, sedangkan `Dual<Number>` menghitung turunan orde dua (Gamma Black-Scholes) dalam satu sapuan maju-lalu-balik.
- **Payoff Diskontinu**: *delta* pathwise pada opsi *digital call* bernilai 0 di hampir semua titik; penghalusan $\Phi((S_T - K)/\varepsilon)$ mengembalikan gradien dengan *tradeoff* bias–varians yang jelas.

## Cara Build, Test & Run

```bash
cd cpp/day002g_aad_greeks

# Build Release + unit test
cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release
cmake --build build_release -j
ctest --test-dir build_release --output-on-failure

# Build AddressSanitizer + UBSan
cmake -S . -B build_asan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZERS=ON
cmake --build build_asan -j
ctest --test-dir build_asan --output-on-failure

# Jalankan benchmark penuh (menulis outputs/day002g/results.json & SVG)
./build_release/day002g_run
```

## Struktur Kode

- `include/aad/tape.hpp` & `include/aad/number.hpp` — `Tape` (*mark/rewind*), `TapeGuard`, dan tipe `Number` dengan operator serta fungsi `exp`, `log`, `sqrt`, `pow`, `abs`, `normal_cdf` (`erfc`), `max`.
- `include/aad/dual.hpp` — `Dual<T>` untuk mode maju dan *forward-over-reverse*.
- `include/aad/cholesky.hpp` — faktorisasi Cholesky terdiferensiasi dari 45 input $\rho_{ij}$ ($i < j$).
- `include/aad/rng.hpp` — `SplitMix64`, `Xoshiro256StarStar`, dan Box–Muller polar buatan sendiri.
- `include/aad/models.hpp` — Black-Scholes dan *geometric basket* tertutup bertemplat.
- `include/aad/mc_engine.hpp` & `src/mc_engine.cpp` — mesin Monte Carlo paralel deterministik (`std::thread`, chunk 4000, reduksi Welford berurut).
- `include/aad/svg.hpp`, `src/svg.cpp`, `app/main.cpp`, `tests/test_aad.cpp`.
