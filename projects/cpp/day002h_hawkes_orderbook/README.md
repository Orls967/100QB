# day002h — Proses Hawkes + Mesin Limit Order Book & Simulasi Pasar (C++20)

Implementasi **proses titik Hawkes**, **mesin limit order book** berkinerja tinggi, dan **simulasi pasar gabungan** dalam **C++20 murni** (hanya *standard library*):

- **Proses Hawkes (Bagian A)**: $\lambda(t) = \mu + \sum_{t_i < t} \alpha e^{-\beta (t - t_i)}$ dengan $\mu = 0.5, \alpha = 0.6, \beta = 1.2$ ($n = 0.5$, $\bar{\lambda} = 1.0$, *Fano factor* $= 4.0$). Dua simulator (*Ogata thinning* & eksak *Dassios–Zhao*), uji penskalaan ulang waktu (KS & Ljung–Box), dan MLE rekursif $O(n)$ dengan *optimizer* BFGS gradien analitik + Hessian numerik ($R = 200$ replikasi paralel `std::thread`).
- **Mesin Limit Order Book (Bagian B)**: prioritas harga-waktu dengan *pool allocator* pra-alokasi, antrean FIFO intrusif per level, peta `id -> slot` $O(1)$, serta **bitmap 2-level 64-bit** (`std::countl_zero` / `std::countr_zero`) untuk pembaruan `best_bid`/`best_ask` $O(1)$ tanpa alokasi *heap* di jalur panas (`0` alokasi setelah *warmup*). Diverifikasi identik bit-for-bit terhadap mesin referensi naif $O(n)$ pada 25+ seed dan aliran `1.000.000` order.
- **Simulasi Pasar Gabungan (Bagian C)**: kedatangan *market order* Hawkes ($T = 20.000$), arah Markov ($p_{\text{same}} = 0.7$), ukuran geometrik, serta *limit/cancel* Poisson untuk mereproduksi fakta stilistik mikrostruktur (*spread*, profil kedalaman, autokorelasi aliran bertanda, *volatility clustering* $|r_t|$, dan regresi *price impact*) sambil menjaga konservasi posisi bersih $= 0$.

## Cara Build, Test & Run

```bash
cd projects/cpp/day002h_hawkes_orderbook

# Build Release + unit test
cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release
cmake --build build_release -j
ctest --test-dir build_release --output-on-failure

# Build UBSan
cmake -S . -B build_ubsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZERS=ON
cmake --build build_ubsan -j
ctest --test-dir build_ubsan --output-on-failure

# Jalankan benchmark penuh (menulis outputs/day002h/results.json & SVG)
./build_release/day002h_run
```

## Struktur Kode

- `include/hawkes_ob/rng.hpp` — `SplitMix64`, `Xoshiro256StarStar`, distribusi eksponensial & geometrik.
- `include/hawkes_ob/stats.hpp` — uji KS 1-sampel & 2-sampel, ACF, Ljung–Box, *Fano factor*, *batch means*.
- `include/hawkes_ob/hawkes.hpp` & `src/hawkes.cpp` — simulator Ogata & Dassios–Zhao, residu penskalaan ulang waktu, MLE $O(n)$ + BFGS, Monte Carlo paralel.
- `include/hawkes_ob/orderbook.hpp` & `src/orderbook.cpp` — `OrderBook` (*pool*, FIFO intrusif, bitmap 2-level).
- `include/hawkes_ob/reference_book.hpp` — mesin referensi $O(n)$ untuk verifikasi kesetaraan eksak.
- `include/hawkes_ob/market_sim.hpp` & `src/market_sim.cpp` — *benchmark* 1M order (*hook* `operator new`) dan simulasi pasar gabungan.
- `include/hawkes_ob/svg.hpp`, `src/svg.cpp`, `app/main.cpp`, `tests/test_hawkes_ob.cpp`.
