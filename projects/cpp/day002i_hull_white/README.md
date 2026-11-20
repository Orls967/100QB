# Day 002i — Hull–White 1-Factor Short-Rate Model (C++20)

Implementasi C++20 murni (hanya standard library, tanpa library eksternal) untuk model short-rate Hull–White 1-faktor:

$$dr_t = (\theta(t) - a\,r_t)\,dt + \sigma\,dW_t, \qquad a = 0.08,\quad \sigma = 0.012$$

di atas kurva zero Nelson–Siegel ($\beta_0 = 0.045$, $\beta_1 = -0.02$, $\beta_2 = 0.015$, $\tau = 2.0$).

## Komponen

1. **Pohon Trinomial Terkalibrasi ke Kurva (`src/trinomial_tree.cpp`)**
   - Konstruksi pohon Hull–White (1994) dengan induksi maju Arrow–Debreu untuk $m \in \{12, 24, 48, 96, 192\}$ langkah/tahun.
   - Reproduksi kurva diskonto awal hingga $< 1.5 \times 10^{-15}$.
   - Harga ZBP, ZBC ($T=1, S=5$, strike forward, 90%, 110%), cap tahunan 5y (strike 4.5%), dan ekstrapolasi Richardson.
2. **Swaption Eropa: Jamshidian vs Pohon vs Monte Carlo (`src/hull_white.cpp`, `src/montecarlo.cpp`)**
   - Dekomposisi Jamshidian (solver akar Brent) untuk swaption payer koterminal 10y ($T_0 \in \{1, 2, 3, 5, 7, 9\}$, $K = 4.5\%$).
   - Monte Carlo antitetik transisi OU eksak (`N = 200,000` path, `dt = 1/52`, integrasi trapesium) dengan *control variate* martingale $D(T_0)P(T_0,T_j) - P(0,T_j)$ dan uji martingale $\mathbb{E}[D(T)] = P(0,T)$ untuk $T \in \{1, 5, 10\}$.
3. **Kalibrasi Parameter (`src/calibration.cpp`)**
   - Kalibrasi $(a, \sigma)$ dalam ruang log via Nelder–Mead dari tebakan awal $(0.2, 0.02)$ dan $(0.02, 0.005)$ (tanpa noise: RMSE $\approx 2 \times 10^{-17}$).
   - 100 replikasi kalibrasi bernoise (noise relatif $0.5\%$) secara paralel (`std::thread` deterministik), melaporkan `mean`, `SD`, `corr(a, sigma)`, dan bilangan kondisi Jacobian.
4. **Swaption Bermudan & LSM (`src/bermudan.cpp`)**
   - Swaption payer Bermudan 10y (eksekusi tahunan $1\dots 9$) pada pohon ($m \in \{12..192\}$ + Richardson) dan batas $\max \text{Euro} \le \text{Bermudan} \le \sum \text{Euro}$.
   - *Longstaff–Schwartz Monte Carlo* dua lintasan (*two-pass*, `N = 100,000` per lintasan, basis polinomial derajat 3 dalam $x$).
5. **DV01 Pergeseran Paralel $\pm 1\text{ bp}$**
   - Perbandingan DV01 swaption Eropa 5x5 (pohon $m=192$ vs Jamshidian) dan Bermudan 10y.

## Cara Build & Menjalankan

```bash
cd projects/cpp/day002i_hull_white
cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release
cmake --build build_release -j
ctest --test-dir build_release --output-on-failure

# Build UBSan
cmake -S . -B build_ubsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZERS=ON
cmake --build build_ubsan -j
ctest --test-dir build_ubsan --output-on-failure

# Jalankan benchmark penuh (dari root repo)
./cpp/day002i_hull_white/build_release/day002i_run
```

## Output

- `outputs/day002i/results.json`
- `outputs/figures/day002i_tree-convergence.svg`
- `outputs/figures/day002i_curve-and-swaptions.svg`
- `outputs/figures/day002i_calib-scatter.svg`
- `outputs/figures/day002i_bermudan-bounds-dv01.svg`
