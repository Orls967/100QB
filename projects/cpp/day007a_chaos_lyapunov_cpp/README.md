# Hari 007a — Kekacauan Deterministik dalam C++20

Implementasi C++20 murni (tanpa pustaka eksternal) untuk:
1. **Peta Logistik** $x_{k+1} = r x_k (1 - x_k)$: sapuan eksponen Lyapunov $\lambda(r)$ pada $r \in [2{,}5, 4{,}0]$ (`1600` titik, `200.000` iterasi pasca-transien) dan uji orbit panjang $10^7$ iterasi di $r = 4$ terhadap $\ln 2$.
2. **Universalitas Feigenbaum**: pencarian akar Newton atas siklus superstabil $f_r^{2^n}(1/2) = 1/2$ ($n = 1\dots 8$), rasio $\delta_8$, dan ekstrapolasi titik akumulasi $r_\infty$.
3. **Sistem Lorenz** ($\sigma = 10, \rho = 28, \beta = 8/3$): integrasi RK4 orde-4 bersama persamaan variasional tangen $3\times 3$ dan ortonormalisasi Gram–Schmidt berkala (algoritma Benettin) untuk spektrum $(\lambda_1, \lambda_2, \lambda_3)$, uji kontraksi volume Liouville $\sum \lambda_i = -(\sigma + 1 + \beta)$, dimensi Kaplan–Yorke $D_{\mathrm{KY}}$, dan pertumbuhan pemisahan dua lintasan $\delta_0 = 10^{-9}$.
4. **Visualisasi SVG Mandiri**: 4 grafik SVG ditulis langsung dari C++ ke `outputs/figures/day007a_*.svg`.

## Cara Menjalankan

```bash
make test
make run
../../../.venv/bin/python py/analyze.py
```

Penjelasan teknis lengkap tersedia di [day007a_sim_chaos-lyapunov-cpp.md](../../../docs/project-explanations/day007a_sim_chaos-lyapunov-cpp.md).
