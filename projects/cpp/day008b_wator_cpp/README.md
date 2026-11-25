# `day008b_wator_cpp` — Simulasi Wa-Tor (Ikan vs Hiu) & Fit Lotka–Volterra (C++20)

Implementasi C++20 murni (`438` baris sumber di `src/`, tanpa pustaka eksternal) untuk model mangsa–pemangsa **Wa-Tor** pada grid torus $W \times H$ dengan aturan Dewdney, RNG deterministik `SplitMix64` $\to$ `xoshiro256++`, fit kuadrat terkecil ke model kontinu **Lotka–Volterra (RK4)**, studi efek ukuran hingga ($W \in \{20, 40, 80, 160\}$, 40 replika, `std::thread`), serta penulisan grafik SVG langsung dari C++.

## Perintah

```bash
make test
make run
.venv/bin/python py/analyze.py
```
