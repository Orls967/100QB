# Day 009c — Metode Lattice Boltzmann (`lbm`) di C++20 Murni

Simulasi dinamika fluida kisi-Boltzmann (D2Q9, D3Q19, serta D3Q27 pembanding) tanpa pustaka eksternal menggunakan C++20 murni dan `std::thread` (`-std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Werror`).

## Perintah
- `make` — kompilasi bersih tanpa peringatan.
- `make test` — jalankan `selftest` (bobot/isotropi, kesetimbangan 1 sel, kekekalan 100 langkah, gelombang geser `N=16`).
- `make run` — jalankan seluruh 65 tahap (`S01`–`S65`), tulis CSV ringkas ke `out/`, SVG ke `outputs/figures/day009c_*.svg`, dan `outputs/day009c/results.json`.
- `make run SCALE=small` — uji cepat sebelum eksekusi penuh.
- `make lines` — verifikasi jumlah baris C++ berada pada `6.000–9.000`.
