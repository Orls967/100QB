# Hari 8i (`day008i`) — Tumpukan Pasir Abelian (C11, ≤ 999 Baris)

Implementasi murni C11 (`-std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off`, hanya `libc` + `libm`, total `Makefile + src/*.c + src/*.h ≤ 999` baris) untuk model tumpukan pasir Abelian Bak–Tang–Wiesenfeld (BTW):

- `selftest`: verifikasi determinan tangan $n=2$ (`192`), identitas $e$ pada seluruh konfigurasi rekuren $n=3$ (`116 640`), dan skenario runtuh tangan $3\times 3$.
- `abelian`: uji sifat Abelian pada konfigurasi acak ($n \in \{8, 16, 24\}$, 100 urutan runtuh acak per konfigurasi vs urutan deterministik, hash FNV-1a 64 atas tinggi + odometer + sink).
- `group`: algoritma pembakaran Dhar, orde grup $|G| = \det(\tilde{\Delta})$ lewat eliminasi Bareiss `__int128` dengan deteksi luapan ($n=2..6$) vs enumerasi brute-force ($n=2,3$), elemen identitas $e = s(2c - s(2c))$ & simetri $D_4$, serta orde elemen acak ($n \le 4$).
- `radius`: penskalaan jari-jari dukungan $R(N)$ untuk $N \in \{10^3, 10^4, 10^5, 10^6\}$ di pusat kisi & kekekalan massa eksak.
- `soc`: kritikalitas terorganisir-sendiri pada kisi $64\times 64$ dan $128\times 128$ ($\ge 2\cdot 10^6$ longsoran pasca-pemanasan).
- `heights`: distribusi tinggi stasioner di $\tfrac{1}{4}$ tengah kisi $128\times 128$ vs konstanta literatur Priezzhev/Majumdar–Dhar.
- `pics`: 6 grafik SVG tulis tangan ke `../../outputs/figures/day008i_*.svg`.
- `all`: menjalankan seluruh tahap dan menulis `out/*.csv` serta `../../outputs/day008i/results.json`.

## Penggunaan

```sh
make
make test
make lines
make run
```
