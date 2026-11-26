# Hari 9a (`day009a`) — Dinamika Molekuler Fluida Lennard-Jones (C11, 60 Tahap)

Implementasi murni C11 (`-std=c11 -O2 -Wall -Wextra -Werror -ffp-contract=off -fno-fast-math`, hanya `libc` + `libm` + `pthread`, target **6.000–9.000 baris** modular di `include/*.h` dan `src/*.c`) untuk simulasi dinamika molekuler (MD) klasik dalam satuan tereduksi ($\sigma = \varepsilon = m = k_B = 1$).

## Ringkasan Cabang & Fitur
- **Potensial (5)**: LJ terpotong-dan-digeser ($r_c=2.5$), LJ terpotong + koreksi ekor analitik ($U_{\text{tail}}, P_{\text{tail}}$), WCA repulsif murni ($r_c=2^{1/6}$), Morse, dan Yukawa; seluruhnya mendukung evaluasi analitik, tabel *lookup* (linear & kubik Hermite), serta verifikasi numerik $F = -dU/dr$.
- **Integrator (4)**: Velocity Verlet (utama, simplektik orde-2), Leapfrog, Euler eksplisit (kontrol negatif), dan r-RESPA dua skala waktu.
- **Ensemble & Termostat/Barostat (6)**: NVE, NVT Langevin BAOAB (momentum pusat massa nol), NVT Nosé–Hoover rantai-1, NVT Berendsen (non-kanonik), NVT Andersen pasangan (kekal momentum), dan NPT Berendsen sederhana.
- **Struktur Tetangga & Tangga Optimasi**: `V0` Brute-force $O(N^2)$, `V1` Cell list $O(N)$, `V2` Verlet list dengan kulit $r_{\text{skin}} \in \{0.2, 0.3, 0.4\}$ + pemicu bangun ulang teraudit, `V3` tata letak SoA, `V4` paralel `pthread` deterministik bit-identik (1/2/4/8 thread), `V5` tata letak terurut-sel SoA.
- **Sistem & Observabel**: 3D kubik (kisi fcc & acak non-tumpang-tindih), 2D segitiga ($N=1024$), fase padat/cair/gas; energi, suhu, tekanan virial, $g(r)$, $S(k)$, bilangan koordinasi, MSD/Einstein $D$, VACF/Green–Kubo $D$, tensor tegangan 5-komponen & viskositas geser $\eta$, parameter orde ikatan $Q_6$, parameter Lindemann $\delta_L$, dan kompresi lintasan delta+RLE *lossless*.
- **60 Tahap (`S01`–`S60`)**: menghasilkan >150 berkas keluaran (CSV ringkas di `out/`, ≥45 grafik SVG tulis tangan ke `../../../outputs/figures/day009a_*.svg`, serta `../../../outputs/day009a/results.json`).

## Penggunaan

```sh
make
make test          # selftest (satu-satunya suite tes)
make lines         # verifikasi 6.000-9.000 baris C
make run           # eksekusi skala penuh 60 tahap
make run SCALE=small  # uji cepat tanpa menimpa klaim resmi
./build/md threads 4  # uji determinisme 1/2/4/8 thread
```
