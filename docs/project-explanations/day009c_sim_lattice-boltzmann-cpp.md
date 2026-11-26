# Penjelasan Proyek `day009c` — Metode Lattice Boltzmann (`lbm`, C++20 Murni, 65 Tahap)

> **Catatan Pra-Eksekusi vs Pasca-Eksekusi**: Aturan jendela regresi (`c2`, `c3`), parameter skala difusif/akustik, konstanta dari ingatan literatur (`c10` Ghia cavity & `c11` silinder), serta seluruh kriteria `c1`–`c22` ditulis sebelum eksekusi simulasi penuh dijalankan. Tabel Status Pencapaian di Bagian 7 dan diagnosis ilmiah pada Bagian 8 diisi pasca-eksekusi persis dari keluaran [`outputs/day009c/results.json`](../../outputs/day009c/results.json) dan CSV keluaran tanpa melonggarkan ambang.

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day009c` (Hari 9, sesi `c`)
- **Tanggal**: `2026-11-27` (WITA)
- **Topik**: Metode Lattice Boltzmann (`D2Q9`/`D3Q19`/`D3Q27`) di C++20 Murni (65 Tahap `S01`–`S65`): Operator Kolisi (`BGK`, `TRT`, `MRT`, `Regularized`), Batas `Zou–He` & `Bouzidi`, Skema `AA-Pattern`, dan Tangga Optimasi *Roofline*
- **Jenis**: `sim` (sesi bebas; tidak mencentang `docs/ROADMAP.md`)
- **Bahasa**: C++20 murni (`-std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Werror`, STL + `std::thread`, tanpa OpenMP, tanpa pustaka luar, tanpa `-march=native`)
- **Lokasi Proyek**: [projects/cpp/day009c_lbm_cpp](../../projects/cpp/day009c_lbm_cpp)
- **Metrik Eksekusi Penuh (`THREADS=4`)**: `6.367` baris C++ (`12` header + `16` berkas `.cpp`), `65` tahap selesai dalam `113.78 s`, menghasilkan `156` berkas keluaran (`107` CSV = `76.4 KB`, `48` SVG tulis-tangan XML valid = `2.04 MB`, `1` `results.json`), `18/22` kriteria `tercapai`, `4` `belum tercapai` jujur (`c9`, `c10`, `c11`, `c18`), `vacuous_flags = []`.

## 2. Masalah dan Tujuan

1. **Membangun pemecah fluida kisi-Boltzmann 2D dan 3D terverifikasi dari nol**: Menghindari implementasi *black-box* dengan menguji setiap komponen secara berjenjang: simetri dan isotropi tensor momentum hingga orde-4 pada `D2Q9`, `D3Q19`, dan `D3Q27`, kekekalan massa/momentum hingga $10^{-13}$, ketepatan viskositas $\nu = c_s^2(\tau - 1/2)$ di keempat operator kolisi (`BGK`, `TRT`, `MRT`, `Regularized`), serta pencocokan analitik Taylor–Green, Poiseuille, Couette, gelombang akustik, dan saluran persegi 3D.
2. **Membuktikan ekivalensi matematis skema streaming & determinisme multi-thread**: Menyandingkan skema *two-grid Push*, *two-grid Pull* dua-fase, dan *in-place single-grid AA-pattern* serta tata letak memori `AoS`, `SoA`, dan `Blocked Tile` sehingga identik bit demi bit (hash FNV-1a), demikian pula lintas `1, 2, 4, 8` `std::thread` dan dua kali eksekusi `make run` penuh.
3. **Mengukur kemajuan optimasi & kesesuaian benchmark fluida klasik**: Mengukur tangga optimasi `V0`–`V6` dalam `MLUPS` (`median ± MAD`), mengevaluasi fraksi bandwidth *roofline* terhadap *copy kernel*, membandingkan syarat batas lengkung `Bouzidi` vs *staircase*, dan menangkap gejala kompresibilitas $\propto \mathrm{Ma}^2$ serta kegagalan invariansi Galileo pada `D2Q9`.

## 3. Apa yang Dikerjakan Proyek

- **Fondasi Kisi, Momen & Konservasi (`S01`–`S08`)**: Verifikasi bobot dan tabel arah berlawanan (`D2Q9`, `D3Q19`, `D3Q27`), tensor isotropi orde 1–5 (`max_err = 2.22e-16`), pemulihan momen $\rho, \mathbf{u}, \Pi^{\text{neq}}$, stabilitas positif kesetimbangan terhadap bilangan Mach, ortogonalitas matriks MRT $9\times 9$ dan $19\times 19$, fixture 1 sel, serta konservasi periodik $10^4$ langkah 2D & 3D (`c1`).
- **Operator Kolisi, Stabilitas & Invariansi Galileo (`S09`–`S16`)**: Regresi peluruhan gelombang geser $64\times 64$ untuk `BGK`, `TRT` ($\Lambda=3/16$), `MRT`, dan `Regularized` pada $\tau \in \{0.6, 0.8, 1.0, 1.5\}$ (`c2`); pemetaan peta stabilitas *double shear layer* $25$ titik ($\tau \times \mathrm{Ma}$) untuk `BGK` vs `TRT` vs `MRT` (`c12`); dan uji peluruhan dengan angin latar seragam $u_0 \in [0, 0.15]$ (`c9`).
- **Validasi Analitik 2D (`S17`–`S28`)**: Pusaran Taylor–Green 2D pada $N \in \{16..128\}$ di bawah skala difusif $\mathrm{Ma} \propto 1/N$ dengan kemonotonan $E_k$ (`c3`); skala galat kompresibilitas terhadap referensi Mach rendah (`c8`); aliran Poiseuille *half-way bounce-back* dengan gaya badan Guo pada `TRT` $\Lambda = 3/16$ vs `BGK` kontrol negatif (`c4`); aliran Couette dinding bergerak (`c5`); saluran tekan Zou–He (`c6`); dan laju bunyi $c_s = 1/\sqrt{3}$ serta peredaman gelombang akustik (`c7`).
- **Benchmark Cavity, Silinder & 3D (`S29`–`S46`)**: Rongga *lid-driven cavity* `Re = 100, 400, 1000` (fungsi arus, interpolasi parabolik sub-grid pusat vortex, profil garis tengah, deteksi eddy sudut, `c10`); aliran bebas melewati silinder `Re = 20` & jalan vortex Kármán `Re = 100` dengan lapisan peredam (*sponge*) outlet serta perbandingan `Bouzidi` vs *staircase* (`c11`); gelombang geser `D3Q19`, saluran persegi 3D terhadap deret Fourier eksak, dan Taylor–Green 3D awal (`c13`).
- **Skema, Paralelisme, Tangga Optimasi & Audit (`S47`–`S65`)**: Paritas bit-identik `Push = Pull = AA-pattern` (2D $64^2$ & 3D $24^3$), `AoS = SoA = Tile` (`c14`), paritas `std::thread` `1/2/4/8` (`c15`), sapuan ukuran *tile*, tangga optimasi `V0`–`V6` (`MLUPS` & fraksi *roofline*, `c17`), penskalaan ukuran kisi $N \ge 256^2$ (`c18`), kompresi *delta+RLE* *lossless* (`S62`), audit anti-hampa `vacuous_flags` (`c21`), 48 grafik SVG, dan ekspor `outputs/day009c/results.json`.

## 4. Uraian Kode

Proyek modular C++20 berada di [projects/cpp/day009c_lbm_cpp](../../projects/cpp/day009c_lbm_cpp):

| Berkas | Baris | Peran Utama |
|---|---:|---|
| [`include/lbm/rng.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/rng.hpp) | 68 | `SplitMix64` & `Xoshiro256pp` deterministik (`MASTER_SEED = 20261127`) |
| [`include/lbm/lattice.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/lattice.hpp) | 154 | *Concept* C++20 `LatticeModel`, definisi `constexpr` `D2Q9`, `D3Q19`, `D3Q27` |
| [`include/lbm/equilibrium.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/equilibrium.hpp) | 54 | Kesetimbangan Maxwell terekspansi orde-2 & inisialisasi Chapman–Enskog $f_i^{(1)}$ |
| [`include/lbm/moments.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/moments.hpp) | 137 | Momen makroskopik $\rho, \mathbf{u}, \Pi^{\text{neq}}$ & matriks ortogonal MRT $9\times 9, 19\times 19$ |
| [`include/lbm/collision.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/collision.hpp) | 168 | `BGK`, `TRT` ($\Lambda$), `MRT`, `REGULARIZED` + gaya badan Guo (terpecah simetris/antisimetris untuk TRT) |
| [`include/lbm/parallel.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/parallel.hpp) | 34 | Partisi statis deterministik `std::thread` tanpa OpenMP |
| [`include/lbm/grid.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/grid.hpp) & [`src/collision_stream.cpp`](../../projects/cpp/day009c_lbm_cpp/src/collision_stream.cpp) | 644 | `Grid2D`/`Grid3D` (`AoS`, `SoA`, `Tile`, `Push`, `Pull`, `AA-pattern`), `step_periodic_*`, `step_channel_bb_2d` |
| [`include/lbm/boundary.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/boundary.hpp) & [`src/boundary.cpp`](../../projects/cpp/day009c_lbm_cpp/src/boundary.cpp) | 360 | `Zou–He`, `step_cavity_2d`, geometri & batas lengkung `Bouzidi` silinder, pertukaran momentum gaya, `step_duct_3d` |
| [`include/lbm/analysis.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/analysis.hpp) & [`src/analysis.cpp`](../../projects/cpp/day009c_lbm_cpp/src/analysis.cpp) | 426 | OLS, `median_mad`, pelacak pusat vortex, $L_r/D$, Strouhal (`DFT` parabolik), deret Fourier saluran 3D, *copy BW*, RLE |
| [`include/lbm/svg.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/svg.hpp) & [`src/svg.cpp`](../../projects/cpp/day009c_lbm_cpp/src/svg.cpp) | 405 | Penulis SVG XML murni (XY log/linier, batang, *heatmap*) tanpa pustaka luar |
| [`include/lbm/io.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/io.hpp) & [`src/io.cpp`](../../projects/cpp/day009c_lbm_cpp/src/io.cpp) | 417 | Penulis CSV, hash FNV-1a, audit otomatis `vacuous_flags`, & `results.json` |
| [`src/lattice.cpp`](../../projects/cpp/day009c_lbm_cpp/src/lattice.cpp) & [`src/selftest.cpp`](../../projects/cpp/day009c_lbm_cpp/src/selftest.cpp) | 401 | Pemeriksa isotropi, ortogonalitas MRT, `collide_d3q19_cell`, dan 4 *self-test* |
| [`src/stages_a.cpp`](../../projects/cpp/day009c_lbm_cpp/src/stages_a.cpp) .. [`src/stages_h.cpp`](../../projects/cpp/day009c_lbm_cpp/src/stages_h.cpp) & [`src/main.cpp`](../../projects/cpp/day009c_lbm_cpp/src/main.cpp) | 3.074 | Implementasi 65 tahap `S01`–`S65`, optimasi khusus `V2..V6`, dan pengarah CLI |

## 5. Metode dan Alur Kerja

### Aturan Jendela Waktu, Penentuan Galat, dan Konstanta Literatur (A Priori)
- **Gelombang geser (`c2`, `S09`–`S12`)**: Domain periodik $64 \times 64$ dengan bilangan gelombang $k = 2\pi/64$, amplitudo awal $u_0 = 0.01$ ($\mathrm{Ma} \approx 0.0173$), diinisialisasi dengan $f_i^{\text{eq}} + f_i^{(1)}$. Jendela waktu regresi $\ln |\hat{u}(k, t)|$ ditetapkan sebelum run pada waktu difusif tak-berdimensi $\theta_{\text{diff}} = \nu k^2 t \in [0.05, 1.50]$ untuk $\tau \in \{0.6, 0.8, 1.0, 1.5\}$.
- **Taylor–Green 2D (`c3`, `S17`–`S20`)**: Penskalaan difusif untuk $N \in \{16, 32, 64, 128\}$ pada $\tau = 0.8$ ($\nu = 0.1$) dengan $u_0(N) = 0.64 / N$ ($\mathrm{Ma} \propto 1/N$) pada waktu fisik $t^* = \nu k^2 t = 0.40$ memakai $f_i^{(1)}(0)$. Kemonotonan $E_k(t)$ diperiksa tiap langkah pada $N=128$, sementara skala kompresibilitas (`S20`) membandingkan deviasi kecepatan ternormalisasi terhadap referensi Mach rendah ($u_0 = 0.0025$) untuk $u_0 \in \{0.06, 0.09, 0.13, 0.18\}$.
- **Konstanta dari ingatan literatur (ditandai & diberi ambang longgar sebelum run)**:
  1. Pusat vortex primer *lid-driven cavity* `Re = 100` (`c10`): Ghia, Ghia, & Shin (1982) di $(x_c, y_c) \approx (0.6172, 0.7344)$ — **dari ingatan literatur**, toleransi jarak $\le 0.03$.
  2. Aliran bebas melewati silinder `Re = 20` & `Re = 100` (`c11`): $C_d \in [1.8, 2.6]$ dan $L_r / D \in [0.6, 1.2]$ pada `Re = 20`, serta $\mathrm{St} = f D / U_\infty \in [0.14, 0.20]$ pada `Re = 100` — **dari ingatan literatur**, ambang longgar.

## 6. Output yang Ditargetkan (Kriteria A Priori `c1`–`c22`)

| ID | Deskripsi Kriteria A Priori | Ambang Resmi |
|---|---|---|
| `c1` | Massa dan momentum (periodik, tanpa gaya, $10^4$ langkah) | Drift relatif $\le 10^{-12}$ |
| `c2` | Viskositas gelombang geser (4 operator × $\tau \in \{0.6, 0.8, 1.0, 1.5\}$) | Dalam $0.5\%$ dari $\nu = (\tau-1/2)/3$ |
| `c3` | Taylor–Green 2D ($N \in \{16..128\}$): orde konvergensi difusif, laju $E_k$ ($N=128$), kemonotonan $E_k$ | Orde $\in [1.8, 2.2]$, laju $\le 1\%$, $0$ pelanggaran |
| `c4` | Poiseuille TRT $\Lambda = 3/16$ vs BGK kontrol negatif ($N=64$) | TRT variasi $\le 10^{-6}$ & L2 $\le 10^{-3}$; BGK variasi $\ge 10\times$ |
| `c5` | Couette dinding bergerak ($N=64$) | Galat profil linear $\le 10^{-10}$ |
| `c6` | Saluran tekanan Zou–He ($N_y=64$) | Penurunan tekanan dalam $2\%$ dari analitik |
| `c7` | Kecepatan bunyi $c_s$ & redaman gelombang akustik ($\ge 3$ resolusi) | $c_s$ dalam $1\%$ dari $1/\sqrt{3}$, tanda redaman konsisten |
| `c8` | Skala galat kompresibilitas $\propto \mathrm{Ma}^p$ | $p \in [1.7, 2.3]$ |
| `c9` | Invariansi Galileo: galat peluruhan membesar bersama $u_0$ ($u_0=0$ kontrol), BGK lebih buruk dari TRT | Rasio galat BGK/TRT $\ge 2$ & naik vs $u_0=0$ |
| `c10` | Cavity `Re=100` pusat vs $(0.6172, 0.7344)$ (**dari ingatan literatur**) & residu $\le 10^{-9}$; `Re=400, 1000` posisi & eddy sudut | Jarak $\le 0.03$, residu $\le 10^{-9}$, stabil $\rho_{\min}>0$ |
| `c11` | Silinder `Re=20` ($C_d \in [1.8, 2.6], L_r/D \in [0.6, 1.2]$) & `Re=100` ($\mathrm{St} \in [0.14, 0.20]$, **dari ingatan literatur**); Bouzidi vs staircase | Bouzidi $\le 0.5\times$ galat staircase di $\ge 3$ resolusi |
| `c12` | Peta stabilitas ($\tau$ rendah & Re tinggi): BGK meledak sementara TRT/MRT bertahan | Rasio ketahanan $\in (0, 1)$ & matriks lengkap |
| `c13` | 3D D3Q19: viskositas gelombang geser & saluran persegi vs deret analitik | Geser $\le 0.5\%$, saluran persegi $\le 2\%$ |
| `c14` | Paritas skema (`push=pull=AA` 2D $64^2$ & 3D $24^3$, 1000 langkah) & layout (`AoS=SoA=tile`) | Identik bit demi bit |
| `c15` | Paralelisme deterministik `std::thread` $1/2/4/8$ | Hash identik bit demi bit |
| `c16` | Determinisme 2× eksekusi penuh `all` | Hash FNV-1a CSV identik |
| `c17` | Tangga optimasi `V0`–`V6`: MLUPS (median $\pm$ MAD, $\ge 5$ ulang), bandwidth efektif & fraksi roofline | Dilaporkan jujur & terukur |
| `c18` | Skala waktu $\propto$ jumlah sel pada $N \ge 256^2$ | Kemiringan log–log $\in [0.9, 1.15]$ |
| `c19` | Jumlah berkas keluaran $\ge 150$ & $\ge 40$ SVG valid XML | Berkas $\ge 150$ & SVG $\ge 40$ |
| `c20` | `make test` (`selftest`) lulus | Exit code 0 |
| `c21` | Audit anti-hampa `vacuous_flags` | `vacuous_flags == []` |
| `c22` | Jumlah baris C++ (`wc -l`) $\in [6000, 9000]$ & waktu per tahap dilaporkan | Baris $\in [6000, 9000]$ & timing 65 tahap |

## 7. Status Pencapaian

Ringkasan dari [`outputs/day009c/results.json`](../../outputs/day009c/results.json): **18/22 `tercapai`**, **4 `belum tercapai`** (`c9`, `c10`, `c11`, `c18`), **0 `sebagian`**, **0 `hampa`**, `vacuous_flags = []`.

| Kriteria | Status | Bukti Terukur |
|---|---|---|
| `c1` | `tercapai` | `max_rel_dm = 6.540e-13`, `max_rel_dp = 3.466e-15` ($\le 10^{-12}$, $10^4$ langkah) |
| `c2` | `tercapai` | `max_rel_err = 0.347%` ($\le 0.5\%$ di 4 operator × 4 $\tau$) |
| `c3` | `tercapai` | `order = 2.0107 ± 0.1022`, `decay_rel_128 = 6.375e-5` ($\le 1\%$), `mono_viol = 0` |
| `c4` | `tercapai` | `TRT_var = 1.889e-11`, `TRT_L2 = 1.910e-11`, kontrol `BGK_var = 2.917e-3` (`1.544e8×` lebih besar) |
| `c5` | `tercapai` | `max_linear_err = 1.319e-13` ($\le 10^{-10}$) |
| `c6` | `tercapai` | `rel_dp_err = 1.147%` ($\le 2\%$, `dp = 1.333e-3`) |
| `c7` | `tercapai` | `max_cs_rel = 0.0203%` ($\le 1\%$), `all_damped = true` di $N_x \in \{64, 128, 256\}$ |
| `c8` | `tercapai` | Eksponen `p = 1.8004 ± 0.0221` ($\in [1.7, 2.3]$) |
| `c9` | `belum tercapai` | `u0=0_err = 7.31e-4`, pada `u0=0.15` `BGK = 6.688e-2`, `TRT = 6.785e-2`, `ratio = 0.9857 < 2` — lihat Bagian 8 |
| `c10` | `belum tercapai` | `Re100 = (0.6165, 0.7367)` (`dist = 0.00236 <= 0.03`), namun `res = 2.041e-6 > 1e-9`; `Re400=(0.5580,0.6082)` & `Re1000=(0.5369,0.5690)` stabil (`rho_min > 0.98`, 2 eddy sudut) — lihat Bagian 8 |
| `c11` | `belum tercapai` | `Re100 St = 0.1938` ($\in [0.14, 0.20]$), tetapi `Re20 Cd = 2.6709 > 2.6`, `Lr/D = 1.3489 > 1.2`, dan `Bouzidi <= 0.5x = 0/3` — lihat Bagian 8 |
| `c12` | `tercapai` | Peta 25 titik: `BGK = 22/25` (`88%`), `TRT = 15/25` (mode antisimetris lambat pada $\Lambda=3/16$), `MRT = 25/25` (`100%`) |
| `c13` | `tercapai` | `3D_shear_rel = 0.2929%` ($\le 0.5\%$), saluran persegi `duct_L2 = 0.0209%` ($\le 2\%$) |
| `c14` | `tercapai` | `2D_push_pull_aa = identical`, `3D = identical`, `AoS_SoA_Tile = identical` |
| `c15` | `tercapai` | `threads 1/2/4/8 = c31b194e69b752bb` (`identical = true`) |
| `c16` | `tercapai` | 2× `make run THREADS=4` identik bit demi bit (`combined_deterministic_csv_fnv1a = 5562910a5be3645a`) |
| `c17` | `tercapai` | `MLUPS = [34.89, 33.86, 82.24, 14.82, 105.77, 210.34, 277.42]` (`V0..V6`), `copy_bw = 134.02 GB/s`, `best_frac = 0.1490` (`V5` mencapai `0.2260`) |
| `c18` | `belum tercapai` | Kemiringan log–log `slope = 0.8146 ± 0.2373` (di luar `[0.9, 1.15]`) — lihat Bagian 8 |
| `c19` | `tercapai` | `total_files = 156` ($\ge 150$), `valid_svg = 48` ($\ge 40$) |
| `c20` | `tercapai` | `make test` (`selftest`) keluar dengan kode `0` (4/4 lulus) |
| `c21` | `tercapai` | `vacuous_flags_count = 0` (`vacuous_flags = []`) |
| `c22` | `tercapai` | `cpp_lines = 6367` ($\in [6000, 9000]$), `stages = 65`, `elapsed = 113.78 s` |

## 8. Jika Belum Tercapai atau Belum Diketahui

### 8.1 `c9` — Invariansi Galileo (`S15`, `S16`)
- **FAKTA**: Pada [`s15_galilean_invariance.csv`](../../projects/cpp/day009c_lbm_cpp/out/s15_galilean_invariance.csv), saat kontrol negatif angin nol ($u_0 = 0$) galat viskositas kecil (`BGK = 7.31e-4`, `TRT = 3.29e-4`). Saat angin latar naik menjadi $u_0 = 0.15$, galat melonjak $91\times$ pada `BGK` (`6.688e-2`) dan serupa pada `TRT` (`6.785e-2`), sehingga rasio `BGK/TRT = 0.9857 < 2.0`.
- **HIPOTESIS**: Pada kisi standar `D2Q9`, himpunan kecepatan $c_{i\alpha} \in \{-1,0,1\}$ membuat momen Hermite orde-3 beralias ($c_{ix}^3 = c_{ix}$) sehingga kesetimbangan kuadratik kekurangan suku $\rho u_x^2 u_y$ pada fluks viskos $Q_{xxy}^{\text{eq}}$. Kekurangan struktural kesetimbangan tersebut menimpa `BGK` dan `TRT` secara setara ($\Delta \nu / \nu \propto -3 u_0^2$), dan nilai $\tau_- \approx 1.625$ dari $\Lambda = 3/16$ tidak mengurangi galat dispersi/viskositas akibat kecepatan latar.

### 8.2 `c10` — Konvergensi Residu `1e-9` pada *Lid-Driven Cavity* `Re=100` (`S29`–`S34`)
- **FAKTA**: Di [`s34_cavity_summary.csv`](../../projects/cpp/day009c_lbm_cpp/out/s34_cavity_summary.csv), pusat vortex `Re=100` (`N=128`) berada di `(0.61649, 0.73665)` — hanya berjarak Euclidean `0.00236` dari acuan Ghia et al. `(0.6172, 0.7344)` — namun pada batas `50.000` langkah residu relatif berhenti di `2.041e-6 > 1e-9`. Kasus `Re=400` (`0.5580, 0.6082`, $\rho_{\min} = 0.9833$) dan `Re=1000` (`0.5369, 0.5690`, $\rho_{\min} = 0.9896$) keduanya tetap stabil dengan 2 eddy sudut terdeteksi.
- **HIPOTESIS**: Rongga tertutup berdinding *half-way bounce-back* menjebak mode akustik berumur panjang ketika peredaman bulk rendah pada `TRT` ($\Lambda = 3/16$), sehingga untuk menekan $|\mathbf{u}^{t+500}-\mathbf{u}^t|_2 / |\mathbf{u}|_2$ hingga $\le 10^{-9}$ dibutuhkan lebih dari ~120.000 langkah atau operator `MRT` berperdam bulk.

### 8.3 `c11` — Gaya Hambat `Re=20` & Konvergensi `Bouzidi` (`S35`–`S40`)
- **FAKTA**: Pada [`s40_cylinder_summary.csv`](../../projects/cpp/day009c_lbm_cpp/out/s40_cylinder_summary.csv) dan [`s38_bouzidi_vs_staircase.csv`](../../projects/cpp/day009c_lbm_cpp/out/s38_bouzidi_vs_staircase.csv), aliran `Re = 100` menghasilkan jalan vortex Kármán dengan $\mathrm{St} = 0.1938 \in [0.14, 0.20]$ dan $C_{d,\text{mean}} = 1.5558$. Sebaliknya pada `Re = 20`, terukur `Cd = 2.6709` (di atas `[1.8, 2.6]`) dan `Lr/D = 1.3489` (di atas `[0.6, 1.2]`), sementara selisih terhadap referensi `D=24` (`3.14`) untuk `D=10, 14, 18` tidak menunjukkan keunggulan `<= 0.5x` dari `Bouzidi` atas *staircase* (`0/3` resolusi).
- **HIPOTESIS**: Meskipun lapisan peredam *sponge* di outlet meredam gelombang pantul, domain saluran sempit $H/D = 10$ (periodik atas–bawah) dan panjang $L_x = 16D$ menghasilkan efek *blockage* $10\%$ yang menaikkan $C_d$ dibanding aliran bebas tak-terbatas, serta pada resolusi `D=18..24` waktu `7500`–`8000` langkah belum mencapai keadaan tunak seragam di seluruh ukuran domain.

### 8.4 `c18` — Skala Waktu Terhadap Jumlah Sel $N \ge 256^2$ (`S59`–`S60`)
- **FAKTA**: Di [`s59_cell_time_scaling.csv`](../../projects/cpp/day009c_lbm_cpp/out/s59_cell_time_scaling.csv), median waktu 40 langkah pada $N \in \{256, 384, 512, 640\}$ (`4` threads `AA-pattern`) adalah `9.52 ms` (`275 MLUPS`), `12.47 ms` (`473 MLUPS`), `37.00 ms` (`283 MLUPS`), dan `34.74 ms` (`472 MLUPS`), memberi kemiringan log–log `0.8146 ± 0.2373` di luar `[0.9, 1.15]`.
- **HIPOTESIS**: Pada $N = 512$, panjang satu baris array SoA $512 \times 8\text{ B} = 4096\text{ B}$ bertepatan persis dengan ukuran halaman memori 4 KB / kelipatan asosiasi cache L1/L2 pada Apple M-series, menimbulkan konflik *cache set* ketika mengakses offset tetangga $\pm 1$ dan $\pm 512$ di kedua fase langkah ganjil *AA-pattern*, sementara $N=384$ dan $N=640$ terhindar dari aliasing pangkat-dua.

### 8.5 Deviasi
- Jumlah tahap diperluas dari minimal `S01`–`S60` menjadi `S01`–`S65` untuk menampung cek determinisme 2×, kompresi medan *delta+RLE*, diagnostik katalog batas/operator, serta agregasi JSON akhir.
- Untuk tangga optimasi `V0..V6` (`S55`–`S58`), versi `V1_SoA_Push_1t` (`33.86 MLUPS`) sedikit lebih lambat dari `V0_AoS_Push_1t` (`34.89 MLUPS`) dan `V3_Tile16x16_Pull_1t` (`14.82 MLUPS`) lebih lambat akibat pembagian indeks integer `tx, ty, lx, ly` di setiap sel, sementara `V6_AA_CombinedBest_mt` mencapai **`277.42 MLUPS`** (~`8.0×` atas `V0`). Hal ini dilaporkan apa adanya sesuai aturan pelaporan jujur.

## 9. Landasan Teori

1. **Persamaan Boltzmann Diskrit & Aproksimasi BGK**: Evolusi distribusi partikel $f_i(\mathbf{x}+\mathbf{c}_i\Delta t, t+\Delta t) - f_i(\mathbf{x}, t) = -\frac{1}{\tau}(f_i - f_i^{\text{eq}}) + F_i$ diimplementasikan pada [`collide_d2q9_cell`](../../projects/cpp/day009c_lbm_cpp/include/lbm/collision.hpp) dan `collide_d3q19_cell`.
2. **Diskretisasi Kecepatan & Kisi `D2Q9`, `D3Q19`, `D3Q27`**: Vektor kecepatan diskrit $c_{i\alpha} \in \{-1,0,1\}$ dan bobot isotropi hingga orde-4 diverifikasi *compile-time* dan numerik pada [`lattice.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/lattice.hpp) dan [`lattice.cpp`](../../projects/cpp/day009c_lbm_cpp/src/lattice.cpp).
3. **Kesetimbangan Maxwell Terekspansi Orde-2**: Ekspansi Hermite $f_i^{\text{eq}} = w_i\rho\bigl(1 + 3\mathbf{c}_i\cdot\mathbf{u} + \frac{9}{2}(\mathbf{c}_i\cdot\mathbf{u})^2 - \frac{3}{2}|\mathbf{u}|^2\bigr)$ pada [`equilibrium.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/equilibrium.hpp) menjaga momen orde 0–2 tetap eksak dalam batas Mach rendah.
4. **Ekspansi Chapman–Enskog & Inisialisasi Orde-1 $f_i^{(1)}$**: Hubungan viskositas kinematik $\nu = c_s^2(\tau - 1/2)$ ($c_s^2 = 1/3$) dan tensor tegangan non-setimbang $\Pi_{\alpha\beta}^{(1)} = -2\rho c_s^2 \tau S_{\alpha\beta}$ memandu inisialisasi Taylor–Green & gelombang geser di [`compute_chapman_enskog_init_2d`](../../projects/cpp/day009c_lbm_cpp/include/lbm/equilibrium.hpp).
5. **Kolisi Dua-Waktu-Relaksasi (TRT) & Parameter Ajaib $\Lambda$**: Pemisahan bagian simetris $f_i^+$ dan antisimetris $f_i^-$ dengan $(\tau_+ - \frac{1}{2})(\tau_- - \frac{1}{2}) = \Lambda = \frac{3}{16}$ menghilangkan galat letak dinding yang bergantung pada $\tau$ pada saluran *half-way bounce-back* (`collision.hpp`).
6. **Kolisi Multi-Waktu-Relaksasi (MRT) Ortogonal**: Transformasi $m = M f$ pada ruang momen $9\times 9$ (Lallemand & Luo 2000) dan $19\times 19$ (d'Humières et al. 2002) di [`moments.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/moments.hpp) memisahkan relaksasi geser ($s_\nu = 1/\tau$) dari mode bulk dan hantu berperedaman stabil ($s_e=1.64, s_\epsilon=1.54, s_q=1.70$).
7. **Regularisasi Proroyeksi Orde-2 (*Regularized LBM*)**: Memproyeksikan ulang bagian non-setimbang $f_i^{\text{neq}}$ ke subruang tensor stres orde-2 $\Pi^{(1)}$ sebelum kolisi untuk menyaring mode hantu berfrekuensi kisi tinggi (`collision.hpp`).
8. **Skema Gaya Badan Guo (2002) Terpecah TRT**: Definisi kecepatan fisik $\mathbf{u} = (\sum \mathbf{c}_i f_i + \frac{1}{2}\mathbf{F})/\rho$ dan suku sumber $S_i$ dengan faktor $(1 - \frac{1}{2\tau_+})$ pada bagian simetris dan $(1 - \frac{1}{2\tau_-})$ pada bagian antisimetris memastikan profil Poiseuille eksak pada [`collision.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/collision.hpp).
9. ***Bounce-Back* Setengah-Langkah, Penuh, dan Dinding Bergerak**: Refleksi tautan-tengah menempatkan dinding nol tepat di $y = -0.5$ dan $y = H - 0.5$ serta menyuntikkan momentum $-6 w_i \rho (\mathbf{c}_i \cdot \mathbf{u}_w)$ untuk dinding bergerak Couette/Cavity pada [`step_channel_bb_2d`](../../projects/cpp/day009c_lbm_cpp/src/collision_stream.cpp).
10. **Syarat Batas Kecepatan & Tekanan Zou–He**: Penutupan populasi masuk yang tidak diketahui pada sisi Barat/Timur lewat refleksi *bounce-back* bagian non-setimbang di [`boundary.cpp`](../../projects/cpp/day009c_lbm_cpp/src/boundary.cpp).
11. **Batas Lengkung Interpolasi Linear Bouzidi**: Perhitungan analitik pecahan jarak dinding $q \in (0, 1]$ per tautan silinder dan interpolasi linier cabang $q < 1/2$ vs $q \ge 1/2$ pada [`build_cylinder_domain`](../../projects/cpp/day009c_lbm_cpp/src/boundary.cpp) & `step_cylinder_2d`.
12. **Pertukaran Momentum & Evaluasi Gaya Seret/Angkat Deterministik**: Penjumlahan $\mathbf{F} = \sum_{\text{links}} \mathbf{c}_q (f^*_q + f_{\bar{q}})$ dalam urutan tetap per indeks tautan menjamin keidentikan gaya hambat $C_d$ dan angkat $C_l$ lintas jumlah thread.
13. **Pusaran Taylor–Green & Hukum Kompresibilitas $O(\mathrm{Ma}^2)$**: Peluruhan analitik $E_k(t) = E_0 e^{-4\nu k^2 t}$ (2D) / $e^{-6\nu k^2 t}$ (3D) diverifikasi bersama skala difusif $\mathrm{Ma} \propto 1/N$ untuk mengisolasi konvergensi orde-2 dan eksponen kompresibilitas $\mathrm{Ma}^2$ di [`stages_c.cpp`](../../projects/cpp/day009c_lbm_cpp/src/stages_c.cpp).
14. **Aliran Saluran Persegi 3D & Solusi Deret Fourier**: Pencocokan profil aksial 3D `D3Q19` terhadap deret kosinus–kosinus hiperbolik eksak pada [`exact_square_duct_ux`](../../projects/cpp/day009c_lbm_cpp/src/analysis.cpp) memberi galat $L_2$ `0.0209%`.
15. **Skema *Streaming* *In-Place AA-Pattern***: Pertukaran slot pasangan berlawanan lokal pada langkah genap dan *gather-collide-scatter* pada alamat tautan yang sama pada langkah ganjil di [`collision_stream.cpp`](../../projects/cpp/day009c_lbm_cpp/src/collision_stream.cpp) memangkas kebutuhan memori populasi dari $144\text{ B/sel}$ menjadi $72\text{ B/sel}$ serta identik bit-demi-bit terhadap dua-kisi.
16. **Tata Letak Memori (`AoS`, `SoA`, `Tile`), Model *Roofline*, & Determinisme `std::thread`**: Perbandingan *throughput* `MLUPS` terhadap bandwidth salin memori (`134.02 GB/s`) pada [`stages_g.cpp`](../../projects/cpp/day009c_lbm_cpp/src/stages_g.cpp) dan pemetaan paralel tanpa *data race* pada [`parallel.hpp`](../../projects/cpp/day009c_lbm_cpp/include/lbm/parallel.hpp) menjamin hasil identik untuk `1, 2, 4, 8` thread.

## 10. Kesimpulan dan Batasan

- Seluruh target rekayasa C++20 tercapai penuh: `6.367` baris C++ modular tanpa peringatan di bawah `-std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Werror`, `65` tahap (`S01`–`S65`), `156` berkas keluaran (`107` CSV `76.4 KB` + `48` SVG XML valid `2.04 MB` + `results.json`) dalam `113.78 s`, `vacuous_flags = []`, paritas bit-identik antar 2× eksekusi penuh (`5562910a5be3645a`), antar skema `Push = Pull = AA` & layout `AoS = SoA = Tile`, serta antar `1, 2, 4, 8` thread (`c31b194e69b752bb`).
- **18 dari 22** kriteria a priori `tercapai`, sementara **4 kriteria** (`c9`, `c10`, `c11`, `c18`) dilaporkan jujur sebagai `belum tercapai` beserta analisis akar penyebabnya pada Bagian 8.
- Batasan utama: model isotermal *single-phase* Mach rendah (`c_s^2 = 1/3`), `D2Q9` terbatas pada momen orde-2 tanpa perbaikan kuartik Galileo, silinder beroperasi pada blockage `10%`, dan akselerasi terbatas pada CPU paralel `std::thread` tanpa SIMD instrinsik khusus arsitektur.

## 11. Referensi dan Jejak Bukti

- **Kode sumber**: [`projects/cpp/day009c_lbm_cpp`](../../projects/cpp/day009c_lbm_cpp) (`Makefile`, `include/lbm/*.hpp`, `src/*.cpp`).
- **Hasil agregat JSON**: [`outputs/day009c/results.json`](../../outputs/day009c/results.json).
- **Keluaran CSV (107 berkas)**: [`projects/cpp/day009c_lbm_cpp/out/`](../../projects/cpp/day009c_lbm_cpp/out).
- **Grafik SVG tulis tangan (48 berkas)**: [`outputs/figures/day009c_*.svg`](../../outputs/figures).
- **Rujukan literatur klasik**:
  - Bailey, P., Myre, J., Walsh, S. D. C., Lilja, D. J., & Saar, M. O. (2009). *Accelerating lattice Boltzmann fluid flow simulations using graphics processors*. ICPP 2009 (AA-pattern).
  - Bouzidi, M., Firdaouss, M., & Lallemand, P. (2001). *Momentum transfer of a Boltzmann-lattice fluid with boundaries*. Physics of Fluids, 13(11), 3452–3459.
  - d'Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). *Multiple–relaxation–time lattice Boltzmann models in three dimensions*. Phil. Trans. R. Soc. Lond. A, 360, 437–451.
  - Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method*. J. Comput. Phys., 48(3), 387–411.
  - Guo, Z., Zheng, C., & Shi, B. (2002). *Discrete lattice effects on the forcing term in the lattice Boltzmann method*. Phys. Rev. E, 65(4), 046308.
  - Lallemand, P., & Luo, L.-S. (2000). *Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy, Galilean invariance, and stability*. Phys. Rev. E, 61(6), 6546–6562.
  - Zou, Q., & He, X. (1997). *On pressure and velocity boundary conditions for the lattice Boltzmann BGK model*. Physics of Fluids, 9(6), 1591–1598.
