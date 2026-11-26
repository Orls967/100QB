# Changelog

Semua perubahan penting pada project **100QB** dicatat di sini, per hari.
Format mengikuti [Keep a Changelog](https://keepachangelog.com/id/1.1.0/).

## [README & Skrip] — 2026-11-27 · Nomor Sesi Berurutan & Baseline Fleksibel

### Diubah
- **`scripts/update_readme.py` & `README.md`**: mengganti kolom `Topik No` menjadi `No. Sesi` berurutan (`001`, `002`, ...) dengan `BASELINE_TOTAL = 100` yang fleksibel serta penghitungan jenis otomatis via `collections.Counter`; menghapus label `bonus` dari keluaran tabel sambil tetap menerima nilai lama pada `meta.json`; menambahkan bagian **"🎯 Tentang nama 100QB dan target"** serta memperbarui kalimat pembuka dan lencana `baseline-100 sesi`.

## [Hari 009b] — 2026-11-27 · Genetika Populasi (Ruby Murni, 65 Tahap)

### Ditambahkan (`day009b`)
- Proyek Ruby stdlib murni `projects/ruby/day009b_popgen_ruby/` (tanpa gem eksternal, lolos bersih `ruby -w` di bawah `ruby 4.0.9`, total **5.156 baris Ruby** di `bin/popgen` + 24 berkas `lib/popgen/*.rb`, **20/23** kriteria `tercapai`, **3** `belum tercapai` jujur (`c8`, `c9`, `c14`), `vacuous_flags = []`):
  - **Sampler Eksak & Rantai Markov `Rational` (`S01`–`S10`)**: SplitMix64 → `Random.new(seed)` deterministik; Binomial dua cabang (inversi $np<10$ + BTPE Kachitvichyanukul–Schmeiser 1988 $np\ge 10$) tervalidasi di 12 pasangan $(n,p)$ (`min chi2 p = 0.0803`, `c18`), Multinomial/Poisson Knuth+PTRS/Eksponensial; matriks transisi WF `Rational` vs `Float` cocok hingga `1.776e-13` (`c1`).
  - **Dinamika Forward WF, Moran & Seleksi (`S11`–`S30`)**: Peluruhan heterozigositas $N \in \{50, 200\}$ (`max z = 1.55`, `c4`); fiksasi netral (`max z = 1.19`, `c2`), waktu fiksasi vs rantai (`max z = 1.38`) & difusi (`8.87%`, `c5`); proses Moran tertutup `1.11e-15`, arah $r>1/r<1$ (`>= 35.4 SE`), kontrol `r=1` (`c3`); seleksi $N=2000$ (`max z = 2.74`, `c11`); stasioner $N=20$ mutasi+seleksi (`TV = 0.01426 <= 0.02932`, `c16`).
  - **Koalesen Kingman, Pulau, Sapuan, LD, ARG-Lite & Inbreeding (`S31`–`S65`)**: $T_{\text{MRCA}}, L, S, \mathrm{Var}(S)$ untuk $n \in \{5,20,100\}$ (`max z = 1.95`, `c6`); SFS Wald $\chi^2$ `p = 0.2875` (`c7`); Tajima's $D$ (`mean_D = -0.1064`, `Var(D) = 0.7949`, `c8` `belum tercapai`); partisi Ewens $n=6$ `p = 0.8164` & $\mathbb{E}[K_n]$ (`max z = 2.11`), namun $\mathbb{E}[H]$ pada $\theta=0.5$ `42.87%` (`c9` `belum tercapai`); model pulau $F_{ST}$ (`min step = 12.4 SE`, kontrol $M=80$ `-0.0001`, `c10`); sapuan selektif (`43.9 SE`, pemulihan `>= 7.9 SE`, `c12`); peluruhan LD & $\sigma_d^2$ vs Ohta–Kimura (`20.11% <= 25%`, `c13`); ARG-lite marginal `max z = 1.49` & ekor $\rho=8 \to 32$ `2.7 SE` (`c14` `belum tercapai`); HWE (`min p = 0.2014`) & *selfing* (`z = 0.45`, `c15`).
  - **Determinisme & Optimasi (`S62`–`S65`)**: 2× `rake run` identik `711328ecfbfdb734` dan `workers = 1, 2, 4` identik bit demi bit `8e12879dfcfa8dc0` (`c17`); tangga optimasi `V0..V4` (`[1.50, 0.05, 0.87, 2.09, 0.87, 2.07, 2.50]` ms, `c23`), kompresi delta+RLE *lossless* hingga `7.50x`, `109` CSV + `45` SVG valid XML `outputs/figures/day009b_*.svg` + `outputs/day009b/results.json` (`155` berkas), serta [day009b_sim_popgen-ruby.md](project-explanations/day009b_sim_popgen-ruby.md).

## [Hari 009a] — 2026-11-27 · Dinamika Molekuler Fluida Lennard-Jones (C11, 60 Tahap)

### Ditambahkan (`day009a`)
- Proyek C11 murni `projects/c/day009a_md_lj_c/` (`-std=c11 -O2 -Wall -Wextra -Werror -ffp-contract=off -fno-fast-math`, hanya `libc` + `libm` + `pthread`, total **8.818 baris C** modular di 13 header `include/*.h` + 21 berkas `src/*.c`, **22/23** kriteria `tercapai`, **1** `belum tercapai` (`c8` sebagai temuan ilmiah jujur), `vacuous_flags = []`):
  - **Fondasi, 5 Potensial & Daftar Tetangga (`S01`–`S18`)**: LJ-shift, LJ+tail, WCA, Morse, dan Yukawa cocok terhadap beda hingga $-dU/dr$ hingga `3.245e-10` (`c2`); kesetaraan gaya brute = cell = Verlet = SoA = threaded = cell-sorted hingga `1.869e-15` (`c1`); skala waktu $N \ge 500$ memberi slope brute `1.801 ± 0.019`, cell `0.973 ± 0.180`, Verlet `1.000 ± 0.011` (`c17`); tangga optimasi `V0`–`V5` (`[1142, 956, 264, 296, 165, 213] µs`, `c18`); `0` luapan daftar dan `0` bangun-ulang Verlet terlambat (`c23`).
  - **Integrator Simplektik & Termostat (`S19`–`S34`)**: reversibilitas Velocity Verlet `6.072e-14` (`c4`), orde $\Delta t$ `1.8708 ± 0.0246` vs Euler `106.740x` lebih buruk (`c5`), kekekalan momentum `6.354e-13` setelah $10^5$ langkah NVE (`c3`); ekuipartisi `max z = 2.33` (`c6`), KS Maxwell `min p = 0.0941` (`c7`); `Var(K)` NVE (`40.9 SE`) dan Berendsen (`99.6 SE`) jauh di bawah kanonik sementara Nosé–Hoover rantai-1 terukur `0.590 ± 0.047` (`c8` `belum tercapai`, temuan non-ergodisitas $M=1$).
  - **Termodinamika, Transport & Pelelehan (`S35`–`S60`)**: isoterm $T=1.5$ monoton naik (`min_step_z = 9.21 SE`) & $T=2.0$ tanpa loop vdW (`c15`); densitas rendah cocok dengan Mayer $B_2(T)\rho$ dalam `0.92 SE` (`c9`); jalur energi $g(r)$ cocok `0.043%` (`c10`); difusi Einstein vs Green–Kubo selisih `1.90%` dan $D(\rho)$ turun monoton (`min_step_z = 22.74 SE`, `c11`); $g(r)$ cair `|g(tail)-1| = 0.0165`, `g(r<0.8)=0` (`c12`); viskositas geser titik tripel $\eta = 3.2379 \pm 0.2997$ (`c13`); pelelehan fcc $Q_6$: `0.5086 → 0.0526` (`256.1 SE`) vs padat dingin `0.5301` (`c14`); determinisme 2× `make run` identik `c94cdec2512e7713` dan 1/2/4/8 thread identik bitwise `ea58c86a4fb8e20f` (`c16`); kompresi delta+RLE *lossless* hingga `7.60x`.
  - **Output**: `117` CSV ringkas (`205 KB`) + `49` SVG tulis tangan valid XML `outputs/figures/day009a_*.svg` + `outputs/day009a/results.json` (`167` berkas total, `c19`), dan [day009a_sim_md-lj-c.md](project-explanations/day009a_sim_md-lj-c.md).

## [Hari 008j] — 2026-11-26 · Diagonalisasi Eksak Rantai Spin-½ Heisenberg/XXZ (Rust ≤ 499 Baris): Lanczos, Keterikatan & Penskalaan CFT

### Ditambahkan (`day008j`)
- Proyek Rust 2021 murni `projects/rust/day008j_heisenberg_ed_rust/` tanpa *crate* eksternal dengan total **494 baris** pasca-`cargo fmt` (`Cargo.toml`: 7, `src/main.rs`: 487; $\le 499$ baris, **12/13** kriteria `tercapai` dan `c7_marshall` + sub-cek `sum_r` bertanda `belum dijalankan` per aturan pangkas batas baris, `vacuous_flags = []`):
  - **Lanczos *On-the-Fly* & Validasi Silang**: konvergensi hingga $N=24$ ($\dim S^z=0 = 2\,704\,156$) dalam maksimal `92` iterasi; selisih maksimum terhadap diagonalisasi Jacobi dense untuk semua $N \le 12, \Delta \in \{0,0.5,1,2\}$ sebesar `4.09e-14` (`c1`); selisih terhadap solusi eksak fermion bebas Jordan–Wigner ($\Delta=0, N=4..24$) sebesar `7.99e-15` (`c2`).
  - **Identitas & CFT**: $|\langle S_{\text{tot}}^2 \rangle| \le 4.47\times 10^{-11}$, selisih $\langle S_0\cdot S_1\rangle - E_0/N \le 7.25\times 10^{-12}$, dan kemonotonan $E_0(\Delta)$ terjaga di seluruh $N$ (`c3`); fit OLS $E_0/N = e_\infty + c_2/N^2$ ($N=10..24, \Delta=1$) memberi $e_\infty = -0.443112 \pm 0.000007$ dan $c_2 = -0.8420 \pm 0.0013$ (`c4`), sementara celah singlet–triplet meluruh dengan kemiringan log–log `-0.9576 ± 0.0012` pada $\Delta=1$ vs tercelah `0.6100` pada $\Delta=2, N=20$ (`c5`).
  - **Entropi Keterikatan & Tolok Ukur**: dekomposisi Schmidt per blok $S^z_A$ pada $N=20, \Delta=1$ menghasilkan kemiringan Calabrese–Cardy `0.3533 ± 0.0030` (galat jejak $\sum\lambda_i - 1 \le 1.32\times 10^{-13}$, `c6`); kecepatan `871.69 ± 2.32` iter/s (`11.22M` state/detik di $N=16$, `c12`); estimasi puncak memori `1908.4 MiB` (`c13`).
  - **Output**: 5 CSV di `out/` + 6 SVG tulis tangan `outputs/figures/day008j_*.svg` + `outputs/day008j/results.json` (identik 0-byte antar dua eksekusi `make run`, `csv_fnv1a = b21b281060c3f526`), dan [day008j_sim_heisenberg-ed-rust.md](project-explanations/day008j_sim_heisenberg-ed-rust.md).

## [Hari 008i] — 2026-11-26 · Tumpukan Pasir Abelian (Sandpile C11 ≤ 999 Baris): Grup Sandpile, Elemen Identitas Fraktal & SOC

### Ditambahkan (`day008i`)
- Proyek C11 murni `projects/c/day008i_sandpile_c/` dengan total **903 baris** (`Makefile + src/*.c + src/*.h ≤ 999` baris, tanpa pustaka luar selain `libc` + `libm`, **12/12** kriteria `c1`–`c12` `tercapai`, `vacuous_flags = []`):
  - **Sifat Abelian & Aljabar Grup Sandpile**: `0` selisih hash FNV-1a atas tinggi, odometer, dan *sink* pada `600 × 100` urutan runtuh acak + gelombang sinkron (`c1`); determinan Bareiss `__int128` cocok persis dengan enumerasi pembakaran Dhar pada $n=2$ (`192`) dan $n=3$ (`100 352`), menghitung eksak $|G_4| = 557\,568\,000$, $|G_5| = 32\,565\,539\,635\,200$, $|G_6| = 19\,872\,369\,301\,840\,986\,112$, serta mendeteksi batas aman `__int128` di `n=6` (`n=7` luapan, `c2`); elemen identitas fraktal $e = s(2c_{\max} - s(2c_{\max}))$ lolos tes rekuren, simetri $D_4$, $s(e+e)=e$, $s(e+r)=r$ (`200/200` pada $n=24$ dan `100 352/100 352` pada $n=3$), dan orde elemen membagi $|G|$ (`15/15`, `c3`).
  - **Geometri Fraktal, SOC & Distribusi Tinggi**: kemiringan jari-jari tumpukan pusat $R(N)$ untuk $N=10^3..10^6$ sebesar `0.493878 ± 0.000505` (`c5`, kekekalan massa `max_abs_diff = 0`, `c4`); eksponen longsoran $2\cdot 10^6$ kejadian $\tau_s(128) = 1.0840 \pm 0.0045$ (`L=64`: `1.0663 ± 0.0075`, rasio ekor $s\ge 4096$ sebesar `7.07x`, `c6`); tinggi rata-rata kuartal-tengah `2.124406 ± 0.000167` vs `2.125` dan `P(0..3)` cocok dengan literatur Priezzhev dalam `max|dP| = 0.000368` (`c7`); kecepatan `114.31M ± 0.57M` runtuhan/detik (`c12`).
  - **Output**: `8` CSV ringkas (`64 KB`), `6` SVG tulis tangan `outputs/figures/day008i_*.svg`, `outputs/day008i/results.json` (identik 0-byte antar dua run, `c8`), dan [day008i_sim_sandpile-c.md](project-explanations/day008i_sim_sandpile-c.md).

## [Hari 008h] — 2026-11-26 · Jaringan Antrean Terminal Bandara Berbasis Agen & Validasi DTMC Eksak Lintas 8 Bahasa (48 Tahap)

### Ditambahkan (`day008h`)
- Proyek oktalingual `projects/c/day008h_airport_queue_octalingual/` dan notebook `days/day008h_sim_airport-queue-octalingual.ipynb` dengan **48 tahap (`s01`–`s48`)** dan **101 berkas keluaran identik per bahasa** (**17/22** kriteria `tercapai`, `1` `sebagian` `c15`, `4` `belum tercapai` `c9`/`c11`/`c12`/`c13` sebagai temuan jujur, `vacuous_flags = []`):
  - **Paritas Bit-Identik & Konservasi**: `0` selisih byte di 8 bahasa pada `S01`, `S05`, `S06`, dan `S40` (`3` hari = `25 920` slot); `0` pelanggaran invarian (`S15`, `little_diff = 0`); *independent queue replay* (`S41`) cocok eksak (`max_diff = 0`).
  - **Validasi DTMC Eksak & Kalibrasi (`C`, `365d`)**: `S09` `rel_err_L = 1.74e-11`; `S10` `9/9` dalam 3 SE; `S11` `4/4` dalam 3 SE di 8/8 bahasa; stasiun terpadat `SECURITY` (`util = 0.6284`, `loss_frac = 0.1133`).
  - **Kecepatan (`slot-pax/detik` median)**: C++20 `61.27M`, C11 `58.49M`, Rust `57.71M`, Go `50.61M`, TS (`node --experimental-strip-types`) `3.42M`, ipynb `566.9k`, Python `542.1k`, Ruby `331.2k`.
  - **Output**: `96` SVG bahasa + `8` pasang grafik gabungan PNG+SVG `outputs/figures/day008h_*.{png,svg}`, `outputs/day008h/results.json`, dan [day008h_sim_airport-queue-octalingual.md](project-explanations/day008h_sim_airport-queue-octalingual.md).

## [Hari 008g] — 2026-11-26 · Simulasi Epidemi Berbasis Agen SEIR Berlapis Lintas 8 Bahasa (39 Tahap)

### Ditambahkan (`day008g`)
- Proyek oktalingual `projects/c/day008g_epidemic_abm_octalingual/` dan notebook `days/day008g_sim_epidemic-abm-octalingual.ipynb` dengan **39 tahap (`s01`–`s39`)** dan **72 berkas keluaran identik per bahasa** (**14/18** kriteria `tercapai`, `1` `sebagian` `c12`, `3` `belum tercapai` `c6`/`c9`/`c11` sebagai temuan jujur, `vacuous_flags = []`):
  - **Paritas Bit-Identik & Konservasi**: `0` selisih byte di 8 bahasa pada `S01`, `S05`, `S06`, dan `S34` (`N=2000, 60d`); `0` pelanggaran invarian (`S12`); *in-process replay* (`S33`) dan *independent count replay* (`S37`) cocok eksak (`max_diff = 0`).
  - **Baseline & Efek Intervensi (`N=100000, 180d`)**: $R_{0,\text{ana}} = 1.7839$, serangan akhir tanpa intervensi `0.6820` vs dengan intervensi `0.0056` ($z = 453.5\text{ SE}$); puncak RS turun dari `250/250` (`396` penolakan) menjadi `7/250` (`0` penolakan).
  - **Kecepatan (`person-days/detik` median)**: C++20 `217.0M`, C11 `211.6M`, Rust `206.8M`, Go `163.1M`, TS (`node --experimental-strip-types`) `23.1M`, ipynb `3.36M`, Python `3.20M`, Ruby `2.03M`.
  - **Output**: `80` SVG per bahasa + `6` pasang grafik gabungan PNG+SVG `outputs/figures/day008g_*.{png,svg}`, `outputs/day008g/results.json`, dan [day008g_sim_epidemic-abm-octalingual.md](project-explanations/day008g_sim_epidemic-abm-octalingual.md).

## [Hari 008f] — 2026-11-26 · Jaringan Rantai Pasok 3 Pabrik/4 Gudang/10 Toko Lintas 8 Bahasa (40 Tahap)

### Diperbaiki (`day008f` — Bagian A Anti-Hampa)
- Mengkalibrasi ulang `projects/c/day008f_supply_chain_octalingual/` di seluruh 8 bahasa (`fill_rate = 0,920166`, `bo_event_frac = 0,073244`, `carrier_util = 0,848972` di `s09_calibration.csv`; store *fill rate* `0,868–0,883` di `s22_service.csv`), mengganti seluruh statistik konstan di `s21`/`s25`/`s31` dengan hitungan nyata dari data, menambahkan detektor `vacuous_flags` dan uji `>= 3 SE` di `compare.py`, serta mengeksekusi ulang seluruh 8 bahasa (paritas bit-identik 0 selisih, `vacuous_flags = []`).

### Diperbaiki (Bagian A — `day008e`)
- Menyelaraskan `outputs/day008e/results.json` dengan konvensi `day008d` (`run_completed_without_errors` terpisah dari `criteria`, tiap kriteria `{status, measured, threshold}`, info lingkungan/toolchain lengkap) serta melacak `rust/Cargo.lock` pada `day008d` dan `day008e` (`d053c0f`).

### Ditambahkan (Bagian C — `day008f`)
- Proyek oktalingual `projects/c/day008f_supply_chain_octalingual/` dan notebook `days/day008f_sim_supply-chain-octalingual.ipynb` dengan **40 tahap (`s01`–`s40`)** dan **103 berkas keluaran identik per bahasa** (**14/16** `tercapai`, `1` `sebagian` cabang langka `c5`, `1` `belum tercapai` saturasi kurir `c10` pada `S=1200`):
  - **Paritas Bit-Identik & Konservasi**: `0` selisih byte di 8 bahasa pada `S01`, `S04`, `S05`, `S06`, dan `S34`; `0` pelanggaran invarian (`S11`); *in-process replay* (`S33`) dan *ledger replay* (`S38`) cocok eksak (`0` selisih).
  - **Kecepatan (`ns/event` median $\pm$ MAD)**: Rust `62,80`, C11 `98,50`, C++20 `104,94`, Go `128,96`, TS `1424,23`, ipynb `6596,50`, Python `6749,13`, Ruby `15048,00`.
  - **Output**: `80` SVG per bahasa + `6` pasang grafik gabungan PNG+SVG `outputs/figures/day008f_*.{png,svg}`, `outputs/day008f/results.json`, dan [day008f_sim_supply-chain-octalingual.md](project-explanations/day008f_sim_supply-chain-octalingual.md).

## [Hari 008e] — 2026-11-26 · Solver Poisson 2D Lintas 8 Bahasa (19 Tahap Identik: C, C++, Rust, Go, Ruby, Python, TS, ipynb)

### Ditambahkan (`day008e`)
- Proyek oktalingual `projects/c/day008e_poisson_solvers_octalingual/` dan notebook `days/day008e_sim_poisson-solvers-octalingual.ipynb` dengan 19 tahap (`s01`–`s19`) yang identik di **8 bahasa** (**11/12** kriteria `c1`–`c12` `tercapai`, `c6` `sebagian` pada `pcg_ssor` $\omega=1$):
  - **Paritas Bit-Identik & Anchor Analitik**: `0` selisih pada `s01` dan `s17` di seluruh 8 bahasa (`8/8`); galat nilai eigen maks `3,165e-13`; orde diskritisasi `2,001445`; $\rho_{\text{emp}}$ Jacobi (`0,9951847`) dan GS (`0,9903927`) pada $N=31$ cocok hingga `< 4e-8` relatif; faktor kontraksi MG V(2,1) `0,07011662` ($N=63$) dan rasio FMG `0,473908`.
  - **Kecepatan (`ns/unknown-iter` median, `rbgs` / `mg` pada `N=63`)**: C++ `0,856 / 5,156`, C `0,954 / 5,493`, Rust `1,412 / 8,723`, Go `1,848 / 13,464`, TS (`node --experimental-strip-types`) `5,718 / 31,756`, ipynb `124,15 / 807,52`, Python `125,10 / 811,68`, Ruby `232,15 / 1299,91`.
  - **Output**: `24` SVG per bahasa + `4` pasang grafik gabungan PNG+SVG `outputs/figures/day008e_*.{png,svg}`, `outputs/day008e/results.json`, dan [day008e_sim_poisson-solvers-octalingual.md](project-explanations/day008e_sim_poisson-solvers-octalingual.md).

## [Hari 008d] — 2026-11-26 · Jalan Acak Kisi Pólya Lintas 8 Bahasa (12 Tahap Identik: C, C++, Rust, Go, Ruby, Python, TS, ipynb)

### Ditambahkan (`day008d`)
- Proyek oktalingual `projects/c/day008d_lattice_walks_octalingual/` dan notebook `days/day008d_sim_lattice-walks-octalingual.ipynb` dengan 12 tahap (`s01`–`s12`) yang identik di **8 bahasa** (**10/10** kriteria `c1`–`c10` `tercapai`):
  - **Paritas Bit-Identik & Anchor Eksak**: `0` selisih pada `s01` dan `s02` di seluruh 8 bahasa; `39/39` MSD, `12/12` peluang kembali, dan `6/6` Gambler's Ruin berada dalam $\le 4\text{ SE}$ di setiap bahasa; renewal eksak $F_3(10^4) = 0{,}337668$ (`0,00287` dari konstanta Watson `0,340537`).
  - **Kecepatan (`ns/langkah` median, `d=2`)**: Rust `1,683`, C++ `1,715`, C `1,719`, Go `23,696`, TS (`node --experimental-strip-types`) `34,213`, ipynb `565,98`, Python `570,43`, Ruby `956,90`.
  - **Output**: `24` SVG per bahasa + `3` pasang grafik gabungan PNG+SVG `outputs/figures/day008d_*.{png,svg}`, `outputs/day008d/results.json`, dan [day008d_sim_lattice-walks-octalingual.md](project-explanations/day008d_sim_lattice-walks-octalingual.md).

## [Hari 008c] — 2026-11-26 · Model Vicsek (Aktif-Materi) di C++20 dengan Cell List O(N) vs Brute Force O(N^2)

### Diperbaiki (Bagian A — `day008a`)
- Menghapus path absolut pada sel keluaran `days/day008a_sim_ising-trilingual.ipynb` menjadi path relatif `outputs/figures` (`35f7a2c`).

### Ditambahkan (Bagian C — `day008c`)
- Proyek C++20 murni `projects/cpp/day008c_vicsek_cpp/` (`427` baris sumber di `src/`, `5/5` tes unit lulus, paralel `std::thread`):
  - **Paritas & Skala**: Cell list cocok bit-identik (`0` selisih) dengan brute force pada `N=1000, 500 langkah` (`c1`), dan lebih cepat `108×` pada `N=16000` dengan kemiringan log–log `1,0745` vs `2,0518` (`c5`, `c6`).
  - **Kurva $\phi(\eta)$ & Temuan Jujur**: Turun monoton (`c4`), titik tengah $\phi=0{,}5$ bergeser `2,887 -> 2,611 -> 2,281` untuk `L=16,32,64`; `c3` berstatus `belum tercapai` karena $L=64$ pada $\eta=0$ hanya mencapai $\phi=0{,}7471$ dalam 2000 langkah termalisasi.
  - **Output**: 4 SVG mandiri `outputs/figures/day008c_*.svg`, `outputs/day008c/results.json`, dan [day008c_sim_vicsek-cpp.md](project-explanations/day008c_sim_vicsek-cpp.md).

## [Hari 008b] — 2026-11-26 · Simulasi Wa-Tor (Ikan vs Hiu) di C++20, Fit Lotka–Volterra RK4 & Efek Ukuran Hingga

### Ditambahkan (`day008b`)
- Proyek C++20 murni `projects/cpp/day008b_wator_cpp/` (`438` baris sumber di `src/`, `4/4` tes unit lulus, paralel `std::thread` untuk 40 replika):
  - **Wa-Tor 200×200 (5000 langkah)**: Pembaruan 2-fase acak Fisher–Yates dengan invarian lulus tiap langkah (`c1`), determinisme bit-identik (`c2`), osilasi periode `99` langkah (`acf = 0,6902`, `c3`), serta fit OLS ke model kontinu Lotka–Volterra RK4 berperiode `88,46` (`rasio = 1,119 <= 2,0`, `c4`).
  - **Efek Ukuran Hingga & SVG Mandiri**: Proporsi kepunahan (40 replika, Wilson 95%) turun dari `1,00` ($W=20, 40$) menjadi `0,65` ($W=80$) dan `0,10` ($W=160$, `c5`), disertai 4 berkas SVG mandiri `outputs/figures/day008b_*.svg` (`c6`), `outputs/day008b/results.json`, dan [day008b_sim_wator-cpp.md](project-explanations/day008b_sim_wator-cpp.md).

## [Hari 008a] — 2026-11-26 · Model Ising 2D Lintas Tiga Bahasa (C11, C++20, Python): Metropolis, Wolff, Swendsen–Wang, FSS & Paritas Bit-Identik

### Diperbaiki (Bagian A — `day006e`)
- Menambahkan eksperimen tindak lanjut SBC `target_accept = 0.95` pada `days/day006e_sim_bayesian-sv-nuts.ipynb` dan mengeksekusi ulang bersih dua kali (`533,2 s <= 600 s`, `SHA-256` identik `f6033b51f2b4dbb8...`): replikasi berdivergensi turun dari `139/200` ke `76/200`, p-nilai $\chi^2_{19}$ menjadi `mu = 0,6838`, `phi = 0,1047`, `sigma_eta = 0,2943` (subset tanpa divergensi `n=124` `mu p=0,2624` dilabeli bias seleksi `HIPOTESIS`, bukan koreksi), serta memisahkan `FAKTA` vs `HIPOTESIS` pada [day006e_sim_bayesian-sv-nuts.md](project-explanations/day006e_sim_bayesian-sv-nuts.md).

### Ditambahkan (Bagian C — `day008a`)
- Proyek tiga bahasa `projects/c/day008a_ising_trilingual/` (`SPEC.md`, `c/`, `cpp/`) dan notebook `days/day008a_sim_ising-trilingual.ipynb` (`14/14` kriteria `c1`–`c14` `tercapai`):
  - **Paritas Bit-Identik & RNG Lintas Bahasa**: `SplitMix64` + `xoshiro256++` + `FNV-1a` menghasilkan `0` selisih pada 3000 langkah dan hash akhir identik di C, C++, dan Python untuk Metropolis (`7458071e0380ee91`), Wolff (`52911df0a1ba8509`), dan Swendsen–Wang (`50d1e256b6df3e61`).
  - **Anchor Eksak & Paralelisme**: Enumerasi $L=4, 5$ cocok hingga `2,92e-16`; hasil `1, 2, 4, 8` *thread* (`pthread` / `std::thread`) identik bit-demi-bit; Onsager (AGM) pada $T \in \{2{,}0, 2{,}6, 3{,}2\}$ dan Yang pada $T < T_c$ semuanya dalam $\le 1\text{ SE}$.
  - **Finite-Size Scaling & Kecepatan**: Perpotongan Binder $U_4$ dalam $\le 0{,}0047$ dari $T_c = 2{,}269185$, $\hat{\gamma}/\hat{\nu} = 1{,}7682 \pm 0{,}0061$, $\hat{\beta}/\hat{\nu} = 0{,}1122 \pm 0{,}0030$, $\hat{z}_{\text{Metropolis}} = 1{,}8625$, $\hat{z}_{\text{Wolff}} = 0{,}4118$, $\hat{z}_{\text{SW}} = 0{,}1998$; kecepatan C/C++/Python terukur median $\pm$ MAD (7 ulang).
- Output `outputs/day008a/results.json`, 8 pasang grafik PNG+SVG `outputs/figures/day008a_*.{png,svg}`, `projects/c/day008a_ising_trilingual/README.md`, dan [day008a_sim_ising-trilingual.md](project-explanations/day008a_sim_ising-trilingual.md).

## [Hari 007a] — 2026-11-25 · Kekacauan Deterministik: Eksponen Lyapunov Peta Logistik & Lorenz, Konstanta Feigenbaum & Dimensi Kaplan–Yorke (C++20, Bonus 10 Menit)

### Ditambahkan (`day007a`)
- Proyek C++20 murni `projects/cpp/day007a_chaos_lyapunov_cpp/` (nol pustaka eksternal, `-O3 -std=c++20 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off`, `4/4` tes unit lulus):
  - **Peta Logistik & Feigenbaum**: Sapuan $\lambda(r)$ pada `1600` titik $r \in [2{,}5, 4{,}0]$ (`200.000` iterasi pasca-transien), pencarian akar Newton siklus superstabil $n = 1\dots 8$ menghasilkan $\delta_8 = 4{,}669191$ dan $r_\infty = 3{,}5699457$ (`c2 = tercapai`), serta dokumentasi jujur perangkap titik tetap $x = 0$ ($\ln 4 = 1{,}386294$) pada orbit *binary64* $10^7$ iterasi di $r = 4$ ($\hat{\lambda}(4) = 1{,}004447$, `c1 = belum tercapai` sebagai temuan artefak *floating-point*).
  - **Spektrum Lyapunov Lorenz & Dimensi Kaplan–Yorke**: Integrasi RK4 persamaan variasional 12-komponen + ortonormalisasi Gram–Schmidt Benettin ($T = 2000, \Delta t = 0{,}005$) menghasilkan $(\lambda_1, \lambda_2, \lambda_3) = (0{,}90968, -0{,}00083, -14{,}57551)$, jumlah kontraksi Liouville $\sum \lambda_i = -13{,}666660$ (`err = 6,2e-6` terhadap $-41/3$, `c3 = tercapai`, `c4 = tercapai`), dimensi Kaplan–Yorke $D_{\mathrm{KY}} = 2{,}06235$ (`c5 = tercapai`), dan pertumbuhan pemisahan $\delta_0 = 10^{-9}$ berkemiringan `0,9359`.
- Output `outputs/day007a/results.json`, CSV ringkasan, 4 berkas SVG mandiri `outputs/figures/day007a_*.svg` (+ PNG bifurkasi), `projects/cpp/day007a_chaos_lyapunov_cpp/README.md`, dan [day007a_sim_chaos-lyapunov-cpp.md](project-explanations/day007a_sim_chaos-lyapunov-cpp.md).

## [Hari 006e] — 2026-11-24 · HMC/NUTS dari Nol untuk Volatilitas Stokastik Bayesian, Kontras Funnel & Kalibrasi SBC (Notebook Python, Bonus)

### Ditambahkan (`day006e`)
- Notebook tereksekusi penuh `days/day006e_sim_bayesian-sv-nuts.ipynb` (NumPy + SciPy murni tanpa pustaka MCMC eksternal; diverifikasi berjalan bersih dua kali berturut-turut dari nol dalam `295,5 s` dan `318,2 s` dengan sidik jari numerik `SHA-256` identik persis):
  - **Model SV Non-Terpusat, Gradien Adjoin Resolven & Multinomial NUTS**: Transformasi tanpa batas + log-Jacobian eksak, *backpropagation* manual lewat operator adjoin resolven AR(1) (`8,88×` lebih cepat via `lfilter` vs loop Python; galat FD `9,58e-9`), integrator *leapfrog* simplektik (`4,44e-16` reversibilitas, kemiringan energi `1,9985`), *dual averaging* + matriks massa diagonal berjendela ala Stan, dan Multinomial NUTS dengan kriteria *generalized U-turn* serta pencabangan progresif terbias (`1819,6 ESS/s` pada Gaussian 50-D).
  - **Diagnostik *Neal's Funnel*, SBC 200 Replikasi & Aplikasi Rezim**: Pada SV sintetis $T = 1000$, parameterisasi non-terpusat mencapai `split-Rhat <= 1,0040`, `ESS bulk >= 1047,8`, liputan jalur $h_t$ 90% sebesar `90,5%`, dan lonjakan efisiensi `136,8×` pada $\sigma_\eta$ atas parameterisasi terpusat (`c4 = sebagian` jujur karena `15/4000` divergensi di ekor $\phi \to 1$); pada `200` replikasi SBC ($T=250, L=99$), seluruh kurva ECDF-diff berada di dalam pita simultan 99% (`<= 0,055 < 0,110`; $\chi^2_{19}$ lulus pada $\phi$ `p=0,6169` dan $\sigma_\eta$ `p=0,7599`, `mu` `p=0,0011` → `c6 = sebagian`), dilengkapi evaluasi *log predictive score* 1-langkah di 30 jendela lonjakan 4-rezim vs EWMA ($\lambda = 0{,}94$) dan tabel ablasi $(\varepsilon, M)$.
- Output `outputs/day006e/results.json` (`8/10` kriteria `tercapai`, `2/10` `sebagian` sebagai temuan jujur), 6 pasang grafik PNG+SVG `outputs/figures/day006e_*.{png,svg}`, dan [day006e_sim_bayesian-sv-nuts.md](project-explanations/day006e_sim_bayesian-sv-nuts.md).

## [Hari 006d] — 2026-11-24 · Rekayasa Performa Solver PDE Heston 2D ADI Craig–Sneyd: Tangga Optimasi Terukur & Roofline (C++20, Bonus)

### Diperbaiki (Bagian A — `day006c`)
- Menjalankan ulang tes penuh (`cargo test --release` 4/4, `go test` 3/3, `unittest` 4/4) dan pipeline `run_all.sh` (`3` detik, exit `0`) pada `projects/rust/day006c_percolation_polyglot/`, menghasilkan 5 pasang grafik `outputs/figures/day006c_perc_*.{png,svg}` dan memperbarui `outputs/day006c/results.json` serta [day006c_sim_percolation-polyglot.md](project-explanations/day006c_sim_percolation-polyglot.md) (`c1`–`c7` `tercapai`, sementara `timing` `789` detik dan `c_wall_clock = belum tercapai` dipertahankan utuh beserta `rerun_note`).

### Ditambahkan (Bagian C — `day006d`)
- Proyek C++20 murni `projects/cpp/day006d_heston_adi_cpp/` (`899` baris C++ di `src/` + `tests/`, nol pustaka eksternal, `-O3 -std=c++20 -fno-fast-math -ffp-contract=off` serta biner pembanding `-mcpu=native`):
  - **Kuadratur Semi-Analitik & Akurasi PDE**: Gauss–Legendre 16-titik komposit (200 panel) untuk Gil–Pelaez (*call* & *put*) dan Lewis ($\mathrm{Im}(z) = -1/2$) dengan fungsi karakteristik Heston *little trap* (paritas `4,55e-13` Set A dan `1,71e-13` Set B; cocok dengan `qb100` Python `8,802660962859`); solver ADI Craig–Sneyd $\theta = 1/2$ mencapai galat `2,54e-3` pada `401×201, M=200` (`6,07e-4` pada `801×401, M=400`) dengan orde konvergensi 3-level `2,28` pada rezim Feller dilanggar (Set A) dan `2,01` pada Set B (`8/8` unit test lulus).
  - **Tangga Optimasi `V0`–`V5` & Analisis Roofline**: Pada `801×401, M=200` di Apple M4 (`STREAM triad` 10T `86,34 GB/s`, puncak register tanpa FMA 10T `90,83 GFLOP/s`), waktu median turun dari `3,2755 ± 0,0282 s` (`50,99 ns`/sel-langkah, `V0`) menjadi `0,2622 ± 0,0214 s` (`4,08 ns`/sel-langkah, `12,98 GFLOP/s`, `V5` `T=10`, **speedup `12,49×`**; `13,33×` pada `401×201`) dengan medan akhir **identik bit-demi-bit** dari `V2` hingga `V5` (`T ∈ {1, 2, 4, 10}`) maupun pada `-mcpu=native`, disertai pelaporan jujur regresi `V2` vs `V1` (`0,856×`) ketika *working set* melampaui cache L2 `16 MiB` (`12,48 MiB` → `17,48 MiB`).
- Output `outputs/day006d/results.json` (`9/9` kriteria `c1`–`c9` `tercapai`), CSV mentah + laporan vektorisasi, 6 pasang grafik PNG+SVG `outputs/figures/day006d_*.{png,svg}`, `projects/cpp/day006d_heston_adi_cpp/README.md`, dan [day006d_sim_heston-adi-cpp.md](project-explanations/day006d_sim_heston-adi-cpp.md).

## [Hari 006c] — 2026-11-24 · Perkolasi Bond Kisi Persegi Poliglot: Newman–Ziff, Dualitas-Diri Eksak & Orkestrator Deterministik (Python + Rust + Go, Bonus 10 Menit)

### Diperbaiki (Tahap 1 — `day006b`)
- `outputs/day006b/results.json`, `docs/project-explanations/day006b_sim_saw-polyglot.md`, `projects/rust/day006b_saw_polyglot/py/saw_analysis.py`, dan `README.md`: mengganti path absolut menjadi relatif, menetapkan `c10_ten_minute_wall_clock = belum tercapai` (`666` detik sesi penuh vs `420` detik pipeline komputasi), menambah poin landasan teori bertautan kode relatif (`12` poin), serta memperbarui tabel progres `README.md`.

### Ditambahkan (Tahap 2 — `day006c`)
- Proyek poliglot `projects/rust/day006c_percolation_polyglot/` (Rust murni `#![forbid(unsafe_code)]` tanpa *crate* `percolate`, orkestrator Go pustaka standar `orchestrate`, serta modul Python `percolation.py`, `finalize.py`, `test_percolation.py`):
  - **Inti Rust & Orkestrator Go Terverifikasi**: Union-Find ganda dengan simpul maya kiri–kanan, Fisher–Yates `Xoshiro256Pp`, serta kolam pekerja pipa JSON deterministik (`cargo clippy`, `cargo test --release` `4/4` termasuk enumerasi eksak $2 \times 3$ `64/128 = 1/2`, `go vet`, dan `go test` `3/3` termasuk invarian `1` vs `4` *worker* lulus pada menit `7:31`).
  - **Disiplin Kunci Waktu (`LOCK` 9:00)**: Eksekusi pipeline penuh `run_all.sh` dihentikan secara patuh saat `dl.sh` mengunci pada menit ke-9; status jujur dicatat pada `outputs/day006c/results.json` dan [day006c_sim_percolation-polyglot.md](project-explanations/day006c_sim_percolation-polyglot.md) (`c7 = tercapai`, `c1`–`c2 = sebagian`, `c3`–`c6 = belum diketahui`, `c_wall_clock = belum tercapai`).

### Diperbaiki (rerun pasca-sesi — day006d Bagian A)
- Menjalankan ulang seluruh tes (`cargo test --release` 4/4, `go test` 3/3, `unittest` 4/4) dan `run_all.sh` penuh (`3` detik, exit `0`), menghasilkan 5 pasang grafik `outputs/figures/day006c_perc_*.{png,svg}` dan memperbarui `outputs/day006c/results.json`: `c1`–`c7` `tercapai` (`\hat{\nu} = 1,3406 \pm 0,0079`, `\hat{D} = 1,8678 \pm 0,0077`, maks `|Z| = 1,4839` di $p=1/2$), sementara `timing` (`789` detik) dan `c_wall_clock = belum tercapai` dipertahankan dari sesi asli beserta `rerun_note`.

## [Hari 006b] — 2026-11-24 · Self-Avoiding Walk (SAW) Poliglot di Kisi Persegi: Enumerasi Eksak, Pivot Monte Carlo & Eksponen Kritis (Python + Rust + Go, Bonus 10 Menit)

### Ditambahkan (`day006b`)
- Proyek poliglot tiga bahasa `projects/rust/day006b_saw_polyglot/` (Rust murni `#![forbid(unsafe_code)]` tanpa *crate* + orkestrator konkuren Go pustaka standar + referensi *brute-force*, regresi statistik & visualisasi Python):
  - **Enumerasi Eksak & Orkestrasi**: DFS bereduksi simetri $D_4$ ($\times 8$) dan peta panas titik ujung $n = 14$ di Rust, dipadukan dengan pemecahan awalan $K = 7$ ($543$ tugas berbobot $4$) pada kolam proses pekerja Rust persisten di Go (`c_12 = 324.932` cocok persis dengan Python *brute-force*; `c_18 = 124.658.732` identik antara Rust 1-proses serta Go `1` vs `10` *worker*, `0` ketidakcocokan terhadap A001411).
  - **Pivot Madras–Sokal & Eksponen Kritis**: PRNG `xoshiro256++` + `StampSet` menghasilkan `100/100` sampel lintasan $N=300$ bebas tabrakan serta taksiran $\hat{\mu} = 2{,}6290 \pm 0{,}0176$, $\hat{\gamma} = 1{,}3834 \pm 0{,}0959$, dan $\hat{\nu} = 0{,}7447 \pm 0{,}0010$ (`9/10` kriteria `c1`–`c9` `tercapai`; `c10_ten_minute_wall_clock = belum tercapai` di `666` detik sesi penuh / `420` detik pipeline komputasi).
- Output `outputs/day006b/results.json`, 4 pasang grafik PNG+SVG `outputs/figures/day006b_*.{png,svg}`, `projects/rust/day006b_saw_polyglot/README.md`, dan [day006b_sim_saw-polyglot.md](project-explanations/day006b_sim_saw-polyglot.md).

## [Hari 005j] — 2026-11-23 · Perbaikan day005i + Laboratorium Model Faktor, PCA BBP, GRS, Shanken, Risk Parity, Black–Litterman, Markowitz QP & Backtest Portofolio (Ruby, Bonus)

### Diperbaiki (Bagian A — `day005i`)
- `projects/ruby/day005i_garch_family/`: meningkatkan replikasi `c6_sandwich_se_restores_t5_coverage` menjadi $R=1000$ ($T \in \{500, 1000, 2500\}$, plus $T=10{.}000$ di `bin/run_extras`: `cov_sw = 0,9380`, `se = 0,0076`, rasio `1,791x`) dengan ambang teoretis *a priori* `|cov_sw - 0,95| <= 3 SE` (`c6` jujur **`belum tercapai`** di $T=2500$: `cov_sw = 0,9130 +- 0,0089`, `z = 4,15`), menaikkan `c19_structural_break_spurious_igarch` dari `25` ke `300` replikasi beserta distribusi lengkap `mean/sd/se/min/p05/p25/median/p75/p95/max`, serta menambahkan benchmark intraday sepadan-informasi (`RV-Naif = 0,02595`, `EWMA-RV = 0,01533`, `GARCH-X = 0,01739` vs `HAR-RV = 0,01175`, seluruh DM $p \le 9{,}33\times 10^{-7}$). `37` berkas RuboCop & `44` tes (`258` asersi) bersih.

### Ditambahkan (Bagian C — `day005j`)
- Project Ruby murni bervolume sangat tinggi `projects/ruby/day005j_factor_portfolio_lab/` ($\ge 3{,}3$, nol gem numerik, `36` modul di `lib/factor_portfolio_lab/`, `83` unit/integrasi tes dengan `712` asersi, `42` kriteria `c1`–`c42`, `26` grafik SVG di `outputs/figures/day005j_*.svg`, `bin/run` selesai dalam `364,21s` dengan 10 pekerja paralel + YJIT ON):
  - **Bagian 1–3 (`SpikedStudy`, `GrsStudy`, `FmStudy`)**: verifikasi BBP $R=500$ ($\ell=0{,}3$ tetap di tepi MP `2,168` vs `2,250`, $\ell=2{,}0$ alignment `0,832` vs teori `0,833`), GRS $R=5000$ (*size* exact `0,0522` vs Wald `0,8482` di $N=50, T=120$), serta Fama–MacBeth + Shanken $R=1000$ (`c_Sh = 0,405`, cakupan intersep `0,848 -> 0,903`, faktor `>= 0,959`).
  - **Bagian 4–6 (`RiskParity`, `BlackLitterman`, `Markowitz`)**: Newton Spinu `6` iterasi (`rc_max_dev <= 5,55e-17`) sepakat dengan CCD (`max_w_diff <= 5,39e-15`), Black–Litterman *importance sampling* $R=500{.}000$ (`max_mu_diff = 5,41e-5`, `ESS = 288.497`), dan QP Markowitz `45` titik (*long-only* `max_kkt = 4,44e-16`, *capped 12%* `9,98e-12`).
  - **Bagian 7–12 (`Resampling`, `BacktestStudy`, `Attribution`, `Stress`, `Robustness`, `Complexity`)**: Michaud $R=500$ (`HHI 0,809 -> 0,221`, `L1 1,702 -> 1,410`), *backtest* rolling dua-rezim `8` strategi $R=300$ pada $c \in \{0, 5, 20\}$ bps, atribusi Euler (`err = 2,08e-17`) & salah spesifikasi $K$, *reverse stress test* Mahalanobis tertutup (`-0,1500` eksak), kekokohan Student-$t(4)$, dan skalabilitas $N \in \{10..200\}$ (`slope_invert = 3,17`).
- Status kriteria `39/42` `tercapai` dan `3` `belum tercapai` (`c30`, `c32`, `c33` dilaporkan jujur); `83` tes (`712` asersi) & `rubocop` (`51` berkas) bersih.
- Output `outputs/day005j/results.json`, `26` grafik SVG `outputs/figures/day005j_*.svg`, `projects/ruby/day005j_factor_portfolio_lab/README.md`, serta [day005j_quant_factor-portfolio-lab.md](project-explanations/day005j_quant_factor-portfolio-lab.md).

## [Hari 005i] — 2026-11-23 · Perbaikan day005h/day005g + Sistem Volatilitas GARCH, QMLE Sandwich, VaR/ES, DCC, HAR-RV & Salah Spesifikasi (Ruby, Bonus)

### Diperbaiki (Bagian A — `day005h` & `day005g`)
- `projects/ruby/day005h_backtest_overfitting/`: mendiagnosis bias negatif sistematis Lo $\hat{\eta}(q)$ pada `c3` (normalisasi keliru terhadap parameter populasi $\mathrm{SR}_1 = 0{,}05$ alih-alih $\widehat{\mathrm{SR}}_1 = 0{,}04864$ pada sampel yang sama, sehingga setelah diperbaiki terukur `0,6157` vs `0,6200` di $q=5$ dan `0,5541` vs `0,5559` di $q=21$ serta menyentuh `[0,6203; 0,5540]` di sapuan sampel $T = 20.000$), mengukur galat standar cakupan *stationary bootstrap* (`0,91375 +- 0,00993`, hanya `+0,38 SE` di atas `0,91` akibat peredaman autokorelasi blok $\tilde{\phi} = 0{,}27$), dan menandai `c3` jujur sebagai **`sebagian`** (`10/11` `tercapai`, `32` tes, `20,66s`).
- `projects/ruby/day005g_mlmc_heston/`: menurunkan batas atas reduksi variansi antitetik Giles–Szpruch khusus skema Euler tanpa koreksi kuadratik Milstein pada $\ell = 1..6$ (`part4_multi_model`: 1D GBM Euler tepat `1,0000`, 1D CIR *Full-Truncation* meluruh `1,0900 -> 1,0611`, 2D Heston $\rho = -0,7$ meluruh `1,2744 -> 1,1869 <= 2/(1+\rho^2) = 1,3423`) dengan tetap mempertahankan ambang CIR `<= 1,04` (`c6` jujur **`belum tercapai`** karena penolakan/pemotongan batas Feller $v_n^+=0$; `22` tes, `90,26s`).

### Ditambahkan (Bagian C — `day005i`)
- Project Ruby murni bervolume besar `projects/ruby/day005i_garch_family/` ($\ge 3{,}3$, nol gem numerik, `26` modul di `lib/garch_family/`, `44` unit/integrasi tes dalam `6` berkas dengan `258` asersi, `22` kriteria `c1`–`c22`, `13` grafik SVG, `bin/run` selesai dalam `70,14s` dan `bin/run_extras` `30,08s`):
  - **Bagian 1 (`MomentsStudy`)**: Verifikasi `15` konfigurasi `Garch11`/`GjrGarch11`/`Egarch11` (`max_var_z = 1,95`, kurtosis `1,95`, ACF kuadrat `2,35`), prediksi multi-horizon (`max_rel_err = 0,57%`), dan paritas skor analitik vs beda hingga `< 1e-5`.
  - **Bagian 2 (`Optim`, `Inference`, `Mle`)**: BFGS ber-skala, kovariansi *Sandwich* Bollerslev–Wooldridge, studi QMLE $T \in \{500, 1000, 2500\}$ (`R=400`, Normal `400/400` & Student-$t(5)$ `400/400`: SE naif gagal `cov = 0,720`, *Sandwich* menaikkan rasio `1,676x` ke `cov = 0,908`; `extras` $T=10\mathrm{k}$ `150/150` `0,700 -> 0,900`).
  - **Bagian 3 (`InnovationStudy`)**: Student-$t$ MLE ($\hat{\nu}_{\mathrm{med}} = 5,048$, efisiensi MSE `2,47x` pada $\alpha$) serta uji batas $1/\nu = 0$ campuran $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_1$ (`zero = 0,577`, ukuran `0,0433`, kuasa `1,000`).
  - **Bagian 4 (`RiskBacktest`)**: Backtest *rolling* OOS 1% & 5% ($T_{\mathrm{oos}} = 3000$): VaR/ES 1% (`garch_n: 1,80%, Z2 = -0,949` vs `garch_t: 0,87%, Z2 = +0,143`, `garch_fhs: 1,03%, Z2 = -0,031`), serta MC `R=500` ukuran `0,046` & kuasa 1% `0,774`.
  - **Bagian 5 (`ModelComparison`)**: Patton QLIKE anchor, seleksi AIC/BIC (`+7,28`), NIC (`1,668x`), Diebold–Mariano (`p = 0,076` pada `1k` OOS `c13 = sebagian`; `0,0075` pada `4k` `extras`), serta ukuran (`0,035`) & kuasa (`0,980`) uji batas $\gamma = 0$.
  - **Bagian 6 (`Dcc`)**: DCC(1,1) 2-tahap (`CCC` `0,0`, `min_det = 0,444`, `corr = 0,963`, `200/200` `a = 0,0512, b = 0,8907`; portofolio `DCC p_cc = 0,812` vs `Uncorr 2,9e-23`).
  - **Bagian 7 (`HarRv`)**: Intraday $M \in \{13..390\}$ *infill* (`max_z = 1,66`, penyusutan MSE `34,59x`), Corsi HAR-RV mengungguli GARCH harian (`QLIKE 0,0117` vs `0,0319`, DM `p = 3,64e-8`).
  - **Bagian 8–9 (`MisspecPortfolio`)**: *Break* $\omega$ menciptakan persistensi semu IGARCH `0,9977` yang pulih ke `0,8435` saat dipisah (CUSUM `100%` di `1511,8`), dan portofolio 5-aset GJR-DCC + stres krisis (`euler_err = 8,67e-19`, `OOS p_uc = 0,0695`, `div_drop = 0,313`, `var1_amp = 3,23x`).
- Status kriteria `21/22` `tercapai` dan `c13 = sebagian`; `44` tes (`258` asersi) & `rubocop` (`37` berkas) bersih.
- Output `outputs/day005i/results.json`, `13` grafik SVG `outputs/figures/day005i_*.svg`, `projects/ruby/day005i_garch_family/README.md`, serta [day005i_quant_garch-family.md](project-explanations/day005i_quant_garch-family.md).

## [Hari 005h] — 2026-11-23 · Perbaikan day005g + Backtest, Statistik Sharpe, Bias Seleksi, Pengujian Ganda & PBO (Ruby, Bonus)

### Diperbaiki (Bagian A — `day005g`)
- `projects/ruby/day005g_mlmc_heston/`: melaporkan fakta terukur *speedup* MLMC Heston (`0,1837x` pada $\varepsilon = 0{,}01$, i.e. `5,44x` lebih lambat dari MC standar; `0,4978x` pada $\varepsilon = 0{,}005$ dalam `bin/run_extras`) bersama perbandingan kemiringan biaya terhadap prediksi teori Giles (asimtotik `2,172`, pra-asimtotik diskret `2,475`, terukur `2,471`) dan titik silang ekstrapolasi ($\varepsilon^* \approx 5{,}05\times 10^{-4}$..$1{,}59\times 10^{-3}$; kontrol Feller OK $\xi=0{,}25$ memberi `speedup = 2,557x`), memperketat ambang `c2`/`c3` ke nilai teori/Feller serta `c6` ke target literatur `>= 2,0x` (status jujur **`belum tercapai`** di `1,253x`/`1,216x` disertai diagnosis skema Euler tanpa koreksi kuadratik Milstein), dan menambah `11` tes unit (`21` tes, `34` asersi).

### Ditambahkan (Bagian C — `day005h`)
- Project Ruby murni mandiri `projects/ruby/day005h_backtest_overfitting/` ($\ge 3{,}3$, nol gem numerik, `SplitMix64`, `Process.fork` + `IO.pipe` + `Marshal` deterministik, JSON penuh, 6 grafik SVG):
  - **Bagian 1**: Generator sintetis, pustaka `300` strategi kausal (`SMA`, `MOM`, `MRZ`), mesin backtest beserta invarian ekuitas/biaya/tanpa *look-ahead*, dan *anchor* eksak `AR(1)` (`z = 0,65`) & `GBM` *buy-and-hold* (`z = 0,55`).
  - **Bagian 2**: Inferensi Sharpe ($R = 20\mathrm{k}$ MC Gaussian, `Student-t(5)`, `Lognormal(0,45)`), faktor koreksi autokorelasi Lo $\eta(q)$ ($q \in \{5, 21\}$), dan *stationary bootstrap* Politis–Romano (`AR(1)` $\phi=0{,}3$: cakupan `0,846 -> 0,914`).
  - **Bagian 3**: $\mathbb{E}[\max Z_i]$ Gauss–Legendre 160-titik vs Euler–Mascheroni (`< 0,018` di $N=300$), penyusutan $N_{\mathrm{eff}} = 231{,}5$ dari korelasi pustaka ($\bar{\rho} = 0{,}229$), serta kalibrasi `DSR > 0,95` (`0,1%` lolos palsu vs `100%` uji naif di $R=2\mathrm{k}$).
  - **Bagian 4**: Pengujian ganda (`Unadjusted`, `Bonferroni`, `Holm`, `Benjamini–Hochberg`, `Romano–Wolf`) di $R = 2\mathrm{k}$ ($N = 300, \pi_0 = 0{,}90, \alpha = 0{,}05$): `Bonf FWER = 0,041`, `BH FDR = 0,0405` (vs teori `0,045`), plus dependensi $\rho = 0{,}4$ & `bin/run_extras` `RW FWER = 0,046 / Power = 0,601`.
  - **Bagian 5**: Enumerasi CSCV eksak seluruh `12.870` pembagian $\binom{16}{8}$ (*anchor* *noise* `0,464 +- 0,036` vs `0,50`, sinyal `0,000`, `0,462 +- 0,035` di $N=300$ `extras`), dan studi biaya $c \in \{0, 5, 20\}$ bps (runtuhnya `MOM(1)` digantikan `SMA(10,40)`).
  - **Bagian 6**: Validasi label tumpang-tindih $h = 10$ (`Naive Shuffled K-Fold` bocor `+0,3741` vs `Purged + Embargoed` `-0,0037 +- 0,0095`, `z = 0,39`).
- Runner utama `day005h` selesai dalam `20,62s` (`Ruby 4.0.9`, YJIT `ON`, `10` *worker*; `11/11` kriteria `tercapai`), `31` tes unit/integrasi (`13.163` asersi) dan `rubocop` (`26` berkas, `0` pelanggaran) lulus bersih.
- Output `outputs/day005h/results.json`, 6 grafik SVG `outputs/figures/day005h_*.svg`, `projects/ruby/day005h_backtest_overfitting/README.md`, serta [day005h_quant_backtest-overfitting.md](project-explanations/day005h_quant_backtest-overfitting.md).

## [Hari 005g] — 2026-11-23 · Perbaikan day005f + Multilevel Monte Carlo (Giles) untuk CIR & Heston (Ruby, Bonus)

### Diperbaiki (Bagian A — `day005f`)
- `projects/ruby/day005f_cdo_tranche/`: meng-commit `base_corr.rb` deterministik (kuadratur 2D eksak dengan substitusi $u = \sqrt{w}$ dan uji konvergensi `< 1e-8`: `1,09e-9`), memperkuat `c4` dengan arah *skew* dan monotonisitas EL *tranche* ternormalisasi serta sensitivitas $(\nu, \rho)$, memperkuat `c5` dengan pemindaian grid $\rho \in [0{,}01, 0{,}94]$ untuk seluruh 5 *tranche* (`0-3%`: 1, `3-7%`: 1, `7-10%`: 0, `10-15%`: 2 akar, `15-30%`: 1), dan menambah 6 unit tes baru (`16` tes, `48` asersi).

### Ditambahkan (Bagian C — `day005g`)
- Project Ruby murni mandiri `projects/ruby/day005g_mlmc_heston/` (topik cadangan sesuai Langkah 0.6 karena Hull–White telah dikerjakan penuh di `day002i`; nol gem numerik, `SplitMix64`, `Process.fork` + `IO.pipe` + `Marshal` deterministik, JSON penuh, 5 grafik SVG):
  - **Bagian 1**: Validasi *anchor* Heston CF (*little trap* + Gil–Pelaez 96-titik: batas BS $\xi=0$ galat `4,26e-14`, *put-call parity* `7,11e-15`, konvergensi kuadratur `8,23e-10`) dan momen eksak CIR.
  - **Bagian 2**: Regresi laju kuat $\beta$ & lemah $\alpha$: GBM Euler (`0,982`) vs Milstein (`1,945`), Heston *Full-Truncation* Feller OK (`1,276`) vs Feller dilanggar (`0,821`).
  - **Bagian 3**: Hukum kompleksitas Giles ($\varepsilon \in \{0{,}08..0{,}01\}$): kemiringan biaya log-log `GBM Milstein = 1,876` (teori `2,000`, *speedup* `14,90x`), `Heston MLMC = 2,471` (pra-asimtotik Giles `2,475`, asimtotik `2,172`) vs `Std MC = 3,000`, dengan fakta terukur Heston `speedup = 0,1837x` (`5,44x` lebih lambat pada $\varepsilon=0{,}01$; *crossover* ekstrapolasi $\varepsilon^* \approx 5{,}05 \times 10^{-4}$..$1{,}59 \times 10^{-3}$), serta kecocokan harga dalam `0,65 SE`.
  - **Bagian 4**: *Antithetic coupling* Giles–Szpruch pada Heston (`level 3` = `1,253x`, `level 4` = `1,216x`, status `belum tercapai` terhadap target literatur `>= 2,0x` karena skema dasar Euler *full-truncation* tanpa suku kuadratik Milstein).
- Runner utama `day005g` selesai dalam `74,50s` (`Ruby 4.0.9`, YJIT `ON`, `10` *worker*; `7/8` `tercapai`, `c6 = belum tercapai`), `21` tes unit/integrasi (`34` asersi) dan `rubocop` lulus bersih.
- Output `outputs/day005g/results.json`, 5 grafik SVG `outputs/figures/day005g_*.svg`, `projects/ruby/day005g_mlmc_heston/README.md`, serta [day005g_quant_mlmc-heston.md](project-explanations/day005g_quant_mlmc-heston.md).

## [Hari 005f] — 2026-11-23 · Perbaikan day005e + Valuasi Tranche CDO: LHP Vasicek, Rekursi ASB & Exponential Twisting (Ruby, Bonus)

### Diperbaiki (Bagian A — `day005e`)
- `projects/ruby/day005e_copula_tail_dependence/`: mendesain ulang `c2_tail_dependence_convergence` dengan $\chi_{\mathrm{Gauss}}(u)$ eksak berhingga via integral ortan bivariat ($\le 3$ SE di $u \in \{0{,}9, 0{,}99, 0{,}999\}$ + monoton turun hingga $1-10^{-10}$) dan SE binomial untuk $t$/Clayton/Gumbel (`11/11` kriteria `tercapai`), mengembalikan `.rubocop.yml` ke `ClassLength: 120` dengan memecah `Runner` menjadi `Criteria`, `Figures`, dan `RunnerParts`, serta menambahkan kriteria `c10_non_gaussian_euler_allocation` dan `c11_aic_selects_t_copula`.

### Ditambahkan (Bagian C — `day005f`)
- Project Ruby murni mandiri `projects/ruby/day005f_cdo_tranche/` (topik cadangan sesuai Langkah 0.6 karena Dupire/SSVI telah dikerjakan penuh di `day003a` dan Heston di `day002e`; nol gem numerik, `SplitMix64`, `Process.fork` + `IO.pipe` + `Marshal` deterministik, JSON penuh, 5 grafik SVG) pada portofolio CDO sintetis 125-nama 5Y:
  - **Bagian 1**: Rekursi DP Andersen–Sidenius–Basu vs binomial eksak (`max_err = 7,18e-16`) dan *brute-force* $2^{12}$, bentuk tertutup LHP Vasicek via ortan bivariat (`max_err = 2,52e-10`), serta studi granularitas $N \in \{25, 50, 125, 500\}$ ($O(1/N)$).
  - **Bagian 2**: Valuasi 6 *tranche* standar homogen vs heterogen 5-sektor dengan konservasi PV proteksi (`6,94e-18`).
  - **Bagian 3**: Kalibrasi *base correlation* terhadap pasar Student-$t(\nu=4)$ via kuadratur 2D eksak (`[0,455; 0,438; 0,432; 0,424; 0,410]`) dan diagnosis kegagalan akar *compound correlation* pada *tranche senior*.
  - **Bagian 4**: *Conditional MC + Exponential Twisting* dua-tahap Glasserman–Li ($R=100\mathrm{k}$: `max |is_z| = 2,64`, penghematan variansi `309,1x` pada *super-senior* `30-100%`; `353,5x` pada `extras` $R=400\mathrm{k}$).
  - **Bagian 5**: *Model risk* *spread* `15-30%` & `30-100%` di bawah Gaussian vs Student-$t(\nu=4)$ vs *random recovery* prosiklis (`27,44` vs `74,31` vs `63,46 bps` pada `15-30%`).
- Runner utama `day005f` selesai dalam `24,59s` (`Ruby 4.0.9`, YJIT `ON`, `10` *worker*; `c1`–`c9 = tercapai` semua), `10` tes unit/integrasi (`31` asersi) dan `rubocop` lulus bersih.
- Output `outputs/day005f/results.json`, 5 grafik SVG `outputs/figures/day005f_*.svg`, `projects/ruby/day005f_cdo_tranche/README.md`, serta [day005f_quant_cdo-tranche.md](project-explanations/day005f_quant_cdo-tranche.md).

## [Hari 005e] — 2026-11-23 · Perbaikan day005d + Kopula, Ketergantungan Ekor & Alokasi Modal Euler (Ruby, Bonus)

### Diperbaiki (Bagian A — `day005d`)
- `projects/rust/day005d_merton_jump_diffusion/`: menambahkan kalibrasi kuotasi IV berderau $N(0, 10\text{ bps})$ ($R=200$, perbandingan 1-maturitas $T=0{,}25$ vs 2-maturitas $\{0{,}25, 1{,}0\}$ dan $\mathrm{Corr}(\hat{\lambda}, \hat{\delta}_J)$), studi multi-start `24` titik acak (`23/24` global), uji kuantitatif skew jangka pendek orde pertama $O(T)$ pada $T \in \{0{,}01, 0{,}02, 0{,}05\}$ (`max_rel_err = 0,082%`), serta pembanding BS variansi-sama ($\sigma_{\mathrm{eq}} = 28{,}37\%$, `ES99 = 9,248` vs `4,017` BS-ATM vs `28,156` Merton).

### Ditambahkan (Bagian C — `day005e`)
- Project Ruby murni mandiri `projects/ruby/day005e_copula_tail_dependence/` ($\ge 3{,}3$, nol gem numerik, `SplitMix64`, `Process.fork` + `IO.pipe` + `Marshal` deterministik, JSON penuh, 7 grafik SVG) pada data sintetis:
  - **Bagian 1**: Kopula Gaussian, Student-$t$, Clayton, Gumbel ($R=200\mathrm{k}$: `max |tau_z| = 1,66`, `spearman_err = 8,83e-5`), konvergensi $\chi(u)$ ($u \in \{0{,}9, 0{,}99, 0{,}999\}$), serta ortan bivariat Gauss–Legendre (`max_err = 5,55e-17`).
  - **Bagian 2**: Alokasi modal Euler ES99 ($d=5$: `sum_err = 0`, `fd_err = 5,50e-12`, `max_mc_z = 0,53`; perbandingan vs kopula $t_4$ dan Clayton).
  - **Bagian 3**: CML vs inversi $\tau$ ($R=500$, $n \in \{250, 1000, 5000\}$), salah spesifikasi $t_4$ vs Gaussian (`LRT = 155,35`, bias ES99 `-7,44%`), dan perbandingan efek marginal vs kopula.
  - **Bagian 4**: *Importance sampling* pergeseran-rerata + *conditional MC* ($R=100\mathrm{k}$: `max |is_z| = 2,45`, penghematan variansi `2.604x`–`5.165.715x`).
- Runner utama `day005e` selesai dalam `9,20s` (`Ruby 4.0.9`, YJIT `ON`, `10` *worker*; `c1, c3–c9 = tercapai`, `c2 = belum tercapai` dilaporkan dengan diagnosis bias berhingga Gaussian), `12` tes unit/integrasi (`62` asersi) dan `rubocop` lulus bersih.
- Output `outputs/day005e/results.json`, 7 grafik SVG `outputs/figures/day005e_*.svg`, `projects/ruby/day005e_copula_tail_dependence/README.md`, serta [day005e_quant_copula-tail-dependence.md](project-explanations/day005e_quant_copula-tail-dependence.md).

## [Hari 005d] — 2026-11-23 · Jump-Diffusion Merton: Seri Eksak, PIDE Crank–Nicolson, Kalibrasi & Risiko Ekor (Rust, Bonus)

### Ditambahkan
- Project Rust stable mandiri `projects/rust/day005d_merton_jump_diffusion/` (topik cadangan `ROADMAP 075` karena topik utama telah dikerjakan penuh di `day002b`; *zero external crates*, `xoshiro256++`, `std::thread::scope` deterministik blok tetap `64`, JSON 17-digit, 4 grafik SVG) pada data sintetis:
  - **Bagian 1**: Seri eksak Merton (1976), *put-call parity*, dan simulasi MC paralel ($R=400\mathrm{k}$: `martingale_z = 0,49`, `max |call_z| = 1,10`, `max |put_z| = 0,75`).
  - **Bagian 2**: Pemecah PIDE Crank–Nicolson + 4 *half-step* Rannacher dengan algoritma Thomas $O(N)$ dan integral lompatan *predictor-corrector* (`401x200`: `max_abs_err = 0,00303`, konvergensi `15,43x`).
  - **Bagian 3**: *Implied volatility smile* ($T \in \{0{,}25; 1; 2\}$: *skew* `35,40%` di $K=80$ vs `22,55%` ATM) dan kalibrasi Nelder–Mead 4D (`price_rmse = 1,15e-8`, `214` iterasi).
  - **Bagian 4**: Risiko ekor *put* OTM (`K=70` rasio `191,14x` BS) dan *short put* 1 bulan ($K=90$, $R=200\mathrm{k}$: `ES99 = 28,156` vs `4,017` BS).
- Runner utama `day005d` selesai dalam `0,49s` pada `10` *worker* (`c1`–`c7 = tercapai` semua), `9` unit test + `2` integration test lolos.
- Output `outputs/day005d/results.json`, 4 grafik SVG `outputs/figures/day005d_*.svg`, `projects/rust/day005d_merton_jump_diffusion/README.md`, serta [day005d_quant_merton-jump-diffusion.md](project-explanations/day005d_quant_merton-jump-diffusion.md).

## [Hari 005c] — 2026-11-23 · Perapihan dokumen day004j + Kriteria Kelly & Risiko Drawdown (Rust, Bonus)

### Diperbaiki (Bagian A — `day004j`)
- `docs/project-explanations/day004j_quant_deep-hedging.md`: mengklarifikasi ketidakabsahan pembandingan delta cabang *entropic* ($\lambda = 10$) terhadap delta BS beserta perbandingan $\rho_{10}(\mathrm{NN}_{\mathrm{ent}}) = 0{,}000329$ vs $\rho_{10}(\Delta_{\mathrm{BS}}) = 0{,}000159$, merinci pembagian waktu runner (`986,90s` vs `1305,62s`: `part2 359s`, `part3 591s`, `part4 35s`) serta saran perbaikan anggaran, dan menandai eksponen skala pita (`1,7131`, rasio `0,1300`) sebagai FAKTA sedangkan dugaan `tanh` sebagai HIPOTESIS.

### Ditambahkan (Bagian C — `day005c`)
- Project Rust stable mandiri `projects/rust/day005c_kelly_drawdown/` (*zero external crates*, `xoshiro256++`, `std::thread::scope` deterministik blok tetap `64`, JSON 17-digit, 7 grafik SVG) pada data sintetis:
  - **Bagian 1**: Kelly biner diskret ($p=0{,}55, b=1$: $f^* = 0{,}10$, $f_{\mathrm{be}} = 0{,}198668$), simulasi $N=1000, R=100{,}000$ (`max_z_log = 2,07`), serta perbandingan terhadap aproksimasi kontinu.
  - **Bagian 2**: Peluang sentuh GBM eksak horizon tak hingga ($P_\infty(0{,}5 \mid f^*) = 0{,}257028$; kasus $r=0$: $P=x$) dan hingga, serta estimator jembatan Brownian tanpa bias pada `72` kombinasi ($R=10^6$: `72/72` dalam `3 SE`, `max |z| = 2,80`).
  - **Bagian 3**: Kelly fraksional $c(2-c)$ (`max_err = 4,44e-16`), peta $(c, T)$, dan akar numerik vs tertutup (`diff = 1,11e-16`).
  - **Bagian 4**: Kerugian galat estimasi $1/(2T_{\mathrm{est}})$ (`max |z| = 1,41`, invarian terhadap $\mu,\sigma$), penyusutan optimal $c^* = 0{,}2381$ (`+0,07561 ± 0,00036` vs plug-in penuh; peluang sentuh 20-thn `5,82%` vs `53,71%`), dan ekor berat Student-$t(\nu=4)$ (`g(2f*) = 0,02408 < 0,03000`, `ruin = 0,387%`).
- Runner utama `day005c` selesai dalam `50,47s` pada `10` *worker* (`c1, c3, c4, c5, c7, c8 = tercapai`, `c2` & `c6 = belum tercapai` dilaporkan dengan diagnosis).
- Output `outputs/day005c/results.json`, 7 grafik SVG `outputs/figures/day005c_*.svg`, `projects/rust/day005c_kelly_drawdown/README.md`, serta [day005c_quant_kelly-drawdown.md](project-explanations/day005c_quant_kelly-drawdown.md).

## [Hari 005b] — 2026-11-23 · Harga Opsi Eropa Black–Scholes & Put-Call Parity (Rust, Bonus)

### Ditambahkan
- Project Rust stable mandiri `projects/rust/day005b_black_scholes/` (*edition 2021*, tanpa crate eksternal, deterministik):
  - Fungsi harga opsi Eropa *call* & *put* Black–Scholes (termasuk *dividend yield* kontinu $q$) dan residual *put-call parity* di `src/lib.rs`.
  - Pendekatan CDF normal standar via rumus rasional Abramowitz & Stegun (7.1.26) dengan galat $\Phi(x) \le 7{,}5 \times 10^{-8}$.
  - Validasi input tanpa *panic* (`BsError`) serta penanganan batas $T = 0$ dan $\sigma = 0$ tanpa pembagian nol.
  - Tiga tes unit (kasus acuan $S=100, K=100, r=0{,}05, q=0, \sigma=0{,}20, T=1$; *put-call parity*; validasi & kasus batas) dan contoh CLI di `src/main.rs`.
- Dokumentasi `projects/rust/day005b_black_scholes/README.md` serta [day005b_quant_black-scholes.md](project-explanations/day005b_quant_black-scholes.md).

## [Hari 005a] — 2026-11-23 · Perbaikan day004h/day004j/day004k + Model SABR (Rust, Bonus)

### Diperbaiki (Bagian A — `day004h`, `day004j`, `day004k`)
- `projects/rust/day004k_amm_lvr`: mengganti `c6` dengan *anchor* kuantitatif LVR lompatan Merton $(\lambda/2)(e^{m+s^2/2}-2e^{m/2+s^2/8}+1) = 0{,}007061$ (`z = 0,22`), memperketat `c4` dengan $\Delta \in \{0{,}25\text{s}; 0{,}5\text{s}\}$ dan ekstrapolasi Richardson ($p_{\mathrm{OLS}} = 0{,}464$, `0,07992 ± 0,00021` vs `0,07988`, `z = 0,21`), serta uji sapuan $\varphi \in \{0{,}0005; 0{,}003; 0{,}01\}$.
- `projects/rust/day004j_deep_hedging`: mengimplementasikan *persistent worker pool* + penyimpanan `StepCache` (`1,47x` lebih cepat pada 20 iterasi dan bit-identik), melatih cabang entropik & biaya hingga `10.000` iterasi, menguji hipotesis aktivasi `tanh` vs `ReLU`/`Softplus` dan fitur $(\delta_{t-1} - \delta_{\mathrm{BS}})$ di `run_extras`, serta menyelaraskan ambang `c10_runtime_under_600s` ke `600s`.
- `projects/rust/day004h_rough_bergomi`: memindahkan eksperimen tambahan ke `run_extras` (`run` kembali `123,25s < 300s`, `c9 = tercapai`), menganalisis 10 *seed* `true_rmse_grid` (median vs IQR), menghitung indeks ekor Hill ($\hat{\alpha}_{S_T} \approx 17$, $\hat{\alpha}_{\mathrm{payoff}^2} \approx 2{,}2..3{,}2$), dan menguji kalibrasi estimator BS bersyarat ($|\mathrm{bias}_H| = 0{,}0108$ dan $0{,}0228 < 0{,}03$).

### Ditambahkan (Bagian C — `day005a`)
- Project Rust stable mandiri `projects/rust/day005a_sabr/` (*zero external crates*, `xoshiro256++`, `std::thread::scope` deterministik blok tetap `64`, JSON 17-digit, 7 grafik SVG) pada data sintetis:
  - **Bagian 1**: CEV $\nu=0$ tiga cara (rumus tertutup $\chi'^2$ non-sentral, PDE Crank–Nicolson 1D, MC transisi eksak & Euler) cocok (`max_pde_rel_err = 1,67e-8`, `max_z <= 1,78`), serta verifikasi massa absorpsi $T=10$ (`0,000514` vs `0,000475`) dan martingal.
  - **Bagian 2**: Ekspansi Hagan (2002) dengan Taylor stabil $z/x(z)$ dan *anchor* $\nu \to 0$ vs CEV eksak (`0,001..1,90 bps`).
  - **Bagian 3**: Referensi $\rho = 0$ kondisional CEV MC ($R = 10^6$), konvergensi Euler *Absorbing* vs *FullTruncation* ($n \in \{100, 200, 400, 800\}$, $R = 2\times 10^6$), *smile* referensi $\rho = -0{,}3$ ($T=1$: `max_abs = 7,00 bps` pada strike likuid), dan peta galat $(\nu, \rho, \beta)$.
  - **Bagian 4**: Kepadatan implisit Hagan Breeden–Litzenberger negatif pada sayap bawah ($T=5$: massa negatif `0,00156`; $T=10$: `0,10195`, `min_q = -88,59`) dan dampaknya pada digital/*butterfly* *floor* `0,5%`.
  - **Bagian 5**: Kalibrasi Nelder–Mead + korelasi Hessian numerik, non-identifiabilitas $\beta \in \{0, 0{,}5, 1\}$, lindung nilai harian CRN ($R = 20.000$) dengan delta Bartlett (`SD = 6,318 bps`) vs naif (`6,598 bps`, `var_diff = -3,62 ± 0,18 bps²`) vs Black ATM vs CEV, serta salah spesifikasi $\beta \in \{0{,}3; 0{,}5; 0{,}8\}$.
- Runner utama `day005a` selesai dalam `222,91s` pada `8` worker (`c1..c5, c7, c8 = tercapai`, `c6 = belum tercapai` pada $T=5$).
- Output `outputs/day005a/results.json`, `outputs/day005a/extras.json`, 7 grafik SVG `outputs/figures/day005a_*.svg`, `projects/rust/day005a_sabr/README.md`, serta [day005a_quant_sabr.md](project-explanations/day005a_quant_sabr.md).

## [Hari 004k] — 2026-11-22 · Perbaikan day004h/day004i/day004j + AMM, LVR & Likuiditas Terkonsentrasi (Rust, Bonus)

### Diperbaiki (Bagian A — `day004h`, `day004i`, `day004j`)
- `projects/rust/day004i_pairs_trading`: menyelaraskan simulasi null Johansen (`q95_h0_drift = 18,000 ± 0,097`) dengan DGP aplikasi ber-*drift* (`0,05*dt`), melaporkan ukuran & *power* baik pada `q95_h0_drift` (`size = 4,87%`, `power = 81,56%`) maupun pada kuantil literatur `15,41` (`size = 10,99%`, `power = 93,46%`), serta menulis ulang §8 dokumen penjelasan.
- `projects/rust/day004h_rough_bergomi`: menguraikan `rmse_at_true` pada grid `(n, R) ∈ {100, 250, 500} × {20k, 80k, 200k}` (membuktikan fluktuasi `0,012..0,076` berasal dari ekor berat $\exp(\eta\widetilde{W})$ alih-alih diskretisasi $n$) dan menambahkan kalibrasi ruang harga berbobot vega (`fit_vega_price`: $\hat{H} = 0,1775$).
- `projects/rust/day004j_deep_hedging`: melatih hingga konvergen (`10.000` iterasi, penskalaan input, bias keluaran `0,5`, *batch* `1024`/`4096`), menambahkan diagnostik residual $\delta = \delta_{\mathrm{BS}} + \mathrm{NN}$, grid pita halus `401` titik dengan SE 3 *seed*, kesepakatan kebijakan DP berbobot frekuensi kunjungan, serta metrik terukur `c1`–`c2`.

### Ditambahkan (Bagian C — `day004k`)
- Project Rust stable mandiri `projects/rust/day004k_amm_lvr/` (*zero external crates*, `xoshiro256++`, `std::thread::scope` deterministik blok tetap `64`, JSON 17-digit, dan 8 grafik SVG) pada data sintetis:
  - **Bagian 1**: CPMM $x\cdot y = L^2$, *impermanent loss* $\mathrm{IL}(r)$, *anchor* martingal GBM ($R = 10^6$, `max_z = 2,31`), LVR *pathwise* Itô (`0,0065501` vs `0,0065516`), dan sapuan $\sigma \in \{0{,}3; 0{,}5; 0{,}8; 1{,}2\}$.
  - **Bagian 2**: Arbitrase per blok (biaya pada input) via kuadratur Gauss–Legendre 64-titik vs MC pada $12$ kombinasi $(\Delta, \varphi)$ (`max_z = 1,65`), limit kontinu gerak Brownian terefleksi ($\Delta = 1\text{s}$ biaya LP `0,07761` vs limit `0,07988`), dan titik impas aliran *noise trader* Poisson.
  - **Bagian 3**: Likuiditas terkonsentrasi $[P_a, P_b]$ ($\pm 5\%, \pm 10\%, \pm 30\%$, penuh), efisiensi modal $\mathrm{CE}$, waktu keluar ODE GBM, serta perbandingan CRN 5 strategi 30 hari (*delta-hedged ±10%* `Sharpe = 2,009`).
  - **Bagian 4**: Salah spesifikasi lompatan Merton (`0,08678` > `0,07973`), volatilitas stokastik CIR (`0,07968` vs `0,07968`, `corr = 0,923`), *drift* $\mu \in \{-0{,}5; +0{,}5\}$ (`LVR/int V dt = 0,08001` & `0,07999`), dan ukuran blok `600s`.
- Seluruh runner `day004k` selesai dalam `86,58s` pada `8` worker (`c1..c9 = tercapai` semua).
- Output `outputs/day004k/results.json`, 8 grafik SVG `outputs/figures/day004k_*.svg`, `projects/rust/day004k_amm_lvr/README.md`, serta [day004k_quant_amm-lvr.md](project-explanations/day004k_quant_amm-lvr.md).

## [Hari 004j] — 2026-11-22 · Perbaikan day004h/day004i + Deep hedging dari nol (Rust, Bonus)

### Diperbaiki (Bagian A — `day004i` & `day004h`)
- `projects/rust/day004i_pairs_trading`: mengaudit kuantil Johansen (`q95 = 17,997`) via fungsional Brownian asimtotik (`T=5000, R=20000`: Case 1 `12,334`, Case 2 `20,266`, Case 3 drift-nol `17,987`, Case 3 drift-non-nol `15,407` vs Osterwald–Lenum `15,41`) dan sapuan $T \in \{250..5000\}$ (`18,21..18,14`), mengganti uji ukuran tautologis dengan simulasi *out-of-sample* independen (`EG = 4,77%`, `Johansen = 4,87%`, `power = 81,57%`), serta mengembalikan toleransi `c1` ke `< 1e-8` dengan polesan Newton dan cek skor/Hessian (`mle_max_diff = 8,17e-9`, `rel_score = 3,02e-9`).
- `projects/rust/day004h_rough_bergomi`: mendiagnosis bias kalibrasi `c8` melalui uji *self-consistency* (`H_hat = 0,100001, rmse = 1,18e-7`), replikasi dengan bank CRN berbeda per replikasi (`rep_diff_crn`: `bias_H = +0,1027 ± 0,0514`), dan uji $R_{\mathrm{calib}} = 80.000$ (`H_hat = 0,1657`), serta memperbarui bagian "Batasan" dokumen penjelasan.

### Ditambahkan (Bagian C — `day004j`)
- Project Rust stable mandiri `projects/rust/day004j_deep_hedging/` (*zero external crates*, `xoshiro256++`, `std::thread::scope` deterministik blok tetap `64`, JSON 17-digit, dan 8 grafik SVG) pada data sintetis:
  - **Bagian 1**: MLP $2 \times 32$ (`tanh`, keluaran $\delta \in \mathbb{R}$), BPTT manual 30 langkah, cek gradien vs Richardson (`< 1e-6`) untuk MSE, risiko entropik, dan CVaR 95%, Adam buatan sendiri (`< 1e-6`), dan determinisme `1` vs `8` *worker* bit-identik.
  - **Bagian 2**: Dunia Black–Scholes tanpa biaya ($R = 100.000$): `p0_MSE = 0,039831` vs `BS = 0,039878`, `SD_BS = 0,006244` vs aproksimasi `0,006307`, namun lantai optimisasi 360 iterasi membuat `c3` (`L2 = 0,0880 > 0,03`) dan `c4` (`slope_NN = -0,141`) berstatus `belum tercapai`.
  - **Bagian 3**: Biaya transaksi proporsional $c \in \{0{,}0005; 0{,}001; 0{,}002\}$, $\lambda \in \{10, 50\}$: pita *no-trade* NN (`eksponen c = 0,2925` vs Whalley–Wilmott `1/3`), *anchor* DP eksak $n=10$ (`DP_indiff = 0,04191` vs `NN_sim = 0,04202`), dan keunggulan NN (`rho10 = 0,002137`) atas BS (`0,002393`).
  - **Bagian 4**: Simulator Heston QE (`z(E[S_T]=1) = 0,51`; `z(MC vs Fourier) = 0,31`) dan perbandingan 5 strategi (`NN+v_t` memangkas `SD` dari `0,01563` menjadi `0,01356`).
  - **Bagian 5**: Uji salah spesifikasi (`BS -> Heston`, `BS -> BS(0,3)`, biaya salah) dan ablasi (lebar `4`, tanpa $\delta_{t-1}$ yang menghilangkan pita *no-trade* dan menaikkan `rho10` ke `0,00519`).
- Seluruh runner penuh selesai dalam `132,64s` pada `8` worker (`c1, c2, c5..c10 = tercapai`; `c3, c4 = belum tercapai`).
- Output `outputs/day004j/results.json`, 8 grafik SVG `outputs/figures/day004j_*.svg`, `projects/rust/day004j_deep_hedging/README.md`, serta [day004j_quant_deep-hedging.md](project-explanations/day004j_quant_deep-hedging.md).

## [Hari 004i] — 2026-11-22 · Perbaikan day004h + Cointegration & trading OU optimal (Rust, Bonus)

### Diperbaiki (Bagian A — `day004h`)
- `projects/rust/day004h_rough_bergomi`: mengganti skema hybrid $\kappa=1$ menjadi titik optimal $L^2$ Bennedsen–Lunde–Pakkanen (`BlpOptimal`) serta membandingkan tiga skema (`RightEndpoint`, `BlpOptimal`, `VarianceMatching`) pada $n \in \{50, 100, 250, 500\}$. Dengan `BlpOptimal`, bias variansi turun dari `+4,94%..+3,23%` menjadi `-0,062%..-0,039%`, sehingga `c3`, `c4`, dan `c6` kini seluruhnya `tercapai`. Menambahkan uji martingal Cholesky + kuantil bootstrap (`[0,9732; 0,9998; 1,0271]`), ambang Bonferroni a priori pada `c2` (`z <= 4,12` untuk `1275` pasangan), serta pemisahan 3 kondisi kalibrasi (`same-scheme-same-n`: `bias_H = +0,0719`, mengonfirmasi `c8 = belum tercapai` akibat *ridge* identifiabilitas datar).

### Ditambahkan (Bagian C — `day004i`)
- Project Rust stable mandiri `projects/rust/day004i_pairs_trading/` (*zero external crates*, `xoshiro256++`, `std::thread::scope` berurutan blok tetap, JSON 17-digit, dan 7 grafik SVG) pada data sintetis:
  - **Bagian 1**: Transisi OU eksak vs Euler, MLE AR(1) tertutup vs numerik (`diff = 2,12e-6`), dan replikasi $R=20.000$ bias Kendall/Marriott–Pope (rasio `1,04..1,09`).
  - **Bagian 2**: Waktu lintasan pertama tertutup vs BVP beda hingga vs MC ($R=200.000$), ambang optimal Bertram $a^* = -0{,}01574$ (`0,445 SD`, `R(a*) = 0,2044`), serta verifikasi MC 200 tahun (`z = 0,64`).
  - **Bagian 3**: Simulasi Dickey–Fuller ($T=1000, R=10^6$: `q05 = -2,8641`), Engle–Granger ($R=20.000$: `q05 = -3,3531`), regresi palsu (`93,26%` $|t|>1{,}96$), dan Johansen trace bivariat (`q95 = 17,997`, *power* 2 thn `81,81%`, `c5 = belum tercapai`).
  - **Bagian 4**: *Backtest* bergulir ($W \in \{1, 2, 5\}$, $R=2.000\times 10$ thn) dengan akuntansi MTM harian, filter EG pada pasangan *spurious* (`-0,00061 ± 0,00038`), *structural break* ($\theta: 16 \to 2$), dan kurva biaya.
- Seluruh runner penuh selesai dalam `3,67s` pada `8` worker (`c1..c4, c6..c8 = tercapai`; `c5 = belum tercapai`).
- Output `outputs/day004i/results.json`, 7 grafik SVG `outputs/figures/day004i_*.svg`, `projects/rust/day004i_pairs_trading/README.md`, serta [day004i_quant_pairs-trading.md](project-explanations/day004i_quant_pairs-trading.md).

## [Hari 004h] — 2026-11-22 · Perbaikan day004g + Volatilitas kasar: model rough Bergomi (Rust, Bonus)

### Diperbaiki (Bagian A — `day004g`)
- `projects/rust/day004g_market_making`: memperbaiki orde ekstrapolasi Richardson di `c4` menggunakan orde empiris $p = 1{,}017$ dengan propagasi SE CRN (`CE_Rich = 66,2684 ± 0,0295`, $z = 0{,}23$); menambahkan simulator waktu-kontinu *Lewis–Shedler thinning* bebas-$\Delta t$ (`CE_cont = 66,2763`, cocok dengan Forward Kolmogorov pada $\mathbb{E}[\mathrm{exec}]$ $z=0{,}59$ dan $\operatorname{Var}(q_T)$ $z=0{,}60$, `c5 = tercapai`) serta audit $\Delta t \in \{10^{-3}; 5\times 10^{-4}; 2{,}5\times 10^{-4}; 1{,}25\times 10^{-4}\}$; memperbaiki `exact_as_hjb` agar memakai tabel penuh $(t, q)$ sehingga reduksi $\varepsilon = 0$ pulih (`66,3191` vs `66,3191`); dan menyelaraskan `c1`–`c9` antara `results.json` dan dokumen penjelasan.

### Ditambahkan (Bagian C — `day004h`)
- Project Rust stable mandiri `projects/rust/day004h_rough_bergomi/` (*zero external crates*, `xoshiro256++`, `std::thread::scope` berurutan blok tetap, JSON 17-digit, dan 7 grafik SVG) yang mengimplementasikan:
  - **Bagian 1**: Kovarians Volterra eksak (${}_2F_1$ vs Gauss–Legendre bebas singularitas, `max_rel_err = 2,12e-10`), referensi Cholesky $2n\times 2n$ ($n=50, R=200.000$), dan bias skema hybrid $\kappa = 1$ titik kanan ($n \in \{50, 100, 250, 500\}$, $H=0{,}1$ dan $H=0{,}5$).
  - **Bagian 2**: Uji martingal variansi/harga & *variance swap* ($T \in \{0{,}25; 0{,}5; 1\}$), serta pembanding MC log-Euler vs estimator bersyarat Black–Scholes (`VRF = 2,16..308,4`, paritas put-call `2,45e-14`).
  - **Bagian 3**: *Smile* $13 \times 5$ ($R=100.000, n=250$), regresi *skew* ATM $\ln\psi(T)$ vs $\ln T$ (`slope(H=0,1) = -0,3850 ± 0,0128`; `slope(H=0,5) = -0,0404 ± 0,0253`), dan pembanding Markov $H=0{,}5$ di $T=0{,}1$.
  - **Bagian 4**: Pasar sintetis ($R=200.000, n=500$), kalibrasi Nelder–Mead CRN ($\hat{H} = 0{,}1629$), profil $H$, korelasi Hessian $H$–$\eta$ (`0,7915`), 30 replikasi (`bias_H = +0,0652`), serta salah spesifikasi horizon ($T \ge 0{,}5$ vs $T \le 0{,}1$).
- Seluruh runner penuh selesai dalam `104,42s` pada `8` worker (`c1, c2, c5, c7, c9 = tercapai`; `c3, c4, c6, c8 = belum tercapai` akibat bias variansi $+3{,}2\%$–$+4{,}9\%$ skema hybrid titik kanan).
- Output `outputs/day004h/results.json`, 7 grafik SVG `outputs/figures/day004h_*.svg`, `projects/rust/day004h_rough_bergomi/README.md`, serta [day004h_quant_rough-bergomi.md](project-explanations/day004h_quant_rough-bergomi.md).

## [Hari 004g] — 2026-11-22 · Market making optimal: Guéant–Lehalle–Tapia, HJB, dan salah spesifikasi (Rust, Bonus)

### Ditambahkan
- Project Rust stable mandiri `projects/rust/day004g_market_making/` (*zero external crates*, `xoshiro256++`, `std::thread::scope`, JSON 17-digit, dan 7 grafik SVG) yang mengimplementasikan:
  - **Bagian 1**: Solver eksak sistem-$v$ (Jacobi siklik, *scaling-and-squaring*, RK4+Richardson dengan galat relatif maks `2,80e-11`), solver HJB FD nonlinier independen bentuk $\theta$ (orde konvergensi $\Delta t$ `1,001–1,003`), verifikasi simetri/monotonisitas kuotasi, celah spektral (`1,277564`) & batas stasioner, hampiran *closed-form*, serta sensitivitas $Q \in \{10, 20, 40\}$.
  - **Bagian 2**: Simulator Monte Carlo CRN ($R = 200.000$ lintasan, $\mathrm{CE}_{\mathrm{Rich}} = 66{,}2880$ vs eksak `66,2616`, $z = 2{,}20$), pembanding ODE Forward Kolmogorov (`c5 = belum tercapai` karena efek diskretisasi eksekusi simultan Bernoulli pada $\Delta t = 5\times 10^{-4}$ yang meredam $\operatorname{Var}(q_T)$ menjadi `60,6586` vs `62,0174`), dan perbandingan 4 strategi (`glt_optimal`, `avellaneda_stoikov`, `symmetric_fixed`, `tight_fixed_0.1`).
  - **Bagian 3**: Studi salah spesifikasi (*adverse selection* $\varepsilon \in \{0; 0{,}02; 0{,}05; 0{,}1\}$, salah $\sigma \in \{0{,}5\sigma, \sigma, 1{,}5\sigma\}$, salah intensitas $A/k \pm 30\%$, dan *drift* $\mu \in \{-0{,}2; -0{,}05; 0; 0{,}05; 0{,}2\}$).
- Seluruh runner penuh selesai dalam `47,05s` pada `8` worker (`c1..c4, c6..c8 = tercapai`, `c5 = belum tercapai`).
- Output `outputs/day004g/results.json`, 7 grafik SVG `outputs/figures/day004g_*.svg`, `projects/rust/day004g_market_making/README.md`, serta [day004g_quant_market-making.md](project-explanations/day004g_quant_market-making.md).

## [Hari 004f] — 2026-11-22 · Audit tail day004e + port lengkap day004a–4e & valuasi proyek baterai AFNS (Rust, Bonus)

### Diperbaiki (Bagian A)
- `projects/go/day004e_power_storage`: menambahkan audit *expanded-tail* `[-1000, 5000] EUR/MWh` (`expanded_tail_38`) dengan *common random numbers* kondisional pada bin interior `fine_33`, uji unit overflow/representasi tail/CRN/determinisme worker, serta menjaga seluruh hasil asli `17/33` bin tidak berubah. Massa tail atas turun dari `1,400%` menjadi `0,000%` dan bawah dari `0,431%` menjadi `0,00616%`, dengan kesepakatan aksi evaluasi `100%` (`48.000/48.000`) terhadap `DP-33` (`selisih = 0,00 EUR/hari`).

### Ditambahkan (Bagian C)
- Project Rust stable mandiri `projects/rust/day004f_energy_transition/` (*zero external crates*, `xoshiro256++`, `std::thread::scope`, JSON 17-digit, dan 7 grafik SVG) yang memport penuh `day004a`–`day004e` serta menambahkan **Panel F** (valuasi proyek baterai 10 tahun didiskontokan kurva AFNS, sensitivitas fit `4b/4c/4d`, durasi, konveksitas, KRD, PV01, dan *hedge* dua obligasi).
- Seluruh runner penuh selesai dalam `587,02s` pada `8` worker dengan `all_technical_passed = true` (`NPV_base = 734.424,80 EUR`, durasi `5,14y`).
- Output `outputs/day004f/results.json`, 7 grafik SVG `outputs/figures/day004f_*.svg`, `projects/rust/day004f_energy_transition/README.md`, serta [day004f_rust-afns-energy-storage.md](project-explanations/day004f_rust-afns-energy-storage.md).

## [Hari 004e] — 2026-11-22 · Harga listrik negatif dan arbitrase baterai (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ mandiri `projects/go/day004e_power_storage/` (*standard library* saja, tanpa impor `day004d`, `MainSeed = 2026112205`, `KernelSeed = 2026112215`, `EvalSeed = 2026112225`, `VerifySeed = 2026112235`, `xoshiro256++`, paralelisme deterministik per chunk, JSON 17-digit, dan 5 grafik SVG):
  - **Verifikasi Teknis (`c1..c8 = tercapai`, `0,70s` pada 8 worker)**: distribusi stasioner eksak 3-rezim (`residual = 1,04e-17`), momen *compound Poisson* & transisi satu langkah cocok terhadap MC `200k` dalam `4 MCSE`, Bellman 4 jam cocok terhadap enumerasi *brute force* `625` lintasan (`diff = 1,14e-13`), residual Bellman nol pada grid 17 & 33 bin, nol pelanggaran SoC/daya, nol pelanggaran dominasi *perfect foresight*, dan determinisme bit-identik `1` vs `8` worker.
  - **Model Pasar Listrik Sintetis (`N = 2.000` hari evaluasi, *burn-in* `60` hari = `48.000` jam)**: `5.608` jam harga negatif (`11,68%` di `1.240/2.000` hari), `2.061` jam `scarcity`, `1.029` jam `oversupply`, `204` jump non-nol (`8` jam multi-jump $N>1$), rentang harga `[-354,69, 1041,92] EUR/MWh`.
  - **Evaluasi Berpasangan CRN 5 Strategi (`4 MWh` baterai, degradasi `2 EUR/MWh`, nilai terminal `45 * 0.9 * soc`)**:
    - `stochastic_dp_33bin`: mean `236,09 ± 4,73 EUR/hari` (median `192,98`, VaR 5% `27,59`, ES 5% `+4,87`, EFC `0,932`, regret `252,30 ± 4,56`).
    - `stochastic_dp_17bin`: mean `231,90 ± 4,75 EUR/hari` (median `189,54`, VaR 5% `22,55`, ES 5% `-1,85`, EFC `1,020`, regret `256,49 ± 4,60`).
    - `threshold_policy` (`P <= 0` charge, `P >= 200` discharge): mean `116,01 ± 4,41 EUR/hari` (median `77,55`, VaR/ES 5% `0,00`, EFC `0,175`, regret `372,38 ± 4,40`).
    - `perfect_foresight_dp` (batas atas bersyarat): mean `488,39 ± 7,16 EUR/hari` (EFC `2,505`); `no_battery` (`0,00 EUR/hari`).
    - Selisih berpasangan 95%: `DP-33 - Threshold = +120,08 ± 3,47` (`[113,29, 126,88]`), `DP-17 - Threshold = +115,89 ± 3,49` (`[109,06, 122,73]`), dan `DP-33 - DP-17 = +4,19 ± 1,26 EUR/hari` (`[1,72, 6,66]`).
- Output `outputs/day004e/results.json`, 5 grafik SVG `outputs/figures/day004e_*.svg`, serta dokumen penjelasan [day004e_quant_power-storage.md](project-explanations/day004e_quant_power-storage.md).

### Diperbarui
- `README.md` (tabel progres baris `4e` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 004d — pemeliharaan A] — 2026-11-22 · Audit penyajian interval & ringkasan robust `day004d`

### Diperbaiki
- `projects/go/day004d_afns_interval_calibration`: menambahkan `physically_valid_covered_count` (memisahkan interval `kp11` yang mencakup nilai benar tetapi berbatas bawah non-positif: `Basic 0/100`, `Studentized 4/100`, `Wald 6/100`), ringkasan lebar robust (`q25_width`, `q75_width`, `iqr_width`, `p95_width`, `max_width`), median lebar per strata `kp11_hat`, `technical_scope_note` pada `acceptance_criteria`, uji unit `TestSummarizeMethodRobustAndPhysicalCounts`, serta klarifikasi pada README dan dokumen penjelasan (mempertahankan seluruh hasil simulasi mentah & cakupan identik bit-demi-bit).

## [Hari 004d] — 2026-11-22 · AFNS — interval bootstrap yang peka bias dan batas parameter (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ mandiri `projects/go/day004d_afns_interval_calibration/` (`zero external modules`, `MainSeed = 2026112204`, `xoshiro256++` + `SplitMix64`, paralelisme goroutine deterministik per *chunk*, `FastLogLikelihood` *steady-state* Riccati, penulis JSON 17 digit & 4 grafik SVG):
  - **Verifikasi Teknis (`c1..c5 = tercapai`, `369,51s` pada 8 worker)**: jangkar Gaussian gabungan $T=5$ (3 tenor) cocok terhadap Kalman + *RTS smoother* (`ll_diff = 5,68e-13`, `fast_ll_diff = 5,54e-13`, `mean_diff = 3,16e-16`, `cov_diff = 6,39e-17`), skor FD vs Richardson ($\le 3,38\cdot 10^{-7}$), determinisme bit-identik `1` vs `8` worker pada fixture cepat, serta `100/100` outer dan `10.000/10.000` refit bootstrap dalam bernilai *finite* (`0` interval invalid).
  - **Perbandingan Berpasangan 5 Metode Interval 95% (`H=20y`, $R=100 \times B=100$ + Hessian $\widehat{\text{SE}}^*$)**:
    - `kp11 (0.05)`: `Percentile` runtuh ke `14,0%` (`[8,5%, 22,1%]`); `Basic` mencapai `95,0%` (`[88,8%, 97,9%]`) tetapi `100/100 (100%)` berbatas bawah $\le 0$ (non-fisik); `Wald` `85,0%` (`79%` non-fisik); `Studentized` `86,0%` (`87%` non-fisik); `Log-Basic` menjamin `0%` non-fisik (`2%` OOB) namun cakupannya hanya `49,0%` (`[39,4%, 58,7%]`).
    - `theta1 (0.045)`: `Studentized` (*bootstrap-t*) menaikkan cakupan dari `52,0%` (`Wald`), `46,0%` (`Percentile`), dan `47,0%` (`Basic`) menjadi **`72,0%`** (`[62,5%, 79,9%]`, median lebar naik $2\times$ menjadi `0,0462`), namun masih sensitif terhadap $\hat{\kappa}_1$ luar (`94,0%` bila $\hat{\kappa}_1 \le 0{,}2079$ vs `50,0%` bila $\hat{\kappa}_1 > 0{,}2079$) dan 2 outlier regularisasi Hessian (`mean_w = 19,12`).
- Output `outputs/day004d/results.json`, 4 grafik SVG `outputs/figures/day004d_*.svg`, serta dokumen penjelasan [day004d_quant_afns-interval-calibration.md](project-explanations/day004d_quant_afns-interval-calibration.md).

### Diperbarui
- `README.md` (tabel progres baris `4d` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 004c — pemeliharaan A] — 2026-11-22 · Diagnostik presisi pencarian garis `day004c`

### Diperbaiki
- `projects/go/day004c_afns_bias_correction`: menambahkan diagnostik penghentian pencarian garis per fit (`min_alpha_tested`, `last_step_norm`, `candidate_bit_identical`, `last_obj_change`, `obj_bit_identical`, `line_search_failure_cause` serta agregasi `step_underflow_plateau` vs `armijo_not_met` vs `invalid_objective_or_gradient`) di `optim.go`, `kalman.go`, `jackknife/study.go`, `bootstrap/study.go`, dan `cmd/run/main.go`, plus uji unit plateau vs gagal Armijo non-plateau di `optim_test.go`.
- Memverifikasi kesetaraan penuh terhadap `results.json` awal (`206,14s`, seluruh parameter, log-likelihood, sampel, dan status identik bit-demi-bit) serta membuktikan seluruh `5.000/5.000` refit bootstrap berhenti pada lantai presisi `float64` (`4.628` `step_underflow_plateau`, `372` `armijo_not_met` di $\alpha = 2^{-34}$ dengan perubahan $\sim 1\text{–}4\text{ ULP}$, `2.840` berobjektif identik bit-demi-bit, `0` tak-valid).

## [Hari 004c] — 2026-11-22 · AFNS — koreksi bias mean-reversion dan inferensi sampel terbatas dengan jackknife serta bootstrap parametrik (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ mandiri `projects/go/day004c_afns_bias_correction/` (`zero external modules`, `MainSeed = 2026112203`, `xoshiro256++` + `SplitMix64`, paralelisme goroutine deterministik per *chunk*, `FastLogLikelihood` *steady-state* Riccati, penulis JSON 17 digit & 4 grafik SVG):
  - **Jangkar Eksak & Verifikasi Teknis (`c1..c5 = tercapai`)**: Gaussian gabungan $T=5$ (3 tenor) cocok terhadap Kalman + *RTS smoother* (`ll_diff = 5,68e-14`, `fast_ll_diff = 5,68e-14`, `mean_diff = 5,55e-16`, `cov_diff = 6,39e-17`), skor FD vs Richardson ($\le 2,75\cdot 10^{-7}$), uji aljabar pembatalan $c/H$ jackknife (`1,39e-17`), dan determinisme bit-identik `1` vs `8` worker pada fixture cepat.
  - **Bagian 1 (Jackknife Temporal Berpasangan $H \in \{20, 80\}$ tahun mingguan, $R=50$, prefiks 80y bersama, tanpa *clamping*)**: pada `20y`, bias `kp11` turun $+0{,}3170 \to +0{,}0994$ (`RMSE` `0,4299 -> 0,3546`), `kp22` $+0{,}2130 \to +0{,}0463$ (`RMSE` `0,4046 -> 0,4165`), `kp33` $+0{,}1821 \to -0{,}0715$ (`RMSE` `0,3059 -> 0,3730`), namun memunculkan `40%` estimasi $\hat{\kappa}_1 \le 0$ (`32%` di `80y`), serta kegagalan batas paruh-sampel pada `20y half1 Rep 28` dan `80y half1/half2 Rep 30` (yang melonjakkan `RMSE` jackknife `theta1` di `80y` menjadi `6,28x`).
  - **Bagian 2 (Bootstrap Parametrik `20y` mingguan, $R=50 \times B=100 = 5{,}000$ refit penuh)**: `5.000/5.000` refit bernilai *finite* (`50/50` interval luar valid memenuhi ambang $\ge 95/100$), memangkas bias sekaligus RMSE seluruh ketiga $K^{\mathbb{P}}_{ii}$ (`kp11`: bias `+0,3170 -> +0,1014`, `RMSE 0,4299 -> 0,3172` dengan `38%` $\le 0$; `kp22`: `+0,2130 -> +0,0034`, `RMSE 0,4046 -> 0,3446`; `kp33`: `+0,1821 -> -0,0282`, `RMSE 0,3059 -> 0,2476`).
  - **Cakupan Interval 95% (Wald Hessian vs Bootstrap Persentil + Wilson 95% CI, denominator $R=50$)**: interval persentil bootstrap naif untuk `kp11` runtuh ke `2,0%` (`1/50`, Wilson `[0,35%, 10,5%]`) dibanding Wald `76,0%` karena bias ke atas ganda pada simulasi bootstrap, dan untuk `theta1` tetap rendah di `42,0%` (`21/50`, Wilson `[29,4%, 55,8%]`) dibanding Wald `46,0%`.
- Output `outputs/day004c/results.json` (`180,13s` pada 8 worker), 4 grafik SVG `outputs/figures/day004c_*.svg`, serta dokumen penjelasan [day004c_quant_afns-bias-correction.md](project-explanations/day004c_quant_afns-bias-correction.md).

### Diperbarui
- `README.md` (tabel progres baris `4c` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 004b — pemeliharaan A] — 2026-11-22 · Diagnostik penghentian optimizer & audit jumlah tes `day004b`

### Diperbaiki
- `projects/go/day004b_afns_inference`: menambahkan alasan penghentian BFGS eksplisit (`gradient_tolerance`, `boundary_stationary`, `line_search_failed`, `max_iterations`, `invalid_objective_or_gradient`) beserta diagnostik `stage1`, `stage2`, `best_stage`, `iterations_best_stage`, `iterations_cumulative`, dan `raw_grad_norm_opt` pada `optim.go`, `kalman.go`, `study.go`, dan `cmd/run/main.go`.
- Menghapus reklasifikasi otomatis `line_search_failed` sebagai `converged = true` dan membatasi pemilihan tahap-2 ke perbaikan ketat `res2.Fun < res1.Fun` (`est_loglik` tetap identik bit-demi-bit pada `max_ll_diff = 0`, sambil mengeliminasi pergeseran koordinat datar $\le 1{,}61\times 10^{-8}$ pada 71 replikasi seri serta membuktikan 4 kasus sebelumnya non-konvergen interior sebenarnya mencapai `max_iterations` 150 di tahap 2).
- Memperbarui [day004b_quant_afns-inference.md](project-explanations/day004b_quant_afns-inference.md) dan [projects/go/day004b_afns_inference/README.md](../projects/go/day004b_afns_inference/README.md) sesuai hasil audit tes terkini (`9` paket, `15` tes utama, `4` subtes, `54,76s`).

## [Hari 004b] — 2026-11-22 · AFNS — bias sampel terbatas, rentang kalender (long-span), dan frekuensi observasi (infill) (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ mandiri `projects/go/day004b_afns_inference/` (`zero external modules`, `MainSeed = 2026112202`, `xoshiro256++`, paralelisme deterministik per replikasi, optimisasi *steady-state* Riccati pada `FastLogLikelihood`, penulis JSON 17 digit & 5 grafik SVG):
  - **Bagian 1 (Jangkar Eksak, Verifikasi Skor & Determinisme)**: jangkar Gaussian gabungan penuh $T=5$ (3 tenor) cocok terhadap Kalman + *RTS smoother* (`ll_diff = 3,41e-13`, `fast_ll_diff = 3,55e-13`, `mean_diff = 2,83e-16`, `cov_diff = 6,39e-17`); skor FD vs Richardson sepanjang studi ($\le 4{,}34\cdot 10^{-7}$); determinisme bit-identik `1` vs `8` worker pada fixture cepat.
  - **Bagian 2a (Replikasi Baseline `1000w`, $R=50$)**: terpisah dari konfigurasi tepat $20$ tahun (`T=1040`), seluruh $11$ parameter konsisten dengan acuan `day004a` dalam $3\text{ MCSE}$ gabungan (`max |z| = 2,82` pada `kp33`, `10/11` dalam $2\text{ MCSE}$).
  - **Bagian 2b (Long-Span Mingguan $\{10, 20, 40, 80\}$ tahun, $R=50$, prefiks 80y bersama)**: bias $K^{\mathbb{P}}_{ii}$ meluruh konsisten $\propto 1/H$ (`kp11`: $+0{,}4657 \to +0{,}2241 \to +0{,}0983 \to +0{,}0451$, OLS slope $+4{,}84 \pm 0{,}94$; `kp22`: $+0{,}4254 \to +0{,}1998 \to +0{,}1007 \to +0{,}0731$, OLS $+4{,}12 \pm 0{,}84$; `kp33`: $+0{,}3590 \to +0{,}2026 \to +0{,}0981 \to +0{,}0452$, OLS $+3{,}56 \pm 0{,}88$), dengan produk $\text{bias} \times H \approx 3{,}6\text{–}5{,}8$.
  - **Bagian 2c (Infill 20 tahun $\Delta \in \{1/12, 1/52, 1/252\}$, $R=50$, memakai ulang `longspan_20y_weekly`)**: peningkatan frekuensi hingga harian ($T=5{,}040$) gagal menurunkan RMSE/bias *drift* (`daily/monthly RMSE`: `kp11 = 0,980`, `kp22 = 1,181`, `kp33 = 0,936`), sementara parameter penampang dan difusi menyusut tajam (`lambda = 0,012`, `sigma_eps = 0,003`, `sigma1 = 0,217`, `sigma2 = 0,100`, `sigma3 = 0,128`).
  - **Ringkasan Utuh & Diagnostik Hessian**: seluruh $50/50$ replikasi dipertahankan tanpa eksklusi kasus gagal konvergen atau batas aktif (`1/50` batas aktif di `1000w`, `10y`, `40y`, `80y`, `20y monthly`; `0/50` di `20y weekly` & `daily`), cakupan Wald + interval binomial Wilson 95% dilaporkan per parameter (`theta1` `34%–48%` di $10\text{–}40\text{y}$, naik ke `72%` di $80\text{y}$).
- Output `outputs/day004b/results.json` (`41,93s` total pada 8 worker), 5 grafik SVG `outputs/figures/day004b_*.svg`, serta dokumen `docs/project-explanations/day004b_quant_afns-inference.md`.

### Diperbarui
- `README.md` (tabel progres baris `4b` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 004a — pemeliharaan B] — 2026-11-22 · Kriteria penerimaan subruang derau PCA `day004a`

### Diperbaiki
- `projects/go/day004a_term_structure`: mengubah penghitungan penerimaan PCA pada `internal/pca/pca.go` dan `cmd/run/main.go` sehingga `c9c_pca_all10_sorted_eigenvalues_within_3se` tetap dicatat sebagai diagnostik efek statistik urutan Marchenko–Pastur pada 7 mode derau identik (`false`), sedangkan status agregat `c9_pca_eigenvalues_cosine_and_signs` (`true`, `"tercapai"`) bertumpu pada PC1–PC3 dalam 3 SE, rata-rata *noise bulk* PC4–PC10 dalam 3 SE, kosinus 3 eigenvektor $> 0{,}99$, dan pola tanda `[0, 1, 2]`. `all_passed` tetap `false` karena status MLE Bagian 3 masih `sebagian`. Memperbarui tes unit serta dokumentasi `day004a`.

## [Hari 004a] — 2026-11-22 · Struktur suku bunga: Nelson–Siegel, AFNS dan filter Kalman (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ murni `projects/go/day004a_term_structure/` (nol modul eksternal, paralelisme *goroutine* deterministik per *chunk*, `xoshiro256++` dengan `MainSeed = 20261122`, aljabar linear LU/Cholesky/Thomas/Jacobi/`MatExpSeries` buatan sendiri, kuadratur Gauss–Legendre, RK4 + ekstrapolasi Richardson, optimizer BFGS, penulis JSON 17 digit & 8 grafik SVG):
  - **Bagian 1 (Nelson–Siegel(–Svensson), Bootstrap, Interpolasi & Risiko)**: evaluasi stabil deret Taylor untuk $\tau$ kecil vs pembatalan naif (`g1RelErr = 5.91e-5` di $\tau=10^{-12}$), verifikasi $y(\tau) = \frac{1}{\tau}\int_0^\tau f(s)ds$ (`max_err = 6.94e-18`) dan limit $\tau \to 0, \infty$, *bootstrap* par semi-tahunan $0.5\text{–}30$ tahun (`max_err = 2.78e-16`), perbandingan stabilitas derau par $R=200$ (`0` *forward* negatif pada keduanya; `TV(f)` `0.300` LogLin vs `0.422` Spline, `yRMSE` `0.823` vs `0.688 bp`), metrik risiko obligasi 10y 4% (`ModDur relErr = 5.84e-11`, `Conv relErr = 1.12e-7`, jumlah KRD = $D_{\text{eff}}$ dengan galat `0`), PV01 swap 10y (`$811.01/1M`), dan lindung nilai netral dua KRD (`5y4% + 10y2%`, residual `< 5e-16`).
  - **Bagian 2 (AFNS, ODE Riccati, Dinamika $\mathbb{P}$ & Premi Jangka Waktu)**: penurunan tertutup $B(\tau)$ (pemuatan NS eksak) dan $-A(\tau)/\tau$ (CDR 2011; `-10.09 bps` di 10y, `-59.98 bps` di 30y) cocok vs RK4/Richardson (`1.16e-15`), transisi segitiga bawah $\Phi(\Delta)$ vs `MatExpSeries` (`1.11e-16`), $Q_\Delta$ Van Loan vs Gauss–Legendre (`4.74e-20`), uji MC $200\text{k}$ transisi 1y (`max_z = 2.68`), harga risiko afinitas $(\lambda_0, \lambda_1)$, dan premi jangka waktu 10y rata-rata jangka panjang tertutup `+173.04 bps` vs MC $100\text{k}$ `+172.21 ± 0.38 bps` (`z = 2.19`).
  - **Bagian 3 (Simulasi, Filter Kalman, RTS Smoother, MLE & Replikasi)**: jangkar Gaussian gabungan penuh $T=5$ (3 tenor) cocok terhadap Kalman & *RTS smoother* (`ll_diff = 2.84e-14`, `mean_diff = 7.27e-16`, `cov_diff = 1.24e-16`), MLE BFGS 11 parameter ($T=1000$, `est_ll = 58231.66 > true_ll = 58225.07`, skor vs Richardson `2.32e-10`), studi replikasi $R=50$ (`50/50` $\ell_{\text{MLE}} \ge \ell_{\text{true}}$; seluruh 5 parameter teridentifikasi kuat $\lambda, \sigma_{1..3}, \sigma_\varepsilon$ memenuhi `|bias| < 2 SE_rep` dan cakupan `90%–96%`, sementara $K^{\mathbb{P}}_{ii}$ menunjukkan bias Kendall $+0.21\text{/thn}$), dan uji kesalahan spesifikasi DNS (`loglik_drop = 767.51`).
  - **Bagian 4 (PCA Perubahan Yield, $T=5000, R=100$)**: tiga faktor pertama menjelaskan `82.90%` varians (`0, 1, 2` perubahan tanda), eigenvalue empiris PC1–PC3 cocok dalam `1.33 SE` dan kosinus eigenvektor `0.99995, 0.99923, 0.99748` ($> 0.99$), rata-rata *bulk* derau PC4–PC10 cocok $2\sigma_\varepsilon^2$ (`z = 0.28`), sementara mode terurut individual terpecah oleh statistik orde Marchenko–Pastur.
- Output `outputs/day004a/results.json`, 8 grafik SVG `outputs/figures/day004a_*.svg`, serta dokumen penjelasan `docs/project-explanations/day004a_quant_term-structure.md`.

### Diperbarui
- `README.md` (tabel progres baris `4a` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 004a — pemeliharaan A] — 2026-11-22 · Pembersihan bukti tes `day003e`/`day003f` & penyelarasan bendera C8 `day003e`

### Diperbaiki
- `day003e` & `day003f`: menjalankan `gofmt -l .`, `go vet ./...`, dan `go test -race -count=1 ./...` pada seluruh paket (`5/5` paket `day003e`, `8/8` paket `day003f`) dan mencatat keluaran nyata pada `docs/project-explanations/day003e_quant_cva-wrong-way.md` dan `day003f_quant_ruin-theory.md`.
- `day003e`: memisahkan bendera penerimaan C8 pada `cmd/run/main.go` menjadi `c8a_collateral_threshold_monotonicity: true`, `c8b_collateral_mpor_clean_schedule_convergence: true`, `c8c_collateral_mpor_raw_coupon_convergence: false`, `c8_collateral_threshold_and_mpor_convergence: false`, `c8_status: "sebagian"`, dan `all_passed: false` konsisten dengan dokumentasi; seluruh angka lain di `outputs/day003e/results.json` terbukti bit-identik.

## [Hari 003f] — 2026-11-21 · Teori ruin Cramér–Lundberg dan simulasi kejadian langka (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ murni `projects/go/day003f_ruin_theory/` (nol modul eksternal, paralelisme *goroutine* deterministik per *chunk*, `xoshiro256++`, sampler Gamma/Erlang/Pareto/Weibull/Lognormal/Poisson buatan sendiri, Brent, Golden Section, LU pivoting parsial, uji KS nol simulasi $O(n)$, penulis JSON 17 digit & 7 grafik SVG):
  - **Bagian 1 (Bentuk Tertutup, Momen & Defisit)**: akar Lundberg Brent & transformasi Laplace rasional LU untuk (A) Eksponensial (`R = 1/6`, `err = 4.16e-16`), (B) Erlang(2) (`R1 = 0.226765`, `err = 1.61e-15`), dan (C) Hipereksponensial-2 ($p_1 = p_2 = 0.5, \beta_1 = 2/3, \beta_2 = 2$, `CV = 1.2247 > 1`, `R1 = 0.130502`, `err = 5.83e-16`); verifikasi momen `200k` di $t=10$ (beserta audit jujur varians sampel Pareto karena $E[X^4] = \infty$) dan defisit Eksponensial(1) lolos uji KS ($p = 0.394, 0.985$).
  - **Bagian 2 (Rekursi Panjer / Beekman)**: pengapitan eksak `psi_low <= psi_exact <= psi_high` pada semua $h \in \{0.1, 0.05, 0.02, 0.01\}$ dan $u \in \{0..50\}$ ber-orde `0.997–0.999`, rasio Cramér–Lundberg $\to 1$, batas Lundberg berlaku, serta evaluasi ekor berat stabil hingga $u = 200$ (`Weibull ratio_vs_asym = 217.02`).
  - **Bagian 3 (MC Biasa & IS Eksponensial)**: MC `1.000.000` cocok `<= 1.47 SE` dan monoton dalam $T$; IS eksponensial menaksir peluang sampai `1.209e-10` ($u=100$) ber-SE relatif `0.045%–0.062%` (VRF waktu `4.82e6x–1.42e11x`); *tilting* salah ($0.7R, 1.3R$) tetap tak-bias bervariansi `3.25x–10.11x` lebih tinggi.
  - **Bagian 4 (Ekor Berat & Asmussen–Kroese)**: pembuktian divergensi $M_X(0.1)$ hingga $x_{\max}=500$ (`1.36e14` Pareto, `4.53e7` Weibull); estimator AK berbobot-ukuran (`200k`) cocok dengan batas Panjer $h=0.01$ sampai $u = 200$ (`relSE = 0.092%` Pareto, VRF waktu hingga `3.49e5x`) vs MC biasa & asimtotik subeksponensial.
  - **Bagian 5 (Reinsurans & Modal)**: kelebihan rugi $\min(X, M)$ ($\xi = 0.4$) dengan solusi tertutup (A), optimisasi retensi $M^*_A = 1.529$ & $M^*_C = 1.760$, serta penghematan modal awal $u^*(\varepsilon)$ `25.0%–33.1%` pada $\varepsilon \in \{10^{-2}, 10^{-3}, 10^{-4}\}$.
- Output `outputs/day003f/results.json`, 7 grafik SVG `outputs/figures/day003f_*.svg`, serta dokumen penjelasan `docs/project-explanations/day003f_quant_ruin-theory.md`.

### Diperbarui
- `README.md` (tabel progres baris `3f` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 003e] — 2026-11-21 · CVA swap suku bunga: Hull–White, intensitas CIR dan wrong-way risk (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ murni `projects/go/day003e_cva_wrong_way/` (nol dependensi eksternal, paralelisme *goroutine* deterministik per *chunk*, `xoshiro256++`, Gamma/Poisson eksak Ahrens–Dieter, kuadratur Golub–Welsch Gauss–Hermite & Gauss–Legendre, penulis JSON 17 digit & 7 grafik SVG):
  - **Bagian 1 (Hull–White 1-Faktor & Eksposur)**: simulasi bersama eksak $(x_t, I_t)$ dengan stabilisasi deret Taylor untuk $\Delta$ kecil; uji martingale (`max_gh_abs_err = 8.88e-16`, `mc_z <= 0.578`); *swaption* `2y × 8y` Jamshidian (`ITM=0.071464`, `ATM=0.026532`, `OTM=0.005768`, selisih vs kuadratur bebas-*kink* `<= 1.11e-16`, `mc_z <= 0.953`); profil eksposur kuartalan `200k` lintasan (`EPE = 0.017726` / `$17,726.10/1M`, puncak `EE = 0.031838` di $t=4$, `max_z = 2.010` di 41 tanggal).
  - **Bagian 2 (Intensitas CIR & Par CDS)**: verifikasi Feller ($2\kappa\theta = 0.025 \ge \xi^2 = 0.01$), kurva survival afin vs ODE Riccati RK4 (`max_err = 1.33e-15`), sampling eksak Poisson–Gamma vs Euler *full-truncation* (`Q(10) = 0.789127`), serta par CDS `1/3/5/10Y` (`126.68–141.75 bps`, `pv_prot_z <= 1.10`).
  - **Bagian 3 (CVA Independen & Sensitivitas)**: `CVA_anal = 0.00191515` (`$1,915.15/1M`), turunan analitik $\partial\text{CVA}/\partial\lambda_0$ vs FD pusat (`abs_err = 8.08e-14`), CS01 (`$1.60 / $5.69 / $7.28` per `$1M`), CVA MC waktu gagal eksplisit $\tau$ & Rao–Blackwellized (`z <= 0.57` di $m \in \{20, 40, 80\}$), dan konvergensi *paired* terhadap $m = 320$ (`9.67e-7 -> 3.81e-7 -> 1.57e-7`, *slope* `-1.31`).
  - **Bagian 4 (*Wrong-Way Risk*, *Netting* & Kolateral)**: kurva CRN $\text{CVA}(\rho)$ naik monoton `$1,080.34 -> $2,962.22/1M` ($\rho: -0.8 \to +0.8$); *netting* `5y payer + 10y receiver` hemat `36.88%` ($\rho=0$) dan identitas berlawanan arah nol eksak; kolateral MPOR $\delta=0.04$ lintas threshold $H \in \{0, 0.01, 0.05, \infty\}$ (`$112.77 -> $1,912.87/1M` di $\rho=0$) serta analisis MPOR $\delta \to 0$ (jadwal bersih $O(\sqrt{\delta})$ *slope* `0.454` vs lantai kupon mentah `$325/1M`).
- Output `outputs/day003e/results.json`, 7 grafik SVG `outputs/figures/day003e_*.svg`, serta dokumen penjelasan `docs/project-explanations/day003e_quant_cva-wrong-way.md`.

### Diperbarui
- `README.md` (tabel progres baris `3e` Go via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 003d] — 2026-11-21 · Risiko ekor: GARCH–EVT, VaR/ES dan backtesting (TypeScript, Bonus)

### Ditambahkan
- Project TypeScript murni `projects/typescript/day003d_tail_risk/` (`"strict": true`, nol dependensi runtime, paralelisme `worker_threads` deterministik, `xoshiro256++` 64-bit limb 32-bit, penulis JSON 17 digit & 7 grafik SVG):
  - **Bagian 1 (Distribusi eksak & EVT)**: Student-$t$ (`regularizedIncompleteBeta`, Newton+bisection, rumus tertutup ES baku vs kuadratur `<= 3.38e-14`), GPD (`<= 1.78e-15`), Pickands–Balkema–de Haan & Hill pada `N = 2.000.000` sampel inversi $t_{\text{std}}(5)$ ($\hat{\xi}_{\text{POT}}$ `0.1031 -> 0.1785` menuju `0.2`).
  - **Bagian 2 (GARCH-t)**: simulasi, rekonstruksi $\sigma_t^2$ (`2.71e-20`), skor analitik vs FD (`2.02e-6`), BFGS vs Nelder–Mead (`1.29e-14`), studi replikasi `R = 200` (`200/200` konvergen, `0` batas aktif, t-MLE unggul `1.53x–2.92x` MSE atas QMLE Gaussian), serta pembuktian aljabar $E[\varepsilon_t^4] = \infty$ ($D=-0.0116$) plus jangkar bermomen-4 berhingga ($t(12), D=0.1625$).
  - **Bagian 3 (Filtered EVT VaR/ES)**: evaluasi 4 metode di `R = 200` ($n=2500$) pada konvensi kerugian $L_t = -r_t$ (`Gaussian VaR 99.9%` bias `-32.41%`, `EVT RMSE 0.1024 < FHS 0.1304` dengan $z=6.21$, plus sensitivitas ambang $q \in \{0.85, 0.90, 0.95\}$).
  - **Bagian 4 (Backtesting)**: ukuran Kupiec (`100k` seri), daya klaster Markov Christoffersen (`20k`: `99.7%` / `81.2%`), distribusi nol simulasian sendiri Acerbi–Székely $Z_2$ (`30k`+`10k`), *rolling backtest* $40 \times 6000$ hari (`true_control rate99 = 1.009% ± 0.025%`, `gaussian_garch` ditolak `97.5%` Kupiec & `100%` $Z_2$), dan perbandingan deteksi perubahan rezim hari ke-3000 (`rolling` adapt $\hat{\alpha}_t \to 0.141$, pulih ke `1.042%` vs `fixed` `0.862%` ditolak `27.5%`).
- Output `outputs/day003d/results.json`, 7 grafik SVG `outputs/figures/day003d_*.svg`, serta dokumen penjelasan `docs/project-explanations/day003d_quant_tail-risk.md`.

### Diperbarui
- `README.md` (tabel progres baris `3d` TypeScript via `scripts/update_readme.py`), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Pemeliharaan] — 2026-11-21 · Rename sesi pertama harian (`1` → `1a`, `3` → `3a`) & sinkronisasi skrip

### Diubah
- Sesi pertama tiap hari kini secara konsisten berlabel huruf `a` (`days/day001a_quant_return-harian-kumulatif.ipynb`, `days/day003a_quant_volatilitas-lokal-dupire-ssvi.ipynb`, `outputs/figures/day001a_return-harian-kumulatif.png`, serta metadata `session: "a"`).
- `scripts/new_day.py`, `scripts/update_readme.py`, `scripts/run_all.py`, dan `tests/test_scripts.py` diselaraskan agar sesi pertama bernama `a` sambil tetap kompatibel membaca metadata lama, serta menyesuaikan struktur `projects/<bahasa>/`.

## [Hari 003c] — 2026-11-21 · Risiko sistemik jaringan: kliring Eisenberg–Noe, kaskade Watts & biaya kebangkrutan (Go, Bonus)

### Ditambahkan
- Project Go 1.22+ murni `projects/go/day003c_systemic_risk/` (hanya standard library, paralelisme goroutine deterministik per chunk, RNG `xoshiro256++`, solver LU buatan sendiri, pembangkit graf CSR Erdős–Rényi Batagelj–Brandes $O(n+m)$, penulis JSON 17 digit & SVG):
  - **C1 (Eisenberg–Noe)**: iterasi Picard, algoritma default fiktif, dan *brute-force* $2^n$ cocok hingga `4.05e-13` pada 210 jaringan acak kecil; anchor rantai & siklus 2 bank (`<= 1.11e-16`); benchmark $n=2.000$ bank (Poisson & Power-Law).
  - **C2 (Kaskade Watts)**: bukti invariansi urutan FIFO vs LIFO pada $n=100.000$; koreksi & verifikasi rumus pendekatan pohon Watts (2002) terhadap simulasi $R=200$ pada 9 nilai $c$ (`z <= 2.30`); batas fase $(c_1, c_2)$ via bisection untuk $\phi \in \{0.10, 0.18, 0.30\}$; distribusi ukuran kaskade 1 benih ($R=2.000$).
  - **C3 (Rogers–Veraart & Uji Stres)**: verifikasi $\alpha=\beta=1 \implies$ EN (`3.55e-15`) & rantai tertutup (`0.0`); uji stres Monte Carlo `20.000` skenario per $s \in \{0.2, 0.4, 0.6, 0.8\}$ (VaR 99%, ES 97.5%, amplifikasi RV/EN `~1.46x`); perbandingan ketahanan jaringan Poisson vs Pangkat-Hukum terhadap guncangan 20 bank acak vs top-20 hub.
- Output `outputs/day003c/results.json`, 5 grafik SVG `outputs/figures/day003c_*.svg`, serta 21 dokumen penjelasan proyek di `docs/project-explanations/` (`day001`–`day003c`).

### Diperbarui
- `README.md` (tabel progres baris `3c` Go dan relabel `2k` Ruby), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 003c — pemeliharaan A] — 2026-11-21 · Perbaikan day003, day003b, relabel 2k, tautan relatif

### Diperbaiki
- `day003` (notebook "3a"): diagnosis gamma SSVI di $S_0$. Selisih gamma grid (`0.0380`) vs beda hingga harga (`0.0234`) bukan bug: di sudut singular $(S_0, t\to0)$ permukaan $\sigma_{loc}$ gamma bergantung pada regularisasi `t_min` dan tidak konvergen; di $S_0 \ne 100$ keduanya cocok (~1e-4 relatif). `pde_delta_gamma` kini memakai grid TETAP untuk beda hingga (opsi `x_bounds` di `solve_backward_european_call`); fungsi `gamma_diagnosis_ssvi`, tes baru, sel notebook, dan `results.json` diperbarui. Bias Euler MC diukur menyusut dengan jumlah langkah (CRN, 50/100/200/400).
- `day003b` (Rust): uji mean AC dengan antitetik dinyatakan vakum; ditambahkan uji mean tanpa antitetik (z O(1)); `sell-all-at-once` dicatat sebagai uji deterministik; dua slope limit kontinu (seluruh rentang dan ekor) plus dekomposisi galat yang membuktikan penjelasan pembatalan; kondisi kernel osilasi dan scan kelas "PD tetapi $x^*$ negatif".
- ID ganda `2g`: Ruby put Amerika dipindah menjadi `2k` (`projects/ruby/day002k_american_put`, `outputs/day002k/`); `outputs/day002g/results.json` tetap hasil AAD C++. Angka kunci Ruby identik dengan riwayat git (hanya waktu yang berbeda).
- Dokumen `day003b`: 14 tautan absolut mesin lokal diganti tautan relatif.

## [Hari 003b] — 2026-11-21 · Eksekusi optimal: Almgren–Chriss, dampak transien, likuiditas stokastik (Rust, Bonus)

### Ditambahkan
- Project Rust 2021 murni `projects/rust/day003b_optimal_execution/` (`#![forbid(unsafe_code)]`, nol crate eksternal runtime): Cholesky, Thomas tridiagonal $O(N)$, Jacobi eigenvalue simetris, RNG `xoshiro256++` + Box–Muller, paralelisme `std::thread` deterministik (chunk tetap), penulis JSON (17 digit) & SVG buatan sendiri.
- Bagian 1–3: Almgren–Chriss (QP KKT vs `sinh` diskrit, limit kontinu, frontier 40 titik, MC `400.000` jalur), Obizhaeva–Wang (solusi umum, anchor eksak, benchmark $N=2000$ $O(N)$ vs $O(N^3)$, power-law, osilasi PD dengan manipulasi harga, non-PD, gabungan), dan DP likuiditas stokastik 2-rezim (rekursi kuadratik, brute-force $N=5$, simulasi `400.000` jalur, nilai adaptivitas `14.09%`, sensitivitas persistensi).
- Output `outputs/day003b/results.json`, 7 grafik SVG `outputs/figures/day003b_*.svg`, serta dokumen penjelasan `docs/project-explanations/day003b_quant_optimal-execution.md`.

### Diperbarui
- `README.md` (tabel progres baris `3b` Rust), `docs/CHANGELOG.md`, dan `docs/project-explanations/README.md`.

## [Hari 003a] — 2026-11-21 · Volatilitas lokal Dupire dari permukaan SSVI (Python, Bonus)

### Ditambahkan
- Paket `src/qb100/day003a_local_vol/` (`bs.py`, `ssvi.py`, `dupire.py`, `grids.py`, `pde_backward.py`, `pde_forward.py`, `mc.py`, `barrier.py`, `greeks.py`, `calibration.py`, `plotting.py`, `runner.py`): permukaan SSVI bebas-arbitrase beserta turunan analitik & uji kontrafaktual, volatilitas lokal Dupire analitik + verifikasi beda hingga orde 4 & anchor $\eta \to 0$, dua solver PDE independen (backward di $x = \ln S$ dan forward Dupire di $K$, grid tak-seragam `sinh`, Crank–Nicolson + 4 langkah Rannacher), Monte Carlo log-Euler antitetik (`N = 400.000`, `200` langkah/tahun, PCG64) & uji *weak order*, *down-and-out call* (Reiner–Rubinstein, PDE absorbing, MC jembatan Brownian), Greeks (delta, gamma, vega-bucket vs pergeseran paralel), dan kalibrasi balik Levenberg–Marquardt buatan sendiri ($R = 50$ replikasi).
- Notebook `days/day003_quant_volatilitas-lokal-dupire-ssvi.ipynb`, unit test `tests/test_day003a_*.py`, `outputs/day003a/results.json` (presisi 17 digit), dan 7 grafik di `outputs/figures/day003a_*.png`.

### Diperbarui
- `README.md` (tabel progres baris `3`) dan `docs/CHANGELOG.md`.

## [Hari 002j] — 2026-11-20 · Longevity risk: Lee–Carter, anuitas stokastik & modal solvabilitas (Ruby, Bonus)

### Ditambahkan
- Project Ruby 3.3+ murni `projects/ruby/day002j_longevity_risk/` (tanpa gem runtime, YJIT aktif, paralelisme `Process.fork` deterministik): SVD Jacobi satu sisi, estimasi Lee–Carter (tanpa & dengan noise $R=200$), fit Gompertz–Makeham (Nelder–Mead), proyeksi stokastik $k_t$ ($N=20.000$) & bootstrap parametrik, anuitas kohor usia 65 (rekursi eksak vs simulasi 200.000 umur), Hukum Total Variansi portofolio, kurva $\mathrm{CV}(N)$, modal VaR 99.5% / ES 99% (simulasi penuh & aproksimasi normal bersyarat), serta stres deterministik 0.5%.
- Unit test `minitest` (`test/test_longevity.rb`) dan `rubocop` lulus bersih.
- Output `outputs/day002j/results.json` dan 6 grafik SVG (`day002j_mortality-surface.svg`, `day002j_life-expectancy-fan.svg`, `day002j_k-moments-qq.svg`, `day002j_cv-vs-n.svg`, `day002j_liability-hist-var.svg`, `day002j_capital-vs-n.svg`) di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2j` Ruby) dan `docs/CHANGELOG.md`.

## [Hari 002i] — 2026-11-20 · Model short-rate Hull–White 1-faktor: pohon trinomial, Jamshidian, kalibrasi & Bermudan LSM (C++, Bonus)

### Ditambahkan
- Project C++20 `projects/cpp/day002i_hull_white/` (hanya standard library): pohon trinomial Hull–White terkalibrasi ke kurva Nelson–Siegel via induksi maju Arrow–Debreu ($m \in \{12..192\}$, galat kurva $< 1.5 \times 10^{-15}$), pricing analitik & pohon ZBP/ZBC dan cap 5y, dekomposisi Jamshidian (Brent) vs pohon vs Monte Carlo antitetik OU + *control variate* martingale (`N = 200.000`, `dt = 1/52`), kalibrasi Nelder–Mead bebas-noise & 100 replikasi bernoise paralel (`std::thread`), swaption Bermudan 10y (pohon + Richardson vs *two-pass* LSM `N = 100.000`), serta DV01 pergeseran paralel $\pm 1\text{ bp}$.
- Unit test `doctest` (`tests/test_hw.cpp`) lulus di build Release dan UBSan (`-fsanitize=undefined`).
- Output `outputs/day002i/results.json` dan grafik `day002i_tree-convergence.svg`, `day002i_curve-and-swaptions.svg`, `day002i_calib-scatter.svg`, `day002i_bermudan-bounds-dv01.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2i` C++) dan `docs/CHANGELOG.md`.

## [Hari 002h] — 2026-11-20 · Hawkes + limit order book: simulasi eksak, MLE, mesin pencocokan & fakta stilistik (C++, Bonus)

### Ditambahkan
- Project C++20 `projects/cpp/day002h_hawkes_orderbook/` (hanya standard library): proses titik Hawkes (*Ogata thinning*, eksak *Dassios–Zhao*, uji penskalaan ulang waktu KS & Ljung–Box, MLE $O(n)$ + BFGS paralel `std::thread` $R=200$), mesin *limit order book* berkinerja tinggi (*pool allocator*, antrean FIFO intrusif, bitmap 2-level 64-bit, `0` alokasi *heap* setelah *warmup*, `23.9M ops/detik`) vs mesin referensi $O(n)$ pada 1.000.000 order, dan simulasi pasar gabungan ($T=20.000$) dengan fakta stilistik lengkap.
- Unit test `doctest` (`tests/test_hawkes_ob.cpp`) lulus di build Release dan UBSan (`-fsanitize=undefined`).
- Output `outputs/day002h/results.json` dan grafik `day002h_hawkes-fano-qq.svg`, `day002h_hawkes-mle-scatter.svg`, `day002h_ob-latency.svg`, `day002h_market-stylized-facts.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2h` C++) dan `docs/CHANGELOG.md`.

## [Hari 002g] — 2026-11-20 · AAD untuk Greeks Monte Carlo: tape mode-balik, pathwise & keranjang 10 aset (C++, Bonus)

### Ditambahkan
- Project C++20 `projects/cpp/day002g_aad_greeks/` tanpa library AD/numerik eksternal: *tape* mode-balik (`mark`/`rewind_to_mark`), `Dual<T>` mode maju, *forward-over-reverse* Gamma (`Dual<Number>`), faktorisasi Cholesky terdiferensiasi, dan mesin Monte Carlo paralel deterministik `std::thread` (`200.000` path) untuk GBM 1-langkah, Opsi Asia 252-langkah, *Digital call* (penghalusan $\Phi((S_T-K)/\varepsilon)$), keranjang aritmetik 10 aset (67 sensitivitas vs beda hingga CRN), dan keranjang geometrik 10 aset (vs gradien rumus tertutup).
- Unit test `doctest` (`tests/test_aad.cpp`) lulus di build Release dan UBSan (`-fsanitize=undefined`).
- Output `outputs/day002g/results.json` dan grafik `day002g_aad-vs-analytic.svg`, `day002g_digital-bias-variance.svg`, `day002g_basket-speedup.svg`, `day002g_tape-growth.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2g` C++) dan `docs/CHANGELOG.md`.

## [Hari 002f] — 2026-11-20 · Inferensi Bayesian state-space: Kalman, particle filter & PMMH (Rust, Bonus)

### Ditambahkan
- Crate Rust `projects/rust/day002f_state_space_pmmh/` (Kalman filter eksak + verifikasi *brute force* $T=10$ via Cholesky, bootstrap particle filter tanpa alokasi di loop panas + resampling sistematik, posterior grid 2D trapesium, Haario adaptive MH & PMMH paralel `rayon` deterministik, split $\hat{R}$, ESS *initial positive sequence*, MCSE, dan Model B Stochastic Volatility).
- 5 unit test `cargo test` (Kalman vs *brute force*, resampling sistematik tak bias, ketidakbiasan PF & normalisasi grid, validasi ESS AR(1) & $\hat{R}$ i.i.d., reproduktibilitas & invariansi jumlah thread).
- Output `outputs/day002f/results.json` dan grafik `day002f_pf-unbiased-and-var.png`, `day002f_model_a_posteriors.png`, `day002f_model_b_sv.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2f`) dan `docs/CHANGELOG.md`.

## [Hari 002e] — 2026-11-20 · Model Heston: pricing Fourier, kalibrasi & Monte Carlo QE (Bonus)

### Ditambahkan
- Paket `src/qb100/heston/` (`cf.py`, `pricing.py`, `calibration.py`, `simulation.py`): fungsi karakteristik *Little Heston Trap* (Albrecher et al. 2007), *pricing* Fourier Gil-Pelaez & COS (Fang–Oosterlee 2008), *implied volatility* (Newton + Brent), kalibrasi *least-squares* dengan *multi-start* & analisis identifiabilitas Jacobian ($100$ replikasi noisy), serta Monte Carlo QE (Andersen 2008) vs Euler *full-truncation* pada kasus pelanggaran kondisi Feller, lengkap dengan *control variate* $S_T$ dan Greeks CRN.
- Unit test `tests/test_heston_pricing.py` & `tests/test_heston_calib_mc.py` (6 test).
- Notebook `days/day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb` beserta grafik `day002e_smile-and-fourier.png`, `day002e_calibration-noisy.png`, dan `day002e_qe-vs-euler.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2e`) dan `docs/CHANGELOG.md`.

## [Hari 002d] — 2026-11-20 · Risiko kredit Vasicek: distribusi kerugian eksak, LHP & importance sampling (Bonus)

### Ditambahkan
- Modul `src/qb100/credit.py` (kopula Gaussian satu faktor Vasicek, DP konvolusi PMF bersyarat + Gauss–Hermite probabilis 200 node, verifikasi varians bivariat $\Phi_2$ dan `scipy.integrate.quad`, VaR/ES diskret, pendekatan LHP, dan 3 estimator Monte Carlo: Naif, Conditional MC via `CubicSpline` 2001-node DP, serta Importance Sampling pada $Z \sim \mathcal{N}(-\theta, 1)$).
- Unit test `tests/test_credit.py` (4 test: PMF/mean/varians bivariat, kesesuaian `quad`, konvergensi LHP pada $N=5000$, dan ketiga estimator tidak bias/reprodusibel/VRF).
- Notebook `days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb` beserta grafik `day002d_pmf-and-tail.png`, `day002d_lhp-vs-exact.png`, `day002d_is-variance-and-vrf.png`, dan `day002d_convergence-loglog.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2d`) dan `docs/CHANGELOG.md`.

## [Hari 002c] — 2026-11-20 · HMM dua rezim: forward-backward, Baum–Welch & Viterbi (Bonus)

### Ditambahkan
- Modul `src/qb100/hmm.py` dengan NumPy murni (simulator HMM Gaussian, forward naif vs forward–backward ber-skala Rabiner, dekode Viterbi ruang log, enumerasi *brute force* untuk verifikasi, Baum–Welch EM *batched* `(S, K)` dengan inisialisasi kuantil $|r|$ deterministik & pengurutan state menurut $\sigma$, seleksi BIC, diagnostik Viterbi, dan diagram reliabilitas).
- Unit test `tests/test_hmm.py` (3 test: verifikasi *brute force* $T=10$ multi-seed, *overflow* versi naif vs stabilitas versi ber-skala $T=2000$, dan monotonisitas/reprodusibilitas EM).
- Notebook `days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb` beserta grafik `day002c_viterbi-and-em-ll.png`, `day002c_param-recovery.png`, dan `day002c_calibration-and-bic.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2c`) dan `docs/CHANGELOG.md`.

## [Hari 002b] — 2026-11-20 · Estimasi kovarians & portofolio varians minimum: Wishart, Marchenko–Pastur, Ledoit–Wolf (Bonus)

### Ditambahkan
- Modul `src/qb100/covariance.py` (model faktor sintetis, kovarians sampel, shrinkage Ledoit–Wolf target identitas buatan sendiri, GMV bentuk tertutup & `SLSQP`, frontier *long-only*, spektrum Marchenko–Pastur, eksperimen Monte Carlo batched & sapuan rasio $q = N/T$).
- Unit test `tests/test_covariance.py` (5 test: SPD model faktor, GMV tertutup vs numerik, sifat & asimtotik Ledoit–Wolf, rasio Wishart & reprodusibilitas, batas eigenvalue Marchenko–Pastur).
- Notebook `days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb` beserta grafik `day002b_wishart-lw.png`, `day002b_marchenko-pastur.png`, `day002b_q-sweep.png`, dan `day002b_long-only-frontier.png` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `2b`) dan `docs/CHANGELOG.md`.

## [Hari 001i] — 2026-11-19 · Opsi Asia (Ruby): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Project Ruby murni `projects/ruby/day001i_asian_options/` (Ruby 3.3+ / YJIT, nol dependensi runtime, `Random.new` + Box–Muller polar buatan sendiri, `Math.erfc`, Black-Scholes, Kemna–Vorst diskret, paralelisasi deterministik `Process.fork` + `IO.pipe` + `Marshal`, buffer pra-alokasi, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator SVG murni Ruby) beserta 7 unit test Minitest (40 assertions) dan RuboCop bersih (0 offense).
- Output hasil `outputs/day001i/results.json` dan grafik `day001i_asian-convergence-loglog.svg`, `day001i_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `scripts/update_readme.py` (menambah `ruby` ke `LANG_DIRS` beserta test di `tests/test_scripts.py`), `README.md` (tabel progres baris `1i`), dan `.gitignore` (`.bundle/`, `projects/ruby/**/vendor/bundle/`).

## [Hari 001h] — 2026-11-19 · Opsi Asia (C): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Project C17 `projects/c/day001h_asian_options/` (CMake, `xoshiro256**` + Polar Marsaglia, `erfc`, Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `pthread`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator SVG murni `fprintf`) beserta 7 unit test.
- Output hasil `outputs/day001h/results.json` dan grafik `day001h_asian-convergence-loglog.svg`, `day001h_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1h` otomatis terbaca tanpa mengubah `update_readme.py`) dan `.gitignore` (`projects/c/**/build*/`).

## [Hari 001g] — 2026-11-19 · Opsi Asia (C++): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Project C++20 `projects/cpp/day001g_asian_options/` (CMake, `Xoshiro256**` + Polar Marsaglia, `std::erfc`, Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `std::thread`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator SVG murni) beserta 6 unit test `doctest`.
- Output hasil `outputs/day001g/results.json` dan grafik `day001g_asian-convergence-loglog.svg`, `day001g_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1g` otomatis terbaca tanpa mengubah `update_readme.py`) dan `.gitignore` (`projects/cpp/**/build*/`, `cmake-build-*/`).

## [Hari 001f] — 2026-11-19 · Opsi Asia (TypeScript): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Package TypeScript `projects/typescript/day001f_asian_options/` (PRNG `Xoshiro128++` + Polar Marsaglia, `erf` Cody, Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `worker_threads`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope log-log, generator grafik SVG murni) beserta 6 unit test Vitest.
- Output hasil `outputs/day001f/results.json` dan grafik `day001f_asian-convergence-loglog.svg`, `day001f_asian-estimators-ci.svg` di `outputs/figures/`.

### Diperbarui
- `scripts/update_readme.py` (`typescript` di `LANG_DIRS`), `README.md` (tabel progres baris `1f`), dan `.gitignore` (`node_modules/`, `projects/typescript/**/dist/`).

## [Hari 001e] — 2026-11-19 · Opsi Asia (Go): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Module Go `projects/go/day001e_asian_options/` (Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel deterministik `sync.WaitGroup` + `math/rand/v2` PCG, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope konvergensi log-log) beserta 5 unit test.
- Output hasil `outputs/day001e/results.json` dan grafik `day001e_asian-convergence-loglog.png`, `day001e_asian-estimators-ci.png` di `outputs/figures/`.

### Diperbarui
- Generalisasi `scripts/update_readme.py` untuk membaca `*/meta.json` lintas bahasa (`rust`, `go`, `ts`, `cpp`, `c`) dan menambah kolom `Bahasa` di tabel progres `README.md`.
- `.gitignore` (`projects/go/**/bin/`, `*.test`).

## [Hari 001d] — 2026-11-19 · Opsi Asia (Rust): Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Crate Rust `projects/rust/day001d_asian_options/` (Black-Scholes, Kemna–Vorst diskret, simulasi GBM paralel `rayon` + `ChaCha8Rng`, estimator Naive/Antithetic/Control Variate, Delta *pathwise* & *bump-and-revalue* CRN, slope konvergensi log-log) beserta 5 unit test.
- Output hasil `outputs/day001d/results.json` dan grafik `day001d_asian-convergence-loglog.png`, `day001d_asian-estimators-ci.png` di `outputs/figures/`.
- Dukungan metadata `projects/rust/*/meta.json` pada `scripts/update_readme.py` beserta unit test.

### Diperbarui
- `README.md` (tabel progres baris `1d` dengan `Topik No = bonus`) dan `.gitignore` (`target/`).

## [Hari 001c] — 2026-11-19 · Opsi Asia: Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Modul `src/qb100/options.py` (Black-Scholes call/put, geometric Asian call diskret Kemna–Vorst, simulasi GBM tervektorisasi, estimator Naive/Antithetic/Control Variate, serta Delta *pathwise* & *bump-and-revalue* CRN) beserta unit test di `tests/test_options.py`.
- Dukungan `--bonus` (`topic_no: null` / tampilan `"bonus"` di tabel progres) pada `scripts/new_day.py`, `scripts/update_readme.py`, dan `tests/test_scripts.py`.
- Notebook off-roadmap `days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb` beserta grafik `day001c_asian-convergence-loglog.png` dan `day001c_asian-estimators-ci.png`.

### Diperbarui
- `README.md` (tabel progres baris `1c` dengan `Topik No = bonus`).

## [Hari 001b] — 2026-11-19 · Estimasi π dengan Monte Carlo

### Ditambahkan
- Helper `running_mean` dan `running_se` di `src/qb100/simulation.py` beserta unit test di `tests/test_simulation.py`.
- Notebook `days/day001b_sim_estimasi-dengan-monte-carlo.ipynb` (topik 002) beserta 4 grafik (`day001b_pi-scatter.png`, `day001b_pi-convergence.png`, `day001b_pi-error-scaling.png`, `day001b_pi-histogram.png`) di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1b`) dan `docs/ROADMAP.md` (centang topik 002).

## [Hari 001] — 2026-11-19 · Return harian & kumulatif

### Ditambahkan
- Aturan dan dukungan multi-sesi per hari (`day001`, `day001b`, ...) berbasis tanggal WITA (`Asia/Makassar`) pada `scripts/new_day.py`, `scripts/update_readme.py`, `scripts/run_all.py`, dan `templates/day_template.ipynb`.
- Fungsi `cumulative_returns` di `src/qb100/finance.py` beserta unit test di `tests/test_finance.py`.
- Notebook `days/day001_quant_return-harian-kumulatif.ipynb` (topik 001) beserta grafik `outputs/figures/day001_return-harian-kumulatif.png`.

### Diperbarui
- `README.md` (tabel progres dengan kolom Hari & Topik No, bagian scripts, dan konvensi penamaan).
- `docs/ROADMAP.md` (centang topik 001).

## [Hari 000] — 2026-11-19 · Setup baseline

### Ditambahkan
- Struktur repo: `days/`, `templates/`, `src/qb100/`, `scripts/`, `tests/`, `data/`, `outputs/figures/`.
- Package `qb100`: `plotting` (style seragam, `save_figure`), `data` (yfinance + cache CSV),
  `finance` (returns, log returns, volatilitas, Sharpe, drawdown), `simulation` (`get_rng`,
  random walk, Monte Carlo, ringkasan MC).
- Template notebook harian `templates/day_template.ipynb` dengan metadata `qb100`.
- Scripts `new_day.py`, `update_readme.py`, `run_all.py`.
- Pytest untuk `finance`, `simulation`, `data`/`plotting`, dan scripts.
- `docs/ROADMAP.md` berisi 100 topik (50 quant, 50 sim) dalam 4 fase.
- `README.md`, `LICENSE` (MIT), `requirements.txt`, `pyproject.toml`, `.env.example`.
