# Hari 004g — Market Making Optimal: Guéant–Lehalle–Tapia, HJB, dan Salah Spesifikasi (Rust)

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day004g`
- **Tanggal**: `2026-11-22`
- **Topik**: Market making optimal: Guéant–Lehalle–Tapia (GLT), HJB, dan studi salah spesifikasi
- **Jenis**: `quant`
- **Bahasa**: Rust (*pure standard library*, `edition = "2021"`, nol *crate* eksternal)
- **Lokasi Proyek**: [projects/rust/day004g_market_making](../../projects/rust/day004g_market_making/README.md)

Proyek ini mengimplementasikan kerangka **Guéant–Lehalle–Tapia (2013)** untuk *market making* optimal dengan batas inventori keras $q \in [-Q, Q]$, memverifikasi reduksi eksponensial terhadap *solver* Hamilton–Jacobi–Bellman (HJB) nonlinier langsung, memvalidasi simulasi Monte Carlo *Common Random Numbers* ($R = 200.000$ lintasan) terhadap fungsi nilai eksak dan persamaan Forward Kolmogorov, membandingkan empat strategi kuotasi, serta mengukur dampak salah spesifikasi (*adverse selection*, volatilitas $\sigma$, intensitas $A/k$, dan *drift* $\mu$).

## 2. Masalah dan Tujuan

1. **Verifikasi tanda dan keakuratan solver eksak**: Memastikan reduksi eksponensial $v_q(t) = \exp(k\,\theta(t, q))$ menghasilkan solusi yang konsisten antara tiga solver independen (dekomposisi eigen Jacobi siklik, eksponensial matriks *scaling-and-squaring*, dan ODE RK4+Richardson) serta terhadap solver beda hingga (FD) nonlinier langsung pada $\theta(t, q)$.
2. **Validasi dinamika stokastik**: Menguji apakah ekstrapolasi $\Delta t \to 0$ dari *certainty equivalent* Monte Carlo ($R = 200.000$) cocok dengan $\theta(0, 0)$ eksak dan apakah $\mathbb{E}[\mathrm{eksekusi}]$ serta $\operatorname{Var}(q_T)$ cocok dengan ODE Forward Kolmogorov.
3. **Evaluasi strategi & ketahanan salah spesifikasi**: Mengukur keunggulan GLT atas Avellaneda–Stoikov (AS), kuotasi simetris tetap, dan kuotasi terlalu rapat (*tight fixed*), serta mengkuantifikasi degradasi performa akibat *adverse selection* ($\varepsilon$), salah $\sigma$, salah $A/k$, dan *drift* tak-nol $\mu$.

## 3. Apa yang Dikerjakan Proyek

Runner biner `run` mengeksekusi empat bagian eksperimen secara deterministik (`main_seed = 20261122`, `8` *worker*) dan menulis:
- Laporan hasil lengkap berpresisi 17 digit ke [outputs/day004g/results.json](../../outputs/day004g/results.json).
- Tujuh grafik SVG diagnostik dan perbandingan ke [outputs/figures/](../../outputs/figures/).

## 4. Uraian Kode

Seluruh modul berada di [projects/rust/day004g_market_making/src/](../../projects/rust/day004g_market_making/src/lib.rs):
- [lib.rs](../../projects/rust/day004g_market_making/src/lib.rs): Konstanta skenario ($T=1{,}0$, $\sigma=0{,}3$, $\gamma=0{,}1$, $k=1{,}5$, $A=140$, $Q=20$), tipe hasil, dan re-export modul.
- [rng.rs](../../projects/rust/day004g_market_making/src/rng.rs) & [parallel.rs](../../projects/rust/day004g_market_making/src/parallel.rs): PRNG `xoshiro256++` di-seed deterministik per indeks lintasan dengan `SplitMix64` dan *worker pool* `std::thread::scope`.
- [linalg.rs](../../projects/rust/day004g_market_making/src/linalg.rs) & [ode.rs](../../projects/rust/day004g_market_making/src/ode.rs): Eigensolver Jacobi siklik, eksponensial matriks *scaling-and-squaring*, dan integrator RK4+Richardson.
- [glt.rs](../../projects/rust/day004g_market_making/src/glt.rs): Solver eksak sistem-$v$, kuotasi optimal $(\delta^b, \delta^a)$, celah spektral $\Lambda_1 - \Lambda_0$, batas stasioner (*ground state*), hampiran *closed-form*, serta sensitivitas $Q \in \{10, 20, 40\}$.
- [hjb_fd.rs](../../projects/rust/day004g_market_making/src/hjb_fd.rs): Solver Euler eksplisit mundur langsung pada persamaan HJB nonlinier bentuk $\theta(t, q)$.
- [kolmogorov.rs](../../projects/rust/day004g_market_making/src/kolmogorov.rs): Integrator ODE Forward Kolmogorov untuk probabilitas keadaan $P_q(t)$, $\mathbb{E}[\mathrm{eksekusi}]$, dan $\operatorname{Var}(q_T)$.
- [strategies.rs](../../projects/rust/day004g_market_making/src/strategies.rs), [sim.rs](../../projects/rust/day004g_market_making/src/sim.rs), & [misspec.rs](../../projects/rust/day004g_market_making/src/misspec.rs): Definisi strategi (`glt_optimal`, `avellaneda_stoikov`, `symmetric_fixed`, `tight_fixed_0.1`, `const_shifted_glt`, `exact_as_hjb`), simulator MC CRN, dan eksperimen salah spesifikasi.
- [jsonout.rs](../../projects/rust/day004g_market_making/src/jsonout.rs), [svg.rs](../../projects/rust/day004g_market_making/src/svg.rs), & [bin/run.rs](../../projects/rust/day004g_market_making/src/bin/run.rs): Serialisasi JSON, pembuatan 7 grafik SVG, dan orkestrasi eksperimen.

## 5. Metode dan Alur Kerja

1. **Alur**: Parameter skenario $\to$ tiga solver eksak sistem-$v$ + solver HJB FD nonlinier $\to$ analisis kuotasi, celah spektral, dan $Q \in \{10, 20, 40\}$ $\to$ Forward Kolmogorov ODE & simulator MC CRN ($R = 200.000$) $\to$ perbandingan 4 strategi & 4 studi salah spesifikasi $\to$ ekspor JSON dan 7 SVG.
2. **Persamaan Kunci**:
   - Transformasi eksponensial $v_q(t) = \exp(k\,\theta(t, q))$ dengan $\alpha = \frac{1}{2}k\gamma\sigma^2 = 0{,}00675$ dan $\eta = A(1 + \gamma/k)^{-(1+k/\gamma)} = 49{,}850378$ mengubah HJB menjadi $\frac{dv}{dt} = M v$, $v(T) = \mathbf{1}$.
   - Kuotasi optimal: $\delta^{b*}(t,q) = \frac{1}{k}\ln(v_q/v_{q+1}) + \frac{1}{\gamma}\ln(1+\gamma/k)$ dan $\delta^{a*}(t,q) = \frac{1}{k}\ln(v_q/v_{q-1}) + \frac{1}{\gamma}\ln(1+\gamma/k)$.

## 6. Output yang Ditargetkan

1. `c1_three_solvers_within_1e10`: Kesepakatan tiga solver eksak sistem-$v$ dengan galat relatif maksimum $< 10^{-10}$.
2. `c2_hjb_nonlinear_order1`: Konvergensi orde $\approx 1$ (`0,85..1,15`) untuk solver HJB FD eksplisit terhadap $\Delta t$.
3. `c3_quote_symmetry_and_monotonicity`: Simetri sempurna ($\delta^b(t,q) = \delta^a(t,-q)$) dan monotonisitas kuotasi terhadap $q$.
4. `c4_mc_utility_matches_exact`: Kesesuaian $\mathrm{CE}_{\mathrm{Rich}}$ (orde empiris $p$) dan $\mathrm{CE}_{\mathrm{cont}}$ (simulator kontinu) terhadap $\theta(0, 0)$ eksak dalam $|z| < 3$.
5. `c5_kolmogorov_matches_mc`: Kesesuaian $\mathbb{E}[\mathrm{eksekusi}]$ dan $\operatorname{Var}(q_T)$ simulator kontinu bebas-$\Delta t$ terhadap Forward Kolmogorov ODE dalam `3 SE` (disertai audit konvergensi $\operatorname{Var}(q_T)$ diskret pada $\Delta t \in \{10^{-3}; 5\times 10^{-4}; 2{,}5\times 10^{-4}; 1{,}25\times 10^{-4}\}$).
6. `c6_glt_dominates_benchmarks_ce`: Strategi `glt_optimal` mengungguli `avellaneda_stoikov`, `symmetric_fixed`, dan `tight_fixed_0.1` pada $\mathrm{CE}$ dan kekayaan berpasangan.
7. `c7_as_adjustment_improves_ce`: Reduksi $\varepsilon = 0$ antara `exact_as_hjb` dan GLT terverifikasi, dan pada $\varepsilon = 0{,}1$ strategi `exact_as_hjb` mengungguli `unadjusted_glt`.
8. `c8_wrong_sigma_ce_degradation`: Degradasi CE akibat salah spesifikasi $\sigma \in \{0{,}5\sigma, \sigma, 1{,}5\sigma\}$ tidak melebihi $0{,}1$ dan tidak melampaui $\mathrm{CE}(\sigma)$ secara signifikan (`< 3 SE`).
9. `c9_runtime_under_300s`: Total waktu runner penuh berada di bawah anggaran `300s`.

## 7. Status Pencapaian

Seluruh runner penuh selesai dalam **`53,53s`** pada `8` worker ([results.json](../../outputs/day004g/results.json)):

| Tujuan / Kriteria | Status | Bukti Terukur |
| --- | --- | --- |
| `c1_three_solvers_within_1e10` | `tercapai` | $\theta(0,0) = 66{,}261588$; galat rel maks Jacobi vs `expm` `4,97e-14`, vs RK4+Rich `2,80e-11 < 1e-10` |
| `c2_hjb_nonlinear_order1` | `tercapai` | Galat `1,48e-4` $\to$ `1,84e-5`; orde empiris `[1,0032; 1,0016; 1,0008]` |
| `c3_quote_symmetry_and_monotonicity` | `tercapai` | `quote_symmetry_max_err = 9,38e-15`; $\delta^b(0,0)=\delta^a(0,0)=0{,}648923$, `monotonicity_ok = true` |
| `c4_mc_utility_matches_exact` | `tercapai` | $p_{\mathrm{emp}} = 1{,}017$; $\mathrm{CE}_{\mathrm{Rich}} = 66{,}2684 \pm 0{,}0295$ (CRN; indep `0,0365`), $z = 0{,}23$; $\mathrm{CE}_{\mathrm{cont}} = 66{,}2763 \pm 0{,}0166$, $z = 0{,}88$ vs `66,2616` |
| `c5_kolmogorov_matches_mc` | `tercapai` | Simulator kontinu: $\mathbb{E}[\mathrm{exec}]$ `106,0240` vs `106,0105` ($z=0{,}59$); $\operatorname{Var}(q_T)$ `62,1197` vs `62,0174` ($z=0{,}60$) |
| `c6_glt_dominates_benchmarks_ce` | `tercapai` | CE GLT `66,3285 ± 0,0164` > AS `66,2892` > Sym `66,0818` > Tight `23,3540`; `z(GLT-AS) = 22,54` |
| `c7_as_adjustment_improves_ce` | `tercapai` | Reduksi $\varepsilon=0$: `unadj` `66,3191` vs `exact_as_hjb` `66,3191` ($\operatorname{Var}(q_T)$ `60,738` vs `60,740`); pada $\varepsilon=0{,}1$: `60,7627` vs `57,5183` |
| `c8_wrong_sigma_ce_degradation` | `tercapai` | CE `0,5σ`: `66,3074 ± 0,0166`; `1,0σ`: `66,3148 ± 0,0165`; `1,5σ`: `66,2944 ± 0,0164`; `max_deg = 0,0203 < 0,1` |
| `c9_runtime_under_300s` | `tercapai` | Total waktu `53,53s < 300s` (`8` worker) |

Hasil tambahan & audit konvergensi $\Delta t$:
- **Konvergensi $\mathrm{CE}_{\mathrm{MC}}(\Delta t)$**: `[66,3901; 66,3285; 66,2981; 66,3092]` pada $\Delta t \in \{10^{-3}; 5\times 10^{-4}; 2{,}5\times 10^{-4}; 1{,}25\times 10^{-4}\}$.
- **Audit $\operatorname{Var}(q_T)$ diskret vs Kolmogorov (`62,0174`)**: `[59,4008; 60,6586; 61,5494; 61,7044]` dengan selisih `[2,6166; 1,3589; 0,4680; 0,3130]`, menyusut mendekati $\times 2$ tiap pembelahan $\Delta t$.
- **Celah spektral & batas stasioner**: $\Lambda_1 - \Lambda_0 = 1{,}277564$; beda kuotasi maksimum vs *ground state* turun dari `0,007579` ($T=1$) menjadi `0` ($T=50$).
- **Salah spesifikasi intensitas & drift**: GLT mengungguli AS di keempat variasi $A/k \pm 30\%$ (`46,4034..88,7541` vs `46,3926..88,6100`); *drift* $\mu = \pm 0{,}2$ menurunkan CE secara simetris dari `66,2920` ke `66,2553` / `66,2549`.

## 8. Catatan Audit Perbaikan (`c4`, `c5`, `c7`, `c8`)

Seluruh 9 kriteria (`c1`–`c9`) kini berstatus `tercapai` setelah tiga perbaikan akar masalah:
1. **Orde Richardson (`c4`)**: Ekstrapolasi awal memakai orde-2 keliru; orde empiris $p = 1{,}017$ dari $\Delta t \in \{10^{-3}; 5\times 10^{-4}; 2{,}5\times 10^{-4}\}$ dan propagasi kovarian CRN memberi $\mathrm{CE}_{\mathrm{Rich}} = 66{,}2684 \pm 0{,}0295$ ($z = 0{,}23$).
2. **Diskretisasi Bernoulli vs Simulator Kontinu (`c5`)**: Uji $\Delta t$ menunjukkan selisih $\operatorname{Var}(q_T)$ diskret terhadap Kolmogorov menyusut (`2,617` → `1,359` → `0,468` → `0,313`), membuktikan efek $O(\Delta t)$ dari eksekusi simultan *bid+ask*; simulator kontinu *Lewis–Shedler thinning* bebas-$\Delta t$ cocok dengan Kolmogorov dalam `< 1 SE`.
3. **Tabel $(t, q)$ untuk `exact_as_hjb` (`c7`)**: Implementasi awal hanya menyimpan irisan $t=0$ sehingga gagal kembali ke GLT saat $\varepsilon = 0$ ($\operatorname{Var}(q_T) = 38{,}18$). Dengan tabel waktu penuh $(t, q)$ di `hjb_fd::solve_nonlinear_hjb_theta_table`, reduksi $\varepsilon = 0$ pulih (`CE = 66,3191` vs `66,3191`, galat kuotasi maks `2,74e-3` pada $\Delta t = 5\times 10^{-4}$ dan menyusut $O(\Delta t)$).

## 9. Landasan Teori (11 Poin)

1. **Model Harga Mid Difusif & Risiko Inventori**: Harga mid berevolusi mengikuti gerak Brown $dS_t = \sigma\,dW_t + \mu\,dt$, sehingga inventori $q_t \neq 0$ mengekspos kekayaan *mark-to-market* terhadap variansi $\sigma^2 q_t^2\,dt$; dipakai sebagai dinamika dasar di [sim.rs](../../projects/rust/day004g_market_making/src/sim.rs).
2. **Fungsi Intensitas Eksekusi Eksponensial**: Laju kedatangan order pasar yang menyentuh kuotasi berjarak $\delta$ berbentuk $\lambda(\delta) = A e^{-k\delta}$ ($A, k > 0$), mencerminkan penurunan probabilitas eksekusi secara eksponensial saat spread melebar; diimplementasikan di [glt.rs](../../projects/rust/day004g_market_making/src/glt.rs) dan [sim.rs](../../projects/rust/day004g_market_making/src/sim.rs).
3. **Utilitas CARA dan Certainty Equivalent**: Agen memaksimalkan $\mathbb{E}[-e^{-\gamma W_T}]$ dengan koefisien penghindaran risiko absolut $\gamma > 0$, yang memisahkan kekayaan tunai dan inventori ke dalam fungsi nilai $u = -e^{-\gamma(x + qs + \theta(t,q))}$ serta menghasilkan $\mathrm{CE} = -\frac{1}{\gamma}\ln\mathbb{E}[e^{-\gamma W_T}]$; dihitung dengan *log-sum-exp* di [sim.rs](../../projects/rust/day004g_market_making/src/sim.rs).
4. **Persamaan Hamilton–Jacobi–Bellman (HJB) Nonlinier**: Kondisi optimalitas kontinu pada ruang keadaan terbatas $q \in [-Q, Q]$ menghasilkan PDE/ODE nonlinier untuk $\theta(t,q)$ dengan kondisi batas refleksi $\mathbf{1}_{q<Q}$ dan $\mathbf{1}_{q>-Q}$; diverifikasi langsung lewat skema Euler eksplisit mundur di [hjb_fd.rs](../../projects/rust/day004g_market_making/src/hjb_fd.rs).
5. **Transformasi Eksponensial Guéant–Lehalle–Tapia (GLT)**: Substitusi $v_q(t) = \exp(k\,\theta(t, q))$ mereduksi persamaan HJB kuasi-linier menjadi sistem ODE linier $\frac{dv}{dt} = M v$ dengan matriks tridiagonal simetris $M_{q,q}=\alpha q^2$, $M_{q,q\pm 1}=-\eta$; diimplementasikan di [glt.rs](../../projects/rust/day004g_market_making/src/glt.rs).
6. **Dekomposisi Spektral Matriks Simetris (Jacobi Siklik)**: Karena $M$ simetris riil berukuran $(2Q+1)\times(2Q+1)$, $v(t) = V \operatorname{diag}(e^{-\Lambda_i(T-t)}) V^\top \mathbf{1}$ dapat dihitung eksak tanpa langkah waktu; diimplementasikan di [linalg.rs](../../projects/rust/day004g_market_making/src/linalg.rs).
7. **Eksponensial Matriks Scaling-and-Squaring & Integrasi RK4-Richardson**: Dua metode pembanding independen ($\exp(-M(T-t))\mathbf{1}$ dan RK4 orde-4 berekstrapolasi) dipakai untuk memastikan galat solver eksak berada di bawah $10^{-10}$; diimplementasikan di [linalg.rs](../../projects/rust/day004g_market_making/src/linalg.rs) dan [ode.rs](../../projects/rust/day004g_market_making/src/ode.rs).
8. **Celah Spektral dan Limit Stasioner Asimtotik**: Selisih dua nilai eigen terkecil $\Lambda_1 - \Lambda_0$ mengatur laju konvergensi $v(t)$ menuju vektor eigen dasar $\psi_0$ saat $T - t \to \infty$, menghasilkan kuotasi stasioner independen waktu; dihitung di [glt.rs](../../projects/rust/day004g_market_making/src/glt.rs).
9. **Persamaan Forward Kolmogorov**: Evolusi deterministik distribusi probabilitas inventori $P_q(t)$ mengikuti sistem ODE linier yang menghubungkan laju transisi $\lambda^b_q, \lambda^a_q$ ke $\mathbb{E}[\mathrm{eksekusi}]$ dan $\operatorname{Var}(q_T)$; diimplementasikan di [kolmogorov.rs](../../projects/rust/day004g_market_making/src/kolmogorov.rs).
10. **Hampiran Kuotasi Avellaneda–Stoikov (AS)**: Strategi tanpa batas keras menggunakan *reservation price* $r = s - q\gamma\sigma^2(T-t)$ dan spread optimal asimtotik yang dipotong di batas $\pm Q$; diimplementasikan sebagai pembanding utama di [strategies.rs](../../projects/rust/day004g_market_making/src/strategies.rs).
11. **Seleksi Merugikan (Adverse Selection) & Common Random Numbers (CRN)**: Dampak informasi privat dimodelkan sebagai lompatan harga pasca-eksekusi $\pm\varepsilon$ yang memodifikasi kuotasi HJB (`exact_as_hjb`), sementara CRN dengan `xoshiro256++` per indeks lintasan memperkecil variansi selisih antar-strategi; diimplementasikan di [strategies.rs](../../projects/rust/day004g_market_making/src/strategies.rs), [misspec.rs](../../projects/rust/day004g_market_making/src/misspec.rs), dan [rng.rs](../../projects/rust/day004g_market_making/src/rng.rs).

## 10. Kesimpulan dan Batasan

1. **Kesimpulan**: Reduksi GLT terverifikasi hingga ketelitian mesin (`2,80e-11` antar-solver eksak, orde konvergensi FD `1,001–1,003`), $\mathrm{CE}_{\mathrm{MC}}$ cocok terhadap $\theta(0,0)$ ($z = 2{,}20$), dan GLT secara konsisten mengungguli Avellaneda–Stoikov ($z = 22{,}54$), kuotasi simetris tetap ($z = 43{,}85$), dan kuotasi rapat. Studi salah spesifikasi menunjukkan *adverse selection* $\varepsilon = 0{,}1$ menimbulkan kerugian CE terbesar (`~8,8` unit), jauh melampaui salah spesifikasi $\sigma$ (`~0,02`) maupun *drift* $\mu = \pm 0{,}2$ (`~0,037`).
2. **Batasan**: Langkah waktu diskret $\Delta t = 5\times 10^{-4}$ menyisakan bias diskretisasi $O(\Delta t)$ pada variansi inventori (`c5 = belum tercapai`); seluruh hasil MC bergantung pada `main_seed = 20261122` dan $R = 200.000$. Catatan: `run_completed_without_errors = true` menandakan eksekusi selesai tanpa galat numerik, bukan bahwa seluruh kriteria ilmiah tercapai.

## 11. Referensi dan Jejak Bukti

- Guéant, O., Lehalle, C.-A., & Fernandez-Tapia, J. (2013). *Dealing with the inventory risk: a solution to the market making problem*. Mathematics and Financial Economics, 7(4), 477–507.
- Avellaneda, M., & Stoikov, S. (2008). *High-frequency trading in a limit order book*. Quantitative Finance, 8(3), 217–224.
- Artefak hasil: [outputs/day004g/results.json](../../outputs/day004g/results.json)
- Proyek & tes: [projects/rust/day004g_market_making/README.md](../../projects/rust/day004g_market_making/README.md), [integration_tests.rs](../../projects/rust/day004g_market_making/tests/integration_tests.rs)
- Grafik SVG (`outputs/figures/`): [quotes-and-spread](../../outputs/figures/day004g_quotes-and-spread.svg), [hjb-dt-convergence](../../outputs/figures/day004g_hjb-dt-convergence.svg), [mc-utility-dt-extrapolation](../../outputs/figures/day004g_mc-utility-dt-extrapolation.svg), [strategy-comparison](../../outputs/figures/day004g_strategy-comparison.svg), [adverse-selection-misspec](../../outputs/figures/day004g_adverse-selection-misspec.svg), [sigma-and-drift-misspec](../../outputs/figures/day004g_sigma-and-drift-misspec.svg), [example-path-s-q-wealth](../../outputs/figures/day004g_example-path-s-q-wealth.svg)
