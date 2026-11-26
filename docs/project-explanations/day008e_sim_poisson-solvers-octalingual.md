# Hari 8e (`day008e`) — Solver Poisson 2D Lintas 8 Bahasa (19 Tahap Identik: C, C++, Rust, Go, Ruby, Python, TS, ipynb)

## 1. Identitas Sesi
- **Kode Sesi**: `day008e` (`bonus` · `sim` · `c+cpp+rust+go+ruby+py+ts+ipynb`)
- **Tanggal**: `2026-11-26` (WITA)
- **Status**: `sebagian` (`11/12` kriteria `c1`–`c12` `tercapai`; `c6_scaling_exponents` berstatus `sebagian` karena `pcg_ssor` dengan $\omega = 1$ tetap memiliki eksponen `0,889428`, di luar rentang `[0,4, 0,7]` — dilaporkan jujur tanpa tuning ulang $\omega$)
- **Artifak Utama**:
  - Proyek oktalingual: `projects/c/day008e_poisson_solvers_octalingual/` (`SPEC.md`, `meta.json`, `Makefile`, `run_all.sh`, `compare.py`, `README.md`, 7 implementasi CLI, serta `out/<lang>/`)
  - Notebook Jupyter: `days/day008e_sim_poisson-solvers-octalingual.ipynb`
  - Ringkasan terstruktur: `outputs/day008e/results.json`
  - 4 pasang grafik gabungan PNG+SVG: `outputs/figures/day008e_{history,scaling,speed,parity}.{png,svg}`

## 2. Masalah & Tujuan
Menyelesaikan persamaan Poisson 2D dengan syarat batas Dirichlet homogen pada domain persegi satuan:
$$-\Delta u(x,y) = f(x,y) \quad \text{pada } (0,1)^2, \qquad u|_{\partial (0,1)^2} = 0$$
menggunakan diskritisasi beda hingga 5-titik orde 2 pada grid $N \times N$ ($h = 1/(N+1)$, $N \in \{15, 31, 63, 127, 255, 511\}$) dengan solusi analitik kombinasi mode eigen dan suku polinomial:
$$u^*(x,y) = \sin(\pi x)\sin(2\pi y) + x(1-x)y(1-y)$$
$$f(x,y) = 5\pi^2 \sin(\pi x)\sin(2\pi y) + 2x(1-x) + 2y(1-y)$$
Tujuan utamanya adalah membangun **alur 19 tahap yang identik (`s01`–`s19`, `all`) secara independen di 8 bahasa** (C11, C++20, Rust, Go, Ruby, Python, TypeScript, dan Jupyter Notebook), memverifikasi paritas **bit-identik IEEE-754** lintas kedelapan bahasa melalui FNV-1a 64-bit, memvalidasi anchor analitik (nilai eigen, orde diskritisasi, radius spektral, faktor kontraksi multigrid, dan FMG), serta membandingkan kecepatan (`ns/unknown-iter`) dan ukuran kode sumber (`<= 700` baris).

## 3. Arsitektur 19 Tahap (`s01`–`s19`) & Matriks Bahasa
Setiap bahasa mengimplementasikan seluruh alur secara mandiri dan menulis berkas keluaran bernama seragam ke `out/<lang>/`:
- `s01_selftest.txt`: uji operator $A u$ pada $N=3$ dan tiga hash FNV-1a 64-bit (`h_u3 = c56de10c6671e6c2`, `h_Au3 = 59ed72ec82ee959c`, `h_test = afc5dc1ca48a2366`).
- `s02_problem.csv`: norma $\|b\|_2$ dan hash `h_b`, `h_ue` untuk setiap $N$.
- `s03_eigen.csv`: verifikasi analitik $A v_{j,k} = \lambda_{j,k} v_{j,k}$ untuk $(j,k) \in \{(1,1),(1,2),(3,5)\}$ pada $N \in \{15, 31, 63\}$.
- `s04`–`s10`: 7 solver (`jacobi`, `gs`, `rbgs`, `sor` dengan $\omega_{\text{opt}}$, `cg`, `pcg_ssor` dengan $\omega=1$, dan `mg` V(2,1)) hingga residual relatif $\|r_k\|_2 / \|b\|_2 \le 10^{-8}$.
- `s11_fmg.csv`: Full Multigrid (1 siklus FMG + V(2,1)) vs galat diskritisasi.
- `s12_order.csv`: regresi log–log galat $\|u - u^*\|_\infty$ terhadap $h$.
- `s13_scaling.csv`: regresi log–log jumlah iterasi terhadap $N$ beserta galat baku kemiringan.
- `s14_history.csv`: riwayat residual relatif $N=63$ (disampling hingga ~200 titik per metode).
- `s15_mg_factor.csv`: faktor kontraksi rata-rata geometrik 5 siklus terakhir MG V(2,1).
- `s16_cost.csv`: tolok ukur kecepatan `ns/unknown-iter` (median $\pm$ MAD dari 7 pengulangan pada $N=63$).
- `s17_parity.csv`: tabel paritas FNV-1a 64-bit (`iters`, `hist_fnv1a`, `sol_fnv1a`) untuk 7 metode $\times$ $N \in \{15, 31\}$.
- `s18` (`figures/{history,scaling,order}.svg`): 3 SVG mandiri per bahasa (total 24 SVG).
- `s19_summary.json`: metadata bahasa, batas $N$, toolchain, dan waktu total.

## 4. Grafik Gabungan (`outputs/figures/`)
- `outputs/figures/day008e_history.{png,svg}`: kurva konvergensi residual relatif $N=63$ untuk 7 metode.
- `outputs/figures/day008e_scaling.{png,svg}`: grafik log–log jumlah iterasi terhadap $N$ beserta kemiringan empiris.
- `outputs/figures/day008e_speed.{png,svg}`: perbandingan kecepatan kernel `rbgs` dan `mg` (`ns/unknown-iter`, median $\pm$ MAD) di 8 bahasa.
- `outputs/figures/day008e_parity.{png,svg}`: matriks paritas bit-identik $8 \times 14$ dan diagram batang jumlah baris sumber (`<= 700`).

## 5. Performa & Ukuran Kode Lintas 8 Bahasa (`N=63`, median $\pm$ MAD)

| Bahasa | Baris Sumber | $N_{\max}$ | `rbgs` (`ns/unk-iter`) | `mg` (`ns/unk-iter`) | Waktu Total (`s`) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **C11** | `500` | `511` | `0.9543 ± 0.0378` | `5.4926 ± 0.1008` | `12.107` |
| **C++20** | `446` | `511` | `0.8564 ± 0.0094` | `5.1556 ± 0.0766` | `8.223` |
| **Rust** | `386` | `511` | `1.4121 ± 0.0273` | `8.7228 ± 0.1911` | `10.812` |
| **Go** | `698` | `511` | `1.8484 ± 0.0034` | `13.4637 ± 0.3517` | `16.006` |
| **TypeScript** | `489` | `255` | `5.7181 ± 0.0556` | `31.7555 ± 0.3895` | `2.464` |
| **Python** | `316` | `127` | `125.1029 ± 0.6593` | `811.6794 ± 3.5589` | `28.548` |
| **Notebook (`ipynb`)** | `356` (sel kode) | `63` | `124.1547 ± 0.1541` | `807.5166 ± 2.0429` | `25.188` |
| **Ruby** | `324` | `127` | `232.1492 ± 8.0373` | `1299.9076 ± 4.1992` | `45.447` |

## 6. Fakta Teramati vs Hipotesis
- **FAKTA**: Seluruh 8 bahasa menghasilkan berkas `s01_selftest.txt` dan `s17_parity.csv` yang **100% bit-identik** (`0` selisih byte pada seluruh 14 kasus $7 \text{ metode} \times \{15, 31\}$).
- **FAKTA**: `PCG-SSOR` dengan $\omega = 1$ tetap mencatat kemiringan scaling iterasi `0,889428 ± 0,012022` (`iters`: `20, 36, 66, 123, 232, 439` untuk $N=15..511$), sehingga tidak masuk rentang `[0,4, 0,7]` pada `c6`.
- **FAKTA**: Notebook `days/day008e_sim_poisson-solvers-octalingual.ipynb` dieksekusi bersih **dua kali** lewat `nbclient` dan menghasilkan SHA-256 identik (`c9f8ee25088aa8d0a452816db5738bb9e5d0a013aebd0e2a73e2306d99a8b337`).
- **HIPOTESIS**: Scaling `PCG-SSOR` mendekati $O(N)$ pada $\omega = 1$ (alih-alih $O(N^{1/2})$ seperti $\omega_{\text{opt}}$ SSOR pada grid 2D) terjadi karena tanpa parameter relaksasi optimal $\omega \approx 2 / (1 + \sin(\pi h))$, bilangan kondisi matriks terprakondisi $\kappa(M_{\text{SSOR}}^{-1} A)$ tetap tumbuh secara kuadratik terhadap $N$ ($\kappa \propto N^2$, bukan $N$), sehingga jumlah iterasi CG tumbuh mendekati $O(N)$ alih-alih $O(N^{1/2})$.

## 7. Ringkasan Hasil & Status Kriteria (`c1`–`c12`)

| Kriteria | Status | Bukti & Hasil Terukur (`outputs/day008e/results.json`) |
| :--- | :---: | :--- |
| `c1_parity_8lang` | `tercapai` | `s01` SHA-256 `e96e6acf6688584c...`, `s17` SHA-256 `0a49eae379bdc868...`, `0` selisih di 8 bahasa |
| `c2_eigen_anchor_le_1e-12` | `tercapai` | `max_eig_rel_err = 3.165e-13 <= 1e-12` (seluruh mode $(1,1),(1,2),(3,5)$ pada $N \le 63$) |
| `c3_res_le_1e-8_and_agree_le_1e-6` | `tercapai` | `max_rel_res = 9.999975e-09 <= 1e-8`, `sol_agree_63_inf = 2.501636e-08 <= 1e-6` |
| `c4_order_in_1.9_2.1` | `tercapai` | `disc_order_slope = 2.001445` ($N \in \{15..511\}$ pada C11/C++/Rust/Go) |
| `c5_rho_n31_within_1pct` | `tercapai` | Jacobi `0.9951847268` vs `0.9951847267` (`1.28e-10`), GS `0.9903926782` vs `0.9903926402` (`3.84e-08`) |
| `c6_scaling_exponents` | `sebagian` | Jac `1.950`, GS `1.962`, RB-GS `1.950`, SOR `1.005`, CG `0.995`, MG `0.043` lulus; **PCG-SSOR `0.889`** di luar `[0,4, 0,7]` pada $\omega=1$ |
| `c7_mg_factor_n63_le_0.2` | `tercapai` | `mg_factor_63 = 0.07011662 <= 0.2` dan `mg_factor_127 = 0.07100588` |
| `c8_fmg_ratio_n63_le_3` | `tercapai` | `fmg_ratio_63 = 0.473908 <= 3.0` |
| `c9_all_finished` | `tercapai` | `8/8` bahasa selesai bersih menulis `s19_summary.json` |
| `c10_selftest_pass` | `tercapai` | `8/8` `selftest` (`s01`) lulus (`ok = 1`) |
| `c11_svgs_valid` | `tercapai` | `24` SVG per bahasa + `4` SVG gabungan non-kosong & lolos `xml.etree.ElementTree` |
| `c12_speed_and_lines_le_700` | `tercapai` | Seluruh sumber `<= 698` baris; tolok ukur `ns/unknown-iter` median $\pm$ MAD tercatat |

## 8. Deviasi & Catatan Rekayasa
1. **Fallback `tsc` ke Node `--experimental-strip-types`**: `tsc` tidak terpasang pada lingkungan ini, sehingga `ts/poisson.ts` dijalankan langsung melalui `node v24.21.0 --experimental-strip-types`.
2. **Pencegahan FMA & `libm` lintas platform**:
   - `sin_spec(x)` memakai polinomial Horner 10-suku dengan reduksi rentang berulang (`while x > PI_VAL`) agar akurat hingga $\le 10^{-12}$ pada mode $(3,5)$ sekalipun dan identik antara `libm` macOS dan V8 `Math.sin` (yang berbeda 1 ULP).
   - Pada Go (`go 1.27.4` arm64), compiler menggabungkan perkalian-penjumlahan menjadi instruksi `FMADDD` meskipun dibungkus konversi `float64(...)`; pembantu `//go:noinline func fm(a, b float64) float64` dan `var piVal` non-konstanta dipakai untuk menjamin paritas bit-identik IEEE-754 dengan C/C++/Rust/Ruby/Python/TS/ipynb.
3. **Beban kerja terukur per kelas bahasa**: C, C++, Rust, dan Go menjalankan hingga $N=511$ ($N_{\text{stat},\max}=127$); TypeScript hingga $N=255$ ($N_{\text{stat},\max}=63$); Ruby dan Python hingga $N=127$ ($N_{\text{stat},\max}=63$); notebook hingga $N=63$.
4. **Status `c6` jujur**: Nilai $\omega = 1$ pada `PCG-SSOR` dipertahankan sesuai spesifikasi dan kegagalan `c6` dicatat sebagai `sebagian` tanpa mengubah parameter.

## 9. Landasan Teori (14 Poin)

1. **Operator Laplacian 2D 5-Titik** — Diskritisasi beda hingga orde 2 pada grid seragam $h = 1/(N+1)$ menghasilkan matriks simetris definit positif $A \in \mathbb{R}^{N^2 \times N^2}$ berbentuk $A = I_N \otimes T_N + T_N \otimes I_N$ dengan $T_N = \text{tridiag}(-1, 2, -1)$, sehingga kondisi batas Dirichlet nol masuk langsung ke suku kanan $b_{i,j} = h^2 f(x_i, y_j)$.
2. **Spektrum Eksak & Mode Eigen** — Vektor eigen diskret $v_{j,k}(r,c) = \sin(j\pi r h)\sin(k\pi c h)$ memenuhi $A v_{j,k} = \lambda_{j,k} v_{j,k}$ dengan $\lambda_{j,k} = 4 - 2\cos(j\pi h) - 2\cos(k\pi h)$, memberi bilangan kondisi $\kappa(A) = \lambda_{N,N}/\lambda_{1,1} \approx 4/(\pi^2 h^2) = O(N^2)$.
3. **Konvergensi Jacobi & Gauss–Seidel** — Untuk matriks 2-siklik konsisten terurut, radius spektral matriks iterasi Jacobi adalah $\rho(G_J) = \cos(\pi h) \approx 1 - \frac{1}{2}\pi^2 h^2$, sedangkan Gauss–Seidel memenuhi $\rho(G_{\text{GS}}) = \rho(G_J)^2 = \cos^2(\pi h) \approx 1 - \pi^2 h^2$, sehingga keduanya memerlukan $O(N^2)$ iterasi.
4. **Red–Black Gauss–Seidel** — Pewarnaan papan catur $(i+j) \bmod 2$ mempartisi titik grid menjadi dua himpunan independen karena stensil 5-titik hanya menghubungkan warna merah dan hitam; pembaruan tetap memiliki $\rho = \cos^2(\pi h)$ namun vektorisasi/paralelisme lebih ramah dan urutan floating-point tetap deterministik.
5. **Teori SOR Young & Parameter Optimal** — Successive Over-Relaxation menggeser nilai eigen dengan faktor relaksasi $\omega \in (0, 2)$; parameter optimal $\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho(G_J)^2}} = \frac{2}{1 + \sin(\pi h)}$ menekan radius spektral menjadi $\rho(G_{\text{SOR}}) = \omega_{\text{opt}} - 1 \approx 1 - 2\pi h$, menurunkan kompleksitas iterasi dari $O(N^2)$ menjadi $O(N)$.
6. **Conjugate Gradient (CG)** — Karena $A$ simetris definit positif, CG meminimalkan norma energi $\|u - u^*\|_A$ pada subruang Krylov $\mathcal{K}_k(A, r_0)$ dengan batas konvergensi $\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1} \approx 1 - \frac{2}{\sqrt{\kappa}} = 1 - O(h)$, menghasilkan penskalaan iterasi $O(N)$ pada 2D tanpa parameter $\omega$.
7. **Prakondisioner SSOR Simetris** — Matriks prakondisioner $M_{\text{SSOR}} = \frac{\omega}{2-\omega}(\frac{1}{\omega}D - L) D^{-1} (\frac{1}{\omega}D - U)$ bersifat SPD dan terfaktorisasi sebagai sapuan maju-mundur segitiga tanpa alokasi matriks tambahan; dengan $\omega = 1$ tetap, pengurangan $\kappa$ bersifat konstan sehingga penskalaan tetap mendekati $O(N)$, sedangkan $\omega \approx \omega_{\text{opt}}$ diperlukan untuk mencapai $O(N^{1/2})$.
8. **Sifat Smoothing & Frekuensi Tinggi** — Smoother lokal seperti RB-GS meredam komponen galat frekuensi tinggi ($\theta_x, \theta_y \in [\pi/2, \pi]$) dengan faktor konstan $\le 0{,}25$ dalam 2–3 sapuan, sementara mode frekuensi rendah ($\theta \approx \pi h$) hampir tidak teredam dan diproyeksikan ke grid kasar.
9. **Operator Restriksi Full-Weighting & Prolongasi Bilinear** — Operator transfer $I_h^{2h}$ (bobot $1/4$ pusat, $1/8$ sisi, $1/16$ sudut dikalikan faktor 4 untuk menyamakan skala operator $A_h$ tanpa $1/h^2$) merupakan adjoin terskala dari interpolasi bilinear $I_{2h}^h$, menjaga simetri dan kestabilan koarse-grid correction.
10. **Siklus V(2,1) Multigrid & Kompleksitas Optimal** — Dengan 2 pra-sapuan RB-GS, koreksi grid kasar rekursif hingga $N=3$ (eliminasi Gauss $9\times 9$), dan 1 pasca-sapuan RB-GS, faktor kontraksi per siklus bebas dari $h$ ($\approx 0{,}070 \le 0{,}2$), menghasilkan kompleksitas total $O(N^2)$ yang optimal secara asimtotik.
11. **Full Multigrid (FMG)** — Dengan membangun tebakan awal pada grid halus melalui interpolasi solusi grid lebih kasar secara rekursif diikuti 1 siklus V(2,1), FMG langsung mencapai tingkat galat diskritisasi $O(h^2)$ dalam satu lintasan bernilai biaya $< 2$ siklus V.
12. **Orde Diskritisasi & Stensil 5-Titik** — Ekspansi Taylor pada $x_i \pm h$ memberi galat pemotongan lokal $\frac{h^2}{12}(\partial_x^4 u + \partial_y^4 u) = O(h^2)$; oleh prinsip maksimum diskret untuk matriks-M, galat global $\|u_h - u^*\|_\infty$ menurun tepat dengan orde 2 ($\text{kemiringan } \approx 2{,}00$).
13. **Determinisme IEEE-754 & Fused Multiply-Add (FMA)** — Instruksi `FMA(a,b,c)` melakukan satu pembulatan akhir alih-alih dua pembulatan pada `round(a*b) + c`, yang cukup memicu selisih 1 ULP dan menyebarkan divergensi iterasi; penonaktifan FMA dan polinomial trigonometri bersama menjamin paritas bit-identik lintas 8 arsitektur bahasa.
14. **Benchmarking Robust dengan Median dan MAD** — Pengukuran latensi `ns/unknown-iter` dari 7 pengulangan diringkas dengan median dan *median absolute deviation* (MAD) karena distribusi waktu eksekusi bersudut kanan akibat *jitter* penjadwal OS dan pengumpul sampah (*garbage collection*).

## 10. Reproduksibilitas & Cara Menjalankan
```bash
make -C projects/c/day008e_poisson_solvers_octalingual test   # build & selftest 7 CLI
make -C projects/c/day008e_poisson_solvers_octalingual run    # jalankan 8 bahasa + compare.py
.venv/bin/python scripts/check_docs.py
.venv/bin/python scripts/update_readme.py --check
```

## 11. Referensi
1. Briggs, W. L., Henson, V. E., & McCormick, S. F. (2000). *A Multigrid Tutorial* (2nd ed.). SIAM.
2. Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.). SIAM.
3. Young, D. M. (1971). *Iterative Solution of Large Linear Systems*. Academic Press.
4. Trottenberg, U., Oosterlee, C. W., & Schüller, A. (2001). *Multigrid*. Academic Press.
5. Goldberg, D. (1991). What every computer scientist should know about floating-point arithmetic. *ACM Computing Surveys*, 23(1), 5–48.
