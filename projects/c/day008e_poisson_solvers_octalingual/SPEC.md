# Spesifikasi Teknis: Solver Poisson 2D Lintas 8 Bahasa (`day008e`)

Dokumen ini adalah kontrak tunggal yang diikuti secara identik oleh **8 bahasa** (**C11**, **C++20**, **Rust**, **Go**, **Ruby**, **Python**, **TypeScript**, dan **Jupyter Notebook `.ipynb`**) untuk menyelesaikan $-\Delta u = f$ pada $(0,1)^2$, $u = 0$ di batas.

---

## 1. Diskretisasi, Solusi Buatan & Urutan Operasi

- **Kisi**: $N \times N$ titik dalam, $h = \frac{1}{N+1}$, indeks baris-utama `idx = i*N + j` ($0 \le i, j < N$), $x_i = (i+1)h, y_j = (j+1)h$.
- **Operator 5-titik (tanpa dibagi $h^2$)**:
  $$(Au)_{i,j} = (((4 c - w) - e) - s) - n$$
  dengan `w = (j > 0) ? u[idx-1] : 0`, `e = (j+1 < N) ? u[idx+1] : 0`, `s = (i > 0) ? u[idx-N] : 0`, `n = (i+1 < N) ? u[idx+N] : 0`.
- **Solusi buatan & sisi kanan**:
  $$u_{\text{exact}}(x,y) = \sin(\pi x)\sin(2\pi y) + x(1-x)y(1-y)$$
  $$f(x,y) = 5\pi^2\sin(\pi x)\sin(2\pi y) + 2x(1-x) + 2y(1-y), \qquad b_{i,j} = h^2 f(x_i, y_j).$$
- **Aritmetika & anti-FMA**:
  - Seluruh perhitungan memakai `float64` tanpa fusi FMA (`-ffp-contract=off` di C/C++, `float64(a*b)` di Go).
  - **Sinus deterministik `sin_spec(x)`**: Reduksi kuadran ke $[0, \pi/2]$ dengan `PI = 3.14159265358979323846`, lalu polinomial Horner 10-suku $\sum_{k=0}^{9} c_k x^{2k+1}$ (`c = [1, -1/6, 1/120, -1/5040, 1/362880, -1/39916800, 1/6227020800, -1/1307674368000, 1/355687428096000, -1/121645100408832000]`) agar bebas perbedaan 1-ulp `libm` antar bahasa/runtime.
  - **Penjumlahan titik & norma**: `s = 0.0; for k = 0..N*N-1: s = s + a[k]*b[k]`, `norm2(a) = sqrt(dot(a,a))`.
  - **FNV-1a 64-bit**: Atas pola bit `uint64` little-endian tiap `float64`, basis `0xCBF29CE484222325`, prima `0x100000001B3`.
- **Konstanta $\omega_{\text{opt}} = \frac{2}{1 + \sin(\pi h)}$ (17 digit desimal)**:
  - `N=15`: `1.673513677715992`
  - `N=31`: `1.8214651907890225`
  - `N=63`: `1.906454701582762`
  - `N=127`: `1.952093233850055`
  - `N=255`: `1.9757544535797149`
  - `N=511`: `1.9878030696593354`
- **Toleransi berhenti**: `||r||_2 / ||b||_2 <= 1e-8`, tebakan awal $u^{(0)} = 0$.

---

## 2. Metode Solver & Anchor Analitik

- **Jacobi**: `u_new[idx] = 0.25 * ((((b[idx] + w) + e) + s) + n)` (`max_iter = 100000`).
- **GS (leksikografis)**: in-place `i=0..N-1, j=0..N-1`, rumus sama (`max_iter = 50000`).
- **RB-GS**: `color = 0` (`(i+j)&1 == 0`) lalu `color = 1` (`max_iter = 50000`).
- **SOR**: in-place `u[idx] = (1.0 - om)*u[idx] + om*(0.25*((((b[idx]+w)+e)+s)+n))` (`max_iter = 10000`).
- **CG**: $r_0 = b, p_0 = b$, $\alpha = (r,r)/(p,Ap)$, $\beta = (r_{\text{new}},r_{\text{new}})/(r,r)$ (`max_iter = 5000`).
- **PCG-SSOR ($\omega=1$)**: sapuan maju $(D-L)y = r$, mundur $(D-U)z = Dy$, $\alpha = (r,z)/(p,Ap)$, berhenti saat $\|r\|_2/\|b\|_2 \le 10^{-8}$ (`max_iter = 3000`).
- **Multigrid V(2,1)**: 2 pre-RB-GS, residual, restriksi *full-weighting* dikali `4.0`, rekursi hingga `N=3` (eliminasi Gauss 9×9 *partial pivoting*), prolongasi bilinear, 1 post-RB-GS (`max_iter = 100`).
- **FMG**: restriksi $b$ hingga `N=3`, selesaikan `N=3`, prolongasi naik + 1 siklus V(2,1) tiap tingkat.
- **Anchor analitik**: $\lambda_{jk} = 4 - 2\cos(j\pi h) - 2\cos(k\pi h)$, $\rho_J = \cos(\pi h)$, $\rho_{GS} = \cos^2(\pi h)$, $\kappa = \cot^2(\pi h/2)$, estimasi batas CG $\approx \frac{1}{2}\sqrt{\kappa}\ln(2/\varepsilon)$.

---

## 3. Beban Kerja & Kriteria A Priori (`c1`–`c12`)

- **Beban $N$**:
  - `c, cpp, rust, go`: `N in {15,31,63,127,255,511}` (SOR, CG, PCG, MG, FMG) dan `N in {15,31,63,127}` (Jacobi, GS, RB-GS).
  - `ts`: `N in {15,31,63,127,255}` (SOR/CG/PCG/MG/FMG) dan `{15,31,63}` (Jacobi/GS/RB-GS).
  - `ruby, py`: `N in {15,31,63,127}` (SOR/CG/PCG/MG/FMG) dan `{15,31,63}` (Jacobi/GS/RB-GS).
  - `ipynb`: `N in {15,31,63}` semua metode.
- **Kriteria (`c1`–`c12`)**: `c1` paritas S01+S17 bit-identik ($N=15,31$ semua metode); `c2` anchor eigen $\le 10^{-12}$; `c3` residu $\le 10^{-8}$ & selisih solusi $N=63 \le 10^{-6}$; `c4` orde diskretisasi $\in [1.9, 2.1]$; `c5` $\rho_J, \rho_{GS}$ ($N=31$) dalam 1% teori; `c6` eksponen iterasi vs $N$: Jacobi/GS $\in [1.8,2.2]$, SOR/CG $\in [0.8,1.2]$, PCG-SSOR $\in [0.4,0.7]$, MG $\in [-0.1,0.1]$; `c7` faktor MG $\le 0.2$ ($N=63$); `c8` FMG $\le 3\times$ galat diskretisasi ($N=63$); `c9` `all` selesai; `c10` `selftest` lulus; `c11` SVG valid XML; `c12` kecepatan & baris kode dilaporkan.
