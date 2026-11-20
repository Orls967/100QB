# Penjelasan Proyek — day003 (alias 3a): Volatilitas Lokal Dupire dari Permukaan SSVI (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day003` / alias `day003a` (Hari 3, Sesi a / pertama)
- **Tanggal**: 2026-11-21 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day003_quant_volatilitas-lokal-dupire-ssvi.ipynb](../../days/day003_quant_volatilitas-lokal-dupire-ssvi.ipynb) · **Paket**: [src/qb100/day003a_local_vol/](../../src/qb100/day003a_local_vol/__init__.py)
- **Ringkasan**: Permukaan total variansi implisit SSVI Gatheral–Jacquier bebas-arbitrase, volatilitas lokal Dupire analitik, dua solver PDE independen (backward di $x=\ln S$ dan forward Dupire di $K$), Monte Carlo log-Euler, opsi barrier Down-and-Out (Reiner–Rubinstein, PDE, MC jembatan Brownian), Greeks & diagnosis lapisan batas gamma di $S_0$, vega-bucket, serta kalibrasi balik Levenberg–Marquardt ($R=50$).

## 2. Masalah dan Tujuan
Bagaimana membangun permukaan *implied volatility* parametrik SSVI yang bebas dari arbitrase kupu-kupu dan kalender, menurunkan volatilitas lokal Dupire $\sigma_{\text{loc}}(S,t)$ secara analitik, dan memverifikasi bahwa dua solver PDE independen (backward di spot dan forward di strike) serta simulasi Monte Carlo menghasilkan harga Eropa, opsi barrier, dan sensitivitas yang konsisten?

## 3. Apa yang Dikerjakan Proyek
Paket `src/qb100/day003a_local_vol/` dan notebook `day003` menjalankan seluruh pipeline Bagian 1–4 beserta diagnosis tambahan dari sesi `day003c` (Bagian A1: diagnosis gamma SSVI di sudut singular $S_0$ dan tabel penyusutan bias Euler MC terhadap jumlah langkah), lalu menyimpan [outputs/day003a/results.json](../../outputs/day003a/results.json) dan 7 grafik PNG.

## 4. Uraian Kode
- [src/qb100/day003a_local_vol/ssvi.py](../../src/qb100/day003a_local_vol/ssvi.py) & [dupire.py](../../src/qb100/day003a_local_vol/dupire.py): rumus SSVI, turunan analitik $\partial_y w, \partial_{yy} w, \partial_T w$, fungsi densitas $g(y,T)$, uji Gatheral–Jacquier, rumus Dupire analitik, dan pembanding FD orde-4.
- [src/qb100/day003a_local_vol/grids.py](../../src/qb100/day003a_local_vol/grids.py), [pde_backward.py](../../src/qb100/day003a_local_vol/pde_backward.py), [pde_forward.py](../../src/qb100/day003a_local_vol/pde_forward.py): generator grid `sinh`, skema Crank–Nicolson + 4 langkah Rannacher, opsi `x_bounds` tetap, dan cell-averaging kink.
- [src/qb100/day003a_local_vol/mc.py](../../src/qb100/day003a_local_vol/mc.py) & [barrier.py](../../src/qb100/day003a_local_vol/barrier.py): MC log-Euler antitetik PCG64, uji weak order CRN, rumus Reiner–Rubinstein, PDE barrier absorbing, dan MC jembatan Brownian.
- [src/qb100/day003a_local_vol/greeks.py](../../src/qb100/day003a_local_vol/greeks.py) & [calibration.py](../../src/qb100/day003a_local_vol/calibration.py): `pde_delta_gamma` (dengan grid tetap untuk FD), `gamma_diagnosis_ssvi`, `vega_bucket_analysis`, dan optimizer Levenberg–Marquardt buatan sendiri.
- [tests/test_day003a_ssvi_dupire.py](../../tests/test_day003a_ssvi_dupire.py), [test_day003a_pde.py](../../tests/test_day003a_pde.py), [test_day003a_mc_barrier.py](../../tests/test_day003a_mc_barrier.py), [test_day003a_greeks_calib.py](../../tests/test_day003a_greeks_calib.py): 17 unit test pytest.

## 5. Metode dan Alur Kerja
Verifikasi turunan analitik SSVI & kondisi bebas-arbitrase $\to$ evaluasi $\sigma_{\text{loc}}^2(y,T) = (\partial_T w)/g(y,T)$ & anchor $\eta\to 0$ / FD4 $\to$ penyelesaian PDE backward & forward pada kasus datar dan SSVI ($5\times 4$ grid) $\to$ simulasi MC Euler & barrier $\to$ ekstraksi Greeks, diagnosis gamma di $S_0$, vega-bucket, dan kalibrasi balik LM $R=50$.

## 6. Output yang Ditargetkan
Turunan analitik vs FD `< 1e-6`; $g > 0, \partial_T w > 0$ dan kontrafaktual $g < 0$; anchor $\eta\to 0$ `< 1e-12` dan Dupire FD4 `< 1e-4`; orde PDE datar $\approx 2$; galat IV PDE SSVI `< 2 bp` dan selisih harga antar-solver `< 1e-5`; MC Euler & barrier dalam 3 SE; barrier flat PDE vs RR `< 1e-5`; Greeks datar `< 1e-5` (delta) dan `< 1e-4` (gamma); vega-bucket konsisten `< 1e-3`; kalibrasi LM tak-bias.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Turunan analitik SSVI, bebas-arbitrase & uji kontrafaktual | `tercapai` | `dw_dy = 9.03e-7`, `d2w_dy2 = 5.87e-7`, `dw_dT = 1.92e-9`; `min_g = 0.2605`, `min_dw_dT = 0.01767`; kontrafaktual `min_g = -1.475 < 0` |
| Anchor Dupire $\eta\to 0$ & konsistensi FD orde-4 | `tercapai` | `eta_zero_max_abs_err = 0.0`, `dupire_fd4_max_rel_err = 2.13e-9 < 1e-4` |
| Konvergensi PDE datar & ekuivalensi 2 solver pada SSVI | `tercapai` | Slope datar backward `-2.0002`, forward `-1.9519`; maks galat IV SSVI `0.127 bp` (bwd) dan `0.119 bp` (fwd); maks selisih harga antar-solver `9.55e-6 < 1e-5` |
| MC Euler, Weak Order & Diagnosis Penyusutan Bias (A1) | `tercapai` | Run utama 200 langkah/thn `max z = 2.344 < 3`, weak order slope `-0.9846`; diagnosis CRN 50/100/200/400 langkah membuktikan bias positif menyusut monoton (mis. ATM `0.1057 -> 0.0573 -> 0.0310 -> 0.0181`) |
| Opsi Barrier (Reiner–Rubinstein, PDE absorbing, MC Brownian Bridge) | `tercapai` | Flat `|PDE - RR| = 6.01e-6 < 1e-5`, paritas `0.0`; SSVI PDE vs MC bridge `z = 0.299 < 3` |
| Greeks datar, Diagnosis Gamma SSVI (A1) & Vega-Bucket | `tercapai` | Flat delta err `1.72e-6`, gamma err `2.17e-7`; di luar sudut ($S_0=90,95,105,110$) gamma grid vs FD cocok `~1e-4` relatif; di sudut $S_0=100$ terdiagnosis lapisan batas $t\to 0$; vega-bucket `rel_diff = 1.67e-5 < 1e-3` |
| Kalibrasi balik Levenberg–Marquardt & studi $R=50$ | `tercapai` | Single-fit RMS vol `3.71e-4` (10 iterasi), $R=50$ mean RMS `4.42e-4`, maks `|bias/SD| = 0.266` |
| Eksekusi Notebook (`run_all.py`) & Unit Test (`pytest -q`) | `tercapai` | `run_all.py --day 3 --inplace` (`1/1 lulus`, 56.1s) dan `pytest -q` dijalankan pada sesi `day003c`: lulus 100% |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Diagnosis jujur untuk dua temuan dari Bagian A1:
- **Gamma SSVI di $S_0 = 100$ (FAKTA TERUKUR vs HIPOTESIS)**:
  - **FAKTA**: Pada $S_0 = 100$, gamma grid (`0.03803` pada $601\times 400$; `0.03682` pada $2401\times 1600$) berbeda dari beda hingga harga (`0.02340` pada $h=0.05$; `0.04293` pada $h=0.5$).
  - **FAKTA**: Hipotesis bahwa pergeseran grid saat $S_0$ di-bump yang menyebabkan selisih telah **dikesampingkan**: dengan `x_bounds` tetap (`pde_delta_gamma` yang telah diperbarui), beda hingga harga memberikan angka yang hampir identik (`0.02340` pada $h=0.05$).
  - **FAKTA**: Di titik spot di luar pusat koordinat SSVI ($S_0 \in \{90, 95, 105, 110\}$), gamma grid dan beda hingga harga cocok sangat rapat (selisih relatif $4\times 10^{-5}$ hingga $2.3\times 10^{-4}$). Demikian pula gamma *strike* $C_{KK}(K=100)$ dari solver forward (`0.016738`) cocok dengan turunan analitik SSVI (`0.016762`).
  - **FAKTA**: Di $S_0 = 100$ tepat, gamma sensitif terhadap pemotongan waktu dekat nol `t_min` di `dupire_local_vol_from_spot`: untuk `t_min` $\in \{10^{-3}, 10^{-4}, 10^{-5}, 10^{-6}\}$, gamma grid berubah drastis `0.02887 -> 0.03692 -> 0.05099 -> 0.06401` (sementara harga hanya berubah `< 5e-4`).
  - **Penjelasan & HIPOTESIS**: Dalam model lokal-vol yang diturunkan dari permukaan SSVI berpusat di $S_0 = 100$, $y = \ln(S/F(t))$ dikalikan dengan $\varphi(\theta(t)) \sim \theta(t)^{-\gamma} \to \infty$ saat $t \to 0$. Ini menciptakan lapisan batas (*boundary layer*) selebar $|S - S_0| \sim O(\sqrt{t})$ di sekitar $(S_0, t=0)$ di mana kelengkungan $\partial_{SS} V(S, 0)$ menjadi singular bila `t_min` $\to 0$.
- **Bias Positif MC Euler (FAKTA TERUKUR)**:
  - **FAKTA**: Kelima z-score MC Euler pada 200 langkah/tahun bernilai positif (`+1.01` hingga `+2.34`). Uji CRN pada `{50, 100, 200, 400}` langkah/tahun di `mc_euler_bias_vs_steps` membuktikan bahwa ini adalah bias diskretisasi Euler orde-1 yang menyusut hampir separuh setiap kali langkah digandakan (misalnya di $K/F=1.0$: `+0.1057 -> +0.0573 -> +0.0310 -> +0.0181`).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Parametrisasi Total Variansi SSVI (Gatheral–Jacquier 2014)**: $w(y,\theta) = \frac{\theta}{2}(1 + \rho\varphi(\theta)y + \sqrt{(\varphi(\theta)y+\rho)^2 + 1-\rho^2})$ dengan fungsi power-law $\varphi(\theta) = \eta/(\theta^\gamma(1+\theta)^{1-\gamma})$.
2. **Fungsi Densitas Kupu-Kupu Gatheral $g(y,T)$**: $g(y,T) = (1 - \frac{y w_y}{2w})^2 - \frac{w_y^2}{4}(\frac{1}{w}+\frac{1}{4}) + \frac{w_{yy}}{2} \ge 0$ menjamin kepositifan densitas risiko-netral $\partial_{KK} C(K,T) \ge 0$ (`ssvi.py`).
3. **Bebas-Arbitrase Kalender & Syarat Cukup Gatheral–Jacquier**: $\partial_T w(y,T)|_y \ge 0$ bersama $\theta\varphi(\theta)(1+|\rho|) \le 4$ dan $\theta\varphi(\theta)^2(1+|\rho|) \le 4$ menjamin permukaan bebas arbitrase statis.
4. **Rumus Volatilitas Lokal Dupire (1994) dalam Koordinat Total Variansi**: $\sigma_{\text{loc}}^2(y,T) = \frac{\partial_T w(y,T)}{g(y,T)}$ menghilangkan ketergantungan eksplisit pada $(r-q)$ karena pembatalan suku drift forward (`dupire.py`).
5. **Dualitas PDE Backward Kolmogorov vs Forward Fokker–Planck / Dupire**: Solver backward menilai satu kontrak $(K,T)$ sebagai fungsi $(S,t)$, sedangkan solver forward Dupire $\partial_T C = \frac{1}{2}\sigma_{\text{loc}}^2 K^2 C_{KK} - (r-q)K C_K - q C$ menghasilkan seluruh permukaan $(K,T)$ dalam satu sapuan dari $C(K,0)=(S_0-K)_+$.
6. **Transformasi Grid Hiperbolik `sinh` & Peredaman Rannacher (1984)**: Pemetaan `sinh` memusatkan simpul di sekitar kink $K$, sementara 4 half-step implisit penuh di awal meredam osilasi frekuensi tinggi Crank–Nicolson akibat ketidakhalusan turunan pertama payoff.
7. **Koreksi Jembatan Brownian untuk Opsi Barrier Diskrit**: Peluang kondisional bahwa proses difusi lokal beku yang berada di $X_k, X_{k+1} > \ln H$ telah menembus $\ln H$ di dalam interval $\Delta t$ adalah $\exp(-2(X_k-\ln H)(X_{k+1}-\ln H)/(\sigma_k^2\Delta t))$ (`barrier.py`).
8. **Rumus Tertutup Barrier Reiner–Rubinstein (1991) & Prinsip Refleksi**: Mengkombinasikan suku Black–Scholes langsung dan bayangan terpantul $(H/S_0)^{2\mu}$ untuk memenuhi syarat batas menyerap $V(H,t)=0$ dan paritas $C_{\text{DO}} + C_{\text{DI}} = C_{\text{vanilla}}$.
9. **Dekomposisi Sensitivitas Vega-Bucket vs Pergeseran Paralel**: Oleh kelinearan variasi pertama operator PDE terhadap kurva $\theta(t)$ sepotong-linear, jumlah sensitivitas empat fungsi topi $\sum_i \partial C/\partial \theta_i$ sama dengan sensitivitas pergeseran paralel hingga galat $O(\delta)$ (`greeks.py`).
10. **Algoritma Levenberg–Marquardt dengan Jacobian Analitik**: Interpolasi adaptif $(J^T J + \lambda \operatorname{diag}(J^T J))\Delta\beta = -J^T r$ antara Gauss–Newton dan penurunan gradien dengan matriks kovarians asimtotik $(J^T J)^{-1}\hat{\sigma}^2$ (`calibration.py`).

## 10. Kesimpulan dan Batasan
Proyek `day003` (`3a`) membangun dan mengaudit secara menyeluruh koherensi antara permukaan SSVI, rumus Dupire, dua solver PDE independen, simulasi Monte Carlo, dan diagnosis singularitas lapisan batas gamma di $S_0$.

## 11. Referensi dan Jejak Bukti
- Gatheral, J., & Jacquier, A. (2014). *Arbitrage-free SVI volatility surfaces*. Quantitative Finance, 14(1), 59–71.
- Dupire, B. (1994). *Pricing with a smile*. Risk, 7(1), 18–20.
- Reiner, E., & Rubinstein, M. (1991). *Breaking down the barriers*. Risk, 4(8), 28–35.
- Artefak: [days/day003_quant_volatilitas-lokal-dupire-ssvi.ipynb](../../days/day003_quant_volatilitas-lokal-dupire-ssvi.ipynb), [outputs/day003a/results.json](../../outputs/day003a/results.json).
