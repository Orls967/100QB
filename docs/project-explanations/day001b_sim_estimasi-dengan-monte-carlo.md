# Penjelasan Proyek — day001b: Estimasi π dengan Monte Carlo (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001b` (Hari 1, Sesi b)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `002` (ROADMAP) · **Jenis**: `sim` · **Bahasa**: Python
- **Notebook**: [days/day001b_sim_estimasi-dengan-monte-carlo.ipynb](../../days/day001b_sim_estimasi-dengan-monte-carlo.ipynb)
- **Ringkasan**: Mengestimasi konstanta $\pi$ melalui metode *hit-or-miss* Monte Carlo pada persegi $[-1,1]^2$, mengukur laju konvergensi $O(N^{-1/2})$, serta memverifikasi distribusi asimtotik normal melalui 200 replikasi independen.

## 2. Masalah dan Tujuan
Banyak integral multidimensi dalam fisika dan keuangan kuantitatif tidak memiliki antiderivatif tertutup sehingga metode kuadratur grid mengalami *curse of dimensionality*. Metode Monte Carlo mengganti integrasi deterministik dengan pengambilan sampel acak yang laju konvergensinya $O(N^{-1/2})$ tidak bergantung pada dimensi ruang.

Proyek ini menggunakan masalah klasik estimasi luas lingkaran satuan ($\pi$) untuk menguji secara empiris Hukum Bilangan Besar, Teorema Limit Pusat, keakuratan selang kepercayaan 95%, dan stabilitas PRNG `PCG64`. Ruang lingkup dibatasi pada sampling seragam 2D tanpa teknik *quasi-Monte Carlo* (Sobol/Halton) maupun *variance reduction*.

## 3. Apa yang Dikerjakan Proyek
Notebook `day001b` membangkitkan $N=20.000$ titik acak $(x_i, y_i) \sim U([-1,1]^2)$ menggunakan `qb100.simulation.get_rng(42)`, menghitung proporsi titik di dalam lingkaran $x_i^2+y_i^2 \le 1$, mengestimasi $\hat{\pi} = 4\hat{p}$ beserta standard error dan CI 95% via `monte_carlo_mean` / `mc_summary`, memetakan konvergensi pada $N \in [100, 100.000]$ serta 200 replikasi independen, dan menyimpan 4 grafik ke `outputs/figures/`.

## 4. Uraian Kode
- [src/qb100/simulation.py](../../src/qb100/simulation.py) — fungsi `get_rng`: membungkus `np.random.default_rng(seed)` berbasis algoritma `PCG64` untuk reprodusibilitas.
- [src/qb100/simulation.py](../../src/qb100/simulation.py) — fungsi `mc_summary` & `monte_carlo`: menghitung mean sampel, standar deviasi (`ddof=1`), standard error $s/\sqrt{N}$, dan selang kepercayaan 95% ($\pm 1.96\,\text{SE}$).
- [src/qb100/simulation.py](../../src/qb100/simulation.py) — fungsi `running_mean` & `running_se`: menghitung rata-rata berjalan dan standard error berjalan untuk memvisualisasikan pita kepercayaan sepanjang bertambahnya sampel.
- [days/day001b_sim_estimasi-dengan-monte-carlo.ipynb](../../days/day001b_sim_estimasi-dengan-monte-carlo.ipynb) — sel 5–6 (sampling 20.000 titik & hitung `15,732/20,000`), sel 7 (CI 95%), sel 8–9 (kurva error vs $N$ hingga `100,000`), sel 10 (200 trial), dan sel 12–15 (4 grafik).
- [tests/test_simulation.py](../../tests/test_simulation.py) — 17 unit test pytest yang menguji determinisme `get_rng`, keakuratan `mc_summary`, dan sifat `running_mean`.

## 5. Metode dan Alur Kerja
Untuk $U_1, U_2 \stackrel{\text{iid}}{\sim} U(-1, 1)$, didefinisikan peubah acak $I_i = 4 \cdot \mathbf{1}\{U_{1,i}^2 + U_{2,i}^2 \le 1\}$. Karena luas lingkaran satuan adalah $\pi$ dan luas persegi $[-1,1]^2$ adalah $4$, maka $E[I_i] = \pi$ dan $\mathrm{Var}(I_i) = 16 \frac{\pi}{4}(1 - \frac{\pi}{4}) \approx 2.6968$ ($\sigma_I \approx 1.6422$).

```mermaid
flowchart LR
    A["Seed=42 -> get_rng (PCG64)"] --> B["Sampel U(-1,1)^2 (N=20.000)"]
    B --> C["Indikator I = 4 * 1{x^2+y^2 <= 1}"]
    C --> D["mc_summary & running_mean"]
    D --> E["Uji N=100..100k & 200 Trial"]
    E --> F["4 Grafik day001b_pi-*.png & pytest"]
```

## 6. Output yang Ditargetkan
Target keberhasilan: (1) notebook lulus eksekusi penuh `scripts/run_all.py`; (2) estimasi $\hat{\pi}$ pada $N=20.000$ memuat $\pi = 3.141593$ di dalam selang kepercayaan 95%; (3) galat absolut menyusut sejalan dengan $O(N^{-1/2})$ hingga $N=100.000$; (4) 4 grafik PNG tersimpan dan seluruh 17 tes di `tests/test_simulation.py` lulus.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Eksekusi penuh notebook `day001b` | `tercapai` | `python scripts/run_all.py --day 1 --inplace` (sesi pemeliharaan `2026-11-21`): `[OK] days/day001b_sim_estimasi-dengan-monte-carlo.ipynb (1.2s)` |
| Estimasi $\hat{\pi}$ ($N=20.000$) dan cakupan CI 95% | `tercapai` | Output sel 6–7: `15,732 / 20,000` di dalam lingkaran, $\hat{\pi} = 3.146400 \pm 0.011589$, CI 95% `[3.123686, 3.169114]` memuat $\pi = 3.141593$ (galat `0.004807`) |
| Konvergensi $N \in [100, 100.000]$ & 200 replikasi independen | `tercapai` | Output sel 9–10: pada $N=100.000$, $\hat{\pi} = 3.14256$ (galat `0.00096`); rata-rata 200 trial = `3.14298` (std `0.02403`) |
| Penyimpanan 4 berkas grafik PNG (nama terverifikasi di `outputs/figures/`) | `tercapai` | Tersedia: [day001b_pi-scatter.png](../../outputs/figures/day001b_pi-scatter.png), [day001b_pi-convergence.png](../../outputs/figures/day001b_pi-convergence.png), [day001b_pi-error-scaling.png](../../outputs/figures/day001b_pi-error-scaling.png), [day001b_pi-histogram.png](../../outputs/figures/day001b_pi-histogram.png) |
| Unit test `tests/test_simulation.py` | `tercapai` | `pytest -q` (sesi pemeliharaan `2026-11-21`): 17/17 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Tidak ada kriteria yang gagal atau belum diketahui. **FAKTA**: Eksekusi notebook `day001b` dijalankan langsung di sesi pemeliharaan `2026-11-21` (`1.2s`), dan nama keempat grafik yang sebelumnya salah di dokumen (`pi-error-loglog` / `pi-trials-hist`) telah diperbaiki ke nama aktual di `outputs/figures/` (`pi-error-scaling` / `pi-histogram`).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Integrasi Monte Carlo Hit-or-Miss**: Metode *hit-or-miss* memandang luas daerah $\Omega \subset [-1,1]^2$ sebagai ekspektasi fungsi indikator $\int_{[-1,1]^2} \mathbf{1}_\Omega(x,y)\,dx\,dy = 4 E[\mathbf{1}_\Omega(U_1,U_2)]$. Pendekatan ini berlaku untuk himpunan terukur apa pun tanpa memerlukan batas analitik yang halus. Pada sel 6 notebook, syarat $x^2+y^2 \le 1$ menghasilkan `15,732` titik dari `20,000` sehingga $\hat{\pi} = 3.146400$.
2. **Distribusi Bernoulli dan Variansi Teoretis Estimator**: Variabel indikator $B_i = \mathbf{1}\{U_{1,i}^2+U_{2,i}^2 \le 1\}$ berdistribusi $\text{Bernoulli}(p)$ dengan $p = \pi/4 \approx 0.785398$. Untuk $I_i = 4 B_i$, simpangan baku populasinya adalah $\sigma_I = 4\sqrt{p(1-p)} \approx 1.6422$. Nilai empiris yang diukur di sel 7 notebook (`Std = 1.638872`) cocok sangat dekat dengan nilai teoretis tersebut, mengonfirmasi kebenaran fungsi `mc_summary`.
3. **Hukum Bilangan Besar Kuat (SLLN)**: Karena barisan $I_1, I_2, \dots$ saling bebas dan terdistribusi identik dengan ekspektasi terhingga $E[|I_1|] = \pi < \infty$, Hukum Bilangan Besar Kuat Kolmogorov menjamin rata-rata berjalan $\bar{I}_N \to \pi$ hampir pasti saat $N \to \infty$. Sel 13 memplot keluaran `running_mean` dan `running_se` dari `src/qb100/simulation.py` pada grafik `day001b_pi-convergence.png`.
4. **Teorema Limit Pusat (CLT)**: Teorema Limit Pusat Lindeberg–Lévy menyatakan bahwa galat terstandarisasi $\sqrt{N}(\bar{I}_N - \pi)/\sigma_I$ konvergen dalam distribusi menuju $\mathcal{N}(0, 1)$. Akibatnya, distribusi sampel $\hat{\pi}$ dari 200 replikasi independen membentuk kurva lonceng Gaussian berpusat di $\pi$ dengan simpangan baku $\approx 1.6422/\sqrt{N}$. Sel 10 dan sel 15 memvisualisasikan distribusi empiris 200 trial ini pada `day001b_pi-histogram.png`.
5. **Laju Konvergensi Universal O(N^-1/2)**: Standard error estimator Monte Carlo selalu menyusut sebagai $\text{SE} = \sigma_I / \sqrt{N}$, terlepas dari dimensi integrasi $d$. Konsekuensinya, mengurangi galat sebesar satu digit desimal ($10\times$ lebih akurat) menuntut $100\times$ lebih banyak sampel acak. Kurva log-log pada sel 14 (`day001b_pi-error-scaling.png`) membandingkan peluruhan galat empiris terhadap garis referensi $N^{-1/2}$.
6. **Selang Kepercayaan Asimtotik Normal 95%**: Berdasarkan CLT, interval $[\hat{\pi} - z_{0.975}\,\text{SE}, \hat{\pi} + z_{0.975}\,\text{SE}]$ dengan $z_{0.975} \approx 1.96$ memiliki peluang cakupan asimtotik $95\%$. Pada sel 7 notebook, fungsi `mc_summary` di `src/qb100/simulation.py` menghasilkan CI `[3.123686, 3.169114]` yang secara tepat mencakup nilai sejati $\pi = 3.141593$.
7. **Pembangkit Bilangan Acak Semu PCG64**: Kualitas simulasi Monte Carlo bergantung pada ketiadaan korelasi kisi (*lattice structure*) pada titik acak berdimensi banyak. Fungsi `get_rng` di `src/qb100/simulation.py` memakai `PCG64` (*Permuted Congruential Generator* 64-bit, O'Neill 2014) yang menjadi bawaan NumPy modern dengan periode $2^{128}$ dan lulus uji statistik TestU01 BigCrush.
8. **Ketakbiasan Estimator Sampel (Unbiasedness)**: Estimator $\hat{\pi}_N = \frac{1}{N}\sum_{i=1}^N I_i$ memenuhi $E[\hat{\pi}_N] = \pi$ untuk setiap $N \ge 1$, artinya tidak ada galat sistematis bahkan pada $N$ kecil. Uji 200 replikasi di sel 10 menghasilkan rata-rata `3.14298` yang selisihnya (`0.00139`) berada dalam $0.82$ standard error rata-rata ($0.02403/\sqrt{200} \approx 0.00170$).
9. **Invariansi Simetri Kuadran Pertama**: Karena persamaan lingkaran $x^2 + y^2 \le 1$ simetris terhadap sumbu-$x$ dan sumbu-$y$, pengambilan sampel di kuadran pertama $[0,1]^2$ secara matematis ekuivalen dengan sampling di $[-1,1]^2$. Namun menggunakan $[-1,1]^2$ memudahkan visualisasi *scatter plot* seluruh lingkaran penuh pada `day001b_pi-scatter.png` (sel 12).
10. **Diagnostik Regresi Log-Log Galat Empiris**: Pada grafik `day001b_pi-error-scaling.png` (sel 14), memplot $\log_{10}|\hat{\pi}_N - \pi|$ terhadap $\log_{10} N$ mengubah hukum pangkat $c N^{-1/2}$ menjadi garis lurus berkemiringan $-0.5$. Fluktuasi di sekitar garis $-0.5$ wajar untuk satu realisasi sampel dan merata saat diulang melintasi banyak benih acak.

## 10. Kesimpulan dan Batasan
Proyek `day001b` memverifikasi secara empiris laju konvergensi $O(N^{-1/2})$ dan sifat asimtotik normal estimator Monte Carlo, sekaligus menyediakan utilitas `qb100.simulation` untuk proyek simulasi lanjutan. Keterbatasannya adalah efisiensi lambat $O(N^{-1/2})$ tanpa metode reduksi variansi atau *low-discrepancy sequence*.

## 11. Referensi dan Jejak Bukti
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*. Springer.
- O'Neill, M. E. (2014). *PCG: A family of simple fast space-efficient statistically good algorithms for random number generation*. Harvey Mudd College Technical Report.
- Artefak: [day001b_sim_estimasi-dengan-monte-carlo.ipynb](../../days/day001b_sim_estimasi-dengan-monte-carlo.ipynb), [src/qb100/simulation.py](../../src/qb100/simulation.py).
