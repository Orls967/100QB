# Penjelasan Proyek — day001b: Estimasi π dengan Monte Carlo (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001b` (Hari 1, Sesi b)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `002` (ROADMAP) · **Jenis**: `sim` · **Bahasa**: Python
- **Notebook**: [days/day001b_sim_estimasi-dengan-monte-carlo.ipynb](../../days/day001b_sim_estimasi-dengan-monte-carlo.ipynb)
- **Ringkasan**: Mengestimasi konstanta $\pi$ melalui pengambilan sampel acak seragam di persegi $[-1,1]^2$ dan menganalisis laju konvergensi $O(N^{-1/2})$ beserta selang kepercayaan 95%.

## 2. Masalah dan Tujuan
Bagaimana cara mengaproksimasi luas daerah lengkung (lingkaran berjari-jari 1) menggunakan simulasi probabilistik, dan seberapa cepat galat estimasi Monte Carlo menyusut ketika jumlah sampel $N$ diperbesar?

## 3. Apa yang Dikerjakan Proyek
Notebook membangkitkan titik acak $(X_i, Y_i) \sim U([-1, 1]^2)$ menggunakan `qb100.simulation`, menghitung proporsi titik di dalam lingkaran $X_i^2 + Y_i^2 \le 1$, mengestimasi $\hat{\pi} = 4 \hat{p}$, mengevaluasi konvergensi hingga $N=100.000$ dan distribusi 200 replikasi independen, serta menyimpan 4 grafik ke `outputs/figures/`.

## 4. Uraian Kode
- [src/qb100/simulation.py](../../src/qb100/simulation.py): fungsi `get_rng` (PCG64), `monte_carlo_mean` (estimasi mean, SE, CI 95%), dan `running_mean`.
- [days/day001b_sim_estimasi-dengan-monte-carlo.ipynb](../../days/day001b_sim_estimasi-dengan-monte-carlo.ipynb): 18 sel notebook yang memuat simulasi geometri, tabel konvergensi, dan 4 visualisasi.
- [tests/test_simulation.py](../../tests/test_simulation.py): 17 unit test pytest untuk RNG, `monte_carlo_mean`, dan `running_mean`.

## 5. Metode dan Alur Kerja
Bangkitkan $(x_i, y_i) \sim U(-1,1)^2$ dengan seed tetap `42` $\to$ indikator $I_i = 4 \cdot \mathbf{1}\{x_i^2+y_i^2 \le 1\}$ $\to$ hitung $\hat{\pi}_N = \frac{1}{N}\sum_{i=1}^N I_i$ dan $\text{SE} = s_I / \sqrt{N}$ $\to$ analisis variasi $N \in [100, 100.000]$ dan 200 trial independen.

## 6. Output yang Ditargetkan
Estimasi $\hat{\pi}$ berada dalam selang kepercayaan 95% dari $\pi = 3.141593$, galat menyusut mengikuti $1/\sqrt{N}$, dan 4 grafik tersimpan di `outputs/figures/`.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Estimasi $\pi$ pada $N=20.000$ dalam CI 95% | `tercapai` | Output sel 6–7: `15,732 / 20,000` titik di dalam lingkaran, $\hat{\pi} = 3.146400 \pm 0.011589$, CI 95% `[3.123686, 3.169114]` memuat $\pi = 3.141593$ (galat `0.004807`) |
| Tabel konvergensi hingga $N=100.000$ & 200 trial | `tercapai` | Output sel 9–10: pada $N=100.000$, $\hat{\pi} = 3.14256$ (galat `0.00096`); rata-rata 200 trial = `3.14298` (std `0.02403`) |
| Penyimpanan 4 grafik visualisasi | `tercapai` | Tersimpan: [day001b_pi-scatter.png](../../outputs/figures/day001b_pi-scatter.png), [day001b_pi-convergence.png](../../outputs/figures/day001b_pi-convergence.png), [day001b_pi-error-loglog.png](../../outputs/figures/day001b_pi-error-loglog.png), [day001b_pi-trials-hist.png](../../outputs/figures/day001b_pi-trials-hist.png) |
| Unit test `tests/test_simulation.py` | `tercapai` | `pytest -q` dijalankan di sesi ini: 17/17 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria tercapai. **FAKTA**: `pytest -q` dijalankan pada sesi ini dan seluruh tes `test_simulation.py` lulus.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Integrasi Monte Carlo via Hit-or-Miss**: Luas daerah $\Omega \subset [-1,1]^2$ sama dengan luas kotak ($4$) dikali peluang $P((X,Y)\in\Omega) = \pi/4$; digunakan di sel 6 notebook.
2. **Percobaan Bernoulli dan Distribusi Binomial**: Setiap titik adalah variabel Bernoulli dengan $p = \pi/4 \approx 0.7854$ dan variansi $p(1-p)$; untuk $I_i = 4 B_i$, $\mathrm{Var}(I_i) = 16 p(1-p) \approx 2.6968$ ($\sigma \approx 1.6422$, cocok dengan empiris `1.6389` di sel 7).
3. **Hukum Bilangan Besar Kuat (SLLN)**: Rata-rata sampel $\hat{\pi}_N \to E[I] = \pi$ hampir pasti saat $N \to \infty$; divisualisasikan dengan `running_mean` pada [day001b_pi-convergence.png](../../outputs/figures/day001b_pi-convergence.png).
4. **Teorema Limit Pusat (CLT)**: Distribusi $\sqrt{N}(\hat{\pi}_N - \pi) \Rightarrow N(0, 16p(1-p))$; diverifikasi dengan histogram 200 replikasi di [day001b_pi-trials-hist.png](../../outputs/figures/day001b_pi-trials-hist.png).
5. **Laju Konvergensi $O(N^{-1/2})$**: Standard error menyusut sebanding dengan $N^{-1/2}$ (kemiringan $-0.5$ pada skala log-log), membutuhkan $100\times$ sampel untuk $10\times$ akurasi ekstra.
6. **Selang Kepercayaan Asimtotik Normal**: Interval $\hat{\pi} \pm z_{1-\alpha/2} \text{SE}$ dengan $z_{0.975} \approx 1.96$ diimplementasikan pada `monte_carlo_mean` di `src/qb100/simulation.py`.
7. **Pembangkit Bilangan Acak Semu PCG64**: Generator `np.random.default_rng(seed)` berbasis *Permuted Congruential Generator* 64-bit di `get_rng` untuk periodisitas panjang dan reprodusibilitas.
8. **Estimator Tak-Bias (Unbiasedness)**: $E[\hat{\pi}_N] = \pi$ untuk setiap $N \ge 1$; terlihat pada rata-rata 200 trial `3.14298` yang sangat dekat dengan $\pi$.
9. **Simetri Kuadran Pertama**: Karena $x^2+y^2$ genap pada kedua koordinat, sampling di $[0,1]^2$ dan $[-1,1]^2$ memiliki distribusi indikator yang identik.
10. **Regresi Log-Log Galat Empiris**: Pemetaan $\log |\hat{\pi}_N - \pi|$ terhadap $\log N$ pada [day001b_pi-error-loglog.png](../../outputs/figures/day001b_pi-error-loglog.png) untuk memeriksa kemiringan $-1/2$.

## 10. Kesimpulan dan Batasan
Proyek `day001b` memvalidasi properti dasar estimator Monte Carlo dan membangun fungsi `monte_carlo_mean` yang dipakai ulang pada sesi-sesi berikutnya.

## 11. Referensi dan Jejak Bukti
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*. Springer.
- Artefak: [day001b_sim_estimasi-dengan-monte-carlo.ipynb](../../days/day001b_sim_estimasi-dengan-monte-carlo.ipynb), [src/qb100/simulation.py](../../src/qb100/simulation.py).
