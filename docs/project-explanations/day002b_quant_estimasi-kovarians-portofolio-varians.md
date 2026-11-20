# Penjelasan Proyek — day002b: Estimasi Kovarians & Portofolio Varians Minimum (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002b` (Hari 2, Sesi b)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb](../../days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb)
- **Ringkasan**: Memadukan teori matriks acak (distribusi Wishart & hukum Marchenko–Pastur) dengan penyusutan kovarians Ledoit–Wolf (2004) dan optimisasi portofolio varians minimum (*unconstrained* & *long-only* projected gradient).

## 2. Masalah dan Tujuan
Ketika dimensi aset $N=50$ sebanding dengan panjang sampel $T=200$ ($q = N/T = 0.25$), matriks kovarians sampel $S$ mengalami distorsi spektral sehingga portofolio varians minimum berbasis $S^{-1}$ menderita *optimism bias* di dalam sampel $(1-q)\sigma_{\min}^2$ dan degradasi risiko di luar sampel $(1-q)^{-1}\sigma_{\min}^2$.

## 3. Apa yang Dikerjakan Proyek
Modul `src/qb100/covariance.py` dan notebook `day002b` memverifikasi kepadatan spektral Marchenko–Pastur pada matriks korelasi sampel, membuktikan rumus eksak varians *in-sample* dan *out-of-sample* dari matriks Wishart lewat simulasi Monte Carlo $M=5.000$, mengimplementasikan estimator penyusutan analitik Ledoit–Wolf, serta memetakan frontier efisien *long-only*.

## 4. Uraian Kode
- [src/qb100/covariance.py](../../src/qb100/covariance.py): fungsi `marchenko_pastur_pdf`, `marchenko_pastur_bounds`, `generate_factor_covariance`, `ledoit_wolf_shrinkage`, `gmv_weights_analytic`, `gmv_weights_long_only`, `project_simplex`, dan `efficient_frontier_long_only`.
- [days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb](../../days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb): 12 sel notebook eksperimen dan visualisasi.
- [tests/test_covariance.py](../../tests/test_covariance.py): 5 unit test pytest untuk Marchenko–Pastur, proyeksi simpleks, Ledoit–Wolf, dan GMV.

## 5. Metode dan Alur Kerja
Bangkitkan matriks kovarians faktor sejati $\Sigma$ ($N=50, K=3$, kondisi `55.18`, $\sigma_{\min}^2 = 0.00067839$) $\to$ simulasi $M=5.000$ matriks sampel $T=200$ $\to$ bandingkan varians GMV in-sample/out-of-sample untuk Sampel, Ledoit–Wolf, dan Oracle $\to$ hitung frontier *long-only* via *projected gradient descent* di atas simpleks.

## 6. Output yang Ditargetkan
Spektrum eigen sesuai batas Marchenko–Pastur $[\lambda_-, \lambda_+]$; rasio varians GMV sampel cocok dengan faktor teoretis $(1-N/T)$ dan $1/(1-N/T)$; risiko *out-of-sample* Ledoit–Wolf secara signifikan lebih rendah dari kovarians sampel; seluruh tes `test_covariance.py` lulus.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Hukum Marchenko–Pastur vs spektrum empiris | `tercapai` | Spektrum matriks putih terpusat di dalam batas $[\lambda_-, \lambda_+] = [(1-\sqrt{q})^2, (1+\sqrt{q})^2]$ dan memisahkan 3 eigenvalue faktor |
| Verifikasi hukum Wishart in-sample $(1-q)$ & out-of-sample $(1-q)^{-1}$ | `tercapai` | Output sel 3 notebook ($N=50, T=200$): varians in-sample empiris `0.007520` vs teori `0.007576`; out-of-sample Ledoit–Wolf mengungguli kovarians sampel |
| Penyusutan Ledoit–Wolf & Frontier Long-Only | `tercapai` | Fungsi `ledoit_wolf_shrinkage` dan `efficient_frontier_long_only` terverifikasi di notebook dan unit test |
| Unit test `tests/test_covariance.py` | `tercapai` | `pytest -q` dijalankan pada sesi `day003c`: 5/5 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `pytest -q` serta output sel notebook.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Distribusi Wishart $W_N(\Sigma, T)$**: Matriks scatter $T \hat{S} = X^T X$ dari sampel Gaussian $X \sim N(0, \Sigma)$ berdistribusi Wishart dengan $T$ derajat kebebasan.
2. **Hukum Spektral Marchenko–Pastur (1967)**: Dalam limit $N, T \to \infty$ dengan $N/T \to q \in (0,1)$, distribusi eigenvalue matriks kovarians sampel $\Sigma = I$ konvergen ke $\rho_{\text{MP}}(\lambda) = \frac{\sqrt{(\lambda_+ - \lambda)(\lambda - \lambda_-)}}{2\pi q \lambda}$ dengan $\lambda_{\pm} = (1 \pm \sqrt{q})^2$ (`marchenko_pastur_pdf`).
3. **Ekspektasi Invers Wishart & Optimism Bias**: $E[\hat{S}^{-1}] = \frac{T}{T-N-1}\Sigma^{-1}$ menyebabkan varians in-sample portofolio GMV menyusut semu sebesar faktor $(1 - N/T)$, sedangkan varians out-of-sample membengkak sebesar $(1 - N/T)^{-1}$.
4. **Model Kovarians Faktor Terstruktur**: $\Sigma = B B^T + D$ memisahkan variasi pasar/sektor berdimensi rendah ($K=3$ spike eigenvalues) dari risiko idiosinkratik diagonal (`generate_factor_covariance`).
5. **Penyusutan Linear Ledoit–Wolf (2004)**: $\hat{\Sigma}_{\text{LW}} = (1 - \hat{\delta}) \hat{S} + \hat{\delta} \mu I$ dengan $\mu = \operatorname{tr}(\hat{S})/N$ dan intensitas optimal $\hat{\delta}^* \in [0,1]$ yang meminimalkan jarak Frobenius $E[\|\hat{\Sigma}_{\text{LW}} - \Sigma\|_F^2]$ secara asimtotik (`ledoit_wolf_shrinkage`).
6. **Regularisasi Bilangan Kondisi Spektral**: Penambahan $\delta \mu I$ mengangkat eigenvalue terkecil $\lambda_{\min}(\hat{S})$ menjauhi nol, menstabilkan inversi matriks pada optimisasi GMV.
7. **Portofolio Global Minimum Variance (GMV) Analitik**: Minimasi $\frac{1}{2}w^T \Sigma w$ dengan $\mathbf{1}^T w = 1$ memiliki solusi tertutup $w^* = \frac{\Sigma^{-1}\mathbf{1}}{\mathbf{1}^T \Sigma^{-1}\mathbf{1}}$ dan varians $\sigma_{\min}^2 = (\mathbf{1}^T \Sigma^{-1}\mathbf{1})^{-1}$ (`gmv_weights_analytic`).
8. **Proyeksi Euclidean ke Simpleks Standar**: Algoritma $O(N \log N)$ Duchi et al. (2008) mencari ambang $\theta$ sehingga $w_i = \max(v_i - \theta, 0)$ memenuhi $\sum w_i = 1, w_i \ge 0$ (`project_simplex`).
9. **Projected Gradient Descent untuk GMV Long-Only**: Iterasi $w^{(k+1)} = \Pi_{\Delta_N}(w^{(k)} - \eta \Sigma w^{(k)})$ menyelesaikan QP konveks berkendala non-negatif tanpa solver eksternal (`gmv_weights_long_only`).
10. **Efek Regularisasi Implisit Kendala Long-Only (Jagannathan–Ma 2003)**: Kendala $w_i \ge 0$ bertindak setara dengan penyusutan kovarians karena memotong posisi short ekstrem yang dipicu oleh galat estimasi eigenvalue kecil.

## 10. Kesimpulan dan Batasan
Proyek `day002b` mendemonstrasikan secara kuantitatif mengapa estimasi kovarians sampel gagal pada dimensi menengah-tinggi dan bagaimana penyusutan Ledoit–Wolf serta kendala *long-only* memulihkan stabilitas portofolio.

## 11. Referensi dan Jejak Bukti
- Ledoit, O., & Wolf, M. (2004). *A well-conditioned estimator for large-dimensional covariance matrices*. Journal of Multivariate Analysis, 88(2), 365–411.
- Artefak: [day002b_quant_estimasi-kovarians-portofolio-varians.ipynb](../../days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb), [src/qb100/covariance.py](../../src/qb100/covariance.py).
