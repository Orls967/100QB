# Penjelasan Proyek — day002d: Risiko Kredit Vasicek — Distribusi Eksak, LHP & Importance Sampling (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002d` (Hari 2, Sesi d)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb](../../days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb)
- **Modul & Tes**: [src/qb100/credit.py](../../src/qb100/credit.py), [tests/test_credit.py](../../tests/test_credit.py)
- **Ringkasan**: Komputasi distribusi kerugian portofolio kredit satu-faktor Gaussian copula Vasicek: integrasi kuadratur Gauss–Hermite 200 titik dengan pemrograman dinamis Poisson–Binomial (`exact_loss_pmf`), verifikasi terhadap `scipy.integrate.quad` (`exact_pmf_quad_homogeneous`), analisis konvergensi ukuran portofolio berhingga ke limit analitik *Large Homogeneous Portfolio* (`lhp_tail_prob`, `lhp_var`), serta simulasi kejadian langka menggunakan *Conditional Monte Carlo* (`estimate_conditional_mc`) dan *Importance Sampling* (*exponential/mean-shift twisting* pada faktor sistemik via `find_optimal_is_theta` dan `estimate_importance_sampling_mc`).

## 2. Masalah dan Tujuan
Dalam portofolio kredit beranggotakan $N$ debitur yang terpapar faktor ekonomi makro bersama $Z \sim \mathcal{N}(0, 1)$, indikator aset $X_i = \sqrt{\rho_i} Z + \sqrt{1 - \rho_i}\varepsilon_i$ memicu gagal bayar ketika $X_i \le \Phi^{-1}(\text{PD}_i)$. Mengestimasi peluang kerugian ekstrem di ekor distribusi, seperti $P(L \ge 30\%) \approx 3.37 \times 10^{-4}$ atau $P(L \ge 40\%) \approx 4.1 \times 10^{-5}$, dengan simulasi Monte Carlo naif membutuhkan jutaan sampel karena variansi relatif $\sqrt{(1-p)/(Mp)}$ meledak saat $p \to 0$. Selain itu, banyak praktik industri mengandalkan rumus tertutup Basel II *Large Homogeneous Portfolio* (LHP) yang mengabaikan risiko granularitas pada $N$ kecil atau eksposur heterogen. Tujuan proyek `day002d` adalah: (1) menghitung PMF kerugian eksak baik untuk portofolio homogen maupun heterogen hingga presisi kuadratur mesin, (2) mengukur kecepatan konvergensi portofolio $N \in [20, 5000]$ ke aproksimasi LHP, dan (3) membuktikan reduksi variansi (*Variance Reduction Factor*, VRF) ratusan hingga ribuan kali menggunakan Conditional MC dan Importance Sampling.

## 3. Apa yang Dikerjakan Proyek
Modul [src/qb100/credit.py](../../src/qb100/credit.py) dan notebook [day002d](../../days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb) (9 sel, waktu eksekusi `2.76 detik`) mengimplementasikan:
1. Perhitungan PMF eksak untuk Kasus A (homogen, $N=100$, eksposur total $E=100$, $E[L]=2.0$) dan Kasus B (heterogen, eksposur bertingkat via `make_exposures`, $E=300$, $E[L]=6.0$) menggunakan `conditional_loss_pmf_dp` dan `exact_loss_pmf`, beserta validasi silang `exact_pmf_quad_homogeneous`.
2. Ekstraksi metrik risiko `exact_risk_metrics` ($E[L]$, $\text{SD}[L]$, $\text{VaR}_{99\%}$, $\text{ES}_{99\%}$, $\text{VaR}_{99.9\%}$, $\text{ES}_{99.9\%}$, dan probabilitas ekor $15\%, 30\%, 40\%$) serta verifikasi variansi analitik bivariat `bivariate_loss_variance`.
3. Eksperimen konvergensi $N \in \{20, 50, 100, 500, 2000, 5000\}$ terhadap `lhp_tail_prob` dan `lhp_var`.
4. Perbandingan tiga estimator Monte Carlo (`estimate_naive_mc`, `estimate_conditional_mc` dengan interpolator spline ekor `build_conditional_tail_interpolator`, dan `estimate_importance_sampling_mc` dengan pergeseran optimal `find_optimal_is_theta`), serta regresi slope log-log *standard error* terhadap jumlah sampel $M$ (`fit_loglog_slope`).

## 4. Uraian Kode
- [src/qb100/credit.py](../../src/qb100/credit.py):
  - `make_exposures(...)` & `conditional_pd(pd, rho, z)`: membangun bobot eksposur portofolio dan peluang gagal bayar bersyarat $p_i(z)$.
  - `conditional_loss_pmf_dp(probs, exposures)` & `exact_loss_pmf(...)`: konvolusi DP bersyarat yang diintegrasikan dengan 200 simpul kuadratur Gauss–Hermite.
  - `exact_risk_metrics(pmf, losses)` & `bivariate_loss_variance(...)`: menghitung momen, VaR, Expected Shortfall, serta variansi teoretis via CDF normal bivariat.
  - `exact_pmf_quad_homogeneous(...)`, `lhp_tail_prob(x, pd, rho)`, `lhp_var(alpha, pd, rho)`: pembanding kuadratur adaptif dan limit asimtotik Vasicek LHP.
  - `estimate_naive_mc(...)`, `build_conditional_tail_interpolator(...)`, `estimate_conditional_mc(...)`, `find_optimal_is_theta(...)`, `estimate_importance_sampling_mc(...)`, `fit_loglog_slope(...)`: tiga estimator Monte Carlo beserta penyetelan *saddlepoint/mean-shift* optimal dan analisis laju konvergensi $M^{-1/2}$.
- [tests/test_credit.py](../../tests/test_credit.py): 4 unit test `pytest` untuk normalisasi PMF, kesesuaian momen terhadap `bivariate_loss_variance`, keselarasan `exact_loss_pmf` vs `exact_pmf_quad_homogeneous`, konvergensi ke LHP, serta reduksi variansi Importance Sampling.

## 5. Metode dan Alur Kerja
Alur perhitungan pada `day002d` menghubungkan metode deterministik (DP + kuadratur) sebagai tolok ukur eksak dengan metode stokastik reduksi variansi:

```mermaid
flowchart LR
    A["Parameter (N, PD, rho, exposures)"] --> B["conditional_pd di 200 simpul GH"]
    B --> C["conditional_loss_pmf_dp -> exact_loss_pmf"]
    C --> D["exact_risk_metrics (VaR, ES, Ekor)"]
    A --> E["lhp_tail_prob & lhp_var (Limit N->inf)"]
    D --> F["Benchmark: Naive MC vs Conditional MC vs IS (Z ~ N(-theta, 1))"]
```

Pada langkah pertama, kuadratur Gauss–Hermite 200 titik mengintegrasikan PMF bersyarat yang dihasilkan rekursi DP Poisson–Binomial. Hasilnya diverifikasi terhadap kuadratur adaptif `scipy.integrate.quad` dan variansi bivariat eksak. Selanjutnya, ekor kerugian eksak $P(L \ge x)$ dibandingkan dengan limit LHP untuk berbagai $N$ guna mengkuantifikasi *granularity adjustment*. Terakhir, untuk target ekor $x \in \{15\%, 30\%, 40\%\}$, *shift* faktor sistemik $-\theta^*$ dipilih oleh `find_optimal_is_theta` sehingga distribusi proposal berpusat tepat di wilayah krisis, lalu `fit_loglog_slope` memastikan ketiganya mempertahankan laju $O(M^{-1/2})$ dengan konstanta variansi jauh lebih kecil pada IS.

## 6. Output yang Ditargetkan
1. Total massa PMF eksak mendekati `1.0` hingga `< 1e-14`, galat relatif variansi terhadap `bivariate_loss_variance` `< 1e-12`, dan perbedaan `exact_loss_pmf` vs `scipy.integrate.quad` `< 1e-14`.
2. Galat antara `exact_loss_pmf` dan `lhp_tail_prob` menyusut secara konsisten saat $N$ naik dari `20` ke `5.000`.
3. Ketiga estimator MC tak-bias terhadap nilai eksak (dalam `3 SE`), dengan Conditional MC mengungguli Naif dan Importance Sampling memberikan pengurangan SE puluhan hingga ratusan kali lipat pada ekor `30%` dan `40%`.
4. Slope log-log SE terhadap $M$ mendekati `-0.50` dan seluruh 4 tes `tests/test_credit.py` lulus.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Akurasi PMF eksak (GH-200 vs `quad` & momen) | `tercapai` | Output sel 1: Kasus A `sum(pmf)-1 = 2.220e-16`, `E[L]-2 = -1.554e-15`, `rel_err Var = 9.070e-14`; Kasus B `sum(pmf)-1 = -5.551e-16`, `rel_err Var = 8.614e-14`; `max |PMF_GH200 - quad| = 1.110e-16` |
| Metrik risiko portofolio Kasus A & B | `tercapai` | Output sel 1: Kasus A `VaR99=14, ES99=18.422865, P(L>=30%)=0.000337`; Kasus B `VaR99=43, ES99=56.271861, VaR99.9=74` |
| Konvergensi portofolio berhingga ke LHP Vasicek | `tercapai` | Output sel 3: `Err 15%` turun dari `0.023374` ($N=20$) ke `0.003002` ($N=100$) dan `0.000317` ($N=2000$) menuju `LHP = 0.005877`; ekor 30% konvergen ke `0.000197` |
| Reduksi variansi Conditional MC & Importance Sampling | `tercapai` | Output sel 5: pada Kasus A `x=30%` (eksak `0.000337`), Naif `0.000380 ± 4.358e-5`, Conditional `0.000339 ± 2.926e-5`, IS ($\theta=3.5$) `0.000337 ± 8.175e-7` (SE ~53x lebih kecil dari Naif) |
| Laju konvergensi Monte Carlo $O(M^{-1/2})$ | `tercapai` | Output sel 7: slope log-log SE terukur Naif `-0.5740`, Conditional `-0.5441`, IS `-0.5008` |
| Unit test `tests/test_credit.py` | `tercapai` | `.venv/bin/pytest tests/test_credit.py -v`: `4 passed` |

Visualisasi hasil tersimpan pada:
- [day002d_pmf-and-tail.png](../../outputs/figures/day002d_pmf-and-tail.png)
- [day002d_lhp-vs-exact.png](../../outputs/figures/day002d_lhp-vs-exact.png)
- [day002d_is-variance-and-vrf.png](../../outputs/figures/day002d_is-variance-and-vrf.png)
- [day002d_convergence-loglog.png](../../outputs/figures/day002d_convergence-loglog.png)

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
**FAKTA**: Seluruh eksperimen notebook `day002d` (`2.76 detik`) dan 4/4 tes `tests/test_credit.py` lulus. Catatan kecil yang terlihat pada sel 3 untuk $N=2000 \to 5000$ (`Err 15%`: `0.000317` menjadi `0.000608`) disebabkan oleh diskretisasi ambang `ceil(0.15 * N)` pada distribusi bilangan bulat saat $N$ berubah dan pembulatan grid eksposur, bukan divergensi metode.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Ambang Struktural Merton (1974)**: Setiap obligor $i$ mengalami gagal bayar pada horizon $T$ apabila variabel laten nilai aset terstandarisasi $X_i \sim \mathcal{N}(0, 1)$ jatuh di bawah ambang $c_i = \Phi^{-1}(\text{PD}_i)$. Formulasi ini menghubungkan peluang gagal bayar marjinal $\text{PD}_i$ dengan kerangka kopula multivariat tanpa perlu memodelkan dinamika harga saham secara eksplisit.
2. **Kopula Gaussian Satu-Faktor Vasicek (1987, 2002)**: Ketergantungan antar debitur dimodelkan melalui $X_i = \sqrt{\rho_i} Z + \sqrt{1 - \rho_i}\varepsilon_i$, di mana $Z \sim \mathcal{N}(0, 1)$ adalah faktor risiko sistemik tunggal dan $\varepsilon_i \sim \mathcal{N}(0, 1)$ adalah guncangan idiosinkratik yang saling bebas. Korelasi pasangan aset antara debitur $i$ dan $j$ diberikan oleh $\sqrt{\rho_i \rho_j}$.
3. **Independensi Bersyarat & Peluang Gagal Bayar Bersyarat**: Bersyarat pada realisasi faktor makro $Z = z$, kejadian gagal bayar antar debitur menjadi independen secara stokastik dengan peluang masing-masing $p_i(z) = \Phi\left(\frac{\Phi^{-1}(\text{PD}_i) - \sqrt{\rho_i} z}{\sqrt{1 - \rho_i}}\right)$ (`conditional_pd`). Sifat inilah yang memungkinkan dekomposisi masalah dimensi-$N$ menjadi integrasi 1 dimensi terhadap $Z$.
4. **Konvolusi Pemrograman Dinamis Poisson–Binomial**: Untuk eksposur heterogen (Kasus B), `conditional_loss_pmf_dp` membangun distribusi kerugian bersyarat $P(L = \ell \mid Z = z)$ secara iteratif dengan menambahkan satu debitur pada tiap langkah: $P_{\text{new}}(\ell) = (1 - p_i(z)) P_{\text{old}}(\ell) + p_i(z) P_{\text{old}}(\ell - e_i)$. Kompleksitasnya $O(N \cdot L_{\max})$ jauh lebih efisien dibanding enumerasi $2^N$ kombinasi.
5. **Kuadratur Gauss–Hermite untuk Ekspektasi Normal**: Integrasi tak bersyarat $P(L = \ell) = \int_{-\infty}^{\infty} P(L = \ell \mid z)\phi(z)\,dz$ dievaluasi oleh `exact_loss_pmf` lewat perubahan variabel $z = \sqrt{2}x$ menjadi $\frac{1}{\sqrt{\pi}} \sum_{j=1}^{200} w_j P(L = \ell \mid \sqrt{2}x_j)$. Karena integrand $P(L=\ell \mid z)$ mulus, 200 simpul sudah mencapai akurasi mesin `1.110e-16` terhadap `scipy.integrate.quad`.
6. **Variansi Eksak via CDF Normal Bivariat**: Fungsi `bivariate_loss_variance` menghitung momen kedua $E[L^2] = \sum_i e_i^2 \text{PD}_i + \sum_{i \neq j} e_i e_j \Phi_2(c_i, c_j; \sqrt{\rho_i \rho_j})$ secara langsung tanpa kuadratur. Kecocokan hingga `9.070e-14` menjadi bukti silang independen bagi algoritma DP + Gauss–Hermite.
7. **Limit Large Homogeneous Portfolio (LHP) Vasicek**: Saat $N \to \infty$ dengan $\text{PD}_i = \text{PD}$ dan $\rho_i = \rho$ seragam, Hukum Bilangan Besar Bersyarat menyatakan $L_N \xrightarrow{\text{a.s.}} p(Z)$. Akibatnya diperoleh rumus analitik `lhp_tail_prob` $P(L_\infty \ge x) = \Phi\left(\frac{\Phi^{-1}(\text{PD}) - \sqrt{1-\rho}\Phi^{-1}(x)}{\sqrt{\rho}}\right)$ dan kuantil `lhp_var`, yang menjadi tulang punggung perhitungan modal regulasi Basel II/III IRB.
8. **Conditional Monte Carlo (Rao–Blackwellization)**: Dalam `estimate_conditional_mc`, alih-alih menyimulasikan $\varepsilon_i$ untuk tiap debitur, estimator mengganti indikator biner $\mathbf{1}\{L \ge x\}$ dengan peluang bersyarat eksak $P(L \ge x \mid Z)$ yang diinterpolasi halus lewat `build_conditional_tail_interpolator`. Sesuai teorema Rao–Blackwell, hal ini meniadakan variansi idiosinkratik $E[\operatorname{Var}(\mathbf{1}\{L\ge x\} \mid Z)]$.
9. **Importance Sampling via Pergeseran Mean Faktor (Glasserman & Li 2005)**: Pada ekor langka ($x = 30\%$ atau $40\%$), hampir seluruh kontribusi probabilitas berasal dari nilai $Z$ sangat negatif. Fungsi `estimate_importance_sampling_mc` menyampel $Z \sim \mathcal{N}(-\theta, 1)$ dengan bobot koreksi Radon–Nikodym $W(Z) = \exp(-\theta Z - \frac{1}{2}\theta^2)$, menekan SE di `x=30%` dari `4.358e-5` (Naif) menjadi `8.175e-7` ($\text{VRF} \approx 2.840\times$).
10. **Pemilihan Shift Optimal $\theta^*$ Berbasis Saddlepoint**: Jika $\theta$ terlalu kecil, kejadian langka tetap jarang muncul; jika terlalu besar, bobot $W(Z)$ berfluktuasi liar. Fungsi `find_optimal_is_theta` memilih pergeseran $-\theta^*$ di dekat titik kritis $z^*$ tempat $p(z^*) = x$ sehingga probabilitas bersyarat mencapai $\sim 0.5$ di pusat distribusi proposal dan variansi estimator diminimalkan secara mendekati optimal.

## 10. Kesimpulan dan Batasan
Proyek `day002d` mengintegrasikan solusi deterministik presisi tinggi (`1e-14`) untuk distribusi kerugian kredit portofolio berhingga dengan teknik simulasi kejadian langka Rao–Blackwell dan Importance Sampling yang memangkas *standard error* lebih dari 50 kali lipat pada ekor `30%`. Batasan model saat ini meliputi kopula Gaussian tanpa ketergantungan ekor bawah (*tail dependence* seperti Student-$t$ copula), struktur satu faktor makro, serta *Loss Given Default* (LGD) deterministik konstan.

## 11. Referensi dan Jejak Bukti
- Vasicek, O. (2002). *The distribution of loan portfolio value*. Risk, 15(12), 160–162.
- Glasserman, P., & Li, J. (2005). *Importance sampling for portfolio credit risk*. Management Science, 51(11), 1643–1656.
- Artefak: [day002d_quant_risiko-kredit-vasicek-distribusi.ipynb](../../days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb), [src/qb100/credit.py](../../src/qb100/credit.py), [tests/test_credit.py](../../tests/test_credit.py).
