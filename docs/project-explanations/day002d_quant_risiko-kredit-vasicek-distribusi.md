# Penjelasan Proyek — day002d: Risiko Kredit Vasicek — Distribusi Eksak, LHP & Importance Sampling (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002d` (Hari 2, Sesi d)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb](../../days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb)
- **Ringkasan**: Distribusi kerugian portofolio kredit satu-faktor Gaussian copula Vasicek melalui pemrograman dinamis + kuadratur Gauss–Hermite, aproksimasi *Large Homogeneous Portfolio* (LHP), serta simulasi kejadian langka via *Conditional Monte Carlo* dan *Importance Sampling* (*exponential twisting*).

## 2. Masalah dan Tujuan
Dalam portofolio kredit berkorelasi ($X_i = \sqrt{\rho} Z + \sqrt{1-\rho}\varepsilon_i$), gagal bayar ekstrem di ekor distribusi ($P(L \ge 30\%)$) memiliki peluang sangat kecil ($\sim 10^{-5}$ hingga $10^{-4}$) sehingga Monte Carlo naif gagal memberikan estimasi yang stabil.

## 3. Apa yang Dikerjakan Proyek
Modul `src/qb100/credit.py` dan notebook `day002d` menghitung distribusi eksak jumlah gagal bayar melalui konvolusi DP bersyarat pada faktor sistemik $Z$ yang diintegrasikan dengan kuadratur Gauss–Hermite 200 titik, menganalisis konvergensi portofolio finite-$N$ ke distribusi tertutup LHP Vasicek, dan membandingkan efisiensi Naive MC, Conditional MC (Rao–Blackwell), serta Importance Sampling (pergeseran mean $Z \sim N(-\theta, 1)$).

## 4. Uraian Kode
- [src/qb100/credit.py](../../src/qb100/credit.py): fungsi `conditional_pd`, `exact_loss_pmf_gh`, `lhp_tail_prob`, `lhp_vasicek_pdf`, `mc_tail_naive`, `mc_tail_conditional`, dan `mc_tail_importance_sampling`.
- [days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb](../../days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb): 9 sel notebook eksperimen dan grafik.
- [tests/test_credit.py](../../tests/test_credit.py): 4 unit test pytest untuk PMF eksak, momen kerugian, konvergensi LHP, dan pengurangan variansi IS.

## 5. Metode dan Alur Kerja
Hitung peluang gagal bayar bersyarat $p(z) = \Phi(\frac{\Phi^{-1}(\text{PD}) - \sqrt{\rho}z}{\sqrt{1-\rho}})$ di simpul Gauss–Hermite $\to$ bangun PMF binomial/Poisson-binomial bersyarat via DP lalu integrasikan terhadap $Z$ $\to$ bandingkan ekor kerugian $P(L \ge x)$ eksak dengan LHP dan tiga estimator MC ($M=200.000$).

## 6. Output yang Ditargetkan
Total massa PMF $= 1$ dan momen cocok hingga `< 1e-13`; kuadratur GH-200 cocok dengan `scipy.integrate.quad` hingga `< 1e-15`; galat LHP menyusut $O(1/N)$; Importance Sampling menghasilkan VRF ratusan hingga ribuan kali pada ekor ekstrem.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Akurasi PMF eksak (Gauss–Hermite 200 vs quad & momen) | `tercapai` | Output sel 1: `|sum(pmf)-1| <= 5.55e-16`, `rel_err Var <= 9.07e-14`, `max |PMF_GH200 - quad| = 1.11e-16` |
| Konvergensi finite-$N$ ke pendekatan LHP Vasicek | `tercapai` | Output sel 3: galat ekor LHP terhadap eksak menurun monoton saat $N$ naik dari 20 ke 1.000 |
| Reduksi variansi Conditional MC & Importance Sampling | `tercapai` | Output sel 5 & 7: ketiganya konsisten dengan nilai eksak dalam 3 SE, VRF IS mencapai ratusan/ribuan kali di ekor 30%, slope log-log SE `-0.5008` |
| Unit test `tests/test_credit.py` | `tercapai` | `pytest -q` dijalankan pada sesi `day003c`: 4/4 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `pytest -q` serta output notebook.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Struktural Merton (1974) & Ambang Gagal Bayar**: Debitur $i$ gagal bayar jika nilai aset terstandarisasi $X_i \le c_i = \Phi^{-1}(\text{PD}_i)$ pada maturitas.
2. **Gaussian Copula Satu-Faktor (Vasicek 1987, 2002)**: Dekomposisi $X_i = \sqrt{\rho_i} Z + \sqrt{1-\rho_i}\varepsilon_i$ dengan faktor makro $Z \sim N(0,1)$ dan komponen idiosinkratik $\varepsilon_i \sim N(0,1)$ iid.
3. **Independensi Bersyarat & Peluang Gagal Bayar Bersyarat**: Bersyarat pada $Z=z$, kejadian gagal bayar antar debitur menjadi independen dengan peluang $p_i(z) = \Phi\left(\frac{\Phi^{-1}(\text{PD}_i) - \sqrt{\rho_i}z}{\sqrt{1-\rho_i}}\right)$ (`conditional_pd`).
4. **Pemrograman Dinamis Distribusi Poisson–Binomial**: Rekursi konvolusi $P(K_m = k | z) = (1-p_m(z)) P(K_{m-1}=k|z) + p_m(z) P(K_{m-1}=k-1|z)$ menghitung PMF kerugian heterogen dalam waktu $O(N^2)$ (`exact_loss_pmf_gh`).
5. **Kuadratur Gauss–Hermite untuk Integrasi Normal**: $\int_{-\infty}^{\infty} f(z)\phi(z)dz = \frac{1}{\sqrt{\pi}}\sum_{j=1}^M w_j f(\sqrt{2}x_j)$ mengintegrasikan polinomial hingga derajat $2M-1$ secara eksak.
6. **Hukum Bilangan Besar Bersyarat & Limit LHP Vasicek**: Saat $N \to \infty$ pada portofolio homogen, fraksi kerugian $L_N \to p(Z)$ hampir pasti, menghasilkan CDF tertutup $F_{\text{LHP}}(x) = \Phi\left(\frac{\sqrt{1-\rho}\Phi^{-1}(x) - \Phi^{-1}(\text{PD})}{\sqrt{\rho}}\right)$ (`lhp_tail_prob`), yang menjadi dasar rumus bobot risiko Basel II IRB.
7. **Hukum Variansi Total untuk Portofolio Kredit**: $\mathrm{Var}(L) = E[\mathrm{Var}(L|Z)] + \mathrm{Var}(E[L|Z])$; suku pertama menyusut $O(1/N)$ sedangkan suku kedua tetap positif akibat korelasi sistemik $\rho > 0$.
8. **Conditional Monte Carlo (Rao–Blackwellization)**: Mengganti indikator simulasi $\mathbf{1}\{L \ge k\}$ dengan peluang ekor binomial eksak bersyarat $P(L \ge k | Z)$ menghilangkan seluruh variansi idiosinkratik (`mc_tail_conditional`).
9. **Importance Sampling via Pergeseran Mean Faktor (Glasserman & Li 2005)**: Mengambil sampel faktor dari $Z \sim N(-\theta, 1)$ mendorong realisasi ke rezim krisis dan dikoreksi dengan rasio likelihood Radon–Nikodym $W(Z) = \exp(-\theta Z - \frac{1}{2}\theta^2)$ (`mc_tail_importance_sampling`).
10. **Pemilihan Parameter Shift Optimal Saddlepoint**: Memilih $-\theta^* \approx z_{\text{kritis}}$ tempat $p(z^*) = x$ memusatkan distribusi sampling tepat pada ambang kerugian langka sehingga rasio likelihood stabil dan VRF maksimal.

## 10. Kesimpulan dan Batasan
Proyek `day002d` mengintegrasikan solusi semi-analitik kuadratur dengan simulasi kejadian langka Rao–Blackwell + Importance Sampling untuk risiko kredit portofolio.

## 11. Referensi dan Jejak Bukti
- Vasicek, O. (2002). *The distribution of loan portfolio value*. Risk, 15(12), 160–162.
- Glasserman, P., & Li, J. (2005). *Importance sampling for portfolio credit risk*. Management Science, 51(11), 1643–1656.
- Artefak: [day002d_quant_risiko-kredit-vasicek-distribusi.ipynb](../../days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb), [src/qb100/credit.py](../../src/qb100/credit.py).
