# Penjelasan Proyek — day002c: HMM Dua Rezim — Forward-Backward, Baum–Welch & Viterbi (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002c` (Hari 2, Sesi c)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb](../../days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb)
- **Ringkasan**: Implementasi NumPy murni untuk *Hidden Markov Model* (HMM) Gaussian: algoritma forward-backward ber-skala vs log-domain, ekspektasi-maksimisasi Baum–Welch, dekoding Viterbi, studi Monte Carlo 200 replikasi, dan seleksi model BIC ($K \in \{1, 2, 3\}$).

## 2. Masalah dan Tujuan
Deret return keuangan mengalami pergantian rezim laten antara fase tenang (*bull/low-vol*) dan bergejolak (*bear/high-vol*). Rekurensi forward naif $\alpha_t(j) = [\sum_i \alpha_{t-1}(i) A_{ij}] b_j(y_t)$ mengalami *floating-point overflow/underflow* pada deret panjang ($T=2.000$).

## 3. Apa yang Dikerjakan Proyek
Modul `src/qb100/hmm.py` dan notebook `day002c` mendemonstrasikan kegagalan forward naif pada $T=2.000$, memverifikasi kesamaan hingga presisi mesin (`< 1e-14`) antara algoritma forward-backward ber-skala Rabiner dan domain-log (`logsumexp`), mengestimasi parameter dengan Baum–Welch EM (monoton naik), mendekode jalur rezim dengan Viterbi (`akurasi 96%`), menjalankan studi MC 200 replikasi, dan memilih jumlah rezim via BIC.

## 4. Uraian Kode
- [src/qb100/hmm.py](../../src/qb100/hmm.py): fungsi `simulate_hmm`, `forward_naive`, `forward_backward_scaled`, `forward_backward_log`, `viterbi_decode`, `baum_welch_fit`, dan `compute_bic`.
- [days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb](../../days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb): 11 sel notebook eksperimen dan grafik.
- [tests/test_hmm.py](../../tests/test_hmm.py): 3 suite unit test pytest untuk overflow/scaling, monotonisitas Baum–Welch, dan seleksi BIC.

## 5. Metode dan Alur Kerja
Simulasi deret $T=2.000$ dari HMM 2-rezim ($\mu_0=0.0005, \sigma_0=0.008, \mu_1=-0.001, \sigma_1=0.025, A_{00}=0.98, A_{11}=0.96$) $\to$ uji `forward_naive` (overflow `inf`) vs `forward_backward_scaled` & `forward_backward_log` $\to$ fit Baum–Welch dengan penyelarasan urutan state (*label switching*) $\to$ dekoding Viterbi & seleksi BIC untuk $K=1,2,3$.

## 6. Output yang Ditargetkan
Forward ber-skala dan log-domain identik hingga `< 1e-13`; log-likelihood Baum–Welch monoton tak-menurun di setiap iterasi; akurasi Viterbi $> 90\%$; BIC memilih $K^* = 2$ dengan benar.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Stabilitas numerik: Scaled vs Log-domain Forward-Backward | `tercapai` | Output sel 3: `forward_naive` overflow (`[inf, inf]`), sementara `|d_LL| = 7.11e-15`, `max|d_gamma| = 3.89e-15`, `max|d_xi| = 3.38e-14` |
| Monotonisitas Baum–Welch EM & Akurasi Viterbi | `tercapai` | Output sel 5: konvergen dalam 20 iterasi (`LL = 5878.034`), `min delta LL = +5.572e-9 >= 0`, akurasi Viterbi = `96.00%` |
| Studi Monte Carlo 200 replikasi | `tercapai` | Output sel 7: `99.50% (199/200)` fit konvergen, parameter tak-bias |
| Seleksi jumlah rezim via BIC ($K \in \{1,2,3\}$) | `tercapai` | Output sel 9: BIC $K=1$: `-10781.64`, **$K=2$: `-11702.86` (minimum)**, $K=3$: `-11665.58` |
| Unit test `tests/test_hmm.py` | `tercapai` | `pytest -q` dijalankan pada sesi `day003c`: 3/3 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `pytest -q` serta output notebook.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Rantai Markov Tersembunyi Waktu-Diskrit**: Proses bivariat $(S_t, Y_t)$ di mana $S_t \in \{0,\dots,K-1\}$ adalah rantai Markov orde-1 dengan matriks transisi $A$ dan $Y_t | (S_t=k) \sim N(\mu_k, \sigma_k^2)$ (`simulate_hmm`).
2. **Distribusi Stasioner Rantai Markov**: Vektor eigen kiri $\pi A = \pi$ untuk 2 rezim adalah $\pi_0 = (1-A_{11})/(2-A_{00}-A_{11})$, menentukan proporsi waktu jangka panjang dan distribusi durasi geometrik $1/(1-A_{kk})$.
3. **Rekurensi Forward-Backward Ber-Skala (Rabiner 1989)**: Menormalkan $\hat{\alpha}_t(k) = \tilde{\alpha}_t(k) / c_t$ dengan $c_t = \sum_k \tilde{\alpha}_t(k)$ menjaga $\sum_k \hat{\alpha}_t(k) = 1$ dan memberikan $\ell(Y) = \sum_{t=1}^T \ln c_t$ tanpa *overflow/underflow* (`forward_backward_scaled`).
4. **Trik Log-Sum-Exp pada Domain Logaritma**: $\ln \sum_i e^{a_i} = m + \ln \sum_i e^{a_i - m}$ dengan $m = \max_i a_i$ memberikan verifikasi independen terhadap skema penskalaan Rabiner (`forward_backward_log`).
5. **Probabilitas Posterior Halus ($\gamma_t$ dan $\xi_t$)**: $\gamma_t(i) = P(S_t=i | Y_{1:T}) \propto \hat{\alpha}_t(i)\hat{\beta}_t(i)$ dan $\xi_t(i,j) = P(S_t=i, S_{t+1}=j | Y_{1:T})$ menjadi statistik cukup pada langkah E algoritma EM.
6. **Algoritma Baum–Welch (EM untuk HMM)**: Memaksimalkan fungsi bantu $Q(\theta, \theta^{(m)}) = E_{\theta^{(m)}}[\ln P(Y, S | \theta) | Y]$ dalam bentuk tertutup untuk $\pi, A, \mu, \sigma^2$ (`baum_welch_fit`).
7. **Monotonisitas Likelihood EM (Dempster–Laird–Rubin 1977)**: Ketaksamaan Jensen menjamin $\ell(\theta^{(m+1)}) \ge \ell(\theta^{(m)})$ pada setiap langkah Baum–Welch, digunakan sebagai invarian tes unit.
8. **Pemrograman Dinamis Viterbi (Max-Sum)**: Rekurensi $\delta_t(j) = \max_i [\delta_{t-1}(i) + \ln A_{ij}] + \ln b_j(y_t)$ dengan *backpointer* menemukan lintasan rezim global paling mungkin $\arg\max_{S_{1:T}} P(S_{1:T} | Y_{1:T})$ (`viterbi_decode`).
9. **Non-Identifiabilitas Permutasi (Label Switching)**: Likelihood HMM invarian terhadap permutasi indeks rezim; diselesaikan dengan mengurutkan state berdasarkan volatilitas $\hat{\sigma}_0 < \hat{\sigma}_1$.
10. **Bayesian Information Criterion (Schwarz 1978)**: $\text{BIC}(K) = -2\ell(\hat{\theta}_K) + p_K \ln T$ dengan $p_K = K^2 + 2K - 1$ menghukum kompleksitas model dan memilih $K^*=2$ (`compute_bic`).

## 10. Kesimpulan dan Batasan
Proyek `day002c` membuktikan kestabilan numerik algoritma HMM ber-skala dan ketepatan estimasi rezim laten pada deret waktu stokastik.

## 11. Referensi dan Jejak Bukti
- Rabiner, L. R. (1989). *A tutorial on hidden Markov models and selected applications in speech recognition*. Proceedings of the IEEE, 77(2), 257–286.
- Hamilton, J. D. (1989). *A new approach to the economic analysis of nonstationary time series and the business cycle*. Econometrica, 57(2), 357–384.
- Artefak: [day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb](../../days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb), [src/qb100/hmm.py](../../src/qb100/hmm.py).
