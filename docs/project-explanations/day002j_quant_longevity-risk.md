# Penjelasan Proyek — day002j: Risiko Umur Panjang (Longevity Risk) — Lee–Carter, Anuitas Stokastik & Modal Solvabilitas (Ruby)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002j` (Hari 2, Sesi j)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Ruby (`ruby`)
- **Folder Proyek**: [ruby/day002j_longevity_risk](../../ruby/day002j_longevity_risk) · **Metadata**: [meta.json](../../ruby/day002j_longevity_risk/meta.json)
- **Ringkasan**: Pemodelan aktuaria mortalitas stokastik dalam Ruby murni tanpa gem runtime: SVD Jacobi satu-sisi buatan sendiri, estimasi model Lee–Carter (1992) & hukum Gompertz–Makeham, proyeksi stokastik RWD $k_t$ & bootstrap parameter, valuasi anuitas kohor usia 65, dekomposisi Hukum Variansi Total (risiko idiosinkratik vs sistematik), serta modal VaR 99.5% / ES 99% Solvency II.

## 2. Masalah dan Tujuan
Dalam portofolio dana pensiun atau anuitas seumur hidup, risiko waktu kematian individu (*idiosyncratic mortality risk*) habis terdiversifikasi saat ukuran kohor $N \to \infty$, tetapi risiko tren peningkatan harapan hidup (*systematic longevity risk*) tidak dapat didiversifikasi dengan memperbesar jumlah peserta.

## 3. Apa yang Dikerjakan Proyek
Proyek Ruby `ruby/day002j_longevity_risk` mengimplementasikan dekomposisi nilai singular (SVD) dari nol untuk mengestimasi permukaan log-mortalitas Lee–Carter, memfit parameter Gompertz–Makeham, memproyeksikan jalur mortalitas masa depan, memverifikasi Hukum Variansi Total dan asimtotik $\mathrm{CV}(N)$ pada portofolio anuitas, serta menghitung modal ekonomi VaR 99.5% dan ES 99%.

## 4. Uraian Kode
- [ruby/day002j_longevity_risk/lib/](../../ruby/day002j_longevity_risk/lib/): modul SVD Jacobi satu-sisi, estimator Lee–Carter dengan normalisasi identifiabilitas, Nelder–Mead Gompertz–Makeham, simulator RWD $k_t$ & anuitas kohor paralel `Process.fork`, serta penulis JSON/SVG.
- [ruby/day002j_longevity_risk/bin/run](../../ruby/day002j_longevity_risk/bin/run) & [test/test_longevity.rb](../../ruby/day002j_longevity_risk/test/test_longevity.rb): runner benchmark YJIT dan unit test `minitest`.

## 5. Metode dan Alur Kerja
Konstruksi matriks log-mortalitas $\ln m_{x,t}$ $\to$ pemusatan baris $a_x$ & SVD peringkat-1 untuk $(b_x, k_t)$ $\to$ estimasi drift RWD $\hat{d}$ dan $\sigma_k$ $\to$ rekursi aktuaria eksak $a_{65}$ vs simulasi umur kematian $\to$ dekomposisi variansi portofolio untuk berbagai ukuran kohor $N$ $\to$ perhitungan VaR 99.5%, ES 99%, dan kejut mortalitas Solvency II $-0.5\%$.

## 6. Output yang Ditargetkan
Rekonstruksi SVD & ortogonalitas `< 1e-13`; pemulihan Lee–Carter tanpa noise `< 1e-14`; anuitas simulasi 200k umur cocok dengan rekursi eksak dalam 3 SE; Hukum Variansi Total cocok; $\mathrm{CV}(N)$ konvergen ke lantai risiko sistematik positif saat $N \to \infty$; ES 99% $\ge$ VaR 99%.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Akurasi SVD Jacobi & Pemulihan Lee–Carter Bebas-Noise | `tercapai` | Galat rekonstruksi SVD `<= 1.005e-14`, pemulihan bersih Lee–Carter `max_a = 8.882e-16` |
| Proyeksi RWD $k_t$ & Bootstrap Parameter | `tercapai` | Momen teoritis dan empiris $k_{T+h}$ konsisten; ketidakpastian parameter memperlebar fan chart harapan hidup |
| Anuitas Eksak vs Simulasi & Hukum Variansi Total | `tercapai` | Simulasi 200.000 umur cocok dengan rekursi aktuaria eksak dalam 3 SE; $\mathrm{CV}(N)$ mendatar pada batas risiko sistematik saat $N$ besar |
| Modal Solvabilitas VaR 99.5%, ES 99% & Stres Deterministik | `tercapai` | VaR 99.5% dan ES 99% simulasi penuh konsisten dengan aproksimasi normal bersyarat di `bagian4` |
| Verifikasi Minitest (`test_longevity.rb`) | `tercapai` | `bundle exec rake test` dijalankan pada sesi `day003c`: 6 runs, 55 assertions, 0 failures, 0 errors |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `bundle exec rake test` serta `outputs/day002j/results.json`.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Mortalitas Faktor Lee–Carter (1992)**: $\ln m_{x,t} = a_x + b_x k_t + \varepsilon_{x,t}$ dengan kendala identifiabilitas $\sum_t k_t = 0$ dan $\sum_x b_x = 1$, sehingga $a_x = \frac{1}{T}\sum_t \ln m_{x,t}$.
2. **Teorema Aproksimasi Peringkat-Rendah Eckart–Young & SVD Jacobi Satu-Sisi**: Suku singular pertama $\sigma_1 u_1 v_1^T$ dari matriks residu $Z_{x,t} = \ln m_{x,t} - a_x$ meminimalkan galat kuadrat terkecil Frobenius; dihitung dengan rotasi ortogonal Jacobi satu-sisi (`bagian1`).
3. **Hukum Mortalitas Parametrik Gompertz–Makeham**: $\mu(x) = A + B c^x$ memisahkan risiko kematian latar belakang independen-usia ($A$) dari penuaan biologis eksponensial ($B c^x$).
4. **Random Walk with Drift (RWD) untuk Indeks Mortalitas $k_t$**: $k_{t} = k_{t-1} + d + \sigma_k e_t$ dengan estimator MLE $\hat{d} = (k_T - k_1)/(T-1)$ menghasilkan variansi ramalan yang tumbuh linear $h \sigma_k^2$.
5. **Jembatan Laju Pusat $m_{x,t}$ ke Probabilitas Kematian $q_{x,t}$**: Di bawah asumsi *uniform distribution of deaths* (UDD) atau gaya mortalitas konstan sepotong, $q_{x,t} = 1 - \exp(-m_{x,t})$ dan peluang bertahan hidup kohor ${}_t p_{x} = \prod_{j=0}^{t-1}(1 - q_{x+j, t_0+j})$.
6. **Rekursi Mundur Nilai Kini Aktuaria Anuitas Seumur Hidup**: $a_x = \sum_{t=1}^{\omega-x} v^t {}_t p_x$ memenuhi rekursi mundur $a_x = v p_x (1 + a_{x+1})$ dengan faktor diskonto $v = e^{-r}$.
7. **Hukum Variansi Total untuk Portofolio Anuitas**: $\mathrm{Var}(\bar{L}_N) = \frac{1}{N} E[\mathrm{Var}(L_i | \mathcal{F}_{\text{mort}})] + \mathrm{Var}(E[L_i | \mathcal{F}_{\text{mort}}])$; suku pertama adalah risiko kematian idiosinkratik $O(1/N)$, suku kedua adalah risiko umur panjang sistematik $O(1)$.
8. **Asimtotik Koefisien Variasi $\mathrm{CV}(N)$**: Saat $N \to \infty$, $\mathrm{CV}(\bar{L}_N) \to \frac{\sqrt{\mathrm{Var}(a_{65}(k))}}{E[a_{65}(k)]} > 0$, membuktikan bahwa diversifikasi ukuran kohor tidak dapat menghilangkan risiko tren umur panjang.
9. **Value-at-Risk (VaR 99.5%) & Expected Shortfall (ES 99%) Solvency II**: Modal ekonomi didefinisikan sebagai $\text{VaR}_{0.995}(L) - E[L]$ (kuantil tahunan 1-dalam-200) dan rata-rata kerugian bersyarat di ekor $\text{ES}_{0.99}(L)$.
10. **Aproksimasi Normal Bersyarat (Conditional CLT Approximation)**: Untuk kohor besar $N$, distribusi $L_N | k$ mendekati Gaussian $\mathcal{N}(N \mu(k), N \sigma^2(k))$, mempercepat evaluasi kuantil ekor ekstrem dibandingkan simulasi individu bersarang.

## 10. Kesimpulan dan Batasan
Proyek `day002j` mendemonstrasikan secara aktuaria dan statistik perbedaan mendasar antara risiko mortalitas idiosinkratik dan risiko umur panjang sistematik.

## 11. Referensi dan Jejak Bukti
- Lee, R. D., & Carter, L. R. (1992). *Modeling and forecasting U.S. mortality*. Journal of the American Statistical Association, 87(419), 659–671.
- Artefak: [ruby/day002j_longevity_risk/README.md](../../ruby/day002j_longevity_risk/README.md), [outputs/day002j/results.json](../../outputs/day002j/results.json).
