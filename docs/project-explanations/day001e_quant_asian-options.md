# Penjelasan Proyek — day001e: Opsi Asia — Monte Carlo, Variance Reduction & Greeks (Go)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001e` (Hari 1, Sesi e)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `bonus` (off-roadmap) · **Jenis**: `quant` · **Bahasa**: Go (`go`)
- **Folder Proyek**: [go/day001e_asian_options](../../go/day001e_asian_options) · **Metadata**: [go/day001e_asian_options/meta.json](../../go/day001e_asian_options/meta.json)
- **Rujukan Silang**: Implementasi kembar dari [day001c (Python)](day001c_quant_opsi-asia-monte-carlo-variance.md) dengan spesifikasi parameter identik ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000, \text{seed}=42$).

## 2. Masalah dan Tujuan
Bagaimana mereplikasi mesin penilaian opsi Asia aritmetik (Naive MC, Antithetic Variates, Control Variate Kemna–Vorst, serta Delta Pathwise & Bump-and-Revalue CRN) di dalam ekosistem **Go** secara deterministik dan tervalidasi silang secara statistik terhadap Python `day001c` (anchor Black–Scholes `10.450583572185565`, Kemna–Vorst `5.5655088313493559`)?

## 3. Apa yang Dikerjakan Proyek
Program biner/runner di `go/day001e_asian_options` menjalankan simulasi $N=200.000$ jalur GBM ($M=252$ langkah), mengevaluasi sanity check Eropa, ketiga estimator harga opsi Asia, kedua estimator Delta, dan kurva konvergensi $N \in [1.000, 200.000]$, lalu menulis hasil ke [outputs/day001e/results.json](../../outputs/day001e/results.json) serta grafik SVG ke `outputs/figures/day001e_*.svg`.

## 4. Uraian Kode
- [go/day001e_asian_options/go.mod](../../go/day001e_asian_options/go.mod): modul Go 1.22+ murni tanpa dependensi eksternal.
- [go/day001e_asian_options/options/options.go](../../go/day001e_asian_options/options/options.go): rumus tertutup BS & Kemna–Vorst (`math.Erfc`), RNG `xoshiro256++` + Box–Muller, worker pool goroutine deterministik per chunk, serta estimator harga & Delta.
- [go/day001e_asian_options/cmd/run/main.go](../../go/day001e_asian_options/cmd/run/main.go) & [options_test.go](../../go/day001e_asian_options/options/options_test.go): runner CLI dan tes `go test`.

## 5. Metode dan Alur Kerja
Arsitektur spesifik **Go**: menggunakan Goroutine worker pool deterministik dengan pembagian chunk tetap dan seed per chunk, serta `math.Erfc` dari standard library Go. Setiap worker menghitung statistik cukup lokal pada jalur yang ditugaskan, lalu direduksi secara deterministik.

## 6. Output yang Ditargetkan
Cocok dengan anchor Black–Scholes (`10.450583572185565`) dan Kemna–Vorst (`5.5655088313493559`), MC Eropa dalam 3 SE dari BS, VRF Antithetic $\approx 2\times$ dan Control Variate $> 1200\times$, harga CV konsisten lintas bahasa di sekitar `~5.782`, Delta pathwise vs bump dalam 3 SE, dan slope log-log SE $\approx -0.50$.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Anchor tertutup BS & Kemna–Vorst vs Python `1c` | `tercapai` | Di [outputs/day001e/results.json](../../outputs/day001e/results.json): `bs_call = 10.450583572185568`, `y_exact = 5.565508831349358` |
| Sanity check MC Eropa (`z < 3`) | `tercapai` | MC Eropa = `10.508448 ± 0.032994`, `z_score = 1.7538 < 3.0` |
| Estimator Naive, Antithetic & Control Variate | `tercapai` | Naive = `5.791336 ± 0.017898`, Antithetic = `5.799159 ± 0.012400` (`VRF = 2.08x`), Control Variate = `5.783085 ± 0.000496` (`VRF = 1304.7x`, `beta = 1.035206`) |
| Konsistensi silang harga CV terhadap Python `1c` (`5.782468`) | `tercapai` | Selisih `+0.000617` konsisten dengan SE `0.000496` |
| Delta Pathwise vs Bump-and-Revalue CRN (`z < 3`) | `tercapai` | Pathwise = `0.589879 ± 0.001179`, Bump = `0.589585 ± 0.001157`, `z_score = 0.1777 < 3.0` |
| Slope konvergensi log-log SE ($\approx -0.50$) | `tercapai` | Slopes: Naive `-0.4958`, Antithetic `-0.4996`, Control Variate `-0.4911` |
| Verifikasi test suite Go | `tercapai` | `go test ./...` dijalankan pada sesi `day003c`: paket `options` lulus bersih (`ok`); waktu inti 200k jalur = `440.7 ms` |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi. **FAKTA**: Test suite `go/day001e_asian_options` telah dijalankan ulang pada sesi `day003c` dan lulus bersih.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Arsitektur Konkurensi & Memori Go**: Konkurensi `sync.WaitGroup` di Go mempartisi $N=200.000$ jalur ke sejumlah chunk tetap yang menulis ke slot berindeks tetap sebelum direduksi berurutan, mencegah *race condition*.
2. **Determinisme PRNG Lintas Worker**: Pembagian stream bilangan acak berbasis indeks jalur/chunk memastikan hasil simulasi tidak bergantung pada penjadwalan thread/proses OS.
3. **Transformasi Box–Muller Normal Standar**: Konversi pasangan seragam $(U_1, U_2) \in (0,1)^2$ menjadi normal standar independen $Z = \sqrt{-2\ln U_1}\cos(2\pi U_2)$.
4. **Komputasi CDF Normal $\Phi(x)$ via `erfc`**: Identitas $\Phi(x) = \frac{1}{2}\operatorname{erfc}(-x/\sqrt{2})$ menjaga akurasi di ekor distribusi.
5. **Simulasi Eksak Log-GBM Diskrit**: Akumulasi log-harga $\ln S_{t_k} = \ln S_{t_{k-1}} + (r - \frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k$ bebas bias diskretisasi waktu.
6. **Rata-Rata Geometrik via Penjumlahan Logaritma**: $\bar{S}_G = \exp(\frac{1}{M}\sum_{k=1}^M \ln S_{t_k})$ mencegah *overflow/underflow* perkalian $M=252$ harga saham.
7. **Koreksi Control Variate Kemna–Vorst**: Memanfaatkan korelasi $\rho > 0.9996$ antara payoff aritmetik dan geometrik pada jalur yang sama untuk mencapai VRF $> 1300\times$.
8. **Antithetic Variates pada Ruang Jalur**: Pembalikan tanda vektor shock $-Z$ pada langkah yang sama memangkas pembangkitan normal menjadi setengahnya.
9. **Estimator Delta Pathwise (IPA)**: $\Delta_{\text{PW}} = e^{-rT} \mathbf{1}\{\bar{S}_A > K\} (\bar{S}_A / S_0)$ dihitung dalam satu sapuan bersama harga opsi.
10. **Common Random Numbers (CRN) pada Beda Hingga**: Penggunaan ulang stream PRNG yang identik untuk $S_0(1+h)$ dan $S_0(1-h)$ menjaga `z_score` selisih terhadap Pathwise di bawah `0.2` SE.

## 10. Kesimpulan dan Batasan
Proyek `day001e` menunjukkan bahwa spesifikasi kuantitatif opsi Asia di `day001c` dapat direplikasi di **Go** dengan kesesuaian rumus tertutup 15 digit dan konsistensi statistik penuh.

## 11. Referensi dan Jejak Bukti
- Rujukan induk: [day001c_quant_opsi-asia-monte-carlo-variance.md](day001c_quant_opsi-asia-monte-carlo-variance.md).
- Artefak: [go/day001e_asian_options/README.md](../../go/day001e_asian_options/README.md), [outputs/day001e/results.json](../../outputs/day001e/results.json), [day001e_asian-se-convergence.svg](../../outputs/figures/day001e_asian-se-convergence.svg), [day001e_asian-ci-comparison.svg](../../outputs/figures/day001e_asian-ci-comparison.svg).
