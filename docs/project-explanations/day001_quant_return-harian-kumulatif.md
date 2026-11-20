# Penjelasan Proyek — day001: Return Harian & Kumulatif (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001` (Hari 1, Sesi pertama)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `001` (ROADMAP) · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day001_quant_return-harian-kumulatif.ipynb](../../days/day001_quant_return-harian-kumulatif.ipynb)
- **Ringkasan**: Menghitung simple return harian, return kumulatif (wealth index), volatilitas tahunan, Sharpe ratio, dan maximum drawdown pada data harga saham (`BBCA.JK`, `^JKSE`, `AAPL`).

## 2. Masalah dan Tujuan
Harga nominal antar aset tidak dapat dibandingkan secara langsung karena perbedaan satuan mata uang dan level harga. Proyek ini menormalkan deret harga penutupan harian menjadi return harian $r_t = P_t/P_{t-1}-1$ dan indeks kekayaan kumulatif $(1+R_t) = \prod_{s\le t}(1+r_s)$ untuk membandingkan kinerja dan risiko lintas aset sepanjang tahun 2024.

## 3. Apa yang Dikerjakan Proyek
Notebook mengambil data harga penutupan harian 3 aset (`BBCA.JK`, `^JKSE`, `AAPL`) sebanyak 258 baris observasi melalui modul `qb100.data`, menghitung deret return harian dan kumulatif beserta tabel ringkasan metrik risiko-imbal hasil melalui `qb100.finance`, serta menyimpan grafik perbandingan ke `outputs/figures/`.

## 4. Uraian Kode
- [src/qb100/data.py](../../src/qb100/data.py): fungsi `fetch_prices` mengambil data `yfinance` dengan fallback cache Parquet/CSV lokal di `data/cache/`.
- [src/qb100/finance.py](../../src/qb100/finance.py): fungsi `simple_returns`, `cumulative_returns`, `annualized_volatility`, `sharpe_ratio`, dan `max_drawdown`.
- [src/qb100/plotting.py](../../src/qb100/plotting.py): fungsi `set_style` dan `save_figure` untuk standarisasi visualisasi.
- [days/day001_quant_return-harian-kumulatif.ipynb](../../days/day001_quant_return-harian-kumulatif.ipynb): 14 sel notebook yang menjalankan pipeline analisis dan menyimpan grafik.
- [tests/test_finance.py](../../tests/test_finance.py): 17 unit test pytest yang memverifikasi seluruh rumus return dan risiko.

## 5. Metode dan Alur Kerja
Input deret harga $P_t$ (`258` hari bursa 2024) $\to$ pembersihan *forward-fill* $\to$ perhitungan $r_t = P_t/P_{t-1}-1$ $\to$ compounding $R_t = \prod_{s=1}^t (1+r_s) - 1$ $\to$ agregasi metrik tahunan ($\sigma_{\text{ann}} = \sigma_d \sqrt{252}$, $\text{SR} = \bar{r}/\sigma_d \sqrt{252}$) $\to$ pembuatan grafik [day001_return-harian-kumulatif.png](../../outputs/figures/day001_return-harian-kumulatif.png).

## 6. Output yang Ditargetkan
Notebook berjalan bersih tanpa error, menghasilkan tabel return kumulatif akhir tahun 2024, tabel Sharpe/volatilitas/drawdown, grafik PNG di `outputs/figures/`, dan lulus seluruh tes di `tests/test_finance.py`.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Perhitungan simple & cumulative return 3 aset | `tercapai` | Output sel 7–8 notebook: return kumulatif akhir 2024 `BBCA.JK = 5.54%`, `^JKSE = -3.33%`, `AAPL = 38.35%` |
| Ringkasan metrik volatilitas & Sharpe ratio | `tercapai` | Output sel 9 notebook: Vol tahunan `BBCA.JK 21.34% (SR 0.35)`, `^JKSE 12.54% (SR -0.20)`, `AAPL 22.14% (SR 1.49)` |
| Visualisasi & penyimpanan grafik | `tercapai` | Berkas [day001_return-harian-kumulatif.png](../../outputs/figures/day001_return-harian-kumulatif.png) tersimpan |
| Unit test modul finance (`pytest -q`) | `tercapai` | `pytest -q` dijalankan di sesi ini: 17/17 tes di `tests/test_finance.py` lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus. **FAKTA**: Notebook tidak dieksekusi ulang via `run_all.py --day 1` pada sesi `day003c` untuk menghemat waktu, namun output sel yang tersimpan di notebook lengkap dan `pytest -q` untuk `tests/test_finance.py` dijalankan dan lulus 100%.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Simple Return Satu Periode**: Definisi $r_t = (P_t - P_{t-1})/P_{t-1}$; aditif lintas aset dalam portofolio, diimplementasikan pada `simple_returns` di `src/qb100/finance.py`.
2. **Return Kumulatif & Compounding Geometrik**: Pertumbuhan modal $(1+R_t) = \prod_{s=1}^t (1+r_s) = P_t/P_0$; diimplementasikan dengan `(1 + r).cumprod() - 1` pada `cumulative_returns`.
3. **Indeks Kekayaan (Wealth Index)**: Normalisasi harga awal ke 1.0 sehingga aset berbeda mata uang (IDR vs USD) dapat diplot pada sumbu yang sama di sel 11 notebook.
4. **Asimetri Keuntungan dan Kerugian Persentase**: Penurunan $-p$ membutuhkan kenaikan $+p/(1-p)$ untuk kembali impas; dijelaskan dan diilustrasikan pada bagian konsep notebook.
5. **Volatilitas Harian dan Annualisasi Akar Waktu**: Standar deviasi sampel $s_r$ diskalakan dengan $\sqrt{252}$ hari bursa di bawah asumsi return tidak berkorelasi serial (`annualized_volatility`).
6. **Sharpe Ratio Tahunan**: Imbal hasil berlebih per unit risiko $\text{SR} = \frac{\bar{r} - r_f/252}{s_r}\sqrt{252}$; diimplementasikan pada `sharpe_ratio` dan menghasilkan `1.49` untuk AAPL.
7. **Maximum Drawdown (MDD)**: Penurunan relatif terdalam dari puncak berjalan $\min_t (W_t / \max_{s\le t} W_s - 1)$; diimplementasikan pada `max_drawdown` di `src/qb100/finance.py`.
8. **Penanganan Missing Data Kalender Bursa Berbeda**: Bursa IDX (`BBCA.JK`, `^JKSE`) dan NASDAQ (`AAPL`) memiliki hari libur berbeda; `fetch_prices` menyelaraskan indeks tanggal dengan *forward-fill*.
9. **Harga Penutupan Tersesuaikan (Adjusted Close)**: Harga yang telah disesuaikan terhadap dividen dan *stock split* agar return mencerminkan total return investor (`auto_adjust=True`).
10. **Reprodusibilitas Cache Data Lokal**: Strategi penyimpanan cache lokal di `src/qb100/data.py` agar analisis tetap deterministik saat dijalankan tanpa koneksi jaringan.

## 10. Kesimpulan dan Batasan
Proyek `day001` meletakkan fondasi modul `qb100.finance` dan `qb100.data`. Keterbatasannya adalah penggunaan simple return belum mengeksplorasi sifat aditif waktu dari log return (yang menjadi topik 003 di ROADMAP).

## 11. Referensi dan Jejak Bukti
- Campbell, J. Y., Lo, A. W., & MacKinlay, A. C. (1997). *The Econometrics of Financial Markets*. Princeton University Press.
- Artefak: [day001_quant_return-harian-kumulatif.ipynb](../../days/day001_quant_return-harian-kumulatif.ipynb), [src/qb100/finance.py](../../src/qb100/finance.py), [day001_return-harian-kumulatif.png](../../outputs/figures/day001_return-harian-kumulatif.png).
