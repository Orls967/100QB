# Hari 9b (`day009b`) — Genetika Populasi (Ruby Murni, 60 Tahap)

Implementasi murni Ruby (`stdlib` saja tanpa gem eksternal, bersih di bawah `ruby -w`, target **5.000–8.000 baris** modular di `bin/popgen` dan `lib/popgen/*.rb`) untuk simulasi dan jangkar eksak genetika populasi klasik serta teori koalesen.

## Fitur & Cabang Utama
- **Sampler Eksak**: SplitMix64 → `Random` deterministik (`20261127`), binomial dua-cabang (inversi & BTPE exact rejection), multinomial berantai, Poisson (inversi & PTRS), eksponensial, geometrik.
- **Rantai Markov Eksak**: Matriks transisi Wright–Fisher dan Moran dalam aritmetika `Rational` eksak ($N \le 12$) vs `Float` ($N \le 40$) dengan eliminasi Gauss, probabilitas fiksasi, waktu fiksasi rata-rata, dan distribusi stasioner.
- **Model Populasi**: Wright–Fisher haploid & diploid, Moran birth-death, seleksi aditif/dominan/overdominan, mutasi dua alel & alel tak-hingga, model pulau 2–8 deme & stepping-stone 1D, Hardy–Weinberg & selfing.
- **Teori Koalesen**: Kingman (`T_MRCA`, `L`, `S`, SFS $\xi_i$, Watterson $\theta_W$, Tajima's $D$, Ewens sampling formula), koalesen terstruktur (sistem linier eksak), sapuan selektif (*hitchhiking*) & peluruhan LD ($D, D', r^2$), serta koalesen rekombinasi dua lokus (**ARG-lite**).
- **Paralelisme & Optimasi**: `Process.fork` + pipa deterministik bit-identik lintas 1/2/4 proses; tangga optimasi `V0`–`V4` (`Hash` objek vs `Array` datar vs binomial inversi/BTPE vs genom bitset `Integer` vs `fork`).
- **60 Tahap (`S01`–`S60`)**: Menghasilkan $\ge 150$ berkas keluaran (CSV ringkas di `out/`, $\ge 40$ SVG tulis tangan valid XML ke `../../../outputs/figures/day009b_*.svg`, dan `../../../outputs/day009b/results.json`).

## Penggunaan

```sh
rake selftest            # satu-satunya suite tes
rake lines               # verifikasi 5.000-8.000 baris Ruby
rake run                 # eksekusi skala penuh 60 tahap
rake run SCALE=small     # uji cepat tanpa menimpa klaim resmi
ruby -w -Ilib bin/popgen workers 4
```
