# `day008c_vicsek_cpp` — Model Vicsek (Aktif-Materi) di C++20 dengan Cell List

Implementasi C++20 murni (`427` baris sumber di `src/`, tanpa pustaka eksternal) untuk model kawanan **Vicsek (1995)** pada kotak periodik $L \times L$ ($\rho = N/L^2 = 4$, $v_0 = 0{,}03$, $r = 1$), membandingkan **Cell List $O(N)$** terhadap **Brute Force $O(N^2)$** hingga paritas bit-identik (`N=1000`, `500` langkah), kurva parameter order $\phi(\eta)$ untuk $L \in \{16, 32, 64\}$, skala waktu hingga $N = 256{.}000$, dan 4 grafik SVG mandiri.

## Perintah

```bash
make test
make run
.venv/bin/python py/analyze.py
```
