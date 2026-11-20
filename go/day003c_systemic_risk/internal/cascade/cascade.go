package cascade

import (
	"math"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
)

// Resolusi scan c yang dinyatakan SEBELUM dijalankan untuk membandingkan lonjakan S simulasi dengan (c1, c2).
const PredeclaredCGridStep = 0.25

// Simulate menjalankan kaskade ambang Watts pada graf CSR tak berarah.
// Jika useLIFO = false memakai antrean FIFO (BFS); jika useLIFO = true memakai tumpukan LIFO (urutan terbalik).
func Simulate(g *graphs.UndirectedCSR, seeds []int, phi float64, useLIFO bool) (int, []bool) {
	n := g.N
	failed := make([]bool, n)
	failNeighbors := make([]int32, n)
	queue := make([]int32, 0, len(seeds)*4+16)

	for _, s := range seeds {
		if s >= 0 && s < n && !failed[s] {
			failed[s] = true
			queue = append(queue, int32(s))
		}
	}

	head := 0
	for {
		var u int32
		if useLIFO {
			if len(queue) == 0 {
				break
			}
			u = queue[len(queue)-1]
			queue = queue[:len(queue)-1]
		} else {
			if head >= len(queue) {
				break
			}
			u = queue[head]
			head++
		}

		start, end := g.Offset[u], g.Offset[u+1]
		for idx := start; idx < end; idx++ {
			v := g.Edges[idx]
			if failed[v] {
				continue
			}
			failNeighbors[v]++
			degV := g.Offset[v+1] - g.Offset[v]
			if degV > 0 && float64(failNeighbors[v])+1e-12 >= phi*float64(degV) {
				failed[v] = true
				queue = append(queue, v)
			}
		}
	}

	total := 0
	for _, f := range failed {
		if f {
			total++
		}
	}
	return total, failed
}

// MeanField menghitung titik tetap terkecil q_inf (mulai dari q=0) dan ukuran kaskade akhir S_teori
// pada pendekatan pohon Erdős-Rényi dengan derajat Poisson(c).
func MeanField(c, rho0, phi float64) (qInf, sFinal float64) {
	if c <= 0 {
		return rho0, rho0
	}
	const kMax = 140
	pk := make([]float64, kMax+1)
	pk[0] = math.Exp(-c)
	for k := 1; k <= kMax; k++ {
		pk[k] = pk[k-1] * c / float64(k)
	}

	// Precompute koefisien binomial C(n, m)
	binom := make([][]float64, kMax+1)
	for n := 0; n <= kMax; n++ {
		binom[n] = make([]float64, n+1)
		binom[n][0] = 1.0
		binom[n][n] = 1.0
		for m := 1; m < n; m++ {
			binom[n][m] = binom[n-1][m-1] + binom[n-1][m]
		}
	}

	q := 0.0
	for it := 0; it < 2000; it++ {
		sumK := 0.0
		for k := 1; k <= kMax; k++ {
			weight := float64(k) * pk[k] / c
			if weight < 1e-18 {
				continue
			}
			// Koreksi literatur Watts (2002) & Gleeson-Cahalane (2007):
			// Simpul v berderajat k yang dicapai dari u (yang belum gagal) membutuhkan m dari k-1
			// tetangga lainnya gagal agar m/k >= phi (bukan (m+1)/k >= phi yang sirkular).
			mMin := int(math.Ceil(phi*float64(k) - 1e-12))
			if mMin < 1 {
				mMin = 1
			}
			if mMin > k-1 {
				continue
			}
			prob := 0.0
			for m := mMin; m <= k-1; m++ {
				prob += binom[k-1][m] * math.Pow(q, float64(m)) * math.Pow(1.0-q, float64(k-1-m))
			}
			sumK += weight * prob
		}
		qNext := rho0 + (1.0-rho0)*sumK
		if math.Abs(qNext-q) < 1e-13 {
			q = qNext
			break
		}
		q = qNext
	}
	qInf = q

	sumS := 0.0
	for k := 1; k <= kMax; k++ {
		if pk[k] < 1e-18 {
			continue
		}
		mMin := int(math.Ceil(phi*float64(k) - 1e-12))
		if mMin < 1 {
			mMin = 1
		}
		if mMin > k {
			continue
		}
		prob := 0.0
		for m := mMin; m <= k; m++ {
			prob += binom[k][m] * math.Pow(qInf, float64(m)) * math.Pow(1.0-qInf, float64(k-m))
		}
		sumS += pk[k] * prob
	}
	sFinal = rho0 + (1.0-rho0)*sumS
	return qInf, sFinal
}

// VulnerableBranching mengevaluasi fungsi kondisi kaskade global Watts:
// F(c) = (sum_{k=1}^{floor(1/phi)} k(k-1) p_k(c)) - c.
// Kaskade global dari benih sangat kecil terjadi bila F(c) >= 0.
func VulnerableBranching(c, phi float64) float64 {
	if c <= 0 {
		return 0
	}
	kVuln := int(math.Floor(1.0 / (phi - 1e-12)))
	if kVuln < 2 {
		return -c
	}
	pk := math.Exp(-c)
	sum := 0.0
	for k := 1; k <= kVuln; k++ {
		pk = pk * c / float64(k)
		if k >= 2 {
			sum += float64(k*(k-1)) * pk
		}
	}
	return sum - c
}

// PhaseBoundaries mencari batas transisi fase (c1, c2) tempat VulnerableBranching(c, phi) = 0 lewat bisection.
// Mengembalikan (c1, c2, true) jika jendela kaskade ada, atau (0, 0, false) jika phi terlalu besar.
func PhaseBoundaries(phi float64) (c1, c2 float64, hasWindow bool) {
	kVuln := int(math.Floor(1.0 / (phi - 1e-12)))
	if kVuln < 2 {
		return 0, 0, false
	}
	cPeak := float64(kVuln - 1)
	if cPeak < 0.5 {
		cPeak = 0.5
	}
	if VulnerableBranching(cPeak, phi) <= 0 {
		return 0, 0, false
	}

	// Bisection c1 pada [1e-4, cPeak]
	lo, hi := 1e-4, cPeak
	for i := 0; i < 80; i++ {
		mid := 0.5 * (lo + hi)
		if VulnerableBranching(mid, phi) >= 0 {
			hi = mid
		} else {
			lo = mid
		}
	}
	c1 = 0.5 * (lo + hi)

	// Bisection c2 pada [cPeak, 40.0]
	lo, hi = cPeak, 40.0
	for VulnerableBranching(hi, phi) > 0 && hi < 200.0 {
		hi *= 2.0
	}
	for i := 0; i < 80; i++ {
		mid := 0.5 * (lo + hi)
		if VulnerableBranching(mid, phi) >= 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	c2 = 0.5 * (lo + hi)
	return c1, c2, true
}
