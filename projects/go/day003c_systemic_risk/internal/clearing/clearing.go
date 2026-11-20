package clearing

import (
	"math"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/linalg"
)

// Result merangkum hasil kliring Eisenberg-Noe.
type Result struct {
	P            []float64
	DefaultCount int
	TotalLoss    float64
	Iterations   int
	HistoryDiff  []float64
}

// Phi mengevaluasi pemetaan kliring Eisenberg-Noe Phi_i(p) = min(p_bar_i, e_i + sum_j Pi_{ji} p_j).
func Phi(net *graphs.InterbankNet, p []float64, out []float64) {
	for i := 0; i < net.N; i++ {
		inflow := net.E[i]
		for _, ed := range net.In[i] {
			inflow += ed.Pi * p[ed.From]
		}
		if inflow < net.PBar[i] {
			out[i] = inflow
		} else {
			out[i] = net.PBar[i]
		}
	}
}

// Picard menjalankan iterasi titik tetap monoton menurun dari p^(0) = p_bar.
func Picard(net *graphs.InterbankNet, tol float64, maxIter int) Result {
	return picardFrom(net, net.PBar, tol, maxIter)
}

// PicardFromZero menjalankan iterasi titik tetap monoton naik dari p^(0) = 0 (untuk uji keunikan).
func PicardFromZero(net *graphs.InterbankNet, tol float64, maxIter int) Result {
	zeros := make([]float64, net.N)
	return picardFrom(net, zeros, tol, maxIter)
}

func picardFrom(net *graphs.InterbankNet, initP []float64, tol float64, maxIter int) Result {
	n := net.N
	p := make([]float64, n)
	copy(p, initP)
	next := make([]float64, n)
	hist := make([]float64, 0, 32)

	it := 0
	for it < maxIter {
		it++
		Phi(net, p, next)
		maxDiff := 0.0
		for i := 0; i < n; i++ {
			d := math.Abs(next[i] - p[i])
			if d > maxDiff {
				maxDiff = d
			}
			p[i] = next[i]
		}
		hist = append(hist, maxDiff)
		if maxDiff <= tol {
			break
		}
	}
	return summarize(net, p, it, hist)
}

// FictitiousDefault menjalankan algoritma default fiktif Eisenberg-Noe (2001) menggunakan solver LU buatan sendiri.
func FictitiousDefault(net *graphs.InterbankNet) Result {
	n := net.N
	p := make([]float64, n)
	copy(p, net.PBar)
	inD := make([]bool, n)

	rounds := 0
	for rounds < n+1 {
		rounds++
		added := false
		for i := 0; i < n; i++ {
			if inD[i] {
				continue
			}
			inflow := net.E[i]
			for _, ed := range net.In[i] {
				inflow += ed.Pi * p[ed.From]
			}
			if inflow < net.PBar[i]-1e-12 {
				inD[i] = true
				added = true
			}
		}
		if !added {
			break
		}

		// Bangun indeks himpunan gagal D
		dIdx := make([]int, 0, n)
		posInD := make([]int, n)
		for i := 0; i < n; i++ {
			if inD[i] {
				posInD[i] = len(dIdx)
				dIdx = append(dIdx, i)
			} else {
				posInD[i] = -1
			}
		}
		m := len(dIdx)
		a := make([][]float64, m)
		b := make([]float64, m)
		for r, bankI := range dIdx {
			row := make([]float64, m)
			row[r] = 1.0
			rhs := net.E[bankI]
			for _, ed := range net.In[bankI] {
				c := posInD[ed.From]
				if c >= 0 {
					row[c] -= ed.Pi
				} else {
					rhs += ed.Pi * net.PBar[ed.From]
				}
			}
			a[r] = row
			b[r] = rhs
		}
		sol, ok := linalg.SolveLU(a, b)
		if !ok {
			// Komponen siklik tanpa arus masuk eksternal (tidak-teratur): jatuh balik ke batas atas Picard
			return Picard(net, 1e-13, 10000)
		}
		for i := 0; i < n; i++ {
			p[i] = net.PBar[i]
		}
		for r, bankI := range dIdx {
			val := sol[r]
			if val < 0 {
				val = 0
			}
			p[bankI] = val
		}
	}
	return summarize(net, p, rounds, nil)
}

// BruteForce mengenumerasi seluruh 2^n himpunan gagal D (untuk n <= 12) dan memilih vektor kliring konsisten terbesar.
func BruteForce(net *graphs.InterbankNet) Result {
	n := net.N
	bestP := make([]float64, n)
	bestSum := -1.0
	posInD := make([]int, n)
	dIdx := make([]int, 0, n)

	totalMasks := 1 << n
	for mask := 0; mask < totalMasks; mask++ {
		dIdx = dIdx[:0]
		for i := 0; i < n; i++ {
			if (mask>>i)&1 == 1 {
				posInD[i] = len(dIdx)
				dIdx = append(dIdx, i)
			} else {
				posInD[i] = -1
			}
		}
		m := len(dIdx)
		candP := make([]float64, n)
		for i := 0; i < n; i++ {
			candP[i] = net.PBar[i]
		}
		if m > 0 {
			a := make([][]float64, m)
			b := make([]float64, m)
			for r, bankI := range dIdx {
				row := make([]float64, m)
				row[r] = 1.0
				rhs := net.E[bankI]
				for _, ed := range net.In[bankI] {
					c := posInD[ed.From]
					if c >= 0 {
						row[c] -= ed.Pi
					} else {
						rhs += ed.Pi * net.PBar[ed.From]
					}
				}
				a[r] = row
				b[r] = rhs
			}
			sol, ok := linalg.SolveLU(a, b)
			if !ok {
				continue
			}
			for r, bankI := range dIdx {
				candP[bankI] = sol[r]
			}
		}

		// Verifikasi konsistensi:
		// 1) Untuk i in D: -1e-11 <= p_i <= p_bar_i + 1e-11 dan p_i == e_i + (Pi^T p)_i
		// 2) Untuk i not in D: e_i + (Pi^T p)_i >= p_bar_i - 1e-11
		valid := true
		sumP := 0.0
		for i := 0; i < n; i++ {
			if candP[i] < -1e-10 || candP[i] > net.PBar[i]+1e-10 {
				valid = false
				break
			}
			inflow := net.E[i]
			for _, ed := range net.In[i] {
				inflow += ed.Pi * candP[ed.From]
			}
			if posInD[i] >= 0 {
				if math.Abs(candP[i]-inflow) > 1e-9 || candP[i] > net.PBar[i]+1e-9 {
					valid = false
					break
				}
			} else {
				if inflow < net.PBar[i]-1e-9 {
					valid = false
					break
				}
			}
			sumP += candP[i]
		}
		if valid && sumP > bestSum {
			bestSum = sumP
			for i := 0; i < n; i++ {
				if candP[i] < 0 {
					candP[i] = 0
				}
				if candP[i] > net.PBar[i] {
					candP[i] = net.PBar[i]
				}
				bestP[i] = candP[i]
			}
		}
	}
	return summarize(net, bestP, totalMasks, nil)
}

func summarize(net *graphs.InterbankNet, p []float64, it int, hist []float64) Result {
	defCnt := 0
	loss := 0.0
	for i := 0; i < net.N; i++ {
		diff := net.PBar[i] - p[i]
		if diff > 1e-9 {
			defCnt++
		}
		if diff > 0 {
			loss += diff
		}
	}
	return Result{
		P:            p,
		DefaultCount: defCnt,
		TotalLoss:    loss,
		Iterations:   it,
		HistoryDiff:  hist,
	}
}
