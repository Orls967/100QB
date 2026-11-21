package stats

import (
	"math"
	"sort"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/parallel"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
)

// Mean menghitung rata-rata aritmetika.
func Mean(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	s := 0.0
	for _, x := range xs {
		s += x
	}
	return s / float64(len(xs))
}

// Variance menghitung varians sampel tak-bias (pembagi n - 1).
func Variance(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return 0
	}
	m := Mean(xs)
	ss := 0.0
	for _, x := range xs {
		d := x - m
		ss += d * d
	}
	return ss / float64(n-1)
}

// SD menghitung simpangan baku sampel.
func SD(xs []float64) float64 {
	return math.Sqrt(Variance(xs))
}

// SE menghitung galat baku rata-rata sampel.
func SE(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	return SD(xs) / math.Sqrt(float64(len(xs)))
}

// SummaryMoments menghitung rata-rata dan galat baku dari jumlah akumulasi (sum, sumSq, n).
func SummaryMoments(sum, sumSq float64, n int) (mean, se float64) {
	if n <= 0 {
		return 0, 0
	}
	fn := float64(n)
	mean = sum / fn
	if n < 2 {
		return mean, 0
	}
	variance := (sumSq - sum*sum/fn) / float64(n-1)
	if variance < 0 {
		variance = 0
	}
	se = math.Sqrt(variance / fn)
	return mean, se
}

// Quantile menghitung kuantil empiris level q dalam [0, 1] pada salinan terurut.
func Quantile(xs []float64, q float64) float64 {
	n := len(xs)
	if n == 0 {
		return 0
	}
	cp := make([]float64, n)
	copy(cp, xs)
	sort.Float64s(cp)
	if q <= 0 {
		return cp[0]
	}
	if q >= 1 {
		return cp[n-1]
	}
	pos := q * float64(n-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	w := pos - float64(lo)
	return (1.0-w)*cp[lo] + w*cp[hi]
}

// LinearSlope menghitung kemiringan regresi linier OLS y = a + b x.
func LinearSlope(xs, ys []float64) float64 {
	n := len(xs)
	if n < 2 || len(ys) != n {
		return 0
	}
	mx := Mean(xs)
	my := Mean(ys)
	num, den := 0.0, 0.0
	for i := 0; i < n; i++ {
		dx := xs[i] - mx
		num += dx * (ys[i] - my)
		den += dx * dx
	}
	if den == 0 {
		return 0
	}
	return num / den
}

// KSOneSampleStat menghitung statistik Kolmogorov-Smirnov D_n satu-sampel terhadap CDF kontinu.
func KSOneSampleStat(samples []float64, cdf func(float64) float64) float64 {
	n := len(samples)
	if n == 0 {
		return 0
	}
	cp := make([]float64, n)
	copy(cp, samples)
	sort.Float64s(cp)
	fn := float64(n)
	dMax := 0.0
	for i, x := range cp {
		fx := cdf(x)
		dPlus := float64(i+1)/fn - fx
		dMinus := fx - float64(i)/fn
		if dPlus > dMax {
			dMax = dPlus
		}
		if dMinus > dMax {
			dMax = dMinus
		}
	}
	return dMax
}

// KSSimulateNullCritical menyimulasikan distribusi nol statistik Kolmogorov-Smirnov D_n
// untuk ukuran sampel n tanpa tabel buku teks, memakai representasi jarak eksponensial eksak
// U_{(i)} = sum_{j<=i} E_j / sum_{j<=n+1} E_j dalam O(n) waktu dan O(1) memori per replikasi.
func KSSimulateNullCritical(n, nRep int, alpha float64, nChunks, nWorkers int, seed uint64) (critVal, empiricalP float64) {
	if n <= 0 || nRep <= 0 {
		return 0, 1
	}
	perChunk := nRep / nChunks
	rem := nRep % nChunks
	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) []float64 {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		out := make([]float64, cnt)
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x2201, uint64(c)))
		fn := float64(n)
		for k := 0; k < cnt; k++ {
			saved := *rg
			total := 0.0
			for j := 0; j <= n; j++ {
				total += rg.ExpFloat64()
			}
			// Putar ulang RNG yang sama untuk melacak statistik orde terurut tanpa pengurutan.
			rg2 := saved
			acc := 0.0
			dMax := 0.0
			for i := 0; i < n; i++ {
				acc += rg2.ExpFloat64()
				uOrd := acc / total
				dPlus := float64(i+1)/fn - uOrd
				dMinus := uOrd - float64(i)/fn
				if dPlus > dMax {
					dMax = dPlus
				}
				if dMinus > dMax {
					dMax = dMinus
				}
			}
			out[k] = dMax
		}
		return out
	})

	all := make([]float64, 0, nRep)
	for _, ch := range chunks {
		all = append(all, ch...)
	}
	critVal = Quantile(all, 1.0-alpha)
	return critVal, 0
}

// KSPValueFromNull menghitung p-value empiris terhadap sampel nol yang disimulasikan sekaligus mengembalikan nilai kritis.
func KSEvalAgainstSimulatedNull(obsStat float64, n, nRep int, alpha float64, nChunks, nWorkers int, seed uint64) (critVal, pValue float64) {
	if n <= 0 || nRep <= 0 {
		return 0, 1
	}
	perChunk := nRep / nChunks
	rem := nRep % nChunks
	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) []float64 {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		out := make([]float64, cnt)
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x2202, uint64(c)))
		fn := float64(n)
		for k := 0; k < cnt; k++ {
			saved := *rg
			total := 0.0
			for j := 0; j <= n; j++ {
				total += rg.ExpFloat64()
			}
			rg2 := saved
			acc := 0.0
			dMax := 0.0
			for i := 0; i < n; i++ {
				acc += rg2.ExpFloat64()
				uOrd := acc / total
				dPlus := float64(i+1)/fn - uOrd
				dMinus := uOrd - float64(i)/fn
				if dPlus > dMax {
					dMax = dPlus
				}
				if dMinus > dMax {
					dMax = dMinus
				}
			}
			out[k] = dMax
		}
		return out
	})

	all := make([]float64, 0, nRep)
	ge := 0
	for _, ch := range chunks {
		for _, v := range ch {
			all = append(all, v)
			if v >= obsStat {
				ge++
			}
		}
	}
	critVal = Quantile(all, 1.0-alpha)
	pValue = float64(ge) / float64(len(all))
	return critVal, pValue
}
