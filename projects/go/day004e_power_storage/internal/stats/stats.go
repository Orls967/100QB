package stats

import (
	"math"
	"sort"
)

// Z95 adalah kuantil normal standar dua sisi untuk tingkat kepercayaan 95%.
const Z95 = 1.959963984540054

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

// SE menghitung galat Monte Carlo rata-rata sampel (SD / sqrt(n)).
func SE(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	return SD(xs) / math.Sqrt(float64(len(xs)))
}

// Quantile menghitung kuantil empiris pada level q di [0, 1] dengan interpolasi linier.
func Quantile(xs []float64, q float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
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

// VaRAndES menghitung VaR ekor bawah pada level alpha (mis. 0.05 = persentil ke-5 profit)
// dan Expected Shortfall (ES) sebagai rata-rata k = ceil(alpha * n) pengamatan terburuk (ekor bawah).
func VaRAndES(xs []float64, alpha float64) (varVal, esVal float64) {
	n := len(xs)
	if n == 0 {
		return math.NaN(), math.NaN()
	}
	cp := make([]float64, n)
	copy(cp, xs)
	sort.Float64s(cp)
	if alpha <= 0 {
		return cp[0], cp[0]
	}
	if alpha >= 1 {
		return cp[n-1], Mean(cp)
	}
	pos := alpha * float64(n-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	w := pos - float64(lo)
	varVal = (1.0-w)*cp[lo] + w*cp[hi]

	k := int(math.Ceil(alpha * float64(n)))
	if k < 1 {
		k = 1
	}
	if k > n {
		k = n
	}
	sum := 0.0
	for i := 0; i < k; i++ {
		sum += cp[i]
	}
	esVal = sum / float64(k)
	return varVal, esVal
}

// PairedDiffSummary menyimpan statistik selisih berpasangan (a_i - b_i) dan interval Monte Carlo 95%.
type PairedDiffSummary struct {
	N            int
	MeanDiff     float64
	SDDiff       float64
	MCSE         float64
	CI95Low      float64
	CI95High     float64
	PositiveFrac float64
	ZStat        float64
}

// PairedDiffCI menghitung rata-rata selisih berpasangan, MCSE, dan interval kepercayaan 95% (normal Z95).
func PairedDiffCI(a, b []float64) PairedDiffSummary {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if n == 0 {
		return PairedDiffSummary{}
	}
	diffs := make([]float64, n)
	pos := 0
	for i := 0; i < n; i++ {
		diffs[i] = a[i] - b[i]
		if diffs[i] > 0 {
			pos++
		}
	}
	m := Mean(diffs)
	sd := SD(diffs)
	se := sd / math.Sqrt(float64(n))
	z := 0.0
	if se > 0 {
		z = m / se
	}
	return PairedDiffSummary{
		N:            n,
		MeanDiff:     m,
		SDDiff:       sd,
		MCSE:         se,
		CI95Low:      m - Z95*se,
		CI95High:     m + Z95*se,
		PositiveFrac: float64(pos) / float64(n),
		ZStat:        z,
	}
}
