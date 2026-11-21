package stats

import (
	"math"
	"sort"
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

// Quantile menghitung kuantil empiris level q dalam [0, 1] dengan interpolasi linier pada salinan terurut.
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
