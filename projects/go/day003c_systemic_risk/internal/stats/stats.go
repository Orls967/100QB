package stats

import (
	"math"
	"sort"
)

// Mean menghitung rata-rata aritmetik.
func Mean(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	sum := 0.0
	for _, v := range xs {
		sum += v
	}
	return sum / float64(len(xs))
}

// Variance menghitung variansi sampel tak-bias (ddof=1).
func Variance(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return 0
	}
	m := Mean(xs)
	ss := 0.0
	for _, v := range xs {
		d := v - m
		ss += d * d
	}
	return ss / float64(n-1)
}

// SE menghitung standard error dari rata-rata sampel.
func SE(xs []float64) float64 {
	if len(xs) < 2 {
		return 0
	}
	return math.Sqrt(Variance(xs) / float64(len(xs)))
}

// VaRAndES menghitung Value-at-Risk pada level alphaVaR (mis. 0.99) dan Expected Shortfall pada alphaES (mis. 0.975).
func VaRAndES(xs []float64, alphaVaR, alphaES float64) (varVal, esVal float64) {
	n := len(xs)
	if n == 0 {
		return 0, 0
	}
	cp := make([]float64, n)
	copy(cp, xs)
	sort.Float64s(cp)

	idxVaR := int(math.Ceil(alphaVaR*float64(n))) - 1
	if idxVaR < 0 {
		idxVaR = 0
	}
	if idxVaR >= n {
		idxVaR = n - 1
	}
	varVal = cp[idxVaR]

	idxES := int(math.Floor(alphaES * float64(n)))
	if idxES < 0 {
		idxES = 0
	}
	if idxES >= n {
		idxES = n - 1
	}
	sum := 0.0
	cnt := 0
	for i := idxES; i < n; i++ {
		sum += cp[i]
		cnt++
	}
	esVal = sum / float64(cnt)
	return varVal, esVal
}
