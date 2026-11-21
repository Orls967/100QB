package stats

import (
	"math"
	"sort"
)

// Z95 adalah kuantil normal standar untuk interval kepercayaan dua arah 95%.
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

// SE menghitung galat baku rata-rata sampel.
func SE(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	return SD(xs) / math.Sqrt(float64(len(xs)))
}

// RMSE menghitung akar rata-rata kuadrat galat terhadap nilai target.
func RMSE(est []float64, target float64) float64 {
	if len(est) == 0 {
		return 0
	}
	ss := 0.0
	for _, x := range est {
		d := x - target
		ss += d * d
	}
	return math.Sqrt(ss / float64(len(est)))
}

// Quantile menghitung kuantil empiris terurut pada level q di [0, 1] dengan interpolasi linier.
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

// WilsonCI menghitung proporsi empiris dan selang kepercayaan skor Wilson untuk proporsi binomial.
func WilsonCI(successes, total int, z float64) (pHat, lower, upper float64) {
	if total <= 0 {
		return 0, 0, 1
	}
	n := float64(total)
	pHat = float64(successes) / n
	z2 := z * z
	denom := 1.0 + z2/n
	center := (pHat + z2/(2.0*n)) / denom
	half := (z / denom) * math.Sqrt(pHat*(1.0-pHat)/n+z2/(4.0*n*n))
	lower = math.Max(0.0, center-half)
	upper = math.Min(1.0, center+half)
	return pHat, lower, upper
}

// CoverageWithValidity menghitung cakupan empiris beserta selang Wilson 95% dengan tetap
// memasukkan seluruh replikasi luar ke dalam denominator (interval tidak valid dihitung tidak tercakup).
func CoverageWithValidity(validFlags, coveredFlags []bool) (validCount, invalidCount, coveredCount, total int, rate, wilsonLo, wilsonHi float64) {
	total = len(coveredFlags)
	for i := 0; i < total; i++ {
		v := i < len(validFlags) && validFlags[i]
		if v {
			validCount++
			if coveredFlags[i] {
				coveredCount++
			}
		} else {
			invalidCount++
		}
	}
	rate, wilsonLo, wilsonHi = WilsonCI(coveredCount, total, Z95)
	return validCount, invalidCount, coveredCount, total, rate, wilsonLo, wilsonHi
}
