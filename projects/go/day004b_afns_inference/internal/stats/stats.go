package stats

import (
	"math"
	"sort"
)

// Z95 adalah kuantil z standar untuk selang kepercayaan dua-arah 95%.
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

// RMSE menghitung akar rata-rata kuadrat deviasi terhadap nilai target atau antar dua vektor.
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

// Quantile menghitung kuantil empiris level q in [0, 1] secara terurut.
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

// WilsonCI menghitung proporsi empiris dan selang kepercayaan skor Wilson untuk data binomial.
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

// LinearTrend merangkum regresi linier OLS dan WLS y = Intercept + Slope * x beserta galat bakunya.
type LinearTrend struct {
	OLSIntercept        float64
	OLSSlope            float64
	OLSSlopeMCSE        float64 // propagasi analitik anggap galat independen: sqrt(sum w_i^2 se_i^2)
	OLSSlopeResidualSE  float64 // SE OLS konvensional dari residual regresi
	WLSIntercept        float64
	WLSSlope            float64
	WLSSlopeMCSE        float64 // propagasi analitik WLS anggap independen: 1 / sqrt(sum tau_i (x_i - xbar_w)^2)
	WLSWeightsAvailable bool
	OLSWeightsX         []float64 // koefisien w_i sehingga OLSSlope = sum_i w_i * y_i
	WLSWeightsX         []float64 // koefisien w_i sehingga WLSSlope = sum_i w_i * y_i
}

// FitLinearTrend menghitung regresi linier OLS dan WLS y = a + b*x.
// Untuk WLS, bobot adalah 1 / se_i^2 bila seluruh se_i > 0.
func FitLinearTrend(xs, ys, se []float64) LinearTrend {
	n := len(xs)
	var res LinearTrend
	if n < 2 || len(ys) != n {
		return res
	}
	mx := Mean(xs)
	my := Mean(ys)
	sxx, sxy := 0.0, 0.0
	for i := 0; i < n; i++ {
		dx := xs[i] - mx
		sxx += dx * dx
		sxy += dx * (ys[i] - my)
	}
	if sxx > 0 {
		res.OLSSlope = sxy / sxx
		res.OLSIntercept = my - res.OLSSlope*mx
		res.OLSWeightsX = make([]float64, n)
		for i := 0; i < n; i++ {
			res.OLSWeightsX[i] = (xs[i] - mx) / sxx
		}
		if len(se) == n {
			varProp := 0.0
			for i := 0; i < n; i++ {
				varProp += res.OLSWeightsX[i] * res.OLSWeightsX[i] * se[i] * se[i]
			}
			res.OLSSlopeMCSE = math.Sqrt(varProp)
		}
		if n > 2 {
			rss := 0.0
			for i := 0; i < n; i++ {
				fit := res.OLSIntercept + res.OLSSlope*xs[i]
				e := ys[i] - fit
				rss += e * e
			}
			s2 := rss / float64(n-2)
			res.OLSSlopeResidualSE = math.Sqrt(s2 / sxx)
		}
	}

	canWLS := len(se) == n
	if canWLS {
		for _, s := range se {
			if s <= 0 || math.IsNaN(s) || math.IsInf(s, 0) {
				canWLS = false
				break
			}
		}
	}
	if canWLS {
		tau := make([]float64, n)
		sumTau, sumTauX, sumTauY := 0.0, 0.0, 0.0
		for i := 0; i < n; i++ {
			tau[i] = 1.0 / (se[i] * se[i])
			sumTau += tau[i]
			sumTauX += tau[i] * xs[i]
			sumTauY += tau[i] * ys[i]
		}
		xBarW := sumTauX / sumTau
		yBarW := sumTauY / sumTau
		sxxW, sxyW := 0.0, 0.0
		for i := 0; i < n; i++ {
			dx := xs[i] - xBarW
			sxxW += tau[i] * dx * dx
			sxyW += tau[i] * dx * (ys[i] - yBarW)
		}
		if sxxW > 0 {
			res.WLSWeightsAvailable = true
			res.WLSSlope = sxyW / sxxW
			res.WLSIntercept = yBarW - res.WLSSlope*xBarW
			res.WLSSlopeMCSE = 1.0 / math.Sqrt(sxxW)
			res.WLSWeightsX = make([]float64, n)
			for i := 0; i < n; i++ {
				res.WLSWeightsX[i] = tau[i] * (xs[i] - xBarW) / sxxW
			}
		}
	}
	if !res.WLSWeightsAvailable {
		res.WLSIntercept = res.OLSIntercept
		res.WLSSlope = res.OLSSlope
		res.WLSSlopeMCSE = res.OLSSlopeMCSE
		res.WLSWeightsX = res.OLSWeightsX
	}
	return res
}
