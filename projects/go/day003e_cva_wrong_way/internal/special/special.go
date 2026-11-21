package special

import (
	"math"
)

const invSqrt2Pi = 0.39894228040143267794

// NormPDF menghitung kerapatan probabilitas normal standar phi(x).
func NormPDF(x float64) float64 {
	return invSqrt2Pi * math.Exp(-0.5*x*x)
}

// NormCDF menghitung fungsi distribusi kumulatif normal standar Phi(x).
func NormCDF(x float64) float64 {
	return 0.5 * math.Erfc(-x/math.Sqrt2)
}

// NormQuantile menghitung kuantil normal standar Phi^{-1}(p) untuk p dalam (0, 1)
// menggunakan pendekatan rasional Acklam diikuti penyempurnaan Halley.
func NormQuantile(p float64) float64 {
	if p <= 0 {
		return math.Inf(-1)
	}
	if p >= 1 {
		return math.Inf(1)
	}
	const (
		a1 = -3.969683028665376e+01
		a2 = 2.209460984245205e+02
		a3 = -2.759285104469687e+02
		a4 = 1.383577518672690e+02
		a5 = -3.066479806614716e+01
		a6 = 2.506628277459239e+00

		b1 = -5.447609879822406e+01
		b2 = 1.615858368580409e+02
		b3 = -1.556989798598866e+02
		b4 = 6.680131188771972e+01
		b5 = -1.328068155288572e+01

		c1 = -7.784894002430293e-03
		c2 = -3.223964580411365e-01
		c3 = -2.400758277161838e+00
		c4 = -2.549732539343734e+00
		c5 = 4.374664141464968e+00
		c6 = 2.938163982698783e+00

		d1 = 7.784695709041462e-03
		d2 = 3.224671290700398e-01
		d3 = 2.445134137142996e+00
		d4 = 3.754408661907416e+00

		pLow  = 0.02425
		pHigh = 1.0 - pLow
	)

	var x float64
	switch {
	case p < pLow:
		q := math.Sqrt(-2.0 * math.Log(p))
		x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q + c6) /
			((((d1*q+d2)*q+d3)*q+d4)*q + 1.0)
	case p <= pHigh:
		q := p - 0.5
		r := q * q
		x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r + a6) * q /
			(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r + 1.0)
	default:
		q := math.Sqrt(-2.0 * math.Log(1.0-p))
		x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q + c6) /
			((((d1*q+d2)*q+d3)*q+d4)*q + 1.0)
	}

	// Satu langkah Halley untuk akurasi penuh float64.
	e := NormCDF(x) - p
	u := e / NormPDF(x)
	x -= u / (1.0 + 0.5*x*u)
	return x
}

// GHNode merepresentasikan simpul dan bobot kuadratur Gauss-Hermite terukur untuk Z ~ N(0, 1),
// sehingga integral f(z) phi(z) dz ~ sum_i w_i f(z_i).
type GHNode struct {
	X float64
	W float64
}

// GaussHermiteStandardNormal menghitung n simpul dan bobot kuadratur Gauss-Hermite
// untuk ekspektasi di bawah N(0, 1) menggunakan algoritma Golub-Welsch (tql2 pada matriks Jacobi).
func GaussHermiteStandardNormal(n int) []GHNode {
	if n <= 0 {
		return nil
	}
	if n == 1 {
		return []GHNode{{X: 0, W: 1}}
	}
	d := make([]float64, n)
	e := make([]float64, n)
	for i := 0; i < n-1; i++ {
		e[i] = math.Sqrt(float64(i + 1))
	}
	return golubWelsch(d, e, 1.0)
}

// GaussLegendre menghitung n simpul dan bobot kuadratur Gauss-Legendre pada selang [-1, 1].
func GaussLegendre(n int) []GHNode {
	if n <= 0 {
		return nil
	}
	if n == 1 {
		return []GHNode{{X: 0, W: 2}}
	}
	d := make([]float64, n)
	e := make([]float64, n)
	for i := 0; i < n-1; i++ {
		k := float64(i + 1)
		e[i] = k / math.Sqrt((2.0*k-1.0)*(2.0*k+1.0))
	}
	return golubWelsch(d, e, 2.0)
}

// GaussLegendreInterval memetakan simpul Gauss-Legendre ke selang [a, b].
func GaussLegendreInterval(n int, a, b float64) []GHNode {
	base := GaussLegendre(n)
	mid := 0.5 * (a + b)
	half := 0.5 * (b - a)
	for i := range base {
		base[i].X = mid + half*base[i].X
		base[i].W *= half
	}
	return base
}

func golubWelsch(d, e []float64, weightSum float64) []GHNode {
	n := len(d)
	v0 := make([]float64, n)
	v0[0] = 1.0

	// Algoritma QL implisit simetris tridiagonal (EISPACK tql2 terfokus pada baris pertama).
	for l := 0; l < n; l++ {
		iter := 0
		for {
			m := l
			for m < n-1 {
				tst := math.Abs(d[m]) + math.Abs(d[m+1])
				if tst+math.Abs(e[m]) == tst {
					break
				}
				m++
			}
			if m == l {
				break
			}
			iter++
			if iter > 60 {
				break
			}
			g := (d[l+1] - d[l]) / (2.0 * e[l])
			r := math.Hypot(g, 1.0)
			if g < 0 {
				g = d[m] - d[l] + e[l]/(g-r)
			} else {
				g = d[m] - d[l] + e[l]/(g+r)
			}
			s := 1.0
			c := 1.0
			p := 0.0
			for i := m - 1; i >= l; i-- {
				f := s * e[i]
				b := c * e[i]
				if math.Abs(f) >= math.Abs(g) {
					c = g / f
					r = math.Sqrt(c*c + 1.0)
					e[i+1] = f * r
					s = 1.0 / r
					c *= s
				} else {
					s = f / g
					r = math.Sqrt(s*s + 1.0)
					e[i+1] = g * r
					c = 1.0 / r
					s *= c
				}
				g = d[i+1] - p
				r = (d[i]-g)*s + 2.0*c*b
				p = s * r
				d[i+1] = g + p
				g = c*r - b

				vi := v0[i]
				vip1 := v0[i+1]
				v0[i+1] = s*vi + c*vip1
				v0[i] = c*vi - s*vip1
			}
			d[l] -= p
			e[l] = g
			e[m] = 0.0
		}
	}

	for i := 0; i < n-1; i++ {
		minIdx := i
		for j := i + 1; j < n; j++ {
			if d[j] < d[minIdx] {
				minIdx = j
			}
		}
		if minIdx != i {
			d[i], d[minIdx] = d[minIdx], d[i]
			v0[i], v0[minIdx] = v0[minIdx], v0[i]
		}
	}

	nodes := make([]GHNode, n)
	wRaw := 0.0
	for i := 0; i < n; i++ {
		w := v0[i] * v0[i]
		nodes[i] = GHNode{X: d[i], W: w}
		wRaw += w
	}
	scale := weightSum / wRaw
	for i := 0; i < n/2; i++ {
		j := n - 1 - i
		xAvg := 0.5 * (nodes[j].X - nodes[i].X)
		wAvg := 0.5 * (nodes[i].W + nodes[j].W) * scale
		nodes[i] = GHNode{X: -xAvg, W: wAvg}
		nodes[j] = GHNode{X: xAvg, W: wAvg}
	}
	if n%2 == 1 {
		mid := n / 2
		nodes[mid].X = 0
		nodes[mid].W *= scale
	}
	return nodes
}
