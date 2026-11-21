package special

import (
	"fmt"
	"math"
)

const invSqrt2Pi = 0.39894228040143267794

// NormPDF menghitung kerapatan probabilitas normal standar phi(x).
func NormPDF(x float64) float64 {
	return invSqrt2Pi * math.Exp(-0.5*x*x)
}

// NormCDF menghitung CDF normal standar Phi(x).
func NormCDF(x float64) float64 {
	return 0.5 * math.Erfc(-x/math.Sqrt2)
}

// NormQuantile menghitung kuantil normal standar Phi^{-1}(p) via pendekatan Acklam + Halley.
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

	e := NormCDF(x) - p
	u := e / NormPDF(x)
	x -= u / (1.0 + 0.5*x*u)
	return x
}

// BrentRoot mencari akar f(x) = 0 di selang pengapit [a, b] memakai metode Brent-Dekker
// hingga toleransi tol.
func BrentRoot(f func(float64) float64, a, b, tol float64) (float64, error) {
	fa := f(a)
	fb := f(b)
	if math.IsNaN(fa) || math.IsNaN(fb) {
		return 0, fmt.Errorf("BrentRoot: NaN endpoint evaluation at [%g, %g]", a, b)
	}
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if (fa > 0 && fb > 0) || (fa < 0 && fb < 0) {
		return 0, fmt.Errorf("BrentRoot: root not bracketed on [%g, %g] (fa=%g, fb=%g)", a, b, fa, fb)
	}

	c := a
	fc := fa
	d := b - a
	e := d
	const epsMach = 2.220446049250313e-16

	for iter := 0; iter < 150; iter++ {
		if (fb > 0 && fc > 0) || (fb < 0 && fc < 0) {
			c = a
			fc = fa
			d = b - a
			e = d
		}
		if math.Abs(fc) < math.Abs(fb) {
			a = b
			b = c
			c = a
			fa = fb
			fb = fc
			fc = fa
		}

		tol1 := 2.0*epsMach*math.Abs(b) + 0.5*tol
		xm := 0.5 * (c - b)
		if math.Abs(xm) <= tol1 || fb == 0 {
			return b, nil
		}

		if math.Abs(e) >= tol1 && math.Abs(fa) > math.Abs(fb) {
			s := fb / fa
			var p, q float64
			if a == c {
				p = 2.0 * xm * s
				q = 1.0 - s
			} else {
				q = fa / fc
				r := fb / fc
				p = s * (2.0*xm*q*(q-r) - (b-a)*(r-1.0))
				q = (q - 1.0) * (r - 1.0) * (s - 1.0)
			}
			if p > 0 {
				q = -q
			} else {
				p = -p
			}
			if 2.0*p < math.Min(3.0*xm*q-math.Abs(tol1*q), math.Abs(e*q)) {
				e = d
				d = p / q
			} else {
				d = xm
				e = d
			}
		} else {
			d = xm
			e = d
		}

		a = b
		fa = fb
		if math.Abs(d) > tol1 {
			b += d
		} else if xm > 0 {
			b += tol1
		} else {
			b -= tol1
		}
		fb = f(b)
	}
	return b, nil
}

// GoldenSectionMin meminimalkan fungsi unimodal 1D f(x) di [a, b] sampai lebar tol.
func GoldenSectionMin(f func(float64) float64, a, b, tol float64) (xMin, fMin float64) {
	const invPhi = 0.6180339887498949
	x1 := b - invPhi*(b-a)
	x2 := a + invPhi*(b-a)
	f1 := f(x1)
	f2 := f(x2)
	for iter := 0; iter < 150 && (b-a) > tol; iter++ {
		if f1 < f2 {
			b = x2
			x2 = x1
			f2 = f1
			x1 = b - invPhi*(b-a)
			f1 = f(x1)
		} else {
			a = x1
			x1 = x2
			f1 = f2
			x2 = a + invPhi*(b-a)
			f2 = f(x2)
		}
	}
	if f1 < f2 {
		return x1, f1
	}
	return x2, f2
}

// SolveLinearLU memecahkan sistem linear Ax = b memakai dekomposisi LU pivoting parsial.
func SolveLinearLU(a [][]float64, b []float64) ([]float64, error) {
	n := len(b)
	if n == 0 || len(a) != n {
		return nil, fmt.Errorf("SolveLinearLU: invalid dimensions")
	}
	mat := make([][]float64, n)
	rhs := make([]float64, n)
	copy(rhs, b)
	for i := 0; i < n; i++ {
		if len(a[i]) != n {
			return nil, fmt.Errorf("SolveLinearLU: non-square matrix")
		}
		mat[i] = make([]float64, n)
		copy(mat[i], a[i])
	}

	for k := 0; k < n; k++ {
		piv := k
		maxVal := math.Abs(mat[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(mat[i][k]); v > maxVal {
				maxVal = v
				piv = i
			}
		}
		if maxVal < 1e-15 {
			return nil, fmt.Errorf("SolveLinearLU: singular matrix at row %d", k)
		}
		if piv != k {
			mat[k], mat[piv] = mat[piv], mat[k]
			rhs[k], rhs[piv] = rhs[piv], rhs[k]
		}
		for i := k + 1; i < n; i++ {
			factor := mat[i][k] / mat[k][k]
			mat[i][k] = factor
			for j := k + 1; j < n; j++ {
				mat[i][j] -= factor * mat[k][j]
			}
			rhs[i] -= factor * rhs[k]
		}
	}

	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := rhs[i]
		for j := i + 1; j < n; j++ {
			sum -= mat[i][j] * x[j]
		}
		x[i] = sum / mat[i][i]
	}
	return x, nil
}
