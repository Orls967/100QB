package kalman

import (
	"fmt"
	"math"
)

func zerosMat(rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for i := range m {
		m[i] = make([]float64, cols)
	}
	return m
}

func dotVec(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

func choleskyMat(a [][]float64) ([][]float64, error) {
	n := len(a)
	l := zerosMat(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			sum := a[i][j]
			for k := 0; k < j; k++ {
				sum -= l[i][k] * l[j][k]
			}
			if i == j {
				if sum <= 0 || math.IsNaN(sum) || math.IsInf(sum, 0) {
					return nil, fmt.Errorf("matrix not positive definite at row %d: diag=%g", i, sum)
				}
				l[i][j] = math.Sqrt(sum)
			} else {
				l[i][j] = sum / l[j][j]
			}
		}
	}
	return l, nil
}

func logDetCholesky(l [][]float64) float64 {
	s := 0.0
	for i := range l {
		s += math.Log(l[i][i])
	}
	return 2.0 * s
}

func solveCholeskyLower(l [][]float64, b []float64) []float64 {
	n := len(l)
	z := make([]float64, n)
	for i := 0; i < n; i++ {
		sum := b[i]
		for j := 0; j < i; j++ {
			sum -= l[i][j] * z[j]
		}
		z[i] = sum / l[i][i]
	}
	return z
}

func solveCholesky(l [][]float64, b []float64) []float64 {
	z := solveCholeskyLower(l, b)
	n := len(l)
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := z[i]
		for j := i + 1; j < n; j++ {
			sum -= l[j][i] * x[j]
		}
		x[i] = sum / l[i][i]
	}
	return x
}

func invertCholesky(l [][]float64) [][]float64 {
	n := len(l)
	inv := zerosMat(n, n)
	e := make([]float64, n)
	for col := 0; col < n; col++ {
		e[col] = 1.0
		sol := solveCholesky(l, e)
		for row := 0; row < n; row++ {
			inv[row][col] = sol[row]
		}
		e[col] = 0.0
	}
	return inv
}

func invertLU(a [][]float64) ([][]float64, error) {
	n := len(a)
	aug := zerosMat(n, 2*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			aug[i][j] = a[i][j]
		}
		aug[i][n+i] = 1.0
	}
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(aug[col][col])
		for row := col + 1; row < n; row++ {
			if v := math.Abs(aug[row][col]); v > best {
				best = v
				piv = row
			}
		}
		if best < 1e-15 {
			return nil, fmt.Errorf("singular matrix at pivot %d", col)
		}
		if piv != col {
			aug[col], aug[piv] = aug[piv], aug[col]
		}
		diag := aug[col][col]
		for j := col; j < 2*n; j++ {
			aug[col][j] /= diag
		}
		for row := 0; row < n; row++ {
			if row == col {
				continue
			}
			fac := aug[row][col]
			if fac == 0 {
				continue
			}
			for j := col; j < 2*n; j++ {
				aug[row][j] -= fac * aug[col][j]
			}
		}
	}
	inv := zerosMat(n, n)
	for i := 0; i < n; i++ {
		copy(inv[i], aug[i][n:])
	}
	return inv, nil
}

func transposeMat(a [][]float64) [][]float64 {
	if len(a) == 0 {
		return nil
	}
	rows, cols := len(a), len(a[0])
	t := zerosMat(cols, rows)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			t[j][i] = a[i][j]
		}
	}
	return t
}

func matMul(a, b [][]float64) [][]float64 {
	r, k, c := len(a), len(b), len(b[0])
	out := zerosMat(r, c)
	for i := 0; i < r; i++ {
		for t := 0; t < k; t++ {
			av := a[i][t]
			if av == 0 {
				continue
			}
			for j := 0; j < c; j++ {
				out[i][j] += av * b[t][j]
			}
		}
	}
	return out
}

func matSub(a, b [][]float64) [][]float64 {
	out := zerosMat(len(a), len(a[0]))
	for i := range a {
		for j := range a[i] {
			out[i][j] = a[i][j] - b[i][j]
		}
	}
	return out
}

// symmetricEigenJacobi menguraikan matriks simetris n x n menjadi nilai eigen terurut menurun
// dan eigenvektor (disimpan sebagai baris vecs[k]).
func symmetricEigenJacobi(a [][]float64, maxSweeps int, tol float64) (vals []float64, vecs [][]float64, sweeps int, err error) {
	n := len(a)
	m := zerosMat(n, n)
	v := zerosMat(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			m[i][j] = 0.5 * (a[i][j] + a[j][i])
		}
		v[i][i] = 1.0
	}
	for sweep := 1; sweep <= maxSweeps; sweep++ {
		maxOff := 0.0
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if av := math.Abs(m[p][q]); av > maxOff {
					maxOff = av
				}
			}
		}
		if maxOff <= tol {
			sweeps = sweep
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				g := m[p][q]
				if math.Abs(g) <= tol {
					continue
				}
				tau := (m[q][q] - m[p][p]) / (2.0 * g)
				var t float64
				if tau >= 0 {
					t = 1.0 / (tau + math.Sqrt(1.0+tau*tau))
				} else {
					t = -1.0 / (-tau + math.Sqrt(1.0+tau*tau))
				}
				c := 1.0 / math.Sqrt(1.0+t*t)
				s := t * c
				app, aqq, apq := m[p][p], m[q][q], m[p][q]
				m[p][p] = c*c*app - 2*s*c*apq + s*s*aqq
				m[q][q] = s*s*app + 2*s*c*apq + c*c*aqq
				m[p][q], m[q][p] = 0, 0
				for r := 0; r < n; r++ {
					if r == p || r == q {
						continue
					}
					mrp, mrq := m[r][p], m[r][q]
					m[r][p] = c*mrp - s*mrq
					m[p][r] = m[r][p]
					m[r][q] = s*mrp + c*mrq
					m[q][r] = m[r][q]
				}
				for r := 0; r < n; r++ {
					vrp, vrq := v[r][p], v[r][q]
					v[r][p] = c*vrp - s*vrq
					v[r][q] = s*vrp + c*vrq
				}
			}
		}
		sweeps = sweep
	}

	vals = make([]float64, n)
	vecs = zerosMat(n, n)
	idx := make([]int, n)
	for i := 0; i < n; i++ {
		vals[i] = m[i][i]
		idx[i] = i
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if vals[j] > vals[i] {
				vals[i], vals[j] = vals[j], vals[i]
				idx[i], idx[j] = idx[j], idx[i]
			}
		}
	}
	for k := 0; k < n; k++ {
		col := idx[k]
		for row := 0; row < n; row++ {
			vecs[k][row] = v[row][col]
		}
	}
	return vals, vecs, sweeps, nil
}
