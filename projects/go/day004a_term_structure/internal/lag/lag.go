package lag

import (
	"fmt"
	"math"
)

// Zeros membuat matriks berukuran rows x cols berisi nol.
func Zeros(rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		m[i] = make([]float64, cols)
	}
	return m
}

// Eye membuat matriks identitas n x n.
func Eye(n int) [][]float64 {
	m := Zeros(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = 1.0
	}
	return m
}

// Clone menyalin matriks 2D agar perubahan tidak memutasi sumber.
func Clone(a [][]float64) [][]float64 {
	m := make([][]float64, len(a))
	for i := range a {
		m[i] = make([]float64, len(a[i]))
		copy(m[i], a[i])
	}
	return m
}

// Transpose mengembalikan transpos A^T.
func Transpose(a [][]float64) [][]float64 {
	if len(a) == 0 {
		return nil
	}
	r, c := len(a), len(a[0])
	t := Zeros(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			t[j][i] = a[i][j]
		}
	}
	return t
}

// MatMul mengalikan dua matriks A (r x k) dan B (k x c).
func MatMul(a, b [][]float64) [][]float64 {
	r := len(a)
	if r == 0 || len(b) == 0 {
		return nil
	}
	k := len(a[0])
	c := len(b[0])
	out := Zeros(r, c)
	for i := 0; i < r; i++ {
		for p := 0; p < k; p++ {
			aip := a[i][p]
			if aip == 0 {
				continue
			}
			for j := 0; j < c; j++ {
				out[i][j] += aip * b[p][j]
			}
		}
	}
	return out
}

// MatVec mengalikan matriks A dengan vektor x.
func MatVec(a [][]float64, x []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		s := 0.0
		for j := range x {
			s += a[i][j] * x[j]
		}
		out[i] = s
	}
	return out
}

// Dot menghitung hasil kali titik dua vektor.
func Dot(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// Add mengembalikan A + B.
func Add(a, b [][]float64) [][]float64 {
	r := len(a)
	out := Zeros(r, len(a[0]))
	for i := 0; i < r; i++ {
		for j := range a[i] {
			out[i][j] = a[i][j] + b[i][j]
		}
	}
	return out
}

// Sub mengembalikan A - B.
func Sub(a, b [][]float64) [][]float64 {
	r := len(a)
	out := Zeros(r, len(a[0]))
	for i := 0; i < r; i++ {
		for j := range a[i] {
			out[i][j] = a[i][j] - b[i][j]
		}
	}
	return out
}

// Scale mengalikan setiap entri matriks dengan skalar s.
func Scale(a [][]float64, s float64) [][]float64 {
	r := len(a)
	out := Zeros(r, len(a[0]))
	for i := 0; i < r; i++ {
		for j := range a[i] {
			out[i][j] = a[i][j] * s
		}
	}
	return out
}

// SolveLU memecahkan Ax = b memakai eliminasi LU dengan pivoting parsial.
func SolveLU(a [][]float64, b []float64) ([]float64, error) {
	n := len(b)
	if n == 0 || len(a) != n {
		return nil, fmt.Errorf("SolveLU: invalid dimensions")
	}
	mat := Clone(a)
	rhs := make([]float64, n)
	copy(rhs, b)

	for k := 0; k < n; k++ {
		if len(mat[k]) != n {
			return nil, fmt.Errorf("SolveLU: non-square matrix")
		}
		piv := k
		maxVal := math.Abs(mat[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(mat[i][k]); v > maxVal {
				maxVal = v
				piv = i
			}
		}
		if maxVal < 1e-15 {
			return nil, fmt.Errorf("SolveLU: singular matrix at row %d", k)
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

// InvertLU menghitung invers matriks A memakai eliminasi LU pivoting parsial.
func InvertLU(a [][]float64) ([][]float64, error) {
	n := len(a)
	if n == 0 {
		return nil, fmt.Errorf("InvertLU: empty matrix")
	}
	inv := Zeros(n, n)
	for j := 0; j < n; j++ {
		e := make([]float64, n)
		e[j] = 1.0
		col, err := SolveLU(a, e)
		if err != nil {
			return nil, err
		}
		for i := 0; i < n; i++ {
			inv[i][j] = col[i]
		}
	}
	return inv, nil
}

// Cholesky menghitung faktor segitiga bawah L sehingga A = L L^T untuk matriks simetris PD.
func Cholesky(a [][]float64) ([][]float64, error) {
	n := len(a)
	if n == 0 {
		return nil, fmt.Errorf("Cholesky: empty matrix")
	}
	l := Zeros(n, n)
	for i := 0; i < n; i++ {
		if len(a[i]) != n {
			return nil, fmt.Errorf("Cholesky: non-square matrix")
		}
		for j := 0; j <= i; j++ {
			sum := a[i][j]
			for k := 0; k < j; k++ {
				sum -= l[i][k] * l[j][k]
			}
			if i == j {
				if sum <= 0 || math.IsNaN(sum) || math.IsInf(sum, 0) {
					return nil, fmt.Errorf("Cholesky: matrix not positive definite at index %d (pivot=%g)", i, sum)
				}
				l[i][j] = math.Sqrt(sum)
			} else {
				l[i][j] = sum / l[j][j]
			}
		}
	}
	return l, nil
}

// SolveCholeskyLower memecahkan L y = b dengan substitusi maju.
func SolveCholeskyLower(l [][]float64, b []float64) []float64 {
	n := len(b)
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		sum := b[i]
		for j := 0; j < i; j++ {
			sum -= l[i][j] * y[j]
		}
		y[i] = sum / l[i][i]
	}
	return y
}

// SolveCholesky memecahkan (L L^T) x = b melalui substitusi maju dan mundur.
func SolveCholesky(l [][]float64, b []float64) []float64 {
	n := len(b)
	y := SolveCholeskyLower(l, b)
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := y[i]
		for j := i + 1; j < n; j++ {
			sum -= l[j][i] * x[j]
		}
		x[i] = sum / l[i][i]
	}
	return x
}

// LogDetCholesky mengembalikan log(|A|) = 2 * sum_i log(L_ii).
func LogDetCholesky(l [][]float64) float64 {
	s := 0.0
	for i := range l {
		s += math.Log(l[i][i])
	}
	return 2.0 * s
}

// InvertCholesky menghitung A^{-1} = (L L^T)^{-1} dari faktor Cholesky L.
func InvertCholesky(l [][]float64) [][]float64 {
	n := len(l)
	inv := Zeros(n, n)
	e := make([]float64, n)
	for j := 0; j < n; j++ {
		if j > 0 {
			e[j-1] = 0
		}
		e[j] = 1.0
		col := SolveCholesky(l, e)
		for i := 0; i < n; i++ {
			inv[i][j] = col[i]
		}
	}
	return inv
}

// SolveTridiagonal memecahkan sistem tridiagonal dalam O(n) via algoritma Thomas.
// sub[i] untuk baris i+1 kolom i (panjang n-1), diag panjang n, sup panjang n-1.
func SolveTridiagonal(sub, diag, sup, rhs []float64) ([]float64, error) {
	n := len(diag)
	if n == 0 || len(rhs) != n || len(sub) != n-1 || len(sup) != n-1 {
		return nil, fmt.Errorf("SolveTridiagonal: invalid dimensions")
	}
	c := make([]float64, n-1)
	d := make([]float64, n)
	if math.Abs(diag[0]) < 1e-15 {
		return nil, fmt.Errorf("SolveTridiagonal: zero pivot at 0")
	}
	if n == 1 {
		return []float64{rhs[0] / diag[0]}, nil
	}
	c[0] = sup[0] / diag[0]
	d[0] = rhs[0] / diag[0]
	for i := 1; i < n; i++ {
		denom := diag[i] - sub[i-1]*c[i-1]
		if math.Abs(denom) < 1e-15 {
			return nil, fmt.Errorf("SolveTridiagonal: zero pivot at %d", i)
		}
		if i < n-1 {
			c[i] = sup[i] / denom
		}
		d[i] = (rhs[i] - sub[i-1]*d[i-1]) / denom
	}
	x := make([]float64, n)
	x[n-1] = d[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = d[i] - c[i]*x[i+1]
	}
	return x, nil
}

// SymmetricEigenJacobi menghitung nilai eigen dan vektor eigen matriks simetris
// memakai algoritma Jacobi siklik. Mengembalikan vals terurut menurun dan vecs[k]
// sebagai vektor eigen satuan ke-k.
func SymmetricEigenJacobi(a [][]float64, maxSweeps int, tol float64) (vals []float64, vecs [][]float64, sweeps int, err error) {
	n := len(a)
	if n == 0 {
		return nil, nil, 0, fmt.Errorf("SymmetricEigenJacobi: empty matrix")
	}
	if maxSweeps <= 0 {
		maxSweeps = 100
	}
	if tol <= 0 {
		tol = 1e-14
	}
	m := Clone(a)
	for i := 0; i < n; i++ {
		if len(m[i]) != n {
			return nil, nil, 0, fmt.Errorf("SymmetricEigenJacobi: non-square matrix")
		}
		for j := i + 1; j < n; j++ {
			avg := 0.5 * (m[i][j] + m[j][i])
			m[i][j] = avg
			m[j][i] = avg
		}
	}
	v := Eye(n)
	normScale := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			normScale += m[i][j] * m[i][j]
		}
	}
	normScale = math.Sqrt(normScale)
	if normScale == 0 {
		normScale = 1.0
	}

	for sw := 0; sw < maxSweeps; sw++ {
		sweeps = sw + 1
		offSq := 0.0
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				offSq += m[p][q] * m[p][q]
			}
		}
		if math.Sqrt(2.0*offSq) <= tol*normScale {
			break
		}
		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				apq := m[p][q]
				if math.Abs(apq) <= 1e-16*normScale {
					continue
				}
				tau := (m[q][q] - m[p][p]) / (2.0 * apq)
				var t float64
				if tau >= 0 {
					t = 1.0 / (tau + math.Sqrt(1.0+tau*tau))
				} else {
					t = -1.0 / (-tau + math.Sqrt(1.0+tau*tau))
				}
				c := 1.0 / math.Sqrt(1.0+t*t)
				s := t * c
				app := m[p][p]
				aqq := m[q][q]
				m[p][p] = app - t*apq
				m[q][q] = aqq + t*apq
				m[p][q] = 0
				m[q][p] = 0

				for r := 0; r < n; r++ {
					if r == p || r == q {
						continue
					}
					mrp := m[r][p]
					mrq := m[r][q]
					nrp := c*mrp - s*mrq
					nrq := s*mrp + c*mrq
					m[r][p], m[p][r] = nrp, nrp
					m[r][q], m[q][r] = nrq, nrq
				}
				for r := 0; r < n; r++ {
					vrp := v[r][p]
					vrq := v[r][q]
					v[r][p] = c*vrp - s*vrq
					v[r][q] = s*vrp + c*vrq
				}
			}
		}
	}

	vals = make([]float64, n)
	vecs = make([][]float64, n)
	for j := 0; j < n; j++ {
		vals[j] = m[j][j]
		col := make([]float64, n)
		for i := 0; i < n; i++ {
			col[i] = v[i][j]
		}
		vecs[j] = col
	}

	// Urutkan nilai eigen menurun.
	for i := 0; i < n-1; i++ {
		maxIdx := i
		for j := i + 1; j < n; j++ {
			if vals[j] > vals[maxIdx] {
				maxIdx = j
			}
		}
		if maxIdx != i {
			vals[i], vals[maxIdx] = vals[maxIdx], vals[i]
			vecs[i], vecs[maxIdx] = vecs[maxIdx], vecs[i]
		}
	}

	// Kanonisasi tanda: komponen tak-nol pertama dibuat positif.
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			if math.Abs(vecs[j][i]) > 1e-12 {
				if vecs[j][i] < 0 {
					for r := 0; r < n; r++ {
						vecs[j][r] = -vecs[j][r]
					}
				}
				break
			}
		}
	}
	return vals, vecs, sweeps, nil
}

// MatExpSeries menghitung exp(A) via scaling-and-squaring dengan deret Taylor terpotong.
func MatExpSeries(a [][]float64, maxTerms int) [][]float64 {
	n := len(a)
	if n == 0 {
		return nil
	}
	if maxTerms <= 0 {
		maxTerms = 30
	}
	normInf := 0.0
	for i := 0; i < n; i++ {
		rowSum := 0.0
		for j := 0; j < n; j++ {
			rowSum += math.Abs(a[i][j])
		}
		if rowSum > normInf {
			normInf = rowSum
		}
	}
	s := 0
	for normInf > 0.5 && s < 30 {
		normInf *= 0.5
		s++
	}
	scaleFactor := math.Ldexp(1.0, -s)
	scaled := Scale(a, scaleFactor)

	res := Eye(n)
	term := Eye(n)
	for k := 1; k <= maxTerms; k++ {
		term = Scale(MatMul(term, scaled), 1.0/float64(k))
		maxEntry := 0.0
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				res[i][j] += term[i][j]
				if v := math.Abs(term[i][j]); v > maxEntry {
					maxEntry = v
				}
			}
		}
		if maxEntry < 1e-17 {
			break
		}
	}
	for i := 0; i < s; i++ {
		res = MatMul(res, res)
	}
	return res
}

// QuadNode menyimpan simpul dan bobot kuadratur.
type QuadNode struct {
	X float64
	W float64
}

// GaussLegendre menghitung n simpul dan bobot kuadratur Gauss-Legendre pada selang [a, b]
// memakai algoritma Golub-Welsch pada matriks Jacobi polinomial Legendre.
func GaussLegendre(n int, a, b float64) []QuadNode {
	if n <= 0 {
		return nil
	}
	if n == 1 {
		return []QuadNode{{X: 0.5 * (a + b), W: b - a}}
	}
	d := make([]float64, n)
	e := make([]float64, n)
	for i := 0; i < n-1; i++ {
		k := float64(i + 1)
		e[i] = k / math.Sqrt((2.0*k-1.0)*(2.0*k+1.0))
	}
	v0 := make([]float64, n)
	v0[0] = 1.0

	const epsMach = 2.220446049250313e-16
	for l := 0; l < n; l++ {
		iter := 0
		for {
			m := l
			for m < n-1 {
				tst := math.Abs(d[m]) + math.Abs(d[m+1])
				if tst+math.Abs(e[m]) == tst || math.Abs(e[m]) <= epsMach*tst {
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
			s, c, p := 1.0, 1.0, 0.0
			for i := m - 1; i >= l; i-- {
				f := s * e[i]
				bv := c * e[i]
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
				r = (d[i]-g)*s + 2.0*c*bv
				p = s * r
				d[i+1] = g + p
				g = c*r - bv

				vi, vip1 := v0[i], v0[i+1]
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

	nodes := make([]QuadNode, n)
	wRaw := 0.0
	for i := 0; i < n; i++ {
		w := v0[i] * v0[i]
		nodes[i] = QuadNode{X: d[i], W: w}
		wRaw += w
	}
	scale := 2.0 / wRaw
	mid := 0.5 * (a + b)
	half := 0.5 * (b - a)
	for i := 0; i < n/2; i++ {
		j := n - 1 - i
		xAvg := 0.5 * (nodes[j].X - nodes[i].X)
		wAvg := 0.5 * (nodes[i].W + nodes[j].W) * scale
		nodes[i] = QuadNode{X: mid - half*xAvg, W: half * wAvg}
		nodes[j] = QuadNode{X: mid + half*xAvg, W: half * wAvg}
	}
	if n%2 == 1 {
		mIdx := n / 2
		nodes[mIdx] = QuadNode{X: mid, W: half * nodes[mIdx].W * scale}
	}
	return nodes
}
