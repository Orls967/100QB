package linalg

import "math"

// SolveLU menyelesaikan sistem linear A * x = b dengan eliminasi Gauss + partial pivoting.
// Mengembalikan (x, true) bila nonsingular, atau (nil, false) bila singular.
func SolveLU(a [][]float64, b []float64) ([]float64, bool) {
	n := len(b)
	if n == 0 {
		return nil, true
	}
	m := make([][]float64, n)
	for i := 0; i < n; i++ {
		row := make([]float64, n+1)
		copy(row, a[i])
		row[n] = b[i]
		m[i] = row
	}

	for k := 0; k < n; k++ {
		pivRow := k
		maxVal := math.Abs(m[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(m[i][k]); v > maxVal {
				maxVal = v
				pivRow = i
			}
		}
		if maxVal < 1e-13 {
			return nil, false
		}
		if pivRow != k {
			m[k], m[pivRow] = m[pivRow], m[k]
		}
		piv := m[k][k]
		for i := k + 1; i < n; i++ {
			factor := m[i][k] / piv
			m[i][k] = 0.0
			for j := k + 1; j <= n; j++ {
				m[i][j] -= factor * m[k][j]
			}
		}
	}

	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := m[i][n]
		for j := i + 1; j < n; j++ {
			sum -= m[i][j] * x[j]
		}
		x[i] = sum / m[i][i]
	}
	return x, true
}
