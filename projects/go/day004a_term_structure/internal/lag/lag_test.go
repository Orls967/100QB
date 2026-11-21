package lag

import (
	"math"
	"testing"
)

func TestLUCholeskyTridiagonalJacobiExpAndGL(t *testing.T) {
	// 1. LU & InvertLU
	a := [][]float64{
		{3, 1, -2},
		{2, 4, 1},
		{-1, 2, 5},
	}
	b := []float64{1, 13, 18}
	x, err := SolveLU(a, b)
	if err != nil {
		t.Fatalf("SolveLU failed: %v", err)
	}
	ax := MatVec(a, x)
	for i := range b {
		if math.Abs(ax[i]-b[i]) > 1e-14 {
			t.Fatalf("LU residual too high at %d: %v", i, ax[i]-b[i])
		}
	}
	invA, err := InvertLU(a)
	if err != nil {
		t.Fatalf("InvertLU failed: %v", err)
	}
	prod := MatMul(a, invA)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			exp := 0.0
			if i == j {
				exp = 1.0
			}
			if math.Abs(prod[i][j]-exp) > 1e-14 {
				t.Fatalf("A*A^-1 mismatch at (%d,%d): %v", i, j, prod[i][j])
			}
		}
	}

	// 2. Cholesky
	spd := [][]float64{
		{4, 2, -1},
		{2, 5, 1},
		{-1, 1, 6},
	}
	l, err := Cholesky(spd)
	if err != nil {
		t.Fatalf("Cholesky failed: %v", err)
	}
	llt := MatMul(l, Transpose(l))
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if math.Abs(llt[i][j]-spd[i][j]) > 1e-14 {
				t.Fatalf("L*L^T mismatch at (%d,%d): %v", i, j, llt[i][j])
			}
		}
	}
	xChol := SolveCholesky(l, b)
	axC := MatVec(spd, xChol)
	for i := range b {
		if math.Abs(axC[i]-b[i]) > 1e-14 {
			t.Fatalf("Cholesky solve mismatch at %d: %v", i, axC[i]-b[i])
		}
	}

	// 3. Thomas tridiagonal
	sub := []float64{-1, -1, -1}
	diag := []float64{4, 4, 4, 4}
	sup := []float64{-1, -1, -1}
	rhs := []float64{5, 5, 10, 23}
	xt, err := SolveTridiagonal(sub, diag, sup, rhs)
	if err != nil {
		t.Fatalf("SolveTridiagonal failed: %v", err)
	}
	tri := [][]float64{
		{4, -1, 0, 0},
		{-1, 4, -1, 0},
		{0, -1, 4, -1},
		{0, 0, -1, 4},
	}
	resTri := MatVec(tri, xt)
	for i := range rhs {
		if math.Abs(resTri[i]-rhs[i]) > 1e-14 {
			t.Fatalf("Tridiagonal mismatch at %d: %v", i, resTri[i]-rhs[i])
		}
	}

	// 4. SymmetricEigenJacobi
	sym := [][]float64{
		{5, 2, 1},
		{2, 4, -1},
		{1, -1, 3},
	}
	vals, vecs, _, err := SymmetricEigenJacobi(sym, 80, 1e-14)
	if err != nil {
		t.Fatalf("Jacobi failed: %v", err)
	}
	for k := 0; k < 3; k++ {
		av := MatVec(sym, vecs[k])
		for i := 0; i < 3; i++ {
			if math.Abs(av[i]-vals[k]*vecs[k][i]) > 1e-13 {
				t.Fatalf("Eigen equation failed for k=%d, i=%d: %v vs %v", k, i, av[i], vals[k]*vecs[k][i])
			}
		}
		for lIdx := k; lIdx < 3; lIdx++ {
			dot := Dot(vecs[k], vecs[lIdx])
			exp := 0.0
			if k == lIdx {
				exp = 1.0
			}
			if math.Abs(dot-exp) > 1e-13 {
				t.Fatalf("Orthonormality failed at (%d,%d): %v", k, lIdx, dot)
			}
		}
	}

	// 5. MatExpSeries untuk matriks segitiga bawah eksak
	kMat := [][]float64{
		{-0.5, 0},
		{0.3, -0.5},
	}
	expK := MatExpSeries(kMat, 30)
	e05 := math.Exp(-0.5)
	if math.Abs(expK[0][0]-e05) > 1e-14 || math.Abs(expK[1][0]-0.3*e05) > 1e-14 || math.Abs(expK[1][1]-e05) > 1e-14 {
		t.Fatalf("MatExpSeries Jordan block failed: %+v", expK)
	}

	// 6. Gauss-Legendre akurasi polinomial orde tinggi
	nodes := GaussLegendre(16, 0.0, 2.0)
	integral := 0.0
	for _, nd := range nodes {
		integral += nd.W * math.Pow(nd.X, 10)
	}
	exact := math.Pow(2.0, 11) / 11.0
	if math.Abs(integral-exact) > 1e-12 {
		t.Fatalf("GaussLegendre failed on x^10: got %v, want %v", integral, exact)
	}
}
