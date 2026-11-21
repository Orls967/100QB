package optim

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
)

func TestBFGSRosenbrockAndGradRichardson(t *testing.T) {
	// Rosenbrock 4D
	rosen := func(x []float64) float64 {
		s := 0.0
		for i := 0; i < len(x)-1; i++ {
			t1 := x[i+1] - x[i]*x[i]
			t2 := 1.0 - x[i]
			s += 100.0*t1*t1 + t2*t2
		}
		return s
	}
	xTest := []float64{-1.2, 1.0, -0.5, 0.8}
	_, _, relErr := VerifyGradVsRichardson(rosen, xTest, nil)
	if relErr >= 1e-5 {
		t.Fatalf("CentralDiff vs Richardson gradient relErr too large: %g", relErr)
	}

	res := MinimizeBFGS(rosen, xTest, nil, nil, 300, 1e-7)
	if !res.Converged || res.Fun > 1e-10 {
		t.Fatalf("BFGS Rosenbrock did not converge: %+v", res)
	}
	for i, xi := range res.X {
		if math.Abs(xi-1.0) > 1e-5 {
			t.Fatalf("Rosenbrock minimum mismatch at %d: %v", i, xi)
		}
	}

	// Uji Hessian pada kuadratik murni f(x) = 0.5 * x^T A x
	hess := CentralDiffHessian(func(x []float64) float64 {
		return 0.5 * (4*x[0]*x[0] + 2*x[0]*x[1] + 6*x[1]*x[1])
	}, []float64{0.5, -0.3}, nil)
	invH, err := lag.InvertLU(hess)
	if err != nil {
		t.Fatalf("Hessian invert failed: %v", err)
	}
	// A = [[4, 1], [1, 6]], A^-1 = (1/23)[[6, -1], [-1, 4]]
	if math.Abs(invH[0][0]-6.0/23.0) > 1e-6 || math.Abs(invH[0][1]+1.0/23.0) > 1e-6 {
		t.Fatalf("Hessian inverse mismatch: %+v", invH)
	}
}
