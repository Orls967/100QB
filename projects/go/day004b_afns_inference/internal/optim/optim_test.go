package optim

import (
	"math"
	"testing"
)

func TestBFGSRosenbrockAndBounds(t *testing.T) {
	rosen := func(x []float64) float64 {
		a := 1.0 - x[0]
		b := x[1] - x[0]*x[0]
		return a*a + 100.0*b*b
	}
	res := MinimizeBFGS(rosen, []float64{-1.2, 1.0}, nil, nil, 250, 1e-7)
	if !res.Converged {
		t.Fatalf("expected convergence on Rosenbrock, grad=%v", res.GradNorm)
	}
	if math.Abs(res.X[0]-1.0) > 1e-5 || math.Abs(res.X[1]-1.0) > 1e-5 {
		t.Fatalf("unexpected Rosenbrock minimum: %v", res.X)
	}

	_, _, relErr := VerifyGradVsRichardson(rosen, []float64{0.3, 0.7}, nil)
	if relErr > 1e-6 {
		t.Fatalf("unexpected central vs Richardson gradient diff: %v", relErr)
	}

	// Uji deteksi batas aktif per parameter
	quad := func(x []float64) float64 {
		return (x[0]+2.0)*(x[0]+2.0) + (x[1]-0.2)*(x[1]-0.2)
	}
	bRes := MinimizeBFGS(quad, []float64{0.5, 0.5}, []float64{0.0, 0.0}, []float64{5.0, 5.0}, 100, 1e-7)
	if !bRes.HitBound || len(bRes.HitBoundByParam) != 2 || !bRes.HitBoundByParam[0] || bRes.HitBoundByParam[1] {
		t.Fatalf("unexpected bound flags: HitBound=%v ByParam=%v X=%v", bRes.HitBound, bRes.HitBoundByParam, bRes.X)
	}
}
