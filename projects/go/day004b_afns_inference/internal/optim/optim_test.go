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
	if bRes.TerminationReason != TermBoundaryStationary || bRes.Converged {
		t.Fatalf("expected boundary_stationary (Converged=false), got reason=%s converged=%v", bRes.TerminationReason, bRes.Converged)
	}
}

func TestBFGSTerminationReasons(t *testing.T) {
	// 1. gradient_tolerance di titik interior
	quad := func(x []float64) float64 {
		return (x[0]-1.0)*(x[0]-1.0) + 4.0*(x[1]+0.5)*(x[1]+0.5)
	}
	rTol := MinimizeBFGS(quad, []float64{3.0, 2.0}, nil, nil, 100, 1e-7)
	if !rTol.Converged || rTol.TerminationReason != TermGradientTol {
		t.Fatalf("expected gradient_tolerance, got %s (conv=%v)", rTol.TerminationReason, rTol.Converged)
	}

	// 2. boundary_stationary pada batas kotak
	rBound := MinimizeBFGS(quad, []float64{0.0, 0.0}, []float64{2.0, -5.0}, []float64{5.0, 5.0}, 100, 1e-7)
	if rBound.Converged || rBound.TerminationReason != TermBoundaryStationary || !rBound.HitBoundByParam[0] {
		t.Fatalf("expected boundary_stationary, got %s (conv=%v hit=%v)", rBound.TerminationReason, rBound.Converged, rBound.HitBoundByParam)
	}

	// 3. line_search_failed ketika gradien FD > 0 menunjuk arah turun, namun seluruh langkah Armijo gagal menurunkan f
	stepDisc := func(x []float64) float64 {
		if x[0] >= 1e-6 {
			return 1.0 + x[0]
		}
		return 1.0
	}
	rLS := MinimizeBFGS(stepDisc, []float64{0.0}, nil, nil, 50, 1e-8)
	if rLS.Converged || rLS.TerminationReason != TermLineSearchFailed {
		t.Fatalf("expected line_search_failed, got %s (conv=%v)", rLS.TerminationReason, rLS.Converged)
	}

	// 4. max_iterations ketika maxIter terlalu kecil untuk konvergen
	rosen := func(x []float64) float64 {
		a := 1.0 - x[0]
		b := x[1] - x[0]*x[0]
		return a*a + 100.0*b*b
	}
	rMax := MinimizeBFGS(rosen, []float64{-1.2, 1.0}, nil, nil, 2, 1e-8)
	if rMax.Converged || rMax.TerminationReason != TermMaxIterations {
		t.Fatalf("expected max_iterations, got %s (conv=%v)", rMax.TerminationReason, rMax.Converged)
	}

	// 5. invalid_objective_or_gradient saat objektif bernilai NaN / penalti
	badObj := func(x []float64) float64 {
		return math.NaN()
	}
	rInv := MinimizeBFGS(badObj, []float64{0.5}, nil, nil, 10, 1e-6)
	if rInv.Converged || rInv.TerminationReason != TermInvalidObjective {
		t.Fatalf("expected invalid_objective_or_gradient, got %s (conv=%v)", rInv.TerminationReason, rInv.Converged)
	}
}
