package optim

import (
	"math"
	"testing"
)

func TestBFGSTerminationReasonsAndRichardson(t *testing.T) {
	rosen := func(x []float64) float64 {
		a := 1.0 - x[0]
		b := x[1] - x[0]*x[0]
		return a*a + 100.0*b*b
	}
	rTol := MinimizeBFGS(rosen, []float64{-1.2, 1.0}, nil, nil, 250, 1e-7)
	if !rTol.Converged || rTol.TerminationReason != TermGradientTol {
		t.Fatalf("expected gradient_tolerance, got %s (conv=%v)", rTol.TerminationReason, rTol.Converged)
	}
	if math.Abs(rTol.X[0]-1.0) > 1e-5 || math.Abs(rTol.X[1]-1.0) > 1e-5 {
		t.Fatalf("unexpected Rosenbrock minimum: %v", rTol.X)
	}
	_, _, relErr := VerifyGradVsRichardson(rosen, []float64{0.3, 0.7}, nil)
	if relErr > 1e-6 {
		t.Fatalf("unexpected gradient Richardson error: %v", relErr)
	}

	quad := func(x []float64) float64 {
		return (x[0]-1.0)*(x[0]-1.0) + 4.0*(x[1]+0.5)*(x[1]+0.5)
	}
	rBound := MinimizeBFGS(quad, []float64{0.0, 0.0}, []float64{2.0, -5.0}, []float64{5.0, 5.0}, 100, 1e-7)
	if rBound.Converged || rBound.TerminationReason != TermBoundaryStationary || !rBound.HitBoundByParam[0] || rBound.HitBoundByParam[1] {
		t.Fatalf("expected boundary_stationary on x[0], got %s (%v)", rBound.TerminationReason, rBound.HitBoundByParam)
	}

	stepDisc := func(x []float64) float64 {
		if x[0] >= 1e-6 {
			return 1.0 + x[0]
		}
		return 1.0
	}
	rLS := MinimizeBFGS(stepDisc, []float64{0.0}, nil, nil, 50, 1e-8)
	if rLS.Converged || rLS.TerminationReason != TermLineSearchFailed {
		t.Fatalf("expected line_search_failed, got %s", rLS.TerminationReason)
	}

	rMax := MinimizeBFGS(rosen, []float64{-1.2, 1.0}, nil, nil, 2, 1e-8)
	if rMax.Converged || rMax.TerminationReason != TermMaxIterations {
		t.Fatalf("expected max_iterations, got %s", rMax.TerminationReason)
	}

	rInv := MinimizeBFGS(func(x []float64) float64 { return math.NaN() }, []float64{0.5}, nil, nil, 10, 1e-6)
	if rInv.Converged || rInv.TerminationReason != TermInvalidObjective {
		t.Fatalf("expected invalid_objective_or_gradient, got %s", rInv.TerminationReason)
	}
}
