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
			return 1.0 + 2.0*x[0]
		}
		if x[0] < 0 {
			return 1.0 - x[0]
		}
		return 1.0
	}
	rLS := MinimizeBFGS(stepDisc, []float64{0.0}, nil, nil, 50, 1e-8)
	if rLS.Converged || rLS.TerminationReason != TermLineSearchFailed || rLS.LineSearchFailureCause != LSCauseArmijoNotMet || rLS.ObjBitIdentical {
		t.Fatalf("expected non-plateau armijo_not_met line_search_failed, got %+v", rLS)
	}
	if rLS.MinAlphaTested <= 0 || rLS.MinAlphaTested > 1e-9 || rLS.LastStepNorm <= 0 {
		t.Fatalf("expected small MinAlphaTested and positive LastStepNorm, got alpha=%g step=%g", rLS.MinAlphaTested, rLS.LastStepNorm)
	}

	// Plateau numerik: nilai objektif identik bit-demi-bit di sekitar skala besar (~1e5) meski gradien beda hingga bernilai non-nol
	plateau := func(x []float64) float64 {
		v := x[0] - 0.5
		if v >= 5e-6 {
			return -1.0e5 + 1e3*v
		}
		return -1.0e5
	}
	rPlat := MinimizeBFGS(plateau, []float64{0.5}, nil, nil, 50, 1e-8)
	if rPlat.Converged || rPlat.TerminationReason != TermLineSearchFailed {
		t.Fatalf("expected plateau line_search_failed, got %+v", rPlat)
	}
	if !rPlat.ObjBitIdentical || rPlat.LastObjChange != 0.0 || rPlat.LineSearchFailureCause != LSCauseArmijoNotMet {
		t.Fatalf("expected bit-identical plateau objective, got objIdent=%v dObj=%g cause=%s", rPlat.ObjBitIdentical, rPlat.LastObjChange, rPlat.LineSearchFailureCause)
	}

	rMax := MinimizeBFGS(rosen, []float64{-1.2, 1.0}, nil, nil, 2, 1e-8)
	if rMax.Converged || rMax.TerminationReason != TermMaxIterations || rMax.LineSearchFailureCause != LSCauseNone {
		t.Fatalf("expected max_iterations with LSCauseNone, got %+v", rMax)
	}

	rInv := MinimizeBFGS(func(x []float64) float64 { return math.NaN() }, []float64{0.5}, nil, nil, 10, 1e-6)
	if rInv.Converged || rInv.TerminationReason != TermInvalidObjective || rInv.LineSearchFailureCause != LSCauseInvalidObjective {
		t.Fatalf("expected invalid_objective_or_gradient, got %+v", rInv)
	}
}
