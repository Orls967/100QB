package special

import (
	"math"
	"testing"
)

func TestBrentGoldenAndLU(t *testing.T) {
	r, err := BrentRoot(func(x float64) float64 { return x*x - 2.0 }, 1.0, 2.0, 1e-14)
	if err != nil || math.Abs(r-math.Sqrt2) > 1e-13 {
		t.Fatalf("BrentRoot failed: r=%v err=%v", r, err)
	}

	xMin, fMin := GoldenSectionMin(func(x float64) float64 { return (x-1.7)*(x-1.7) + 2.5 }, 0.0, 5.0, 1e-10)
	if math.Abs(xMin-1.7) > 5e-8 || math.Abs(fMin-2.5) > 1e-12 {
		t.Fatalf("GoldenSectionMin failed: x=%v f=%v", xMin, fMin)
	}

	a := [][]float64{{2, 1}, {1, -3}}
	b := []float64{5, -8}
	x, err := SolveLinearLU(a, b)
	if err != nil || math.Abs(x[0]-1.0) > 1e-14 || math.Abs(x[1]-3.0) > 1e-14 {
		t.Fatalf("SolveLinearLU failed: x=%v err=%v", x, err)
	}
}
