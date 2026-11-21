package ode

import (
	"math"
	"testing"
)

func TestRK4AndRichardsonOrder(t *testing.T) {
	// Sistem osilator harmonik: y1' = y2, y2' = -y1, y(0) = (1, 0) -> y(T) = (cos T, -sin T)
	deriv := func(_ float64, y, dy []float64) {
		dy[0] = y[1]
		dy[1] = -y[0]
	}
	y0 := []float64{1, 0}
	tEnd := 5.0
	exact0 := math.Cos(tEnd)
	exact1 := -math.Sin(tEnd)

	yCoarse := RK4System(0, tEnd, y0, 200, deriv)
	yFine := RK4System(0, tEnd, y0, 400, deriv)
	yRich := RK4Richardson(0, tEnd, y0, 200, deriv)

	errCoarse := math.Abs(yCoarse[0] - exact0)
	errFine := math.Abs(yFine[0] - exact0)
	errRich := math.Hypot(yRich[0]-exact0, yRich[1]-exact1)

	ratio := errCoarse / errFine
	if ratio < 15.0 || ratio > 17.5 {
		t.Fatalf("RK4 order ratio expected ~16, got %v", ratio)
	}
	if errRich > 1e-10 {
		t.Fatalf("Richardson error too large: %v", errRich)
	}
}
