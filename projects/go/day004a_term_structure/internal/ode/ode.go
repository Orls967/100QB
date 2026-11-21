package ode

import "math"

// RK4System mengintegrasikan dy/dt = deriv(t, y) dari t0 ke t1 dengan nSteps langkah seragam.
func RK4System(t0, t1 float64, y0 []float64, nSteps int, deriv func(t float64, y, dydt []float64)) []float64 {
	dim := len(y0)
	y := make([]float64, dim)
	copy(y, y0)
	if nSteps <= 0 || t1 == t0 {
		return y
	}
	h := (t1 - t0) / float64(nSteps)
	k1 := make([]float64, dim)
	k2 := make([]float64, dim)
	k3 := make([]float64, dim)
	k4 := make([]float64, dim)
	yt := make([]float64, dim)

	for step := 0; step < nSteps; step++ {
		t := t0 + float64(step)*h
		deriv(t, y, k1)
		for i := 0; i < dim; i++ {
			yt[i] = y[i] + 0.5*h*k1[i]
		}
		deriv(t+0.5*h, yt, k2)
		for i := 0; i < dim; i++ {
			yt[i] = y[i] + 0.5*h*k2[i]
		}
		deriv(t+0.5*h, yt, k3)
		for i := 0; i < dim; i++ {
			yt[i] = y[i] + h*k3[i]
		}
		deriv(t+h, yt, k4)
		for i := 0; i < dim; i++ {
			y[i] += (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])
		}
	}
	return y
}

// RK4Grid mengintegrasikan secara berurutan melintasi titik grid meningkat mulai dari t0.
func RK4Grid(t0 float64, grid []float64, y0 []float64, subStepsPerUnit float64, deriv func(t float64, y, dydt []float64)) [][]float64 {
	out := make([][]float64, len(grid))
	curT := t0
	curY := make([]float64, len(y0))
	copy(curY, y0)
	for i, targetT := range grid {
		dt := targetT - curT
		if dt > 0 {
			steps := int(math.Ceil(dt * subStepsPerUnit))
			if steps < 4 {
				steps = 4
			}
			curY = RK4System(curT, targetT, curY, steps, deriv)
			curT = targetT
		}
		cp := make([]float64, len(curY))
		copy(cp, curY)
		out[i] = cp
	}
	return out
}

// RK4Richardson meningkatkan akurasi RK4 lewat ekstrapolasi Richardson (16*y_{2N} - y_N)/15.
func RK4Richardson(t0, t1 float64, y0 []float64, nStepsCoarse int, deriv func(t float64, y, dydt []float64)) []float64 {
	yc := RK4System(t0, t1, y0, nStepsCoarse, deriv)
	yf := RK4System(t0, t1, y0, 2*nStepsCoarse, deriv)
	out := make([]float64, len(y0))
	for i := range out {
		out[i] = (16.0*yf[i] - yc[i]) / 15.0
	}
	return out
}
