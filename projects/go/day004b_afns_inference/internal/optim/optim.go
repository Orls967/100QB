package optim

import (
	"math"
)

// CentralDiffGrad menghitung gradien beda hingga pusat:
// g_i = (f(x + h_i e_i) - f(x - h_i e_i)) / (2 h_i).
func CentralDiffGrad(f func([]float64) float64, x, steps []float64) []float64 {
	n := len(x)
	g := make([]float64, n)
	xp := make([]float64, n)
	copy(xp, x)
	for i := 0; i < n; i++ {
		h := 1e-5 * math.Max(1.0, math.Abs(x[i]))
		if len(steps) == n && steps[i] > 0 {
			h = steps[i]
		}
		orig := xp[i]
		xp[i] = orig + h
		fp := f(xp)
		xp[i] = orig - h
		fm := f(xp)
		xp[i] = orig
		g[i] = (fp - fm) / (2.0 * h)
	}
	return g
}

// RichardsonGrad mengevaluasi gradien beda hingga pusat orde-4 melalui ekstrapolasi Richardson
// g_R = (4 * g(h/2) - g(h)) / 3.
func RichardsonGrad(f func([]float64) float64, x, steps []float64) []float64 {
	n := len(x)
	hCoarse := make([]float64, n)
	hFine := make([]float64, n)
	for i := 0; i < n; i++ {
		h := 2e-5 * math.Max(1.0, math.Abs(x[i]))
		if len(steps) == n && steps[i] > 0 {
			h = steps[i]
		}
		hCoarse[i] = h
		hFine[i] = 0.5 * h
	}
	g1 := CentralDiffGrad(f, x, hCoarse)
	g2 := CentralDiffGrad(f, x, hFine)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = (4.0*g2[i] - g1[i]) / 3.0
	}
	return out
}

// VerifyGradVsRichardson membandingkan gradien beda hingga pusat terhadap gradien Richardson
// dan mengembalikan galat relatif L2.
func VerifyGradVsRichardson(f func([]float64) float64, x, steps []float64) (cdGrad, richGrad []float64, relErr float64) {
	cdGrad = CentralDiffGrad(f, x, steps)
	richGrad = RichardsonGrad(f, x, steps)
	diffSq, normSq := 0.0, 0.0
	for i := range x {
		d := cdGrad[i] - richGrad[i]
		diffSq += d * d
		normSq += richGrad[i] * richGrad[i]
	}
	relErr = math.Sqrt(diffSq) / math.Max(math.Sqrt(normSq), 1e-12)
	return cdGrad, richGrad, relErr
}

// CentralDiffHessian menghitung matriks Hessian simetris f di x dengan beda hingga pusat.
func CentralDiffHessian(f func([]float64) float64, x, steps []float64) [][]float64 {
	n := len(x)
	h := make([]float64, n)
	for i := 0; i < n; i++ {
		if len(steps) == n && steps[i] > 0 {
			h[i] = steps[i]
		} else {
			h[i] = 1e-4 * math.Max(1.0, math.Abs(x[i]))
		}
	}
	f0 := f(x)
	hess := make([][]float64, n)
	for i := 0; i < n; i++ {
		hess[i] = make([]float64, n)
	}
	xp := make([]float64, n)
	copy(xp, x)

	for i := 0; i < n; i++ {
		orig := xp[i]
		xp[i] = orig + h[i]
		fp := f(xp)
		xp[i] = orig - h[i]
		fm := f(xp)
		xp[i] = orig
		hess[i][i] = (fp - 2.0*f0 + fm) / (h[i] * h[i])
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			oi, oj := xp[i], xp[j]
			xp[i], xp[j] = oi+h[i], oj+h[j]
			fpp := f(xp)
			xp[i], xp[j] = oi+h[i], oj-h[j]
			fpm := f(xp)
			xp[i], xp[j] = oi-h[i], oj+h[j]
			fmp := f(xp)
			xp[i], xp[j] = oi-h[i], oj-h[j]
			fmm := f(xp)
			xp[i], xp[j] = oi, oj
			val := (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j])
			hess[i][j] = val
			hess[j][i] = val
		}
	}
	return hess
}

const (
	// TermGradientTol menandakan konvergensi ketat di titik interior karena norma gradien <= gradTol.
	TermGradientTol = "gradient_tolerance"
	// TermBoundaryStationary menandakan solusi berhenti pada batas kotak dengan gradien terproyeksi <= gradTol.
	TermBoundaryStationary = "boundary_stationary"
	// TermLineSearchFailed menandakan pencarian garis Armijo gagal menemukan langkah turun atau stagnan sebelum gradTol.
	TermLineSearchFailed = "line_search_failed"
	// TermMaxIterations menandakan penghentian karena mencapai batas iterasi sebelum gradTol terpenuhi.
	TermMaxIterations = "max_iterations"
	// TermInvalidObjective menandakan nilai objektif atau gradien non-finite atau melebihi penalti tidak valid.
	TermInvalidObjective = "invalid_objective_or_gradient"
)

// BFGSResult merangkum keluaran minimisasi quasi-Newton BFGS dengan kendala kotak opsional.
type BFGSResult struct {
	X                 []float64
	Fun               float64
	RawGradNorm       float64
	GradNorm          float64
	ScaledGradNorm    float64
	Iterations        int
	FuncEvals         int
	Converged         bool
	HitBound          bool
	HitBoundByParam   []bool
	TerminationReason string
}

func clampVec(x, lo, hi []float64) {
	if len(lo) == len(x) && len(hi) == len(x) {
		for i := range x {
			if x[i] < lo[i] {
				x[i] = lo[i]
			} else if x[i] > hi[i] {
				x[i] = hi[i]
			}
		}
	}
}

func isInvalidState(val float64, g []float64) bool {
	if math.IsNaN(val) || math.IsInf(val, 0) || val >= 1e29 {
		return true
	}
	for _, gi := range g {
		if math.IsNaN(gi) || math.IsInf(gi, 0) {
			return true
		}
	}
	return false
}

// MinimizeBFGS meminimalkan fungsi f(x) memakai BFGS dengan backtracking Armijo
// dan penjepit kotak [lo, hi] pengaman. Konvergensi ketat hanya tercapai bila
// norma gradien memenuhi gradTol di titik interior; alasan penghentian lain dicatat eksplisit.
func MinimizeBFGS(f func([]float64) float64, x0, lo, hi []float64, maxIter int, gradTol float64) BFGSResult {
	n := len(x0)
	x := make([]float64, n)
	copy(x, x0)
	clampVec(x, lo, hi)

	evals := 0
	evalF := func(v []float64) float64 {
		evals++
		return f(v)
	}
	gradF := func(v []float64) []float64 {
		return CentralDiffGrad(evalF, v, nil)
	}

	rawInfNorm := func(v []float64) float64 {
		m := 0.0
		for _, a := range v {
			if math.IsNaN(a) {
				return math.NaN()
			}
			if math.Abs(a) > m {
				m = math.Abs(a)
			}
		}
		return m
	}

	checkBounds := func(vec []float64) (bool, []bool) {
		hit := false
		byParam := make([]bool, n)
		if len(lo) == n && len(hi) == n {
			for i := 0; i < n; i++ {
				if vec[i] <= lo[i]+1e-9 || vec[i] >= hi[i]-1e-9 {
					hit = true
					byParam[i] = true
				}
			}
		}
		return hit, byParam
	}

	infNorm := func(v []float64) float64 {
		m := 0.0
		for i, a := range v {
			if math.IsNaN(a) {
				return math.NaN()
			}
			// Abaikan komponen gradien yang mendorong ke luar batas aktif
			if len(lo) == len(v) && len(hi) == len(v) {
				if x[i] <= lo[i]+1e-12 && a > 0 {
					continue
				}
				if x[i] >= hi[i]-1e-12 && a < 0 {
					continue
				}
			}
			if math.Abs(a) > m {
				m = math.Abs(a)
			}
		}
		return m
	}

	fx := evalF(x)
	g := gradF(x)

	hMat := make([][]float64, n)
	for i := 0; i < n; i++ {
		hMat[i] = make([]float64, n)
		hMat[i][i] = 1.0
	}

	termReason := ""
	iters := 0
	xCand := make([]float64, n)
	d := make([]float64, n)
	s := make([]float64, n)
	y := make([]float64, n)
	hy := make([]float64, n)

	if isInvalidState(fx, g) {
		termReason = TermInvalidObjective
	}

	for iter := 0; termReason == "" && iter < maxIter; iter++ {
		gn := infNorm(g)
		if gn <= gradTol {
			hb, _ := checkBounds(x)
			if hb {
				termReason = TermBoundaryStationary
			} else {
				termReason = TermGradientTol
			}
			break
		}
		iters = iter + 1

		slope := 0.0
		for i := 0; i < n; i++ {
			sum := 0.0
			for j := 0; j < n; j++ {
				sum -= hMat[i][j] * g[j]
			}
			d[i] = sum
			slope += g[i] * d[i]
		}
		if slope >= 0 || math.IsNaN(slope) {
			for i := 0; i < n; i++ {
				for j := 0; j < n; j++ {
					hMat[i][j] = 0
				}
				hMat[i][i] = 1.0
				d[i] = -g[i]
			}
		}

		alpha := 1.0
		accepted := false
		fCand := fx
		const c1 = 1e-4
		for ls := 0; ls < 35; ls++ {
			for i := 0; i < n; i++ {
				xCand[i] = x[i] + alpha*d[i]
			}
			clampVec(xCand, lo, hi)
			dirDeriv := 0.0
			for i := 0; i < n; i++ {
				dirDeriv += g[i] * (xCand[i] - x[i])
			}
			fCand = evalF(xCand)
			if !math.IsNaN(fCand) && !math.IsInf(fCand, 0) && fCand < 1e29 && fCand <= fx+c1*dirDeriv {
				accepted = true
				break
			}
			alpha *= 0.5
		}
		if !accepted {
			termReason = TermLineSearchFailed
			break
		}

		stepNorm := 0.0
		for i := 0; i < n; i++ {
			s[i] = xCand[i] - x[i]
			x[i] = xCand[i]
			stepNorm += s[i] * s[i]
		}
		stepNorm = math.Sqrt(stepNorm)
		fPrev := fx
		fx = fCand
		gNew := gradF(x)
		if isInvalidState(fx, gNew) {
			g = gNew
			termReason = TermInvalidObjective
			break
		}

		sy := 0.0
		sNormSq := 0.0
		yNormSq := 0.0
		for i := 0; i < n; i++ {
			y[i] = gNew[i] - g[i]
			sy += s[i] * y[i]
			sNormSq += s[i] * s[i]
			yNormSq += y[i] * y[i]
		}
		g = gNew

		if stepNorm < 1e-13 && math.Abs(fx-fPrev) < 1e-12*math.Max(1.0, math.Abs(fx)) {
			if infNorm(g) <= gradTol {
				hb, _ := checkBounds(x)
				if hb {
					termReason = TermBoundaryStationary
				} else {
					termReason = TermGradientTol
				}
			} else {
				termReason = TermLineSearchFailed
			}
			break
		}

		if sy > 1e-12*math.Sqrt(sNormSq*yNormSq) {
			rho := 1.0 / sy
			if iter == 0 {
				scale := sy / yNormSq
				for i := 0; i < n; i++ {
					hMat[i][i] = scale
				}
			}
			yHy := 0.0
			for i := 0; i < n; i++ {
				sum := 0.0
				for j := 0; j < n; j++ {
					sum += hMat[i][j] * y[j]
				}
				hy[i] = sum
				yHy += y[i] * sum
			}
			coef := (1.0 + rho*yHy) * rho
			for i := 0; i < n; i++ {
				for j := i; j < n; j++ {
					val := hMat[i][j] + coef*s[i]*s[j] - rho*(hy[i]*s[j]+s[i]*hy[j])
					hMat[i][j] = val
					hMat[j][i] = val
				}
			}
		}
	}

	hitBound, hitByParam := checkBounds(x)
	gnFinal := infNorm(g)
	rawFinal := rawInfNorm(g)
	if termReason == "" {
		if gnFinal <= gradTol {
			if hitBound {
				termReason = TermBoundaryStationary
			} else {
				termReason = TermGradientTol
			}
		} else {
			termReason = TermMaxIterations
		}
	}
	scaledFinal := gnFinal / math.Max(1.0, math.Abs(fx))
	return BFGSResult{
		X:                 x,
		Fun:               fx,
		RawGradNorm:       rawFinal,
		GradNorm:          gnFinal,
		ScaledGradNorm:    scaledFinal,
		Iterations:        iters,
		FuncEvals:         evals,
		Converged:         termReason == TermGradientTol,
		HitBound:          hitBound,
		HitBoundByParam:   hitByParam,
		TerminationReason: termReason,
	}
}
