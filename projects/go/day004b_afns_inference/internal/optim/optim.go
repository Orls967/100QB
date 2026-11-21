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

// BFGSResult merangkum keluaran minimisasi quasi-Newton BFGS dengan kendala kotak opsional.
type BFGSResult struct {
	X               []float64
	Fun             float64
	GradNorm        float64
	ScaledGradNorm  float64
	Iterations      int
	FuncEvals       int
	Converged       bool
	HitBound        bool
	HitBoundByParam []bool
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

// MinimizeBFGS meminimalkan fungsi f(x) memakai BFGS dengan backtracking Armijo
// dan penjepit kotak [lo, hi] pengaman. Konvergensi memeriksa norma inf gradien tak-diskalakan
// maupun gradien relatif terhadap |f(x)| untuk mengakomodasi skala log-likelihood sampel besar.
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

	fx := evalF(x)
	g := gradF(x)

	hMat := make([][]float64, n)
	for i := 0; i < n; i++ {
		hMat[i] = make([]float64, n)
		hMat[i][i] = 1.0
	}

	infNorm := func(v []float64) float64 {
		m := 0.0
		for i, a := range v {
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

	converged := false
	iters := 0
	xCand := make([]float64, n)
	d := make([]float64, n)
	s := make([]float64, n)
	y := make([]float64, n)
	hy := make([]float64, n)

	for iter := 0; iter < maxIter; iter++ {
		gn := infNorm(g)
		if gn <= gradTol {
			converged = true
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
			if !math.IsNaN(fCand) && !math.IsInf(fCand, 0) && fCand <= fx+c1*dirDeriv {
				accepted = true
				break
			}
			alpha *= 0.5
		}
		if !accepted {
			scaledGrad := infNorm(g) / math.Max(1.0, math.Abs(fx))
			if scaledGrad <= 1e-6 || infNorm(g) <= math.Max(gradTol*100.0, 5e-4) {
				converged = true
			}
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

		scaledGrad := infNorm(g) / math.Max(1.0, math.Abs(fx))
		if stepNorm < 1e-13 && math.Abs(fx-fPrev) < 1e-12*math.Max(1.0, math.Abs(fx)) {
			converged = infNorm(g) <= math.Max(gradTol*100.0, 1e-4) || scaledGrad <= 1e-6
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

	hitBound := false
	hitByParam := make([]bool, n)
	if len(lo) == n && len(hi) == n {
		for i := 0; i < n; i++ {
			if x[i] <= lo[i]+1e-9 || x[i] >= hi[i]-1e-9 {
				hitBound = true
				hitByParam[i] = true
			}
		}
	}
	gnFinal := infNorm(g)
	scaledFinal := gnFinal / math.Max(1.0, math.Abs(fx))
	return BFGSResult{
		X:               x,
		Fun:             fx,
		GradNorm:        gnFinal,
		ScaledGradNorm:  scaledFinal,
		Iterations:      iters,
		FuncEvals:       evals,
		Converged:       converged || gnFinal <= 1e-4 || (scaledFinal <= 1e-6 && iters < maxIter),
		HitBound:        hitBound,
		HitBoundByParam: hitByParam,
	}
}
