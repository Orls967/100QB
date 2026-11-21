package kalman

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/afns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/optim"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/parallel"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/stats"
)

// StandardTenors adalah 10 jatuh tempo standar (dalam tahun) untuk Bagian 3 & 4.
var StandardTenors = []float64{0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30}

const (
	WeeklyDt        = 1.0 / 52.0
	DefaultSigmaEps = 0.0005 // 5 bp
	log2Pi          = 1.83787706640934548356
)

// ModelSpec menggabungkan parameter AFNS (atau DNS bila UseAFNSAdj=false) dan deviasi baku derau pengukuran.
type ModelSpec struct {
	AFNS       afns.Params
	SigmaEps   float64
	UseAFNSAdj bool
}

// DefaultSpec mengembalikan spesifikasi model benar untuk Bagian 3 dan 4.
func DefaultSpec() ModelSpec {
	return ModelSpec{
		AFNS:       afns.DefaultParams(),
		SigmaEps:   DefaultSigmaEps,
		UseAFNSAdj: true,
	}
}

// FilterOutput menyimpan keluaran lengkap filter Kalman dan opsional RTS smoother.
type FilterOutput struct {
	LogLikelihood float64
	XPred         [][3]float64
	PPred         [][3][3]float64
	XFilt         [][3]float64
	PFilt         [][3][3]float64
	XSmooth       [][3]float64
	PSmooth       [][3][3]float64
}

func buildObsModel(spec ModelSpec, taus []float64) (aVec []float64, hMat [][3]float64) {
	m := len(taus)
	aVec = make([]float64, m)
	hMat = make([][3]float64, m)
	for j, tau := range taus {
		hMat[j] = spec.AFNS.YieldLoadings(tau)
		if spec.UseAFNSAdj {
			aVec[j], _ = spec.AFNS.YieldAdjustment(tau)
		}
	}
	return aVec, hMat
}

// FastLogLikelihood menghitung log-likelihood eksak filter Kalman tanpa alokasi heap di loop waktu
// untuk mempercepat optimisasi BFGS dan studi replikasi.
func FastLogLikelihood(spec ModelSpec, taus []float64, dt float64, yields [][]float64) (float64, error) {
	m := len(taus)
	tSteps := len(yields)
	if m == 0 || m > 16 || tSteps == 0 || spec.SigmaEps <= 0 {
		return math.Inf(-1), fmt.Errorf("FastLogLikelihood: invalid inputs")
	}

	var aVec [16]float64
	var hMat [16][3]float64
	for j, tau := range taus {
		hMat[j] = spec.AFNS.YieldLoadings(tau)
		if spec.UseAFNSAdj {
			aVec[j], _ = spec.AFNS.YieldAdjustment(tau)
		}
	}

	phi := spec.AFNS.TransitionMatrixClosed(dt)
	qDt := spec.AFNS.TransitionCovClosed(dt)
	qInf := spec.AFNS.StationaryCov()
	theta := spec.AFNS.ThetaP
	rVar := spec.SigmaEps * spec.SigmaEps

	xPred := theta
	pPred := qInf
	var xFilt [3]float64
	var pFilt [3][3]float64

	var fMat [16][16]float64
	var lMat [16][16]float64
	var hp [16][3]float64
	var wMat [16][3]float64
	var v [16]float64
	var z [16]float64

	totalLL := 0.0
	constHalfM := float64(m) * log2Pi

	for t := 0; t < tSteps; t++ {
		if t > 0 {
			dx0 := xFilt[0] - theta[0]
			dx1 := xFilt[1] - theta[1]
			dx2 := xFilt[2] - theta[2]
			for i := 0; i < 3; i++ {
				xPred[i] = theta[i] + phi[i][0]*dx0 + phi[i][1]*dx1 + phi[i][2]*dx2
			}
			var tmp [3][3]float64
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					tmp[i][j] = phi[i][0]*pFilt[0][j] + phi[i][1]*pFilt[1][j] + phi[i][2]*pFilt[2][j]
				}
			}
			for i := 0; i < 3; i++ {
				for j := i; j < 3; j++ {
					val := qDt[i][j] + tmp[i][0]*phi[j][0] + tmp[i][1]*phi[j][1] + tmp[i][2]*phi[j][2]
					pPred[i][j], pPred[j][i] = val, val
				}
			}
		}

		yt := yields[t]
		for j := 0; j < m; j++ {
			predY := aVec[j] + hMat[j][0]*xPred[0] + hMat[j][1]*xPred[1] + hMat[j][2]*xPred[2]
			v[j] = yt[j] - predY
			hp[j][0] = hMat[j][0]*pPred[0][0] + hMat[j][1]*pPred[1][0] + hMat[j][2]*pPred[2][0]
			hp[j][1] = hMat[j][0]*pPred[0][1] + hMat[j][1]*pPred[1][1] + hMat[j][2]*pPred[2][1]
			hp[j][2] = hMat[j][0]*pPred[0][2] + hMat[j][1]*pPred[1][2] + hMat[j][2]*pPred[2][2]
		}

		for i := 0; i < m; i++ {
			for j := 0; j <= i; j++ {
				s := hp[i][0]*hMat[j][0] + hp[i][1]*hMat[j][1] + hp[i][2]*hMat[j][2]
				if i == j {
					s += rVar
				}
				fMat[i][j], fMat[j][i] = s, s
			}
		}

		// Cholesky F = L L^T
		logDet := 0.0
		for i := 0; i < m; i++ {
			for j := 0; j <= i; j++ {
				sum := fMat[i][j]
				for k := 0; k < j; k++ {
					sum -= lMat[i][k] * lMat[j][k]
				}
				if i == j {
					if sum <= 0 || math.IsNaN(sum) || math.IsInf(sum, 0) {
						return math.Inf(-1), fmt.Errorf("non-PD innovation covariance at t=%d", t)
					}
					diag := math.Sqrt(sum)
					lMat[i][j] = diag
					logDet += math.Log(diag)
				} else {
					lMat[i][j] = sum / lMat[j][j]
				}
			}
		}
		logDet *= 2.0

		quad := 0.0
		for i := 0; i < m; i++ {
			sum := v[i]
			for j := 0; j < i; j++ {
				sum -= lMat[i][j] * z[j]
			}
			zi := sum / lMat[i][i]
			z[i] = zi
			quad += zi * zi
		}
		totalLL += -0.5 * (constHalfM + logDet + quad)

		if t+1 == tSteps {
			break
		}

		// Hitung W = F^{-1} (H P) untuk pembaruan keadaan dan Joseph kovarians
		for col := 0; col < 3; col++ {
			for i := 0; i < m; i++ {
				sum := hp[i][col]
				for j := 0; j < i; j++ {
					sum -= lMat[i][j] * wMat[j][col]
				}
				wMat[i][col] = sum / lMat[i][i]
			}
			for i := m - 1; i >= 0; i-- {
				sum := wMat[i][col]
				for j := i + 1; j < m; j++ {
					sum -= lMat[j][i] * wMat[j][col]
				}
				wMat[i][col] = sum / lMat[i][i]
			}
		}

		for i := 0; i < 3; i++ {
			upd := 0.0
			for j := 0; j < m; j++ {
				upd += wMat[j][i] * v[j]
			}
			xFilt[i] = xPred[i] + upd
		}

		// Joseph form: I - K H = I - W^T H (3x3)
		var ikh [3][3]float64
		for i := 0; i < 3; i++ {
			for k := 0; k < 3; k++ {
				kh := 0.0
				for j := 0; j < m; j++ {
					kh += wMat[j][i] * hMat[j][k]
				}
				if i == k {
					ikh[i][k] = 1.0 - kh
				} else {
					ikh[i][k] = -kh
				}
			}
		}
		var tmpP [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				tmpP[i][j] = ikh[i][0]*pPred[0][j] + ikh[i][1]*pPred[1][j] + ikh[i][2]*pPred[2][j]
			}
		}
		for i := 0; i < 3; i++ {
			for j := i; j < 3; j++ {
				val := tmpP[i][0]*ikh[j][0] + tmpP[i][1]*ikh[j][1] + tmpP[i][2]*ikh[j][2]
				kk := 0.0
				for r := 0; r < m; r++ {
					kk += wMat[r][i] * wMat[r][j]
				}
				val += rVar * kk
				pFilt[i][j], pFilt[j][i] = val, val
			}
		}
	}
	return totalLL, nil
}

// FilterAndSmooth menjalankan filter Kalman maju serta RTS smoother mundur penuh.
func FilterAndSmooth(spec ModelSpec, taus []float64, dt float64, yields [][]float64) (FilterOutput, error) {
	m := len(taus)
	tSteps := len(yields)
	if m == 0 || tSteps == 0 {
		return FilterOutput{}, fmt.Errorf("FilterAndSmooth: empty inputs")
	}
	aVec, hMat := buildObsModel(spec, taus)
	phi := spec.AFNS.TransitionMatrixClosed(dt)
	qDt := spec.AFNS.TransitionCovClosed(dt)
	qInf := spec.AFNS.StationaryCov()
	theta := spec.AFNS.ThetaP
	rVar := spec.SigmaEps * spec.SigmaEps

	out := FilterOutput{
		XPred:   make([][3]float64, tSteps),
		PPred:   make([][3][3]float64, tSteps),
		XFilt:   make([][3]float64, tSteps),
		PFilt:   make([][3][3]float64, tSteps),
		XSmooth: make([][3]float64, tSteps),
		PSmooth: make([][3][3]float64, tSteps),
	}

	fMat := lag.Zeros(m, m)
	v := make([]float64, m)
	hpCol := make([]float64, m)

	for t := 0; t < tSteps; t++ {
		var xPred [3]float64
		var pPred [3][3]float64
		if t == 0 {
			xPred = theta
			pPred = qInf
		} else {
			xpPrev := out.XFilt[t-1]
			ppPrev := out.PFilt[t-1]
			for i := 0; i < 3; i++ {
				xPred[i] = theta[i]
				for k := 0; k < 3; k++ {
					xPred[i] += phi[i][k] * (xpPrev[k] - theta[k])
				}
			}
			var tmp [3][3]float64
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					tmp[i][j] = phi[i][0]*ppPrev[0][j] + phi[i][1]*ppPrev[1][j] + phi[i][2]*ppPrev[2][j]
				}
			}
			for i := 0; i < 3; i++ {
				for j := i; j < 3; j++ {
					val := qDt[i][j] + tmp[i][0]*phi[j][0] + tmp[i][1]*phi[j][1] + tmp[i][2]*phi[j][2]
					pPred[i][j], pPred[j][i] = val, val
				}
			}
		}
		out.XPred[t] = xPred
		out.PPred[t] = pPred

		hp := make([][3]float64, m)
		for j := 0; j < m; j++ {
			v[j] = yields[t][j] - aVec[j] - (hMat[j][0]*xPred[0] + hMat[j][1]*xPred[1] + hMat[j][2]*xPred[2])
			for k := 0; k < 3; k++ {
				hp[j][k] = hMat[j][0]*pPred[0][k] + hMat[j][1]*pPred[1][k] + hMat[j][2]*pPred[2][k]
			}
		}
		for i := 0; i < m; i++ {
			for j := i; j < m; j++ {
				s := hp[i][0]*hMat[j][0] + hp[i][1]*hMat[j][1] + hp[i][2]*hMat[j][2]
				if i == j {
					s += rVar
				}
				fMat[i][j], fMat[j][i] = s, s
			}
		}
		cholF, err := lag.Cholesky(fMat)
		if err != nil {
			return out, err
		}
		z := lag.SolveCholeskyLower(cholF, v)
		quad := lag.Dot(z, z)
		out.LogLikelihood += -0.5 * (float64(m)*log2Pi + lag.LogDetCholesky(cholF) + quad)

		// Gain K = (F^{-1} H P)^T -> W = F^{-1} (H P)
		wMat := make([][3]float64, m)
		for col := 0; col < 3; col++ {
			for i := 0; i < m; i++ {
				hpCol[i] = hp[i][col]
			}
			sol := lag.SolveCholesky(cholF, hpCol)
			for i := 0; i < m; i++ {
				wMat[i][col] = sol[i]
			}
		}

		var xFilt [3]float64
		for i := 0; i < 3; i++ {
			upd := 0.0
			for j := 0; j < m; j++ {
				upd += wMat[j][i] * v[j]
			}
			xFilt[i] = xPred[i] + upd
		}

		var ikh [3][3]float64
		for i := 0; i < 3; i++ {
			for k := 0; k < 3; k++ {
				kh := 0.0
				for j := 0; j < m; j++ {
					kh += wMat[j][i] * hMat[j][k]
				}
				if i == k {
					ikh[i][k] = 1.0 - kh
				} else {
					ikh[i][k] = -kh
				}
			}
		}
		var tmpP, pFilt [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				tmpP[i][j] = ikh[i][0]*pPred[0][j] + ikh[i][1]*pPred[1][j] + ikh[i][2]*pPred[2][j]
			}
		}
		for i := 0; i < 3; i++ {
			for j := i; j < 3; j++ {
				val := tmpP[i][0]*ikh[j][0] + tmpP[i][1]*ikh[j][1] + tmpP[i][2]*ikh[j][2]
				kk := 0.0
				for r := 0; r < m; r++ {
					kk += wMat[r][i] * wMat[r][j]
				}
				val += rVar * kk
				pFilt[i][j], pFilt[j][i] = val, val
			}
		}
		out.XFilt[t] = xFilt
		out.PFilt[t] = pFilt
	}

	// Pass mundur Rauch-Tung-Striebel (RTS)
	out.XSmooth[tSteps-1] = out.XFilt[tSteps-1]
	out.PSmooth[tSteps-1] = out.PFilt[tSteps-1]
	pNextSlice := lag.Zeros(3, 3)
	for t := tSteps - 2; t >= 0; t-- {
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				pNextSlice[i][j] = out.PPred[t+1][i][j]
			}
		}
		invPNext, err := lag.InvertLU(pNextSlice)
		if err != nil {
			return out, err
		}
		// G_t = P_{t|t} * Phi^T * P_{t+1|t}^{-1}
		var pPhiT, gMat [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				pPhiT[i][j] = out.PFilt[t][i][0]*phi[j][0] + out.PFilt[t][i][1]*phi[j][1] + out.PFilt[t][i][2]*phi[j][2]
			}
		}
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				gMat[i][j] = pPhiT[i][0]*invPNext[0][j] + pPhiT[i][1]*invPNext[1][j] + pPhiT[i][2]*invPNext[2][j]
			}
		}

		var dx [3]float64
		for i := 0; i < 3; i++ {
			dx[i] = out.XSmooth[t+1][i] - out.XPred[t+1][i]
		}
		var xs [3]float64
		for i := 0; i < 3; i++ {
			xs[i] = out.XFilt[t][i] + gMat[i][0]*dx[0] + gMat[i][1]*dx[1] + gMat[i][2]*dx[2]
		}
		out.XSmooth[t] = xs

		var dP, tmpG, ps [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				dP[i][j] = out.PSmooth[t+1][i][j] - out.PPred[t+1][i][j]
			}
		}
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				tmpG[i][j] = gMat[i][0]*dP[0][j] + gMat[i][1]*dP[1][j] + gMat[i][2]*dP[2][j]
			}
		}
		for i := 0; i < 3; i++ {
			for j := i; j < 3; j++ {
				val := out.PFilt[t][i][j] + tmpG[i][0]*gMat[j][0] + tmpG[i][1]*gMat[j][1] + tmpG[i][2]*gMat[j][2]
				ps[i][j], ps[j][i] = val, val
			}
		}
		out.PSmooth[t] = ps
	}
	return out, nil
}

// JointGaussianAnchor memverifikasi log-likelihood Kalman dan smoother RTS terhadap
// distribusi Gaussian gabungan penuh eksak untuk sampel pendek (mis. T=5, 3 tenor).
type JointGaussianAnchor struct {
	T                    int
	NumTenors            int
	KalmanLogLik         float64
	JointLogLik          float64
	LogLikAbsDiff        float64
	MaxSmoothMeanAbsDiff float64
	MaxSmoothCovAbsDiff  float64
}

func matPow3(phi [3][3]float64, k int) [3][3]float64 {
	res := [3][3]float64{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
	for step := 0; step < k; step++ {
		var nxt [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				nxt[i][j] = res[i][0]*phi[0][j] + res[i][1]*phi[1][j] + res[i][2]*phi[2][j]
			}
		}
		res = nxt
	}
	return res
}

// VerifyJointGaussianAnchor membangun matriks kovarians gabungan penuh atas seluruh observasi Y
// (dimensi T*M) dan keadaan X (dimensi 3*T), menghitung log-likelihood serta E[X|Y], Cov(X|Y)
// lewat Cholesky penuh, dan membandingkan terhadap Kalman Filter + RTS Smoother.
func VerifyJointGaussianAnchor(spec ModelSpec, taus []float64, dt float64, yields [][]float64) (JointGaussianAnchor, error) {
	tSteps := len(yields)
	m := len(taus)
	dimY := tSteps * m
	dimX := tSteps * 3

	aVec, hMat := buildObsModel(spec, taus)
	phi := spec.AFNS.TransitionMatrixClosed(dt)
	qInf := spec.AFNS.StationaryCov()
	theta := spec.AFNS.ThetaP
	rVar := spec.SigmaEps * spec.SigmaEps

	muX := make([]float64, dimX)
	muY := make([]float64, dimY)
	yVec := make([]float64, dimY)
	for t := 0; t < tSteps; t++ {
		for i := 0; i < 3; i++ {
			muX[3*t+i] = theta[i]
		}
		for j := 0; j < m; j++ {
			muY[m*t+j] = aVec[j] + hMat[j][0]*theta[0] + hMat[j][1]*theta[1] + hMat[j][2]*theta[2]
			yVec[m*t+j] = yields[t][j]
		}
	}

	sigmaXX := lag.Zeros(dimX, dimX)
	for t := 0; t < tSteps; t++ {
		for s := 0; s <= t; s++ {
			phiPow := matPow3(phi, t-s)
			var blk [3][3]float64
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					blk[i][j] = phiPow[i][0]*qInf[0][j] + phiPow[i][1]*qInf[1][j] + phiPow[i][2]*qInf[2][j]
				}
			}
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					sigmaXX[3*t+i][3*s+j] = blk[i][j]
					if t != s {
						sigmaXX[3*s+j][3*t+i] = blk[i][j]
					}
				}
			}
		}
	}

	// Sigma_XY = Sigma_XX * H_block^T (dimX x dimY), Sigma_YY = H_block * Sigma_XX * H_block^T + R_block
	sigmaXY := lag.Zeros(dimX, dimY)
	for row := 0; row < dimX; row++ {
		for s := 0; s < tSteps; s++ {
			for j := 0; j < m; j++ {
				sigmaXY[row][m*s+j] = sigmaXX[row][3*s+0]*hMat[j][0] +
					sigmaXX[row][3*s+1]*hMat[j][1] +
					sigmaXX[row][3*s+2]*hMat[j][2]
			}
		}
	}

	sigmaYY := lag.Zeros(dimY, dimY)
	for t := 0; t < tSteps; t++ {
		for i := 0; i < m; i++ {
			for s := 0; s < tSteps; s++ {
				for j := 0; j < m; j++ {
					val := hMat[i][0]*sigmaXY[3*t+0][m*s+j] +
						hMat[i][1]*sigmaXY[3*t+1][m*s+j] +
						hMat[i][2]*sigmaXY[3*t+2][m*s+j]
					if t == s && i == j {
						val += rVar
					}
					sigmaYY[m*t+i][m*s+j] = val
				}
			}
		}
	}

	cholYY, err := lag.Cholesky(sigmaYY)
	if err != nil {
		return JointGaussianAnchor{}, err
	}
	diffY := make([]float64, dimY)
	for i := 0; i < dimY; i++ {
		diffY[i] = yVec[i] - muY[i]
	}
	zY := lag.SolveCholeskyLower(cholYY, diffY)
	jointLL := -0.5 * (float64(dimY)*log2Pi + lag.LogDetCholesky(cholYY) + lag.Dot(zY, zY))

	// Kondisional Gaussian gabungan:
	// E[X | Y] = mu_X + Sigma_XY * Sigma_YY^{-1} (Y - mu_Y)
	// Cov(X | Y) = Sigma_XX - Sigma_XY * Sigma_YY^{-1} * Sigma_XY^T
	alpha := lag.SolveCholesky(cholYY, diffY)
	condMean := make([]float64, dimX)
	for i := 0; i < dimX; i++ {
		condMean[i] = muX[i] + lag.Dot(sigmaXY[i], alpha)
	}
	invYY := lag.InvertCholesky(cholYY)
	covPost := lag.Sub(sigmaXX, lag.MatMul(sigmaXY, lag.MatMul(invYY, lag.Transpose(sigmaXY))))

	filtOut, err := FilterAndSmooth(spec, taus, dt, yields)
	if err != nil {
		return JointGaussianAnchor{}, err
	}

	maxMeanDiff, maxCovDiff := 0.0, 0.0
	for t := 0; t < tSteps; t++ {
		for i := 0; i < 3; i++ {
			if d := math.Abs(filtOut.XSmooth[t][i] - condMean[3*t+i]); d > maxMeanDiff {
				maxMeanDiff = d
			}
			for j := 0; j < 3; j++ {
				if d := math.Abs(filtOut.PSmooth[t][i][j] - covPost[3*t+i][3*t+j]); d > maxCovDiff {
					maxCovDiff = d
				}
			}
		}
	}

	return JointGaussianAnchor{
		T:                    tSteps,
		NumTenors:            m,
		KalmanLogLik:         filtOut.LogLikelihood,
		JointLogLik:          jointLL,
		LogLikAbsDiff:        math.Abs(filtOut.LogLikelihood - jointLL),
		MaxSmoothMeanAbsDiff: maxMeanDiff,
		MaxSmoothCovAbsDiff:  maxCovDiff,
	}, nil
}

// ParamNames mengembalikan 11 nama parameter bebas untuk pelaporan MLE.
var ParamNames = [11]string{
	"lambda",
	"kp11", "kp22", "kp33",
	"theta1", "theta2", "theta3",
	"sigma1", "sigma2", "sigma3",
	"sigma_eps",
}

// SpecToVec mengonversi ModelSpec ke vektor parameter asli panjang 11.
func SpecToVec(spec ModelSpec) [11]float64 {
	return [11]float64{
		spec.AFNS.Lambda,
		spec.AFNS.KP[0][0], spec.AFNS.KP[1][1], spec.AFNS.KP[2][2],
		spec.AFNS.ThetaP[0], spec.AFNS.ThetaP[1], spec.AFNS.ThetaP[2],
		spec.AFNS.Sigma[0], spec.AFNS.Sigma[1], spec.AFNS.Sigma[2],
		spec.SigmaEps,
	}
}

// VecToSpec membangun ModelSpec dari vektor parameter asli panjang 11.
func VecToSpec(v [11]float64, useAFNSAdj bool) ModelSpec {
	return ModelSpec{
		AFNS: afns.Params{
			Lambda: v[0],
			KP: [3][3]float64{
				{v[1], 0, 0},
				{0, v[2], 0},
				{0, 0, v[3]},
			},
			ThetaP: [3]float64{v[4], v[5], v[6]},
			Sigma:  [3]float64{v[7], v[8], v[9]},
		},
		SigmaEps:   v[10],
		UseAFNSAdj: useAFNSAdj,
	}
}

func packUnconstrained(v [11]float64) []float64 {
	return []float64{
		math.Log(v[0]),
		math.Log(v[1]), math.Log(v[2]), math.Log(v[3]),
		v[4] / 0.02, v[5] / 0.02, v[6] / 0.02,
		math.Log(v[7] / 0.01), math.Log(v[8] / 0.01), math.Log(v[9] / 0.01),
		math.Log(v[10] / 0.0005),
	}
}

func unpackUnconstrained(u []float64) [11]float64 {
	return [11]float64{
		math.Exp(u[0]),
		math.Exp(u[1]), math.Exp(u[2]), math.Exp(u[3]),
		u[4] * 0.02, u[5] * 0.02, u[6] * 0.02,
		0.01 * math.Exp(u[7]), 0.01 * math.Exp(u[8]), 0.01 * math.Exp(u[9]),
		0.0005 * math.Exp(u[10]),
	}
}

func uBounds() (lo, hi []float64) {
	loV := [11]float64{0.05, 5e-4, 5e-4, 5e-4, -0.5, -0.5, -0.5, 1e-5, 1e-5, 1e-5, 1e-6}
	hiV := [11]float64{5.0, 5.0, 5.0, 5.0, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.05}
	lo = packUnconstrained(loV)
	hi = packUnconstrained(hiV)
	lo[4], lo[5], lo[6] = -0.5/0.02, -0.5/0.02, -0.5/0.02
	hi[4], hi[5], hi[6] = 0.5/0.02, 0.5/0.02, 0.5/0.02
	return lo, hi
}

func jacobianDiag(u []float64) [11]float64 {
	v := unpackUnconstrained(u)
	return [11]float64{
		v[0],
		v[1], v[2], v[3],
		0.02, 0.02, 0.02,
		v[7], v[8], v[9],
		v[10],
	}
}

// PerturbInitialGuess menghasilkan tebakan awal acak yang bukan nilai benar memakai pengali [0.5, 1.8].
func PerturbInitialGuess(trueSpec ModelSpec, rg *rng.RNG) ModelSpec {
	v := SpecToVec(trueSpec)
	var initV [11]float64
	for i := 0; i < 11; i++ {
		mult := rg.UniformRange(0.5, 1.8)
		if math.Abs(mult-1.0) < 0.15 {
			if mult < 1.0 {
				mult = 0.72
			} else {
				mult = 1.35
			}
		}
		initV[i] = v[i] * mult
	}
	return VecToSpec(initV, trueSpec.UseAFNSAdj)
}

// MLEResult menyimpan hasil optimisasi MLE BFGS, verifikasi skor, SE Hessian, dan diagnostik.
type MLEResult struct {
	InitParams         [11]float64
	EstParams          [11]float64
	TrueParams         [11]float64
	InitLogLik         float64
	TrueLogLik         float64
	EstLogLik          float64
	GradRelErrInit     float64
	GradRelErrTrue     float64
	GradAbsDiffOpt     float64
	GradNormOpt        float64
	HessianSE          [11]float64
	HessianMinEig      float64
	HessianRegularized bool
	Iterations         int
	Converged          bool
	HitBound           bool
}

// EstimateMLE mengestimasi 11 parameter model via BFGS dari tebakan awal initSpec,
// memverifikasi gradien beda hingga vs Richardson, dan menghitung SE dari invers Hessian.
func EstimateMLE(initSpec, trueSpec ModelSpec, taus []float64, dt float64, yields [][]float64) MLEResult {
	useAFNS := initSpec.UseAFNSAdj
	negLL := func(u []float64) float64 {
		v := unpackUnconstrained(u)
		sp := VecToSpec(v, useAFNS)
		ll, err := FastLogLikelihood(sp, taus, dt, yields)
		if err != nil || math.IsNaN(ll) || math.IsInf(ll, 0) {
			return 1e30
		}
		return -ll
	}

	uInit := packUnconstrained(SpecToVec(initSpec))
	uTrue := packUnconstrained(SpecToVec(trueSpec))
	lo, hi := uBounds()
	_, _, gradErrInit := optim.VerifyGradVsRichardson(negLL, uInit, nil)
	_, _, gradErrTrue := optim.VerifyGradVsRichardson(negLL, uTrue, nil)

	initLL, _ := FastLogLikelihood(initSpec, taus, dt, yields)
	trueLL, _ := FastLogLikelihood(trueSpec, taus, dt, yields)

	res1 := optim.MinimizeBFGS(negLL, uInit, lo, hi, 200, 1e-6)
	res2 := optim.MinimizeBFGS(negLL, res1.X, lo, hi, 150, 1e-6)
	best := res1
	if res2.Fun <= res1.Fun {
		best = res2
	}
	totalIters := res1.Iterations + res2.Iterations

	cdOpt, richOpt, _ := optim.VerifyGradVsRichardson(negLL, best.X, nil)
	gradAbsDiffOpt := 0.0
	for i := range cdOpt {
		d := cdOpt[i] - richOpt[i]
		gradAbsDiffOpt += d * d
	}
	gradAbsDiffOpt = math.Sqrt(gradAbsDiffOpt)
	estV := unpackUnconstrained(best.X)

	// Hitung matriks informasi observasi (-nabla^2 log L = nabla^2 negLL) dalam koordinat u,
	// lalu transformasi kovarians via Jacobian diagonal J: Cov(psi) = J * H_u^{-1} * J^T.
	hessU := optim.CentralDiffHessian(negLL, best.X, nil)
	vals, vecs, _, err := lag.SymmetricEigenJacobi(hessU, 80, 1e-13)
	var ses [11]float64
	minEig := 0.0
	reg := false
	if err == nil && len(vals) == 11 {
		minEig = vals[10]
		invHessU := lag.Zeros(11, 11)
		const eigFloor = 1e-6
		for k := 0; k < 11; k++ {
			ev := vals[k]
			if ev < eigFloor {
				ev = eigFloor
				reg = true
			}
			invEv := 1.0 / ev
			for i := 0; i < 11; i++ {
				for j := 0; j < 11; j++ {
					invHessU[i][j] += invEv * vecs[k][i] * vecs[k][j]
				}
			}
		}
		jac := jacobianDiag(best.X)
		for i := 0; i < 11; i++ {
			varPsi := jac[i] * jac[i] * invHessU[i][i]
			if varPsi > 0 {
				ses[i] = math.Sqrt(varPsi)
			}
		}
	}

	return MLEResult{
		InitParams:         SpecToVec(initSpec),
		EstParams:          estV,
		TrueParams:         SpecToVec(trueSpec),
		InitLogLik:         initLL,
		TrueLogLik:         trueLL,
		EstLogLik:          -best.Fun,
		GradRelErrInit:     gradErrInit,
		GradRelErrTrue:     gradErrTrue,
		GradAbsDiffOpt:     gradAbsDiffOpt,
		GradNormOpt:        best.GradNorm,
		HessianSE:          ses,
		HessianMinEig:      minEig,
		HessianRegularized: reg,
		Iterations:         totalIters,
		Converged:          best.Converged,
		HitBound:           best.HitBound,
	}
}

// ParamReplicationSummary merangkum bias, SD, RMSE, SE Hessian, dan cakupan 95% satu parameter.
type ParamReplicationSummary struct {
	Name          string
	TrueValue     float64
	MeanEstimate  float64
	Bias          float64
	RepSE         float64
	BiasInRepSE   float64
	EmpiricalSD   float64
	RMSE          float64
	MeanHessianSE float64
	Coverage95    float64
	CoveredCount  int
	TotalReps     int
}

// ReplicationStudySummary merangkum seluruh studi replikasi MLE R simulasi.
type ReplicationStudySummary struct {
	NumReps        int
	TSteps         int
	ConvergedCount int
	LLGeTrueCount  int
	BoundHitCount  int
	Params         [11]ParamReplicationSummary
	RepEstMatrix   [][]float64
}

// RunReplicationStudy menjalankan R simulasi independen panjang T, mengestimasi MLE 11 parameter
// secara paralel deterministik, dan menghitung statistik replikasi serta cakupan 95%.
func RunReplicationStudy(trueSpec ModelSpec, taus []float64, dt float64, tSteps, nRep, nChunks, nWorkers int, seed uint64) ReplicationStudySummary {
	trueV := SpecToVec(trueSpec)
	perChunk := nRep / nChunks
	rem := nRep % nChunks

	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) []MLEResult {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		out := make([]MLEResult, cnt)
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x3401, uint64(c)))
		for r := 0; r < cnt; r++ {
			_, yields := afns.SimulateStationaryPath(trueSpec.AFNS, taus, tSteps, dt, trueSpec.SigmaEps, rg)
			initSpec := PerturbInitialGuess(trueSpec, rg)
			out[r] = EstimateMLE(initSpec, trueSpec, taus, dt, yields)
		}
		return out
	})

	all := make([]MLEResult, 0, nRep)
	for _, ch := range chunks {
		all = append(all, ch...)
	}

	convCnt, llGeCnt, hitCnt := 0, 0, 0
	estMatrix := make([][]float64, nRep)
	for r, res := range all {
		if res.Converged {
			convCnt++
		}
		if res.EstLogLik >= res.TrueLogLik-1e-6 {
			llGeCnt++
		}
		if res.HitBound {
			hitCnt++
		}
		row := make([]float64, 11)
		copy(row, res.EstParams[:])
		estMatrix[r] = row
	}

	var paramSummaries [11]ParamReplicationSummary
	for pIdx := 0; pIdx < 11; pIdx++ {
		ests := make([]float64, nRep)
		ses := make([]float64, nRep)
		cov := 0
		sqErr := 0.0
		for r, res := range all {
			e := res.EstParams[pIdx]
			s := res.HessianSE[pIdx]
			ests[r] = e
			ses[r] = s
			d := e - trueV[pIdx]
			sqErr += d * d
			if s > 0 && math.Abs(d) <= 1.96*s {
				cov++
			}
		}
		meanEst := stats.Mean(ests)
		empSD := stats.SD(ests)
		repSE := empSD / math.Sqrt(float64(nRep))
		bias := meanEst - trueV[pIdx]
		paramSummaries[pIdx] = ParamReplicationSummary{
			Name:          ParamNames[pIdx],
			TrueValue:     trueV[pIdx],
			MeanEstimate:  meanEst,
			Bias:          bias,
			RepSE:         repSE,
			BiasInRepSE:   math.Abs(bias) / math.Max(repSE, 1e-30),
			EmpiricalSD:   empSD,
			RMSE:          math.Sqrt(sqErr / float64(nRep)),
			MeanHessianSE: stats.Mean(ses),
			Coverage95:    float64(cov) / float64(nRep),
			CoveredCount:  cov,
			TotalReps:     nRep,
		}
	}

	return ReplicationStudySummary{
		NumReps:        nRep,
		TSteps:         tSteps,
		ConvergedCount: convCnt,
		LLGeTrueCount:  llGeCnt,
		BoundHitCount:  hitCnt,
		Params:         paramSummaries,
		RepEstMatrix:   estMatrix,
	}
}

// MisspecificationReport membandingkan estimasi DNS (benar) vs AFNS (salah spesifikasi)
// ketika data dibangkitkan dari Dynamic Nelson-Siegel (tanpa suku penyesuaian -A/tau).
type MisspecificationReport struct {
	DNSLogLik      float64
	AFNSLogLik     float64
	LogLikDrop     float64
	DNSEst         [11]float64
	AFNSEst        [11]float64
	LambdaBiasDNS  float64
	LambdaBiasAFNS float64
	Theta1BiasDNS  float64
	Theta1BiasAFNS float64
}

// RunDNSMisspecification membangkitkan data dari Dynamic NS (tanpa -A/tau) dan mengestimasi DNS vs AFNS.
func RunDNSMisspecification(trueSpec ModelSpec, taus []float64, dt float64, tSteps int, seed uint64) MisspecificationReport {
	rg := rng.New(rng.DeriveChunkSeed(seed, 0x3501, 0))
	_, yields := afns.SimulateDynamicNSPath(trueSpec.AFNS, taus, tSteps, dt, trueSpec.SigmaEps, rg)

	initBase := PerturbInitialGuess(trueSpec, rg)
	dnsTrue := trueSpec
	dnsTrue.UseAFNSAdj = false
	dnsInit := initBase
	dnsInit.UseAFNSAdj = false

	afnsTrue := trueSpec
	afnsTrue.UseAFNSAdj = true
	afnsInit := initBase
	afnsInit.UseAFNSAdj = true

	resDNS := EstimateMLE(dnsInit, dnsTrue, taus, dt, yields)
	resAFNS := EstimateMLE(afnsInit, afnsTrue, taus, dt, yields)

	return MisspecificationReport{
		DNSLogLik:      resDNS.EstLogLik,
		AFNSLogLik:     resAFNS.EstLogLik,
		LogLikDrop:     resDNS.EstLogLik - resAFNS.EstLogLik,
		DNSEst:         resDNS.EstParams,
		AFNSEst:        resAFNS.EstParams,
		LambdaBiasDNS:  resDNS.EstParams[0] - trueSpec.AFNS.Lambda,
		LambdaBiasAFNS: resAFNS.EstParams[0] - trueSpec.AFNS.Lambda,
		Theta1BiasDNS:  resDNS.EstParams[4] - trueSpec.AFNS.ThetaP[0],
		Theta1BiasAFNS: resAFNS.EstParams[4] - trueSpec.AFNS.ThetaP[0],
	}
}
