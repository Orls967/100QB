package kalman

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/afns"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/optim"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/rng"
)

// StandardTenors adalah 10 jatuh tempo standar (tahun).
var StandardTenors = []float64{0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30}

const (
	WeeklyDt        = 1.0 / 52.0
	DefaultSigmaEps = 0.0005 // 5 bp
	log2Pi          = 1.83787706640934548356
)

// ModelSpec menggabungkan parameter AFNS dan deviasi baku derau pengukuran.
type ModelSpec struct {
	AFNS       afns.Params
	SigmaEps   float64
	UseAFNSAdj bool
}

// DefaultSpec mengembalikan spesifikasi AFNS acuan yang sama dengan day004b.
func DefaultSpec() ModelSpec {
	return ModelSpec{
		AFNS:       afns.DefaultParams(),
		SigmaEps:   DefaultSigmaEps,
		UseAFNSAdj: true,
	}
}

// FilterOutput menyimpan state prediksi, terfilter, smoother RTS, dan log-likelihood.
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

// FastLogLikelihood mengevaluasi log-likelihood inovasi filter Kalman bebas alokasi heap
// dan mengaktifkan pembekuan steady-state Riccati saat max|dP| <= 1e-15*max|P| (setelah t >= 5).
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
	var prevPPred [3][3]float64
	var xFilt [3]float64
	var pFilt [3][3]float64

	var fMat [16][16]float64
	var lMat [16][16]float64
	var invDiagL [16]float64
	var hp [16][3]float64
	var wMat [16][3]float64
	var v [16]float64
	var z [16]float64

	totalLL := 0.0
	constHalfM := float64(m) * log2Pi
	var logDet float64
	steady := false

	for t := 0; t < tSteps; t++ {
		if t > 0 {
			dx0 := xFilt[0] - theta[0]
			dx1 := xFilt[1] - theta[1]
			dx2 := xFilt[2] - theta[2]
			for i := 0; i < 3; i++ {
				xPred[i] = theta[i] + phi[i][0]*dx0 + phi[i][1]*dx1 + phi[i][2]*dx2
			}
			if !steady {
				prevPPred = pPred
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
		}

		yt := yields[t]
		for j := 0; j < m; j++ {
			val := yt[j]
			if math.IsNaN(val) || math.IsInf(val, 0) {
				return math.NaN(), fmt.Errorf("non-finite observation at (%d,%d)", t, j)
			}
			predY := aVec[j] + hMat[j][0]*xPred[0] + hMat[j][1]*xPred[1] + hMat[j][2]*xPred[2]
			v[j] = val - predY
		}

		if !steady {
			for j := 0; j < m; j++ {
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

			logDet = 0.0
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
						invDiagL[i] = 1.0 / diag
						logDet += math.Log(diag)
					} else {
						lMat[i][j] = sum * invDiagL[j]
					}
				}
			}
			logDet *= 2.0
		}

		quad := 0.0
		for i := 0; i < m; i++ {
			sum := v[i]
			for j := 0; j < i; j++ {
				sum -= lMat[i][j] * z[j]
			}
			zi := sum * invDiagL[i]
			z[i] = zi
			quad += zi * zi
		}
		totalLL += -0.5 * (constHalfM + logDet + quad)

		if t+1 == tSteps {
			break
		}

		if !steady {
			for col := 0; col < 3; col++ {
				for i := 0; i < m; i++ {
					sum := hp[i][col]
					for j := 0; j < i; j++ {
						sum -= lMat[i][j] * wMat[j][col]
					}
					wMat[i][col] = sum * invDiagL[i]
				}
				for i := m - 1; i >= 0; i-- {
					sum := wMat[i][col]
					for j := i + 1; j < m; j++ {
						sum -= lMat[j][i] * wMat[j][col]
					}
					wMat[i][col] = sum * invDiagL[i]
				}
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

			if t >= 5 {
				maxDiff, maxScale := 0.0, 0.0
				for i := 0; i < 3; i++ {
					if d := math.Abs(pPred[i][i]); d > maxScale {
						maxScale = d
					}
					for j := 0; j < 3; j++ {
						if d := math.Abs(pPred[i][j] - prevPPred[i][j]); d > maxDiff {
							maxDiff = d
						}
					}
				}
				if maxScale > 0 && maxDiff <= 1e-15*maxScale {
					steady = true
				}
			}
		}

		for i := 0; i < 3; i++ {
			upd := 0.0
			for j := 0; j < m; j++ {
				upd += wMat[j][i] * v[j]
			}
			xFilt[i] = xPred[i] + upd
		}
	}
	return totalLL, nil
}

// FilterAndSmooth menjalankan filter Kalman Joseph penuh dan smoother RTS mundur.
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

	fMat := zerosMat(m, m)
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
		cholF, err := choleskyMat(fMat)
		if err != nil {
			return out, err
		}
		z := solveCholeskyLower(cholF, v)
		quad := dotVec(z, z)
		out.LogLikelihood += -0.5 * (float64(m)*log2Pi + logDetCholesky(cholF) + quad)

		wMat := make([][3]float64, m)
		for col := 0; col < 3; col++ {
			for i := 0; i < m; i++ {
				hpCol[i] = hp[i][col]
			}
			sol := solveCholesky(cholF, hpCol)
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

	out.XSmooth[tSteps-1] = out.XFilt[tSteps-1]
	out.PSmooth[tSteps-1] = out.PFilt[tSteps-1]
	pNextSlice := zerosMat(3, 3)
	for t := tSteps - 2; t >= 0; t-- {
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				pNextSlice[i][j] = out.PPred[t+1][i][j]
			}
		}
		invPNext, err := invertLU(pNextSlice)
		if err != nil {
			return out, err
		}
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

// JointGaussianAnchor memverifikasi log-likelihood Kalman dan smoother RTS terhadap Gaussian gabungan penuh.
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

// VerifyJointGaussianAnchor mencocokkan log-likelihood dan momen kondisional Kalman/RTS
// dengan rumus Gaussian gabungan penuh pada sampel pendek.
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

	sigmaXX := zerosMat(dimX, dimX)
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

	sigmaXY := zerosMat(dimX, dimY)
	for row := 0; row < dimX; row++ {
		for s := 0; s < tSteps; s++ {
			for j := 0; j < m; j++ {
				sigmaXY[row][m*s+j] = sigmaXX[row][3*s+0]*hMat[j][0] +
					sigmaXX[row][3*s+1]*hMat[j][1] +
					sigmaXX[row][3*s+2]*hMat[j][2]
			}
		}
	}

	sigmaYY := zerosMat(dimY, dimY)
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

	cholYY, err := choleskyMat(sigmaYY)
	if err != nil {
		return JointGaussianAnchor{}, err
	}
	diffY := make([]float64, dimY)
	for i := 0; i < dimY; i++ {
		diffY[i] = yVec[i] - muY[i]
	}
	zY := solveCholeskyLower(cholYY, diffY)
	jointLL := -0.5 * (float64(dimY)*log2Pi + logDetCholesky(cholYY) + dotVec(zY, zY))

	alpha := solveCholesky(cholYY, diffY)
	condMean := make([]float64, dimX)
	for i := 0; i < dimX; i++ {
		condMean[i] = muX[i] + dotVec(sigmaXY[i], alpha)
	}
	invYY := invertCholesky(cholYY)
	covPost := matSub(sigmaXX, matMul(sigmaXY, matMul(invYY, transposeMat(sigmaXY))))

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

// ParamNames adalah nama 11 parameter bebas AFNS MLE dalam urutan standar.
var ParamNames = [11]string{
	"lambda",
	"kp11", "kp22", "kp33",
	"theta1", "theta2", "theta3",
	"sigma1", "sigma2", "sigma3",
	"sigma_eps",
}

// IsPositiveParam menandai apakah parameter pada indeks i wajib positif secara fisik (> 0).
var IsPositiveParam = [11]bool{
	true,
	true, true, true,
	false, false, false,
	true, true, true,
	true,
}

// LowerBoundsOriginal adalah batas pengaman bawah dalam skala asli.
var LowerBoundsOriginal = [11]float64{0.05, 5e-4, 5e-4, 5e-4, -0.5, -0.5, -0.5, 1e-5, 1e-5, 1e-5, 1e-6}

// UpperBoundsOriginal adalah batas pengaman atas dalam skala asli.
var UpperBoundsOriginal = [11]float64{5.0, 5.0, 5.0, 5.0, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.05}

// SpecToVec mengonversi ModelSpec ke vektor 11 parameter skala asli.
func SpecToVec(spec ModelSpec) [11]float64 {
	return [11]float64{
		spec.AFNS.Lambda,
		spec.AFNS.KP[0][0], spec.AFNS.KP[1][1], spec.AFNS.KP[2][2],
		spec.AFNS.ThetaP[0], spec.AFNS.ThetaP[1], spec.AFNS.ThetaP[2],
		spec.AFNS.Sigma[0], spec.AFNS.Sigma[1], spec.AFNS.Sigma[2],
		spec.SigmaEps,
	}
}

// VecToSpec membangun ModelSpec dari vektor 11 parameter skala asli.
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

// PackUnconstrained mentransformasikan parameter asli ke koordinat terdiskalakan u.
func PackUnconstrained(v [11]float64) []float64 {
	return []float64{
		math.Log(v[0]),
		math.Log(v[1]), math.Log(v[2]), math.Log(v[3]),
		v[4] / 0.02, v[5] / 0.02, v[6] / 0.02,
		math.Log(v[7] / 0.01), math.Log(v[8] / 0.01), math.Log(v[9] / 0.01),
		math.Log(v[10] / 0.0005),
	}
}

// UnpackUnconstrained mengembalikan parameter ke skala asli.
func UnpackUnconstrained(u []float64) [11]float64 {
	return [11]float64{
		math.Exp(u[0]),
		math.Exp(u[1]), math.Exp(u[2]), math.Exp(u[3]),
		u[4] * 0.02, u[5] * 0.02, u[6] * 0.02,
		0.01 * math.Exp(u[7]), 0.01 * math.Exp(u[8]), 0.01 * math.Exp(u[9]),
		0.0005 * math.Exp(u[10]),
	}
}

func uBounds() (lo, hi []float64) {
	lo = PackUnconstrained(LowerBoundsOriginal)
	hi = PackUnconstrained(UpperBoundsOriginal)
	lo[4], lo[5], lo[6] = LowerBoundsOriginal[4]/0.02, LowerBoundsOriginal[5]/0.02, LowerBoundsOriginal[6]/0.02
	hi[4], hi[5], hi[6] = UpperBoundsOriginal[4]/0.02, UpperBoundsOriginal[5]/0.02, UpperBoundsOriginal[6]/0.02
	return lo, hi
}

func jacobianDiag(u []float64) [11]float64 {
	v := UnpackUnconstrained(u)
	return [11]float64{
		v[0],
		v[1], v[2], v[3],
		0.02, 0.02, 0.02,
		v[7], v[8], v[9],
		v[10],
	}
}

// PerturbInitialGuess membangkitkan tebakan awal acak yang terpisah dari nilai acuan dengan pengali [0.5, 1.8].
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

// StageDiagnostic menyimpan diagnostik penghentian dan pencarian garis untuk satu tahap BFGS.
type StageDiagnostic struct {
	Stage                  int
	Iterations             int
	FuncEvals              int
	LogLik                 float64
	RawGradNorm            float64
	GradNorm               float64
	ScaledGradNorm         float64
	Converged              bool
	HitBound               bool
	TerminationReason      string
	MinAlphaTested         float64
	LastStepNorm           float64
	CandidateBitIdentical  bool
	LastObjChange          float64
	ObjBitIdentical        bool
	LineSearchFailureCause string
}

// MLEResult menyimpan hasil estimasi MLE 11 parameter beserta status finite dan diagnostik penghentian.
type MLEResult struct {
	Finite                 bool
	InitParams             [11]float64
	EstParams              [11]float64
	TrueParams             [11]float64
	InitLogLik             float64
	TrueLogLik             float64
	EstLogLik              float64
	GradRelErrInit         float64
	GradRelErrTrue         float64
	GradAbsDiffOpt         float64
	RawGradNormOpt         float64
	GradNormOpt            float64
	ScaledGradNormOpt      float64
	HessianSE              [11]float64
	HessianMinEig          float64
	HessianRegularized     bool
	Iterations             int
	BestStageIterations    int
	CumulativeIterations   int
	BestStage              int
	Stage1                 StageDiagnostic
	Stage2                 StageDiagnostic
	Converged              bool
	TerminationReason      string
	HitBound               bool
	HitBoundByParam        [11]bool
	MinAlphaTested         float64
	LastStepNorm           float64
	CandidateBitIdentical  bool
	LastObjChange          float64
	ObjBitIdentical        bool
	LineSearchFailureCause string
}

func makeStageDiag(stage int, res optim.BFGSResult) StageDiagnostic {
	return StageDiagnostic{
		Stage:                  stage,
		Iterations:             res.Iterations,
		FuncEvals:              res.FuncEvals,
		LogLik:                 -res.Fun,
		RawGradNorm:            res.RawGradNorm,
		GradNorm:               res.GradNorm,
		ScaledGradNorm:         res.ScaledGradNorm,
		Converged:              res.Converged,
		HitBound:               res.HitBound,
		TerminationReason:      res.TerminationReason,
		MinAlphaTested:         res.MinAlphaTested,
		LastStepNorm:           res.LastStepNorm,
		CandidateBitIdentical:  res.CandidateBitIdentical,
		LastObjChange:          res.LastObjChange,
		ObjBitIdentical:        res.ObjBitIdentical,
		LineSearchFailureCause: res.LineSearchFailureCause,
	}
}

func isAllFinite(v [11]float64, ll float64) bool {
	if math.IsNaN(ll) || math.IsInf(ll, 0) {
		return false
	}
	for _, x := range v {
		if math.IsNaN(x) || math.IsInf(x, 0) {
			return false
		}
	}
	return true
}

// EstimateMLE mengestimasi MLE penuh dengan verifikasi skor dan Hessian SE.
func EstimateMLE(initSpec, trueSpec ModelSpec, taus []float64, dt float64, yields [][]float64) MLEResult {
	return EstimateMLEWithOptions(initSpec, trueSpec, taus, dt, yields, true)
}

// EstimateMLEWithOptions menjalankan dua tahap BFGS pada koordinat u dengan opsi penghitungan Hessian SE.
func EstimateMLEWithOptions(initSpec, trueSpec ModelSpec, taus []float64, dt float64, yields [][]float64, computeHessian bool) MLEResult {
	useAFNS := initSpec.UseAFNSAdj
	negLL := func(u []float64) float64 {
		v := UnpackUnconstrained(u)
		sp := VecToSpec(v, useAFNS)
		ll, err := FastLogLikelihood(sp, taus, dt, yields)
		if err != nil || math.IsNaN(ll) || math.IsInf(ll, 0) {
			return 1e30
		}
		return -ll
	}

	uInit := PackUnconstrained(SpecToVec(initSpec))
	uTrue := PackUnconstrained(SpecToVec(trueSpec))
	lo, hi := uBounds()

	var gradErrInit, gradErrTrue float64
	if computeHessian {
		_, _, gradErrInit = optim.VerifyGradVsRichardson(negLL, uInit, nil)
		_, _, gradErrTrue = optim.VerifyGradVsRichardson(negLL, uTrue, nil)
	}

	initLL, initErr := FastLogLikelihood(initSpec, taus, dt, yields)
	trueLL, trueErr := FastLogLikelihood(trueSpec, taus, dt, yields)
	if initErr != nil {
		initLL = math.NaN()
	}
	if trueErr != nil {
		trueLL = math.NaN()
	}

	res1 := optim.MinimizeBFGS(negLL, uInit, lo, hi, 200, 1e-6)
	res2 := optim.MinimizeBFGS(negLL, res1.X, lo, hi, 150, 1e-6)
	best := res1
	bestStage := 1
	if res2.Fun < res1.Fun {
		best = res2
		bestStage = 2
	}
	cumulativeIters := res1.Iterations + res2.Iterations
	bestStageIters := best.Iterations

	estV := UnpackUnconstrained(best.X)
	estLL := -best.Fun
	finite := best.TerminationReason != optim.TermInvalidObjective && best.Fun < 1e29 && isAllFinite(estV, estLL)

	var gradAbsDiffOpt float64
	var ses [11]float64
	var minEig float64
	var reg bool

	if finite && computeHessian {
		cdOpt, richOpt, _ := optim.VerifyGradVsRichardson(negLL, best.X, nil)
		for i := range cdOpt {
			d := cdOpt[i] - richOpt[i]
			gradAbsDiffOpt += d * d
		}
		gradAbsDiffOpt = math.Sqrt(gradAbsDiffOpt)

		hessU := optim.CentralDiffHessian(negLL, best.X, nil)
		vals, vecs, _, err := symmetricEigenJacobi(hessU, 80, 1e-13)
		if err == nil && len(vals) == 11 {
			minEig = vals[10]
			invHessU := zerosMat(11, 11)
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
	}

	var hitByParam [11]bool
	for i := 0; i < 11 && i < len(best.HitBoundByParam); i++ {
		hitByParam[i] = best.HitBoundByParam[i]
	}

	return MLEResult{
		Finite:                 finite,
		InitParams:             SpecToVec(initSpec),
		EstParams:              estV,
		TrueParams:             SpecToVec(trueSpec),
		InitLogLik:             initLL,
		TrueLogLik:             trueLL,
		EstLogLik:              estLL,
		GradRelErrInit:         gradErrInit,
		GradRelErrTrue:         gradErrTrue,
		GradAbsDiffOpt:         gradAbsDiffOpt,
		RawGradNormOpt:         best.RawGradNorm,
		GradNormOpt:            best.GradNorm,
		ScaledGradNormOpt:      best.ScaledGradNorm,
		HessianSE:              ses,
		HessianMinEig:          minEig,
		HessianRegularized:     reg,
		Iterations:             bestStageIters,
		BestStageIterations:    bestStageIters,
		CumulativeIterations:   cumulativeIters,
		BestStage:              bestStage,
		Stage1:                 makeStageDiag(1, res1),
		Stage2:                 makeStageDiag(2, res2),
		Converged:              best.Converged,
		TerminationReason:      best.TerminationReason,
		HitBound:               best.HitBound,
		HitBoundByParam:        hitByParam,
		MinAlphaTested:         best.MinAlphaTested,
		LastStepNorm:           best.LastStepNorm,
		CandidateBitIdentical:  best.CandidateBitIdentical,
		LastObjChange:          best.LastObjChange,
		ObjBitIdentical:        best.ObjBitIdentical,
		LineSearchFailureCause: best.LineSearchFailureCause,
	}
}
