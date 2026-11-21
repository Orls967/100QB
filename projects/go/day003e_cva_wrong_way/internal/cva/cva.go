package cva

import (
	"math"
	"sort"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/cir"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/hw"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/parallel"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/special"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/stats"
)

// QuarterlyGrid membangun kisi waktu kuartalan seragam t_j = j * 0.25 untuk horizon T tahun.
func QuarterlyGrid(maturityYears int) []float64 {
	n := maturityYears*4 + 1
	grid := make([]float64, n)
	for j := 0; j < n; j++ {
		grid[j] = float64(j) * 0.25
	}
	return grid
}

// AnalyticResult menyimpan CVA semi-analitik independen dan sensitivitas analitik vs beda hingga.
type AnalyticResult struct {
	CVAAnalytic        float64
	DLambda0Analytic   float64
	DLambda0CentralFD  float64
	DLambda0AbsErr     float64
	DLambda0RelErr     float64
	CS01Lambda0_1bp    float64
	CS01Theta_1bp      float64
	CS01Parallel_1bp   float64
	DSigmaHWFixedK     float64
	SigmaBump10bpFixed float64
	DStrikeK           float64
	StrikeBump1bp      float64
}

// AnalyticCVAFromEEDisc mengevaluasi rumus semi-analitik CVA = (1-R) sum_j EE_disc(t_j) (Q(t_{j-1}) - Q(t_j)).
func AnalyticCVAFromEEDisc(cirP cir.Params, tGrid, eeDisc []float64) (cvaVal, dLambda0 float64) {
	lgd := cirP.LGD()
	qPrev := cirP.SurvivalProb(tGrid[0])
	dqPrev := cirP.SurvivalDLambda0(tGrid[0])
	for j := 1; j < len(tGrid); j++ {
		qCurr := cirP.SurvivalProb(tGrid[j])
		dqCurr := cirP.SurvivalDLambda0(tGrid[j])
		cvaVal += lgd * eeDisc[j] * (qPrev - qCurr)
		dLambda0 += lgd * eeDisc[j] * (dqPrev - dqCurr)
		qPrev = qCurr
		dqPrev = dqCurr
	}
	return cvaVal, dLambda0
}

func eeDiscProfileGH(hwP hw.Params, swap hw.SwapSpec, tGrid []float64, nQuad int) []float64 {
	out := make([]float64, len(tGrid))
	for j, t := range tGrid {
		_, disc, _ := swap.ExpectedExposureSmooth(hwP, t, nQuad)
		out[j] = disc
	}
	return out
}

// ComputeAnalyticCVAAndSensitivities menghitung CVA independen semi-analitik dan seluruh sensitivitas (lambda0, theta, sigma, K).
func ComputeAnalyticCVAAndSensitivities(hwP hw.Params, cirP cir.Params, swap hw.SwapSpec, tGrid []float64, nGH int) AnalyticResult {
	eeDisc := eeDiscProfileGH(hwP, swap, tGrid, nGH)
	cvaBase, dLamAnal := AnalyticCVAFromEEDisc(cirP, tGrid, eeDisc)

	epsLam := 1e-6
	cirUp := cirP
	cirUp.Lambda0 += epsLam
	cirDn := cirP
	cirDn.Lambda0 -= epsLam
	cvaUp, _ := AnalyticCVAFromEEDisc(cirUp, tGrid, eeDisc)
	cvaDn, _ := AnalyticCVAFromEEDisc(cirDn, tGrid, eeDisc)
	dLamFD := (cvaUp - cvaDn) / (2.0 * epsLam)
	absErr := math.Abs(dLamAnal - dLamFD)
	relErr := absErr / math.Max(math.Abs(dLamAnal), 1e-16)

	const oneBp = 1e-4
	cirLam1 := cirP
	cirLam1.Lambda0 += oneBp
	cvaLam1, _ := AnalyticCVAFromEEDisc(cirLam1, tGrid, eeDisc)

	cirTh1 := cirP
	cirTh1.Theta += oneBp
	cvaTh1, _ := AnalyticCVAFromEEDisc(cirTh1, tGrid, eeDisc)

	cirPar1 := cirP
	cirPar1.Lambda0 += oneBp
	cirPar1.Theta += oneBp
	cvaPar1, _ := AnalyticCVAFromEEDisc(cirPar1, tGrid, eeDisc)

	// Sensitivitas sigma_HW menjaga strike swap K tetap (mengisolasi efek opsionalitas CVA).
	epsSig := 1e-4
	hwSigUp := hwP
	hwSigUp.Sigma += epsSig
	hwSigDn := hwP
	hwSigDn.Sigma -= epsSig
	cvaSigUp, _ := AnalyticCVAFromEEDisc(cirP, tGrid, eeDiscProfileGH(hwSigUp, swap, tGrid, nGH))
	cvaSigDn, _ := AnalyticCVAFromEEDisc(cirP, tGrid, eeDiscProfileGH(hwSigDn, swap, tGrid, nGH))
	dSigma := (cvaSigUp - cvaSigDn) / (2.0 * epsSig)

	hwSig10 := hwP
	hwSig10.Sigma += 0.001 // +10 bps
	cvaSig10, _ := AnalyticCVAFromEEDisc(cirP, tGrid, eeDiscProfileGH(hwSig10, swap, tGrid, nGH))

	// Sensitivitas terhadap tingkat tetap swap K (+- 1 bp).
	swapUp := swap
	swapUp.StrikeK += oneBp
	swapDn := swap
	swapDn.StrikeK -= oneBp
	cvaKUp, _ := AnalyticCVAFromEEDisc(cirP, tGrid, eeDiscProfileGH(hwP, swapUp, tGrid, nGH))
	cvaKDn, _ := AnalyticCVAFromEEDisc(cirP, tGrid, eeDiscProfileGH(hwP, swapDn, tGrid, nGH))
	dK := (cvaKUp - cvaKDn) / (2.0 * oneBp)

	return AnalyticResult{
		CVAAnalytic:        cvaBase,
		DLambda0Analytic:   dLamAnal,
		DLambda0CentralFD:  dLamFD,
		DLambda0AbsErr:     absErr,
		DLambda0RelErr:     relErr,
		CS01Lambda0_1bp:    cvaLam1 - cvaBase,
		CS01Theta_1bp:      cvaTh1 - cvaBase,
		CS01Parallel_1bp:   cvaPar1 - cvaBase,
		DSigmaHWFixedK:     dSigma,
		SigmaBump10bpFixed: cvaSig10 - cvaBase,
		DStrikeK:           dK,
		StrikeBump1bp:      cvaKUp - cvaBase,
	}
}

// DiscretizationRow merangkum estimasi CVA MC pada resolusi sub-langkah/tahun berbeda (rho = 0).
type DiscretizationRow struct {
	SubStepsPerYear int
	CVATauBucket    float64
	SETauBucket     float64
	ZScoreTauBucket float64
	CVATauInterp    float64
	SETauInterp     float64
	CVARaoBlackwell float64
	SERaoBlackwell  float64
	ZScoreRB        float64
	AbsErrRBVsAnal  float64
	PairedBiasVs320 float64
	PairedSEVs320   float64
	Q10YMC          float64
	Q10YBias        float64
	Q10PairedBias   float64
	IntRMSEVs320    float64
}

// VerifyIndependentCVAMC menghitung CVA MC dengan sampling waktu gagal eksplisit (u ~ Exp(1) dan integral trapesium)
// serta estimator Rao-Blackwellized pada sub-langkah/tahun [20, 40, 80] memakai increment Brown bersarang terhadap referensi halus 320/tahun.
func VerifyIndependentCVAMC(hwP hw.Params, cirP cir.Params, swap hw.SwapSpec, tGrid []float64, cvaAnal float64, subStepsList []int, nPaths, nChunks, nWorkers int, seed uint64) []DiscretizationRow {
	const refM = 320
	hwMaxM := 80
	for _, m := range subStepsList {
		if m > hwMaxM {
			hwMaxM = m
		}
	}
	maxT := tGrid[len(tGrid)-1]
	nHWSteps := int(math.Round(maxT * float64(hwMaxM)))
	nRefSteps := int(math.Round(maxT * float64(refM)))
	q10Anal := cirP.SurvivalProb(maxT)
	lgd := cirP.LGD()

	rows := make([]DiscretizationRow, len(subStepsList))
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks

	type resSums struct {
		sumTauB, sumSqTauB       []float64
		sumTauI, sumSqTauI       []float64
		sumRB, sumSqRB           []float64
		sumDiffRB, sumSqDiffRB   []float64
		sumDiffQ10, sumSqIntDiff []float64
		sumQ10                   []float64
		count                    int
	}
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) resSums {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		nRes := len(subStepsList)
		rs := resSums{
			sumTauB:      make([]float64, nRes),
			sumSqTauB:    make([]float64, nRes),
			sumTauI:      make([]float64, nRes),
			sumSqTauI:    make([]float64, nRes),
			sumRB:        make([]float64, nRes),
			sumSqRB:      make([]float64, nRes),
			sumDiffRB:    make([]float64, nRes),
			sumSqDiffRB:  make([]float64, nRes),
			sumDiffQ10:   make([]float64, nRes),
			sumSqIntDiff: make([]float64, nRes),
			sumQ10:       make([]float64, nRes),
			count:        cnt,
		}
		dtHW := 1.0 / float64(hwMaxM)
		spHW := hwP.NewStepParams(dtHW)
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x701, uint64(c)))

		zLRef := make([]float64, nRefSteps)
		nQ := len(tGrid) - 1
		stepsPerQuarterHW := hwMaxM / 4
		xAll := make([]float64, nHWSteps+1)
		iAll := make([]float64, nHWSteps+1)
		xQ := make([]float64, nQ+1)
		dQ := make([]float64, nQ+1)
		vPosQ := make([]float64, nQ+1)
		dQ[0] = 1.0

		dtRef := 1.0 / float64(refM)
		sqrtDtRef := math.Sqrt(dtRef)
		stepsPerQRef := refM / 4

		for k := 0; k < cnt; k++ {
			uExp := rg.ExpFloat64()
			xt, it := 0.0, 0.0
			for s := 0; s < nHWSteps; s++ {
				z1 := rg.NormFloat64()
				z2 := rg.NormFloat64()
				xt, it = spHW.Advance(xt, it, z1, z2)
				xAll[s+1] = xt
				iAll[s+1] = it
			}
			for j := 1; j <= nQ; j++ {
				idx := j * stepsPerQuarterHW
				xQ[j] = xAll[idx]
				dQ[j] = hwP.DiscountFactor(tGrid[j], iAll[idx])
				vPosQ[j] = math.Max(swap.Value(hwP, tGrid[j], xQ[j]), 0)
			}

			// Hitung lintasan referensi halus m = 320/tahun memakai inovasi Z_ref.
			lamRef := cirP.Lambda0
			intRef := 0.0
			sPrevRef := 1.0
			rbRef := 0.0
			for s := 0; s < nRefSteps; s++ {
				zLRef[s] = rg.NormFloat64()
				nextLam := cirP.StepFullTruncationEuler(lamRef, dtRef, sqrtDtRef, zLRef[s])
				intRef += 0.5 * (lamRef + nextLam) * dtRef
				lamRef = nextLam
				if (s+1)%stepsPerQRef == 0 {
					qIdx := (s + 1) / stepsPerQRef
					sCurr := math.Exp(-intRef)
					rbRef += lgd * dQ[qIdx] * vPosQ[qIdx] * (sPrevRef - sCurr)
					sPrevRef = sCurr
				}
			}
			q10Ref := math.Exp(-intRef)

			for ir, m := range subStepsList {
				group := refM / m
				dtM := 1.0 / float64(m)
				sqrtDtM := math.Sqrt(dtM)
				invSqrtGroup := 1.0 / math.Sqrt(float64(group))
				stepsM := int(math.Round(maxT * float64(m)))
				stepsPerQ := m / 4
				hwGroup := hwMaxM / m

				lam := cirP.Lambda0
				intLam := 0.0
				sPrevQ := 1.0
				rbSum := 0.0
				tauBucketLoss := 0.0
				tauInterpLoss := 0.0
				defaulted := false

				for sm := 0; sm < stepsM; sm++ {
					zSum := 0.0
					base := sm * group
					for g := 0; g < group; g++ {
						zSum += zLRef[base+g]
					}
					zM := zSum * invSqrtGroup
					nextLam := cirP.StepFullTruncationEuler(lam, dtM, sqrtDtM, zM)
					dInt := 0.5 * (lam + nextLam) * dtM
					intBefore := intLam
					intLam += dInt

					if !defaulted && intLam >= uExp {
						defaulted = true
						frac := 0.5
						if dInt > 0 {
							frac = (uExp - intBefore) / dInt
							if frac < 0 {
								frac = 0
							} else if frac > 1 {
								frac = 1
							}
						}
						tTau := (float64(sm) + frac) * dtM
						fineEnd := (sm + 1) * hwGroup
						fineStart := sm * hwGroup
						if fineEnd >= len(xAll) {
							fineEnd = len(xAll) - 1
						}
						if fineStart >= len(xAll) {
							fineStart = len(xAll) - 1
						}
						xTau := (1.0-frac)*xAll[fineStart] + frac*xAll[fineEnd]
						iTau := (1.0-frac)*iAll[fineStart] + frac*iAll[fineEnd]
						tauInterpLoss = lgd * hwP.DiscountFactor(tTau, iTau) * math.Max(swap.Value(hwP, tTau, xTau), 0)

						qIdx := (sm / stepsPerQ) + 1
						if qIdx <= nQ {
							tauBucketLoss = lgd * dQ[qIdx] * vPosQ[qIdx]
						}
					}

					lam = nextLam
					if (sm+1)%stepsPerQ == 0 {
						qIdx := (sm + 1) / stepsPerQ
						sCurrQ := math.Exp(-intLam)
						rbSum += lgd * dQ[qIdx] * vPosQ[qIdx] * (sPrevQ - sCurrQ)
						sPrevQ = sCurrQ
					}
				}

				q10M := math.Exp(-intLam)
				diffRB := rbSum - rbRef
				diffInt := intLam - intRef
				rs.sumTauB[ir] += tauBucketLoss
				rs.sumSqTauB[ir] += tauBucketLoss * tauBucketLoss
				rs.sumTauI[ir] += tauInterpLoss
				rs.sumSqTauI[ir] += tauInterpLoss * tauInterpLoss
				rs.sumRB[ir] += rbSum
				rs.sumSqRB[ir] += rbSum * rbSum
				rs.sumDiffRB[ir] += diffRB
				rs.sumSqDiffRB[ir] += diffRB * diffRB
				rs.sumDiffQ10[ir] += q10M - q10Ref
				rs.sumSqIntDiff[ir] += diffInt * diffInt
				rs.sumQ10[ir] += q10M
			}
		}
		return rs
	})

	for ir, m := range subStepsList {
		sTB, sqTB, sTI, sqTI, sRB, sqRB, sDRB, sqDRB, sDQ, sSqI, sQ10 := 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
		for _, r := range res {
			sTB += r.sumTauB[ir]
			sqTB += r.sumSqTauB[ir]
			sTI += r.sumTauI[ir]
			sqTI += r.sumSqTauI[ir]
			sRB += r.sumRB[ir]
			sqRB += r.sumSqRB[ir]
			sDRB += r.sumDiffRB[ir]
			sqDRB += r.sumSqDiffRB[ir]
			sDQ += r.sumDiffQ10[ir]
			sSqI += r.sumSqIntDiff[ir]
			sQ10 += r.sumQ10[ir]
		}
		mTB, seTB := stats.SummaryMoments(sTB, sqTB, nPaths)
		mTI, seTI := stats.SummaryMoments(sTI, sqTI, nPaths)
		mRB, seRB := stats.SummaryMoments(sRB, sqRB, nPaths)
		mDiffRB, seDiffRB := stats.SummaryMoments(sDRB, sqDRB, nPaths)
		mQ10 := sQ10 / float64(nPaths)
		zTB := 0.0
		if seTB > 0 {
			zTB = math.Abs(mTB-cvaAnal) / seTB
		}
		zRB := 0.0
		if seRB > 0 {
			zRB = math.Abs(mRB-cvaAnal) / seRB
		}
		rows[ir] = DiscretizationRow{
			SubStepsPerYear: m,
			CVATauBucket:    mTB,
			SETauBucket:     seTB,
			ZScoreTauBucket: zTB,
			CVATauInterp:    mTI,
			SETauInterp:     seTI,
			CVARaoBlackwell: mRB,
			SERaoBlackwell:  seRB,
			ZScoreRB:        zRB,
			AbsErrRBVsAnal:  math.Abs(mRB - cvaAnal),
			PairedBiasVs320: math.Abs(mDiffRB),
			PairedSEVs320:   seDiffRB,
			Q10YMC:          mQ10,
			Q10YBias:        mQ10 - q10Anal,
			Q10PairedBias:   math.Abs(sDQ / float64(nPaths)),
			IntRMSEVs320:    math.Sqrt(sSqI / float64(nPaths)),
		}
	}
	return rows
}

// WWRRow merangkum CVA di bawah korelasi dW^r dan dW^lambda sebesar rho dengan Common Random Numbers (CRN).
type WWRRow struct {
	Rho               float64
	CVARaoBlackwell   float64
	SERaoBlackwell    float64
	CVATauBucket      float64
	SETauBucket       float64
	RatioVsRhoZeroRB  float64
	RatioVsRhoZeroTau float64
}

// ComputeWWRCurve mengevaluasi CVA(rho) untuk daftar rho memakai CRN (lintasan suku bunga & inovasi ortogonal sama).
func ComputeWWRCurve(hwP hw.Params, cirP cir.Params, swap hw.SwapSpec, tGrid, rhos []float64, subStepsPerYear, nPaths, nChunks, nWorkers int, seed uint64) []WWRRow {
	maxT := tGrid[len(tGrid)-1]
	nSteps := int(math.Round(maxT * float64(subStepsPerYear)))
	dt := 1.0 / float64(subStepsPerYear)
	sqrtDt := math.Sqrt(dt)
	stepsPerQ := subStepsPerYear / 4
	sp := hwP.NewStepParams(dt)
	nQ := len(tGrid) - 1
	nRho := len(rhos)
	lgd := cirP.LGD()

	orth := make([]float64, nRho)
	for i, rho := range rhos {
		orth[i] = math.Sqrt(math.Max(0, 1.0-rho*rho))
	}

	type chunkWWR struct {
		sumRB, sumSqRB   []float64
		sumTau, sumSqTau []float64
		count            int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkWWR {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		cw := chunkWWR{
			sumRB:    make([]float64, nRho),
			sumSqRB:  make([]float64, nRho),
			sumTau:   make([]float64, nRho),
			sumSqTau: make([]float64, nRho),
			count:    cnt,
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x801, uint64(c)))
		lam := make([]float64, nRho)
		intLam := make([]float64, nRho)
		sPrev := make([]float64, nRho)
		rbLoss := make([]float64, nRho)
		tauLoss := make([]float64, nRho)
		defaulted := make([]bool, nRho)

		for k := 0; k < cnt; k++ {
			uExp := rg.ExpFloat64()
			xt, it := 0.0, 0.0
			for ir := 0; ir < nRho; ir++ {
				lam[ir] = cirP.Lambda0
				intLam[ir] = 0
				sPrev[ir] = 1.0
				rbLoss[ir] = 0
				tauLoss[ir] = 0
				defaulted[ir] = false
			}

			for s := 1; s <= nSteps; s++ {
				z1 := rg.NormFloat64()
				z2 := rg.NormFloat64()
				zPerp := rg.NormFloat64()
				xt, it = sp.Advance(xt, it, z1, z2)

				isQuarter := (s % stepsPerQ) == 0
				var qPay float64
				if isQuarter {
					qIdx := s / stepsPerQ
					if qIdx <= nQ {
						qPay = lgd * hwP.DiscountFactor(tGrid[qIdx], it) * math.Max(swap.Value(hwP, tGrid[qIdx], xt), 0)
					}
				}

				for ir, rho := range rhos {
					zL := rho*z1 + orth[ir]*zPerp
					nextLam := cirP.StepFullTruncationEuler(lam[ir], dt, sqrtDt, zL)
					intLam[ir] += 0.5 * (lam[ir] + nextLam) * dt
					lam[ir] = nextLam

					if isQuarter {
						sCurr := math.Exp(-intLam[ir])
						rbLoss[ir] += qPay * (sPrev[ir] - sCurr)
						sPrev[ir] = sCurr
						if !defaulted[ir] && intLam[ir] >= uExp {
							defaulted[ir] = true
							tauLoss[ir] = qPay
						}
					}
				}
			}

			for ir := 0; ir < nRho; ir++ {
				cw.sumRB[ir] += rbLoss[ir]
				cw.sumSqRB[ir] += rbLoss[ir] * rbLoss[ir]
				cw.sumTau[ir] += tauLoss[ir]
				cw.sumSqTau[ir] += tauLoss[ir] * tauLoss[ir]
			}
		}
		return cw
	})

	rows := make([]WWRRow, nRho)
	cvaZeroRB, cvaZeroTau := 0.0, 0.0
	for ir, rho := range rhos {
		sRB, sqRB, sT, sqT := 0.0, 0.0, 0.0, 0.0
		for _, r := range res {
			sRB += r.sumRB[ir]
			sqRB += r.sumSqRB[ir]
			sT += r.sumTau[ir]
			sqT += r.sumSqTau[ir]
		}
		mRB, seRB := stats.SummaryMoments(sRB, sqRB, nPaths)
		mT, seT := stats.SummaryMoments(sT, sqT, nPaths)
		rows[ir] = WWRRow{
			Rho:             rho,
			CVARaoBlackwell: mRB,
			SERaoBlackwell:  seRB,
			CVATauBucket:    mT,
			SETauBucket:     seT,
		}
		if math.Abs(rho) < 1e-12 {
			cvaZeroRB = mRB
			cvaZeroTau = mT
		}
	}
	for ir := range rows {
		if cvaZeroRB > 0 {
			rows[ir].RatioVsRhoZeroRB = rows[ir].CVARaoBlackwell / cvaZeroRB
		}
		if cvaZeroTau > 0 {
			rows[ir].RatioVsRhoZeroTau = rows[ir].CVATauBucket / cvaZeroTau
		}
	}
	return rows
}

// NettingProfile menyimpan profil EE & CVA untuk swap 1, swap 2, jumlah mandiri, dan portofolio ter-netting.
type NettingResult struct {
	Swap1EE           []float64
	Swap2EE           []float64
	StandaloneSumEE   []float64
	NettedEE          []float64
	Swap1CVAAnalytic  float64
	Swap2CVAAnalytic  float64
	NettedCVAAnalytic float64
	MaxPointwiseViol  float64
	OppositeNetMaxEE  float64
	OppositeNetCVA    float64
	ByRho             []NettingRhoRow
}

// NettingRhoRow menyimpan CVA portofolio netting pada rho tertentu.
type NettingRhoRow struct {
	Rho            float64
	Swap1CVA       float64
	Swap2CVA       float64
	StandaloneSum  float64
	NettedCVA      float64
	NettedCVASE    float64
	BenefitAbs     float64
	BenefitPercent float64
}

// ComputeNetting mengevaluasi benefit netting portofolio swap 5y payer + 10y receiver (serta identitas 10y payer + 10y receiver = 0).
func ComputeNetting(hwP hw.Params, cirP cir.Params, s1, s2, s10Rec hw.SwapSpec, tGrid, rhos []float64, nGH, subStepsPerYear, nPaths, nChunks, nWorkers int, seed uint64) NettingResult {
	nT := len(tGrid)
	nodes := special.GaussHermiteStandardNormal(nGH)
	ee1 := make([]float64, nT)
	ee2 := make([]float64, nT)
	eeSum := make([]float64, nT)
	eeNet := make([]float64, nT)
	eeOpp := make([]float64, nT)
	d1 := make([]float64, nT)
	d2 := make([]float64, nT)
	dNet := make([]float64, nT)
	dOpp := make([]float64, nT)
	maxViol := 0.0

	s10Pay := hwP.NewAnnualSwap(10, s10Rec.StrikeK, true)

	for j, t := range tGrid {
		if t <= 1e-12 || t >= math.Max(s1.Maturity, s2.Maturity)-1e-12 {
			continue
		}
		sp := hwP.NewStepParams(t)
		for _, nx := range nodes {
			xt := sp.L11 * nx.X
			condD := hwP.CondDiscountGivenX(t, xt)
			v1 := s1.Value(hwP, t, xt)
			v2 := s2.Value(hwP, t, xt)
			p1 := math.Max(v1, 0)
			p2 := math.Max(v2, 0)
			pNet := math.Max(v1+v2, 0)
			if pNet-(p1+p2) > maxViol {
				maxViol = pNet - (p1 + p2)
			}
			vOpp := s10Pay.Value(hwP, t, xt) + s10Rec.Value(hwP, t, xt)
			pOpp := math.Max(vOpp, 0)

			ee1[j] += nx.W * p1
			ee2[j] += nx.W * p2
			eeNet[j] += nx.W * pNet
			eeOpp[j] += nx.W * pOpp
			d1[j] += nx.W * condD * p1
			d2[j] += nx.W * condD * p2
			dNet[j] += nx.W * condD * pNet
			dOpp[j] += nx.W * condD * pOpp
		}
		eeSum[j] = ee1[j] + ee2[j]
	}

	cva1, _ := AnalyticCVAFromEEDisc(cirP, tGrid, d1)
	cva2, _ := AnalyticCVAFromEEDisc(cirP, tGrid, d2)
	cvaNet, _ := AnalyticCVAFromEEDisc(cirP, tGrid, dNet)
	cvaOpp, _ := AnalyticCVAFromEEDisc(cirP, tGrid, dOpp)
	maxOppEE := 0.0
	for _, v := range eeOpp {
		if v > maxOppEE {
			maxOppEE = v
		}
	}

	// Hitung CVA netting across rhos memakai MC CRN.
	maxT := tGrid[nT-1]
	nSteps := int(math.Round(maxT * float64(subStepsPerYear)))
	dt := 1.0 / float64(subStepsPerYear)
	sqrtDt := math.Sqrt(dt)
	stepsPerQ := subStepsPerYear / 4
	spStep := hwP.NewStepParams(dt)
	nRho := len(rhos)
	orth := make([]float64, nRho)
	for i, rho := range rhos {
		orth[i] = math.Sqrt(math.Max(0, 1.0-rho*rho))
	}
	lgd := cirP.LGD()

	type chunkNet struct {
		s1, s2, sNet, sqNet []float64
		count               int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkNet {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		cn := chunkNet{
			s1:    make([]float64, nRho),
			s2:    make([]float64, nRho),
			sNet:  make([]float64, nRho),
			sqNet: make([]float64, nRho),
			count: cnt,
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x901, uint64(c)))
		lam := make([]float64, nRho)
		intLam := make([]float64, nRho)
		sPrev := make([]float64, nRho)
		l1 := make([]float64, nRho)
		l2 := make([]float64, nRho)
		lNet := make([]float64, nRho)

		for k := 0; k < cnt; k++ {
			xt, it := 0.0, 0.0
			for ir := 0; ir < nRho; ir++ {
				lam[ir] = cirP.Lambda0
				intLam[ir] = 0
				sPrev[ir] = 1.0
				l1[ir], l2[ir], lNet[ir] = 0, 0, 0
			}
			for s := 1; s <= nSteps; s++ {
				z1 := rg.NormFloat64()
				z2 := rg.NormFloat64()
				zP := rg.NormFloat64()
				xt, it = spStep.Advance(xt, it, z1, z2)
				isQ := (s % stepsPerQ) == 0
				var p1Q, p2Q, pNetQ float64
				if isQ {
					qIdx := s / stepsPerQ
					tQ := tGrid[qIdx]
					df := hwP.DiscountFactor(tQ, it)
					v1 := s1.Value(hwP, tQ, xt)
					v2 := s2.Value(hwP, tQ, xt)
					p1Q = lgd * df * math.Max(v1, 0)
					p2Q = lgd * df * math.Max(v2, 0)
					pNetQ = lgd * df * math.Max(v1+v2, 0)
				}
				for ir, rho := range rhos {
					zL := rho*z1 + orth[ir]*zP
					nextLam := cirP.StepFullTruncationEuler(lam[ir], dt, sqrtDt, zL)
					intLam[ir] += 0.5 * (lam[ir] + nextLam) * dt
					lam[ir] = nextLam
					if isQ {
						sCurr := math.Exp(-intLam[ir])
						dpd := sPrev[ir] - sCurr
						sPrev[ir] = sCurr
						l1[ir] += p1Q * dpd
						l2[ir] += p2Q * dpd
						lNet[ir] += pNetQ * dpd
					}
				}
			}
			for ir := 0; ir < nRho; ir++ {
				cn.s1[ir] += l1[ir]
				cn.s2[ir] += l2[ir]
				cn.sNet[ir] += lNet[ir]
				cn.sqNet[ir] += lNet[ir] * lNet[ir]
			}
		}
		return cn
	})

	byRho := make([]NettingRhoRow, nRho)
	for ir, rho := range rhos {
		acc1, acc2, accN, accSqN := 0.0, 0.0, 0.0, 0.0
		for _, r := range res {
			acc1 += r.s1[ir]
			acc2 += r.s2[ir]
			accN += r.sNet[ir]
			accSqN += r.sqNet[ir]
		}
		m1 := acc1 / float64(nPaths)
		m2 := acc2 / float64(nPaths)
		mNet, seNet := stats.SummaryMoments(accN, accSqN, nPaths)
		sumStd := m1 + m2
		ben := sumStd - mNet
		pct := 0.0
		if sumStd > 0 {
			pct = 100.0 * ben / sumStd
		}
		byRho[ir] = NettingRhoRow{
			Rho:            rho,
			Swap1CVA:       m1,
			Swap2CVA:       m2,
			StandaloneSum:  sumStd,
			NettedCVA:      mNet,
			NettedCVASE:    seNet,
			BenefitAbs:     ben,
			BenefitPercent: pct,
		}
	}

	return NettingResult{
		Swap1EE:           ee1,
		Swap2EE:           ee2,
		StandaloneSumEE:   eeSum,
		NettedEE:          eeNet,
		Swap1CVAAnalytic:  cva1,
		Swap2CVAAnalytic:  cva2,
		NettedCVAAnalytic: cvaNet,
		MaxPointwiseViol:  maxViol,
		OppositeNetMaxEE:  maxOppEE,
		OppositeNetCVA:    cvaOpp,
		ByRho:             byRho,
	}
}

// CollateralCVARow menyimpan CVA terkolateralisasi untuk kombinasi (delta, H, rho)
// baik untuk definisi jadwal bersih konsisten maupun definisi formula kupon-melompat mentah.
type CollateralCVARow struct {
	Delta        float64
	Threshold    float64
	Rho          float64
	CVA          float64
	SE           float64
	EPE          float64
	CVARawCoupon float64
	SERawCoupon  float64
	EPERawCoupon float64
}

type eventNode struct {
	t          float64
	dt         float64
	sqrtDt     float64
	isQuarter  bool
	quarterIdx int
	isLagMark  bool
	lagForQIdx int
	sp         hw.StepParams
}

func buildCollateralSchedule(tGrid []float64, delta float64, subStepsPerYear int) []eventNode {
	type timeTag struct {
		t       float64
		isQ     bool
		qIdx    int
		isLag   bool
		lagForQ int
	}
	maxT := tGrid[len(tGrid)-1]
	nSub := int(math.Round(maxT * float64(subStepsPerYear)))
	dtSub := 1.0 / float64(subStepsPerYear)
	stepsPerQ := subStepsPerYear / 4

	tags := make([]timeTag, 0, nSub+len(tGrid))
	for s := 1; s <= nSub; s++ {
		t := float64(s) * dtSub
		isQ := (s % stepsPerQ) == 0
		qIdx := 0
		if isQ {
			qIdx = s / stepsPerQ
		}
		tags = append(tags, timeTag{t: t, isQ: isQ, qIdx: qIdx})
	}
	for j := 1; j < len(tGrid); j++ {
		tLag := tGrid[j] - delta
		if tLag <= 1e-12 {
			continue
		}
		// Gabungkan ke simpul yang sudah ada bila sama dalam batas 1e-12.
		merged := false
		for i := range tags {
			if math.Abs(tags[i].t-tLag) < 1e-12 {
				tags[i].isLag = true
				tags[i].lagForQ = j
				merged = true
				break
			}
		}
		if !merged {
			tags = append(tags, timeTag{t: tLag, isLag: true, lagForQ: j})
		}
	}
	sort.Slice(tags, func(i, j int) bool {
		return tags[i].t < tags[j].t
	})

	hwP := hw.DefaultParams()
	nodes := make([]eventNode, len(tags))
	prevT := 0.0
	for i, tg := range tags {
		dt := tg.t - prevT
		nodes[i] = eventNode{
			t:          tg.t,
			dt:         dt,
			sqrtDt:     math.Sqrt(dt),
			isQuarter:  tg.isQ,
			quarterIdx: tg.qIdx,
			isLagMark:  tg.isLag,
			lagForQIdx: tg.lagForQ,
			sp:         hwP.NewStepParams(dt),
		}
		prevT = tg.t
	}
	return nodes
}

// ComputeCollateralCVA mengevaluasi CVA terkolateralisasi untuk pasangan (delta, threshold) pada rho in [0, 0.8].
func ComputeCollateralCVA(hwP hw.Params, cirP cir.Params, swap hw.SwapSpec, tGrid []float64, delta float64, thresholds, rhos []float64, subStepsPerYear, nPaths, nChunks, nWorkers int, seed uint64) []CollateralCVARow {
	sched := buildCollateralSchedule(tGrid, delta, subStepsPerYear)
	nH := len(thresholds)
	nRho := len(rhos)
	nQ := len(tGrid) - 1
	lgd := cirP.LGD()

	orth := make([]float64, nRho)
	for i, rho := range rhos {
		orth[i] = math.Sqrt(math.Max(0, 1.0-rho*rho))
	}

	type chunkColCVA struct {
		sumCVA, sumSqCVA       [][]float64 // [nRho][nH] jadwal konsisten
		sumRawCVA, sumSqRawCVA [][]float64 // [nRho][nH] formula kupon mentah
		sumEE, sumRawEE        [][]float64 // [nH][nQ+1]
		count                  int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkColCVA {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		cc := chunkColCVA{
			sumCVA:      make([][]float64, nRho),
			sumSqCVA:    make([][]float64, nRho),
			sumRawCVA:   make([][]float64, nRho),
			sumSqRawCVA: make([][]float64, nRho),
			sumEE:       make([][]float64, nH),
			sumRawEE:    make([][]float64, nH),
			count:       cnt,
		}
		for ir := 0; ir < nRho; ir++ {
			cc.sumCVA[ir] = make([]float64, nH)
			cc.sumSqCVA[ir] = make([]float64, nH)
			cc.sumRawCVA[ir] = make([]float64, nH)
			cc.sumSqRawCVA[ir] = make([]float64, nH)
		}
		for ih := 0; ih < nH; ih++ {
			cc.sumEE[ih] = make([]float64, nQ+1)
			cc.sumRawEE[ih] = make([]float64, nQ+1)
		}

		// Garam chunk diikat pada salt tetap 0xA01 (plus kuantisasi delta agar tidak bentrok) untuk stabilitas CRN.
		deltaSalt := uint64(math.Round(delta * 10000.0))
		rg := rng.New(rng.DeriveChunkSeed(seed, 0xA01^deltaSalt, uint64(c)))
		lam := make([]float64, nRho)
		intLam := make([]float64, nRho)
		sPrev := make([]float64, nRho)
		loss := make([][]float64, nRho)
		lossRaw := make([][]float64, nRho)
		for ir := 0; ir < nRho; ir++ {
			loss[ir] = make([]float64, nH)
			lossRaw[ir] = make([]float64, nH)
		}
		vLagCleanByQ := make([]float64, nQ+1)
		vLagRawByQ := make([]float64, nQ+1)

		for k := 0; k < cnt; k++ {
			xt, it := 0.0, 0.0
			for q := 0; q <= nQ; q++ {
				vLagCleanByQ[q] = 0
				vLagRawByQ[q] = 0
			}
			for ir := 0; ir < nRho; ir++ {
				lam[ir] = cirP.Lambda0
				intLam[ir] = 0
				sPrev[ir] = 1.0
				for ih := 0; ih < nH; ih++ {
					loss[ir][ih] = 0
					lossRaw[ir][ih] = 0
				}
			}

			for _, ev := range sched {
				z1 := rg.NormFloat64()
				z2 := rg.NormFloat64()
				zP := rg.NormFloat64()
				xt, it = ev.sp.Advance(xt, it, z1, z2)

				if ev.isLagMark {
					tQ := tGrid[ev.lagForQIdx]
					vLagCleanByQ[ev.lagForQIdx] = swap.ValueWithScheduleAfter(hwP, ev.t, xt, tQ)
					vLagRawByQ[ev.lagForQIdx] = swap.Value(hwP, ev.t, xt)
				}

				for ir, rho := range rhos {
					zL := rho*z1 + orth[ir]*zP
					nextLam := cirP.StepFullTruncationEuler(lam[ir], ev.dt, ev.sqrtDt, zL)
					intLam[ir] += 0.5 * (lam[ir] + nextLam) * ev.dt
					lam[ir] = nextLam
				}

				if ev.isQuarter {
					qIdx := ev.quarterIdx
					vt := swap.Value(hwP, ev.t, xt)
					df := hwP.DiscountFactor(ev.t, it)
					vLagC := vLagCleanByQ[qIdx]
					vLagR := vLagRawByQ[qIdx]

					for ir := 0; ir < nRho; ir++ {
						sCurr := math.Exp(-intLam[ir])
						dpd := sPrev[ir] - sCurr
						sPrev[ir] = sCurr
						for ih, h := range thresholds {
							cClean, cRaw := 0.0, 0.0
							if !math.IsInf(h, 1) {
								cClean = math.Max(vLagC-h, 0)
								cRaw = math.Max(vLagR-h, 0)
							}
							expC := math.Max(vt-cClean, 0)
							expR := math.Max(vt-cRaw, 0)
							if ir == 0 {
								cc.sumEE[ih][qIdx] += expC
								cc.sumRawEE[ih][qIdx] += expR
							}
							loss[ir][ih] += lgd * df * expC * dpd
							lossRaw[ir][ih] += lgd * df * expR * dpd
						}
					}
				}
			}
			for ir := 0; ir < nRho; ir++ {
				for ih := 0; ih < nH; ih++ {
					lc := loss[ir][ih]
					cc.sumCVA[ir][ih] += lc
					cc.sumSqCVA[ir][ih] += lc * lc
					lr := lossRaw[ir][ih]
					cc.sumRawCVA[ir][ih] += lr
					cc.sumSqRawCVA[ir][ih] += lr * lr
				}
			}
		}
		return cc
	})

	epeByH := make([]float64, nH)
	epeRawByH := make([]float64, nH)
	totalT := tGrid[nQ] - tGrid[0]
	for ih := 0; ih < nH; ih++ {
		eeCurve := make([]float64, nQ+1)
		eeRawCurve := make([]float64, nQ+1)
		for q := 0; q <= nQ; q++ {
			s, sr := 0.0, 0.0
			for _, r := range res {
				s += r.sumEE[ih][q]
				sr += r.sumRawEE[ih][q]
			}
			eeCurve[q] = s / float64(nPaths)
			eeRawCurve[q] = sr / float64(nPaths)
		}
		area, areaR := 0.0, 0.0
		for q := 1; q <= nQ; q++ {
			dt := tGrid[q] - tGrid[q-1]
			area += 0.5 * (eeCurve[q-1] + eeCurve[q]) * dt
			areaR += 0.5 * (eeRawCurve[q-1] + eeRawCurve[q]) * dt
		}
		if totalT > 0 {
			epeByH[ih] = area / totalT
			epeRawByH[ih] = areaR / totalT
		}
	}

	out := make([]CollateralCVARow, 0, nRho*nH)
	for ir, rho := range rhos {
		for ih, h := range thresholds {
			s, sq, sr, sqr := 0.0, 0.0, 0.0, 0.0
			for _, r := range res {
				s += r.sumCVA[ir][ih]
				sq += r.sumSqCVA[ir][ih]
				sr += r.sumRawCVA[ir][ih]
				sqr += r.sumSqRawCVA[ir][ih]
			}
			m, se := stats.SummaryMoments(s, sq, nPaths)
			mr, ser := stats.SummaryMoments(sr, sqr, nPaths)
			out = append(out, CollateralCVARow{
				Delta:        delta,
				Threshold:    h,
				Rho:          rho,
				CVA:          m,
				SE:           se,
				EPE:          epeByH[ih],
				CVARawCoupon: mr,
				SERawCoupon:  ser,
				EPERawCoupon: epeRawByH[ih],
			})
		}
	}
	return out
}
