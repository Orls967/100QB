package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/cir"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/cva"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/hw"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/jsonout"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/parallel"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/stats"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/svg"
)

const (
	nChunks        = 32
	notionalScale  = 1_000_000.0
	nPathsExposure = 200_000
	nPathsCore     = 100_000
	subStepsStd    = 40
	nQuadOrder     = 32
)

func findRepoRoot() string {
	wd, _ := os.Getwd()
	for cur := wd; cur != "/" && cur != "."; cur = filepath.Dir(cur) {
		if _, err := os.Stat(filepath.Join(cur, "outputs")); err == nil {
			return cur
		}
	}
	return wd
}

func logLogSlope(xs, ys []float64) float64 {
	lx := make([]float64, len(xs))
	ly := make([]float64, len(ys))
	for i := range xs {
		lx[i] = math.Log(xs[i])
		ly[i] = math.Log(ys[i])
	}
	return stats.LinearSlope(lx, ly)
}

func formatThresholdLabel(h float64) string {
	if math.IsInf(h, 1) {
		return "Inf"
	}
	return fmt.Sprintf("%.4g", h)
}

func main() {
	tStart := time.Now()
	root := findRepoRoot()
	outPath := filepath.Join(root, "outputs", "day003e", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	nWorkers := parallel.DefaultWorkers()

	hwP := hw.DefaultParams()
	cirP := cir.DefaultParams()
	tGrid := cva.QuarterlyGrid(10)

	// =========================================================================
	// BAGIAN 1: Hull-White 1-Faktor, Swap Berjalan, Profil Eksposur & Swaption
	// =========================================================================
	swap10Pay := hwP.NewAnnualSwap(10, 0, true)
	swap10Rec := hwP.NewAnnualSwap(10, swap10Pay.StrikeK, false)
	swap5Pay := hwP.NewAnnualSwap(5, 0, true)

	martTimes := [][2]float64{{2.0, 10.0}, {5.0, 10.0}, {8.0, 10.0}}
	martObjs := make([]jsonout.Obj, 0, len(martTimes))
	martPass := true
	maxMartGHErr := 0.0
	maxMartMCZ := 0.0
	for _, pair := range martTimes {
		lh1, lh2, p0T, err1, err2 := hwP.MartingaleCheckGH(pair[0], pair[1], nQuadOrder)
		mcM, mcSE, _, zMC := hwP.MartingaleCheckMC(pair[0], pair[1], nPathsExposure, nChunks, nWorkers, rng.MainSeed)
		if err1 > maxMartGHErr {
			maxMartGHErr = err1
		}
		if err2 > maxMartGHErr {
			maxMartGHErr = err2
		}
		if zMC > maxMartMCZ {
			maxMartMCZ = zMC
		}
		if err1 >= 1e-12 || err2 >= 1e-12 || math.Abs(mcM-p0T) >= 5e-4 {
			martPass = false
		}
		martObjs = append(martObjs, jsonout.Obj{
			{Key: "t", Val: pair[0]},
			{Key: "T", Val: pair[1]},
			{Key: "p0_T", Val: p0T},
			{Key: "gh_1d", Val: lh1},
			{Key: "gh_2d", Val: lh2},
			{Key: "gh_1d_abs_err", Val: err1},
			{Key: "gh_2d_abs_err", Val: err2},
			{Key: "mc_mean", Val: mcM},
			{Key: "mc_se", Val: mcSE},
			{Key: "mc_abs_err", Val: math.Abs(mcM - p0T)},
			{Key: "mc_z_score", Val: zMC},
		})
	}

	// Swaption 2y x 8y Jamshidian vs kuadratur bebas-kink & MC
	payDates2x8 := []float64{3, 4, 5, 6, 7, 8, 9, 10}
	kFwd := hwP.ForwardSwapRate(2.0, payDates2x8)
	swaptionStrikes := []float64{kFwd - 0.01, kFwd, kFwd + 0.01}
	swaptionLabels := []string{"ITM", "ATM", "OTM"}
	swaptionRes := hwP.VerifySwaption(2.0, payDates2x8, swaptionStrikes, swaptionLabels, nQuadOrder, nPathsExposure, nChunks, nWorkers, rng.MainSeed)

	swaptionPass := true
	maxJamErr := 0.0
	maxSwapMCZ := 0.0
	swaptionObjs := make([]jsonout.Obj, 0, len(swaptionRes))
	for _, r := range swaptionRes {
		if r.AbsErrJamVsGH > maxJamErr {
			maxJamErr = r.AbsErrJamVsGH
		}
		if r.ZScoreMC > maxSwapMCZ {
			maxSwapMCZ = r.ZScoreMC
		}
		if r.AbsErrJamVsGH >= 1e-8 || r.ZScoreMC > 3.0 {
			swaptionPass = false
		}
		swaptionObjs = append(swaptionObjs, jsonout.Obj{
			{Key: "label", Val: r.Label},
			{Key: "strike", Val: r.Strike},
			{Key: "x_star", Val: r.XStar},
			{Key: "sum_ci_xi", Val: r.CheckSumCiXi},
			{Key: "jamshidian", Val: r.Jamshidian},
			{Key: "quad_1d", Val: r.QuadGH1D},
			{Key: "quad_2d", Val: r.QuadGH2D},
			{Key: "abs_err_jam_vs_quad2d", Val: r.AbsErrJamVsGH},
			{Key: "mc_mean", Val: r.MCMean},
			{Key: "mc_se", Val: r.MCSE},
			{Key: "mc_z_score", Val: r.ZScoreMC},
		})
	}

	// Profil eksposur swap 10y ATM (200.000 lintasan) + uji determinisme 1 vs nWorkers
	expPts, epeGH, epeMC := hwP.ComputeSwapExposureProfile(swap10Pay, tGrid, nQuadOrder, nPathsExposure, nChunks, nWorkers, rng.MainSeed)
	_, _, epeMC1W := hwP.ComputeSwapExposureProfile(swap10Pay, tGrid, nQuadOrder, 16000, nChunks, 1, rng.MainSeed)
	_, _, epeMCNW := hwP.ComputeSwapExposureProfile(swap10Pay, tGrid, nQuadOrder, 16000, nChunks, nWorkers, rng.MainSeed)
	workerBitIdentical := math.Float64bits(epeMC1W) == math.Float64bits(epeMCNW)

	expPass := true
	maxEEZ := 0.0
	peakEE, peakEET := 0.0, 0.0
	peakPFE99, peakPFE99T := 0.0, 0.0
	tSeries := make([]float64, len(expPts))
	eeGHSeries := make([]float64, len(expPts))
	eeMCSeries := make([]float64, len(expPts))
	eeDiscSeries := make([]float64, len(expPts))
	pfe95Series := make([]float64, len(expPts))
	pfe99Series := make([]float64, len(expPts))
	expObjs := make([]jsonout.Obj, 0, len(expPts))

	for j, ep := range expPts {
		tSeries[j] = ep.T
		eeGHSeries[j] = ep.EEGH
		eeMCSeries[j] = ep.EEMC
		eeDiscSeries[j] = ep.EEDiscGH
		pfe95Series[j] = ep.PFE95MC
		pfe99Series[j] = ep.PFE99MC
		if ep.EEGH > peakEE {
			peakEE, peakEET = ep.EEGH, ep.T
		}
		if ep.PFE99MC > peakPFE99 {
			peakPFE99, peakPFE99T = ep.PFE99MC, ep.T
		}
		z := 0.0
		if ep.EESEMC > 0 {
			z = math.Abs(ep.EEMC-ep.EEGH) / ep.EESEMC
			if z > maxEEZ {
				maxEEZ = z
			}
			if z > 3.5 {
				expPass = false
			}
		}
		expObjs = append(expObjs, jsonout.Obj{
			{Key: "t", Val: ep.T},
			{Key: "ee_analytic", Val: ep.EEGH},
			{Key: "ee_disc_analytic", Val: ep.EEDiscGH},
			{Key: "ee_disc_2d", Val: ep.EEDiscGH2D},
			{Key: "ee_mc", Val: ep.EEMC},
			{Key: "ee_mc_se", Val: ep.EESEMC},
			{Key: "ee_disc_mc", Val: ep.EEDiscMC},
			{Key: "ee_disc_mc_se", Val: ep.EEDiscSEMC},
			{Key: "ee_z_score", Val: z},
			{Key: "pfe_95_mc", Val: ep.PFE95MC},
			{Key: "pfe_99_mc", Val: ep.PFE99MC},
		})
	}

	// =========================================================================
	// BAGIAN 2: Intensitas Gagal CIR, RK4 Riccati, Exact Sampling & Par CDS
	// =========================================================================
	fellerLHS, fellerRHS, fellerOK := cirP.FellerCheck()
	cirHorizons := []float64{0.5, 1, 2, 3, 5, 7, 10}
	rk4Pass := true
	maxRK4Err := 0.0
	cirCurveObjs := make([]jsonout.Obj, 0, len(cirHorizons))
	for _, T := range cirHorizons {
		qAnal := cirP.SurvivalProb(T)
		qRK4 := cirP.SurvivalRiccatiRK4(T, 2000)
		diff := math.Abs(qAnal - qRK4)
		if diff > maxRK4Err {
			maxRK4Err = diff
		}
		if diff >= 1e-10 {
			rk4Pass = false
		}
		aCoeff, bCoeff := cirP.AffineCoeffs(T)
		cirCurveObjs = append(cirCurveObjs, jsonout.Obj{
			{Key: "T", Val: T},
			{Key: "A", Val: aCoeff},
			{Key: "B", Val: bCoeff},
			{Key: "q_analytic", Val: qAnal},
			{Key: "q_rk4", Val: qRK4},
			{Key: "abs_err_rk4", Val: diff},
		})
	}

	cirSubSteps := []int{20, 40, 80}
	survRows := cirP.VerifySurvivalMC(cirHorizons, cirSubSteps, nPathsCore, nChunks, nWorkers, rng.MainSeed)
	survPass := true
	survMCObjs := make([]jsonout.Obj, 0, len(survRows))
	qExact40 := make([]float64, 0, len(cirHorizons))
	qEuler40 := make([]float64, 0, len(cirHorizons))
	qAnalCurve := make([]float64, 0, len(cirHorizons))
	qRK4Curve := make([]float64, 0, len(cirHorizons))
	for _, T := range cirHorizons {
		qAnalCurve = append(qAnalCurve, cirP.SurvivalProb(T))
		qRK4Curve = append(qRK4Curve, cirP.SurvivalRiccatiRK4(T, 2000))
	}
	for _, r := range survRows {
		zEx := math.Abs(r.ExactCIRMC-r.AnalyticQ) / math.Max(r.ExactCIRSE, 1e-16)
		zEu := math.Abs(r.EulerFTMC-r.AnalyticQ) / math.Max(r.EulerFTSE, 1e-16)
		if zEx > 3.5 || zEu > 3.5 {
			survPass = false
		}
		if r.SubStepsPerYear == 40 {
			qExact40 = append(qExact40, r.ExactCIRMC)
			qEuler40 = append(qEuler40, r.EulerFTMC)
		}
		survMCObjs = append(survMCObjs, jsonout.Obj{
			{Key: "horizon", Val: r.Horizon},
			{Key: "sub_steps_per_year", Val: r.SubStepsPerYear},
			{Key: "analytic_q", Val: r.AnalyticQ},
			{Key: "exact_cir_mc", Val: r.ExactCIRMC},
			{Key: "exact_cir_se", Val: r.ExactCIRSE},
			{Key: "exact_cir_z", Val: zEx},
			{Key: "euler_ft_mc", Val: r.EulerFTMC},
			{Key: "euler_ft_se", Val: r.EulerFTSE},
			{Key: "euler_ft_z", Val: zEu},
		})
	}

	cdsMats := []float64{1, 3, 5, 10}
	cdsRes := cirP.VerifyParCDS(cdsMats, hwP.R0, nPathsCore, nChunks, nWorkers, rng.MainSeed)
	cdsPass := true
	cdsObjs := make([]jsonout.Obj, 0, len(cdsRes))
	for _, c := range cdsRes {
		if c.ZScoreProt > 3.0 {
			cdsPass = false
		}
		cdsObjs = append(cdsObjs, jsonout.Obj{
			{Key: "maturity", Val: c.Maturity},
			{Key: "analytic_pv_prot", Val: c.AnalPVProt},
			{Key: "analytic_rpv01", Val: c.AnalRPV01},
			{Key: "analytic_par_spread", Val: c.AnalSpread},
			{Key: "analytic_par_spread_bps", Val: c.AnalSpread * 10_000.0},
			{Key: "mc_pv_prot", Val: c.MCPVProt},
			{Key: "mc_pv_prot_se", Val: c.SEPVProt},
			{Key: "mc_rpv01", Val: c.MCRPV01},
			{Key: "mc_rpv01_se", Val: c.SERPV01},
			{Key: "mc_par_spread_bps", Val: c.MCSpread * 10_000.0},
			{Key: "pv_prot_z_score", Val: c.ZScoreProt},
		})
	}

	// =========================================================================
	// BAGIAN 3: CVA Independen, Sensitivitas & Konvergensi Diskretisasi
	// =========================================================================
	analCVA := cva.ComputeAnalyticCVAAndSensitivities(hwP, cirP, swap10Pay, tGrid, nQuadOrder)
	discRows := cva.VerifyIndependentCVAMC(hwP, cirP, swap10Pay, tGrid, analCVA.CVAAnalytic, cirSubSteps, nPathsCore, nChunks, nWorkers, rng.MainSeed)

	discPass := true
	discObjs := make([]jsonout.Obj, 0, len(discRows))
	mXS := make([]float64, len(discRows))
	cvaBiasYS := make([]float64, len(discRows))
	intRMSEYS := make([]float64, len(discRows))
	for i, r := range discRows {
		mXS[i] = float64(r.SubStepsPerYear)
		cvaBiasYS[i] = math.Max(r.PairedBiasVs320, 1e-12)
		intRMSEYS[i] = math.Max(r.IntRMSEVs320, 1e-12)
		if r.ZScoreTauBucket > 3.0 || r.ZScoreRB > 3.0 {
			discPass = false
		}
		if i > 0 && (r.PairedBiasVs320 >= discRows[i-1].PairedBiasVs320 || r.IntRMSEVs320 >= discRows[i-1].IntRMSEVs320) {
			discPass = false
		}
		discObjs = append(discObjs, jsonout.Obj{
			{Key: "sub_steps_per_year", Val: r.SubStepsPerYear},
			{Key: "cva_tau_bucket_mc", Val: r.CVATauBucket},
			{Key: "cva_tau_bucket_se", Val: r.SETauBucket},
			{Key: "cva_tau_bucket_z", Val: r.ZScoreTauBucket},
			{Key: "cva_tau_interp_mc", Val: r.CVATauInterp},
			{Key: "cva_tau_interp_se", Val: r.SETauInterp},
			{Key: "cva_rb_mc", Val: r.CVARaoBlackwell},
			{Key: "cva_rb_se", Val: r.SERaoBlackwell},
			{Key: "cva_rb_z", Val: r.ZScoreRB},
			{Key: "abs_err_rb_vs_analytic", Val: r.AbsErrRBVsAnal},
			{Key: "paired_cva_bias_vs_320", Val: r.PairedBiasVs320},
			{Key: "paired_cva_se_vs_320", Val: r.PairedSEVs320},
			{Key: "q10_mc", Val: r.Q10YMC},
			{Key: "q10_bias_vs_analytic", Val: r.Q10YBias},
			{Key: "q10_paired_bias_vs_320", Val: r.Q10PairedBias},
			{Key: "integral_rmse_vs_320", Val: r.IntRMSEVs320},
		})
	}
	slopeCVABias := logLogSlope(mXS, cvaBiasYS)
	slopeIntRMSE := logLogSlope(mXS, intRMSEYS)

	// =========================================================================
	// BAGIAN 4.1: Wrong-Way Risk CVA(rho) dengan CRN
	// =========================================================================
	rhos := []float64{-0.8, -0.4, 0, 0.4, 0.8}
	wwrRows := cva.ComputeWWRCurve(hwP, cirP, swap10Pay, tGrid, rhos, subStepsStd, nPathsCore, nChunks, nWorkers, rng.MainSeed)
	wwrPass := true
	wwrObjs := make([]jsonout.Obj, 0, len(wwrRows))
	wwrRhos := make([]float64, len(wwrRows))
	wwrRB := make([]float64, len(wwrRows))
	wwrRBLo := make([]float64, len(wwrRows))
	wwrRBHi := make([]float64, len(wwrRows))
	wwrTau := make([]float64, len(wwrRows))
	wwrTauLo := make([]float64, len(wwrRows))
	wwrTauHi := make([]float64, len(wwrRows))
	for i, r := range wwrRows {
		wwrRhos[i] = r.Rho
		wwrRB[i] = r.CVARaoBlackwell
		wwrRBLo[i] = r.CVARaoBlackwell - 1.96*r.SERaoBlackwell
		wwrRBHi[i] = r.CVARaoBlackwell + 1.96*r.SERaoBlackwell
		wwrTau[i] = r.CVATauBucket
		wwrTauLo[i] = r.CVATauBucket - 1.96*r.SETauBucket
		wwrTauHi[i] = r.CVATauBucket + 1.96*r.SETauBucket
		if i > 0 && r.CVARaoBlackwell <= wwrRows[i-1].CVARaoBlackwell {
			wwrPass = false
		}
		wwrObjs = append(wwrObjs, jsonout.Obj{
			{Key: "rho", Val: r.Rho},
			{Key: "cva_rb", Val: r.CVARaoBlackwell},
			{Key: "cva_rb_se", Val: r.SERaoBlackwell},
			{Key: "cva_rb_1m", Val: r.CVARaoBlackwell * notionalScale},
			{Key: "cva_tau_bucket", Val: r.CVATauBucket},
			{Key: "cva_tau_bucket_se", Val: r.SETauBucket},
			{Key: "ratio_vs_rho0_rb", Val: r.RatioVsRhoZeroRB},
			{Key: "ratio_vs_rho0_tau", Val: r.RatioVsRhoZeroTau},
		})
	}

	// =========================================================================
	// BAGIAN 4.2: Netting Portofolio (5y Payer + 10y Receiver & Identitas Nol)
	// =========================================================================
	netRhos := []float64{-0.8, 0, 0.8}
	netRes := cva.ComputeNetting(hwP, cirP, swap5Pay, swap10Rec, swap10Rec, tGrid, netRhos, nQuadOrder, subStepsStd, nPathsCore, nChunks, nWorkers, rng.MainSeed)
	netPass := netRes.MaxPointwiseViol <= 1e-14 && netRes.OppositeNetMaxEE <= 1e-14 && math.Abs(netRes.OppositeNetCVA) <= 1e-14
	netRhoObjs := make([]jsonout.Obj, 0, len(netRes.ByRho))
	for _, r := range netRes.ByRho {
		if r.NettedCVA > r.StandaloneSum {
			netPass = false
		}
		netRhoObjs = append(netRhoObjs, jsonout.Obj{
			{Key: "rho", Val: r.Rho},
			{Key: "swap1_5y_payer_cva", Val: r.Swap1CVA},
			{Key: "swap2_10y_receiver_cva", Val: r.Swap2CVA},
			{Key: "standalone_sum_cva", Val: r.StandaloneSum},
			{Key: "netted_cva", Val: r.NettedCVA},
			{Key: "netted_cva_se", Val: r.NettedCVASE},
			{Key: "benefit_abs", Val: r.BenefitAbs},
			{Key: "benefit_percent", Val: r.BenefitPercent},
		})
	}

	// =========================================================================
	// BAGIAN 4.3: Kolateral Variasi Margin (H in {0, 0.01, 0.05, Inf} & MPOR)
	// =========================================================================
	thresholds := []float64{0, 0.01, 0.05, math.Inf(1)}
	colRhos := []float64{0, 0.8}
	colEEByH, _, colEPEByH := hwP.ComputeCollateralExposureProfile(swap10Pay, tGrid, 0.04, thresholds, nPathsExposure, nChunks, nWorkers, rng.MainSeed)
	colThreshRows := cva.ComputeCollateralCVA(hwP, cirP, swap10Pay, tGrid, 0.04, thresholds, colRhos, subStepsStd, nPathsCore, nChunks, nWorkers, rng.MainSeed)

	mporDeltas := []float64{0.04, 0.01, 0.0025}
	mporRows := make([]cva.CollateralCVARow, 0, len(mporDeltas)*len(colRhos))
	for _, d := range mporDeltas {
		rows := cva.ComputeCollateralCVA(hwP, cirP, swap10Pay, tGrid, d, []float64{0}, colRhos, subStepsStd, nPathsCore, nChunks, nWorkers, rng.MainSeed)
		mporRows = append(mporRows, rows...)
	}

	c8aThreshMonotone := true
	colThreshObjs := make([]jsonout.Obj, 0, len(colThreshRows))
	for _, rho := range colRhos {
		prev := -1.0
		for _, r := range colThreshRows {
			if r.Rho != rho {
				continue
			}
			if prev >= 0 && r.CVA <= prev {
				c8aThreshMonotone = false
			}
			prev = r.CVA
		}
	}
	for _, r := range colThreshRows {
		colThreshObjs = append(colThreshObjs, jsonout.Obj{
			{Key: "delta", Val: r.Delta},
			{Key: "threshold", Val: formatThresholdLabel(r.Threshold)},
			{Key: "rho", Val: r.Rho},
			{Key: "cva", Val: r.CVA},
			{Key: "cva_se", Val: r.SE},
			{Key: "cva_1m", Val: r.CVA * notionalScale},
			{Key: "epe_exact_two_time_h", Val: colEPEByH[r.Threshold]},
			{Key: "epe_path_grid", Val: r.EPE},
			{Key: "cva_raw_coupon_1m", Val: r.CVARawCoupon * notionalScale},
			{Key: "epe_raw_coupon", Val: r.EPERawCoupon},
		})
	}

	mporObjs := make([]jsonout.Obj, 0, len(mporRows))
	mporRho0X := make([]float64, 0, len(mporDeltas))
	mporEPE := make([]float64, 0, len(mporDeltas))
	mporCVARho0 := make([]float64, 0, len(mporDeltas))
	mporCVARho08 := make([]float64, 0, len(mporDeltas))
	mporCVARawRho0 := make([]float64, 0, len(mporDeltas))
	c8bMPORCleanConv := true
	for _, rho := range colRhos {
		prev := math.Inf(1)
		for _, r := range mporRows {
			if r.Rho != rho {
				continue
			}
			if r.CVA >= prev {
				c8bMPORCleanConv = false
			}
			prev = r.CVA
		}
	}
	for _, r := range mporRows {
		if r.Rho == 0 {
			mporRho0X = append(mporRho0X, r.Delta)
			mporEPE = append(mporEPE, r.EPE)
			mporCVARho0 = append(mporCVARho0, r.CVA)
			mporCVARawRho0 = append(mporCVARawRho0, r.CVARawCoupon)
		} else {
			mporCVARho08 = append(mporCVARho08, r.CVA)
		}
		mporObjs = append(mporObjs, jsonout.Obj{
			{Key: "delta", Val: r.Delta},
			{Key: "threshold", Val: "0"},
			{Key: "rho", Val: r.Rho},
			{Key: "cva", Val: r.CVA},
			{Key: "cva_se", Val: r.SE},
			{Key: "cva_1m", Val: r.CVA * notionalScale},
			{Key: "epe", Val: r.EPE},
			{Key: "cva_raw_coupon", Val: r.CVARawCoupon},
			{Key: "cva_raw_coupon_1m", Val: r.CVARawCoupon * notionalScale},
			{Key: "epe_raw_coupon", Val: r.EPERawCoupon},
		})
	}
	mporSlopeRho0 := logLogSlope(mporRho0X, mporCVARho0)
	mporSlopeEPE := logLogSlope(mporRho0X, mporEPE)
	if mporSlopeRho0 <= 0.40 || len(mporCVARho0) == 0 || mporCVARho0[len(mporCVARho0)-1] >= 1e-4 {
		c8bMPORCleanConv = false
	}
	c8cMPORRawCouponConv := len(mporCVARawRho0) > 0 && mporCVARawRho0[len(mporCVARawRho0)-1] < 1e-4

	// =========================================================================
	// GENERATE 7 SVG FIGURES
	// =========================================================================
	cvaAnalLine := make([]float64, len(wwrRhos))
	for i := range cvaAnalLine {
		cvaAnalLine[i] = analCVA.CVAAnalytic
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_ee-pfe-profile.svg"), svg.ChartOptions{
		Title:  "Profil Eksposur Swap 10Y ATM (Hull-White 1-Faktor)",
		XLabel: "Waktu (Tahun)",
		YLabel: "Eksposur (per 1 nosional)",
	}, []svg.Series{
		{Label: "PFE 99% (MC)", X: tSeries, Y: pfe99Series, Color: svg.Palette[1]},
		{Label: "PFE 95% (MC)", X: tSeries, Y: pfe95Series, Color: svg.Palette[4]},
		{Label: "EE Kuadratur", X: tSeries, Y: eeGHSeries, Color: svg.Palette[0]},
		{Label: "EE MC (200k)", X: tSeries, Y: eeMCSeries, Color: svg.Palette[2], Dashed: true},
		{Label: "EE Terdiskonto", X: tSeries, Y: eeDiscSeries, Color: svg.Palette[3]},
	}); err != nil {
		panic(err)
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_cir-survival.svg"), svg.ChartOptions{
		Title:  "Kurva Probabilitas Survival CIR: Analitik vs RK4 vs MC",
		XLabel: "Jatuh Tempo T (Tahun)",
		YLabel: "Probabilitas Survival Q(T)",
	}, []svg.Series{
		{Label: "Analitik Afin", X: cirHorizons, Y: qAnalCurve, Color: svg.Palette[0]},
		{Label: "ODE Riccati RK4", X: cirHorizons, Y: qRK4Curve, Color: svg.Palette[1], Dashed: true},
		{Label: "MC Poisson-Gamma (40/thn)", X: cirHorizons, Y: qExact40, Color: svg.Palette[2]},
		{Label: "MC Euler FT (40/thn)", X: cirHorizons, Y: qEuler40, Color: svg.Palette[4], Dashed: true},
	}); err != nil {
		panic(err)
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_cva-vs-rho.svg"), svg.ChartOptions{
		Title:         "Wrong-Way Risk: CVA terhadap Korelasi Suku Bunga-Kredit (rho)",
		XLabel:        "Korelasi rho = Corr(dW^r, dW^lambda)",
		YLabel:        "CVA (per 1 nosional)",
		LegendTopLeft: true,
	}, []svg.Series{
		{Label: "CVA Rao-Blackwell (±1.96 SE)", X: wwrRhos, Y: wwrRB, Lower: wwrRBLo, Upper: wwrRBHi, Color: svg.Palette[0]},
		{Label: "CVA Tau Eksplisit u~Exp(1)", X: wwrRhos, Y: wwrTau, Lower: wwrTauLo, Upper: wwrTauHi, Color: svg.Palette[1], Dashed: true},
		{Label: "CVA Analitik Independen", X: wwrRhos, Y: cvaAnalLine, Color: svg.Palette[2], Dashed: true},
	}); err != nil {
		panic(err)
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_trapezoid-discretization-bias.svg"), svg.ChartOptions{
		Title:  "Konvergensi Diskretisasi Integral Intensitas (Log-Log)",
		XLabel: "Sub-langkah per Tahun (m)",
		YLabel: "Galat terhadap Referensi m=320",
		LogX:   true,
		LogY:   true,
	}, []svg.Series{
		{Label: fmt.Sprintf("RMSE Integral (slope=%.2f)", slopeIntRMSE), X: mXS, Y: intRMSEYS, Color: svg.Palette[0]},
		{Label: fmt.Sprintf("|Bias CVA Terpasang| (slope=%.2f)", slopeCVABias), X: mXS, Y: cvaBiasYS, Color: svg.Palette[1]},
	}); err != nil {
		panic(err)
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_netting-exposure.svg"), svg.ChartOptions{
		Title:  "Manfaat Netting: 5Y Payer + 10Y Receiver vs Jumlah Mandiri",
		XLabel: "Waktu (Tahun)",
		YLabel: "Expected Exposure (per 1 nosional)",
	}, []svg.Series{
		{Label: "Jumlah EE Mandiri (EE1 + EE2)", X: tGrid, Y: netRes.StandaloneSumEE, Color: svg.Palette[1]},
		{Label: "EE Portofolio Ter-netting", X: tGrid, Y: netRes.NettedEE, Color: svg.Palette[2]},
		{Label: "EE 10Y Receiver", X: tGrid, Y: netRes.Swap2EE, Color: svg.Palette[0], Dashed: true},
		{Label: "EE 5Y Payer", X: tGrid, Y: netRes.Swap1EE, Color: svg.Palette[4], Dashed: true},
	}); err != nil {
		panic(err)
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_collateral-thresholds.svg"), svg.ChartOptions{
		Title:  "Profil Eksposur Terkolateralisasi (MPOR delta = 0.04 thn)",
		XLabel: "Waktu (Tahun)",
		YLabel: "Expected Exposure Terkolateralisasi",
	}, []svg.Series{
		{Label: "H = Inf (Tanpa Kolateral)", X: tGrid, Y: colEEByH[math.Inf(1)], Color: svg.Palette[1]},
		{Label: "H = 0.05", X: tGrid, Y: colEEByH[0.05], Color: svg.Palette[4]},
		{Label: "H = 0.01", X: tGrid, Y: colEEByH[0.01], Color: svg.Palette[0]},
		{Label: "H = 0 (Full VM)", X: tGrid, Y: colEEByH[0], Color: svg.Palette[2]},
	}); err != nil {
		panic(err)
	}

	if err := svg.WriteLineChart(filepath.Join(figDir, "day003e_mpor-convergence.svg"), svg.ChartOptions{
		Title:         "Konvergensi Eksposur & CVA Kolateral H=0 saat delta -> 0",
		XLabel:        "Margin Period of Risk delta (Tahun)",
		YLabel:        "EPE / CVA (per 1 nosional)",
		LogX:          true,
		LogY:          true,
		LegendTopLeft: true,
	}, []svg.Series{
		{Label: fmt.Sprintf("EPE Jadwal Bersih (slope=%.2f)", mporSlopeEPE), X: mporRho0X, Y: mporEPE, Color: svg.Palette[0]},
		{Label: "CVA Kolateral rho=0.8 (bersih)", X: mporRho0X, Y: mporCVARho08, Color: svg.Palette[1]},
		{Label: fmt.Sprintf("CVA rho=0 bersih (slope=%.2f)", mporSlopeRho0), X: mporRho0X, Y: mporCVARho0, Color: svg.Palette[2]},
		{Label: "CVA rho=0 formula kupon mentah", X: mporRho0X, Y: mporCVARawRho0, Color: svg.Palette[4], Dashed: true},
	}); err != nil {
		panic(err)
	}

	elapsed := time.Since(tStart).Seconds()
	crit1Pass := martPass
	crit2Pass := swaptionPass && expPass
	crit3Pass := analCVA.DLambda0AbsErr < 1e-6 && analCVA.CS01Lambda0_1bp > 0
	crit4Pass := fellerOK && rk4Pass && survPass && cdsPass
	crit5Pass := discPass
	crit6Pass := wwrPass
	crit7Pass := netPass
	crit8Pass := c8aThreshMonotone && c8bMPORCleanConv && c8cMPORRawCouponConv
	c8Status := "belum tercapai"
	if c8aThreshMonotone && c8bMPORCleanConv && c8cMPORRawCouponConv {
		c8Status = "tercapai"
	} else if c8aThreshMonotone && c8bMPORCleanConv {
		c8Status = "sebagian"
	}

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day003e"},
			{Key: "date", Val: "2026-11-21"},
			{Key: "language", Val: "go"},
			{Key: "main_seed", Val: rng.MainSeed},
			{Key: "n_workers", Val: nWorkers},
			{Key: "n_chunks", Val: nChunks},
			{Key: "worker_invariance_1_vs_n_bit_identical", Val: workerBitIdentical},
			{Key: "n_paths_exposure", Val: nPathsExposure},
			{Key: "n_paths_core_mc", Val: nPathsCore},
			{Key: "notional_unit", Val: 1.0},
			{Key: "notional_reference_scale", Val: notionalScale},
			{Key: "elapsed_seconds", Val: elapsed},
		}},
		{Key: "parameters", Val: jsonout.Obj{
			{Key: "hull_white", Val: jsonout.Obj{
				{Key: "r0", Val: hwP.R0},
				{Key: "a", Val: hwP.A},
				{Key: "sigma", Val: hwP.Sigma},
				{Key: "par_swap_rate_10y", Val: swap10Pay.StrikeK},
				{Key: "par_swap_rate_5y", Val: swap5Pay.StrikeK},
				{Key: "forward_swap_rate_2y8y", Val: kFwd},
			}},
			{Key: "cir_credit", Val: jsonout.Obj{
				{Key: "lambda0", Val: cirP.Lambda0},
				{Key: "kappa", Val: cirP.Kappa},
				{Key: "theta", Val: cirP.Theta},
				{Key: "xi", Val: cirP.Xi},
				{Key: "recovery", Val: cirP.Recovery},
				{Key: "lgd", Val: cirP.LGD()},
				{Key: "feller_2kappa_theta", Val: fellerLHS},
				{Key: "feller_xi_sq", Val: fellerRHS},
				{Key: "feller_satisfied", Val: fellerOK},
			}},
		}},
		{Key: "part1_hull_white", Val: jsonout.Obj{
			{Key: "martingale_checks", Val: martObjs},
			{Key: "max_gh_martingale_abs_err", Val: maxMartGHErr},
			{Key: "max_mc_martingale_z", Val: maxMartMCZ},
			{Key: "swaption_2y8y", Val: swaptionObjs},
			{Key: "max_swaption_jamshidian_vs_quad2d_err", Val: maxJamErr},
			{Key: "max_swaption_mc_z", Val: maxSwapMCZ},
			{Key: "epe_analytic", Val: epeGH},
			{Key: "epe_analytic_1m", Val: epeGH * notionalScale},
			{Key: "epe_mc", Val: epeMC},
			{Key: "epe_mc_1m", Val: epeMC * notionalScale},
			{Key: "peak_ee", Val: peakEE},
			{Key: "peak_ee_time", Val: peakEET},
			{Key: "peak_pfe99", Val: peakPFE99},
			{Key: "peak_pfe99_time", Val: peakPFE99T},
			{Key: "max_ee_mc_vs_analytic_z", Val: maxEEZ},
			{Key: "exposure_profile", Val: expObjs},
		}},
		{Key: "part2_cir_credit", Val: jsonout.Obj{
			{Key: "survival_curve_analytic_vs_rk4", Val: cirCurveObjs},
			{Key: "max_survival_rk4_abs_err", Val: maxRK4Err},
			{Key: "survival_mc_checks", Val: survMCObjs},
			{Key: "par_cds_spreads", Val: cdsObjs},
		}},
		{Key: "part3_independent_cva", Val: jsonout.Obj{
			{Key: "cva_analytic", Val: analCVA.CVAAnalytic},
			{Key: "cva_analytic_1m", Val: analCVA.CVAAnalytic * notionalScale},
			{Key: "d_cva_d_lambda0_analytic", Val: analCVA.DLambda0Analytic},
			{Key: "d_cva_d_lambda0_central_fd", Val: analCVA.DLambda0CentralFD},
			{Key: "d_cva_d_lambda0_abs_err", Val: analCVA.DLambda0AbsErr},
			{Key: "d_cva_d_lambda0_rel_err", Val: analCVA.DLambda0RelErr},
			{Key: "cs01_lambda0_1bp", Val: analCVA.CS01Lambda0_1bp},
			{Key: "cs01_lambda0_1bp_1m", Val: analCVA.CS01Lambda0_1bp * notionalScale},
			{Key: "cs01_theta_1bp", Val: analCVA.CS01Theta_1bp},
			{Key: "cs01_theta_1bp_1m", Val: analCVA.CS01Theta_1bp * notionalScale},
			{Key: "cs01_parallel_1bp", Val: analCVA.CS01Parallel_1bp},
			{Key: "cs01_parallel_1bp_1m", Val: analCVA.CS01Parallel_1bp * notionalScale},
			{Key: "d_cva_d_sigma_hw_fixed_k", Val: analCVA.DSigmaHWFixedK},
			{Key: "cva_sigma_bump_10bp_fixed_k_1m", Val: analCVA.SigmaBump10bpFixed * notionalScale},
			{Key: "d_cva_d_strike_k", Val: analCVA.DStrikeK},
			{Key: "cva_strike_bump_1bp_1m", Val: analCVA.StrikeBump1bp * notionalScale},
			{Key: "discretization_mc_convergence", Val: discObjs},
			{Key: "loglog_slope_paired_cva_bias", Val: slopeCVABias},
			{Key: "loglog_slope_integral_rmse", Val: slopeIntRMSE},
		}},
		{Key: "part4_wwr_netting_collateral", Val: jsonout.Obj{
			{Key: "wwr_curve", Val: wwrObjs},
			{Key: "netting", Val: jsonout.Obj{
				{Key: "swap1_5y_payer_cva_analytic", Val: netRes.Swap1CVAAnalytic},
				{Key: "swap2_10y_receiver_cva_analytic", Val: netRes.Swap2CVAAnalytic},
				{Key: "netted_cva_analytic", Val: netRes.NettedCVAAnalytic},
				{Key: "max_pointwise_subadditivity_violation", Val: netRes.MaxPointwiseViol},
				{Key: "opposite_sign_10y_netted_max_ee", Val: netRes.OppositeNetMaxEE},
				{Key: "opposite_sign_10y_netted_cva", Val: netRes.OppositeNetCVA},
				{Key: "by_rho", Val: netRhoObjs},
			}},
			{Key: "collateral_by_threshold_delta004", Val: colThreshObjs},
			{Key: "collateral_mpor_convergence_h0", Val: mporObjs},
			{Key: "mpor_loglog_slope_cva_clean_rho0", Val: mporSlopeRho0},
			{Key: "mpor_loglog_slope_epe_clean", Val: mporSlopeEPE},
		}},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_hull_white_martingale", Val: crit1Pass},
			{Key: "c2_jamshidian_and_ee_profile", Val: crit2Pass},
			{Key: "c3_cva_d_lambda0_analytic_vs_fd", Val: crit3Pass},
			{Key: "c4_cir_feller_rk4_exact_cds", Val: crit4Pass},
			{Key: "c5_independent_cva_mc_and_discretization", Val: crit5Pass},
			{Key: "c6_wwr_monotonic_cva_rho", Val: crit6Pass},
			{Key: "c7_netting_subadditivity_and_zero_identity", Val: crit7Pass},
			{Key: "c8a_collateral_threshold_monotonicity", Val: c8aThreshMonotone},
			{Key: "c8b_collateral_mpor_clean_schedule_convergence", Val: c8bMPORCleanConv},
			{Key: "c8c_collateral_mpor_raw_coupon_convergence", Val: c8cMPORRawCouponConv},
			{Key: "c8_collateral_threshold_and_mpor_convergence", Val: crit8Pass},
			{Key: "c8_status", Val: c8Status},
			{Key: "all_passed", Val: crit1Pass && crit2Pass && crit3Pass && crit4Pass && crit5Pass && crit6Pass && crit7Pass && crit8Pass},
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		panic(err)
	}
	fmt.Printf("day003e selesai dalam %.2fs | CVA anal=%.6f ($%.2f/1M) | dCVA/dlam0 err=%.3g | all_passed=%v\n",
		elapsed, analCVA.CVAAnalytic, analCVA.CVAAnalytic*notionalScale, analCVA.DLambda0AbsErr,
		crit1Pass && crit2Pass && crit3Pass && crit4Pass && crit5Pass && crit6Pass && crit7Pass && crit8Pass)
}
