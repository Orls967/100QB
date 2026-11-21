package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/afns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/interp"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/jsonout"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/kalman"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/parallel"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/pca"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/risk"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/svg"
)

const (
	nChunks = 32
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

func mat3ToSlice(m [3][3]float64) [][]float64 {
	out := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		out[i] = []float64{m[i][0], m[i][1], m[i][2]}
	}
	return out
}

func main() {
	tTotalStart := time.Now()
	root := findRepoRoot()
	outPath := filepath.Join(root, "outputs", "day004a", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	nWorkers := parallel.DefaultWorkers()

	// =========================================================================
	// BAGIAN 1: Nelson-Siegel(-Svensson), Bootstrap, Interpolasi & Risiko
	// =========================================================================
	tP1Start := time.Now()
	nsP := ns.DefaultNS()
	nssP := ns.DefaultNSS()

	idTaus := []float64{0.1, 1.0, 5.0, 30.0}
	idObjs := make([]jsonout.Obj, 0, len(idTaus)*2)
	maxIdentityErr := 0.0
	for _, tc := range []struct {
		name string
		p    ns.Params
	}{{"NS", nsP}, {"NSS", nssP}} {
		for _, tau := range idTaus {
			y := tc.p.Yield(tau)
			avgF := tc.p.ForwardIntegralAvg(tau, 32)
			err := math.Abs(y - avgF)
			if err > maxIdentityErr {
				maxIdentityErr = err
			}
			idObjs = append(idObjs, jsonout.Obj{
				{Key: "curve", Val: tc.name},
				{Key: "tau", Val: tau},
				{Key: "yield", Val: y},
				{Key: "gl32_avg_forward", Val: avgF},
				{Key: "abs_err", Val: err},
			})
		}
	}

	nsLim0 := nsP.Yield(1e-12)
	nsLimInf := nsP.Yield(1e9)
	nssLim0 := nssP.Yield(1e-12)
	nssLimInf := nssP.Yield(1e9)
	limitsOK := math.Abs(nsLim0-(nsP.L+nsP.S)) < 1e-11 && math.Abs(nsLimInf-nsP.L) < 1e-8 &&
		math.Abs(nssLim0-(nssP.L+nssP.S)) < 1e-11 && math.Abs(nssLimInf-nssP.L) < 1e-8

	smallTaus := []float64{1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16}
	smallRows := nsP.CompareSmallTauStability(smallTaus)
	smallObjs := make([]jsonout.Obj, 0, len(smallRows))
	for _, r := range smallRows {
		smallObjs = append(smallObjs, jsonout.Obj{
			{Key: "tau", Val: r.Tau},
			{Key: "g1_stable", Val: r.G1Stable},
			{Key: "g1_naive", Val: r.G1Naive},
			{Key: "g1_rel_err", Val: r.G1RelErr},
			{Key: "g2_stable", Val: r.G2Stable},
			{Key: "g2_naive", Val: r.G2Naive},
			{Key: "g2_rel_err", Val: r.G2RelErr},
			{Key: "yield_stable", Val: r.YieldStable},
			{Key: "yield_naive", Val: r.YieldNaive},
			{Key: "yield_abs_err", Val: r.YieldAbsErr},
		})
	}

	// Bootstrap dari kupon par semi-tahunan 0.5..30 tahun
	mats := bootstrap.SemiAnnualGrid(30.0)
	bootSummary := make([]jsonout.Obj, 0, 2)
	maxBootErr := 0.0
	for _, tc := range []struct {
		name string
		p    ns.Params
	}{{"NS", nsP}, {"NSS", nssP}} {
		par, trueDF, trueY := bootstrap.ParRatesFromCurve(tc.p, mats, 0.5)
		recDF, recY, err := bootstrap.BootstrapFromParRates(mats, par, 0.5)
		if err != nil {
			panic(err)
		}
		mDF, mY := 0.0, 0.0
		for i := range mats {
			if d := math.Abs(recDF[i] - trueDF[i]); d > mDF {
				mDF = d
			}
			if d := math.Abs(recY[i] - trueY[i]); d > mY {
				mY = d
			}
		}
		if mDF > maxBootErr {
			maxBootErr = mDF
		}
		if mY > maxBootErr {
			maxBootErr = mY
		}
		samplePts := make([]jsonout.Obj, 0, 7)
		for _, idx := range []int{0, 1, 3, 9, 19, 39, 59} {
			samplePts = append(samplePts, jsonout.Obj{
				{Key: "maturity", Val: mats[idx]},
				{Key: "par_rate", Val: par[idx]},
				{Key: "true_df", Val: trueDF[idx]},
				{Key: "boot_df", Val: recDF[idx]},
				{Key: "df_abs_err", Val: math.Abs(recDF[idx] - trueDF[idx])},
				{Key: "zero_yield", Val: recY[idx]},
			})
		}
		bootSummary = append(bootSummary, jsonout.Obj{
			{Key: "curve", Val: tc.name},
			{Key: "n_maturities", Val: len(mats)},
			{Key: "max_df_abs_err", Val: mDF},
			{Key: "max_yield_abs_err", Val: mY},
			{Key: "sample_points", Val: samplePts},
		})
	}

	// Interpolasi & derau par R = 200 replikasi (sigma = 0.5 bp)
	interpRes := interp.RunNoisyBootstrapComparison(nsP, 0.5e-4, 200, nChunks, nWorkers, rng.MainSeed)
	interp1W := interp.RunNoisyBootstrapComparison(nsP, 0.5e-4, 32, nChunks, 1, rng.MainSeed)
	interpNW := interp.RunNoisyBootstrapComparison(nsP, 0.5e-4, 32, nChunks, nWorkers, rng.MainSeed)
	interpWorkerEq := math.Float64bits(interp1W.LogLinearP.ForwardTVMean) == math.Float64bits(interpNW.LogLinearP.ForwardTVMean) &&
		math.Float64bits(interp1W.CubicSpline.ForwardTVMean) == math.Float64bits(interpNW.CubicSpline.ForwardTVMean)

	// Risiko obligasi 10y 4%, KRD, PV01 swap 10y, lindung nilai 2 obligasi
	const bumpDur = 1e-6
	const bumpConv = 1e-4
	bondRisk := risk.AnalyzeBondRisk(nsP, 10.0, 0.04, 2.0, bumpDur, bumpConv)
	tBond, cfBond := risk.BondCashFlows(10.0, 0.04, 2.0)
	krdReport := risk.ComputeKeyRateDurations(nsP, tBond, cfBond, risk.StandardKeyRateNodes, bumpDur)
	swapPV01 := risk.AnalyzeSwapPV01(nsP, 10.0, 2.0)
	hedgeRep, err := risk.HedgeTwoKRDs(nsP, risk.StandardKeyRateNodes, 3, 5, bumpDur)
	if err != nil {
		panic(err)
	}

	krdObjs := make([]jsonout.Obj, len(krdReport.Nodes))
	krdCats := make([]string, len(krdReport.Nodes))
	for i, node := range krdReport.Nodes {
		krdCats[i] = fmt.Sprintf("%.0fy", node)
		krdObjs[i] = jsonout.Obj{
			{Key: "node_year", Val: node},
			{Key: "krd_fd", Val: krdReport.KRDFD[i]},
			{Key: "krd_closed", Val: krdReport.KRDClosed[i]},
			{Key: "dollar_krd", Val: krdReport.DollarKRDClosed[i]},
			{Key: "hedged_dollar_krd", Val: hedgeRep.HedgedDollarKRDs[i]},
		}
	}
	elapsedP1 := time.Since(tP1Start).Seconds()

	// =========================================================================
	// BAGIAN 2: AFNS, ODE Riccati, Dinamika P & Premi Jangka Waktu
	// =========================================================================
	tP2Start := time.Now()
	afnsP := afns.DefaultParams()

	// 2.1 Loadings & Riccati tertutup vs RK4 + Richardson
	riccatiTaus := []float64{0.05, 0.25, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 25, 30}
	maxLoadErr, maxRicBErr, maxRicAErr, maxRicAdjErr := 0.0, 0.0, 0.0, 0.0
	adjNegMono := true
	prevAdj := 0.0
	ricObjs := make([]jsonout.Obj, 0, len(riccatiTaus))
	var adj10, adj30 [4]float64
	for _, tau := range riccatiTaus {
		h := afnsP.YieldLoadings(tau)
		dLoad := math.Max(math.Abs(h[0]-1.0), math.Max(math.Abs(h[1]-ns.G1(tau, afnsP.Lambda)), math.Abs(h[2]-ns.G2(tau, afnsP.Lambda))))
		if dLoad > maxLoadErr {
			maxLoadErr = dLoad
		}
		bCl := afnsP.BVec(tau)
		a1, a2, a3, adjTot := afnsP.YieldAdjustmentComponents(tau)
		_, aCl := afnsP.YieldAdjustment(tau)
		if adjTot >= 0 || adjTot >= prevAdj {
			adjNegMono = false
		}
		prevAdj = adjTot
		if tau == 10.0 {
			adj10 = [4]float64{a1, a2, a3, adjTot}
		}
		if tau == 30.0 {
			adj30 = [4]float64{a1, a2, a3, adjTot}
		}

		steps := int(math.Ceil(tau * 400))
		if steps < 400 {
			steps = 400
		}
		bRK, aRK, bRich, aRich := afnsP.SolveRiccatiRK4(tau, steps)
		for i := 0; i < 3; i++ {
			if d := math.Abs(bCl[i] - bRK[i]); d > maxRicBErr {
				maxRicBErr = d
			}
			if d := math.Abs(bCl[i] - bRich[i]); d > maxRicBErr {
				maxRicBErr = d
			}
		}
		eARK := math.Abs(aCl - aRK)
		eARich := math.Abs(aCl - aRich)
		if eARK > maxRicAErr {
			maxRicAErr = eARK
		}
		if eARich > maxRicAErr {
			maxRicAErr = eARich
		}
		if d := math.Abs(adjTot - (-aRich / tau)); d > maxRicAdjErr {
			maxRicAdjErr = d
		}
		ricObjs = append(ricObjs, jsonout.Obj{
			{Key: "tau", Val: tau},
			{Key: "adj_level", Val: a1},
			{Key: "adj_slope", Val: a2},
			{Key: "adj_curvature", Val: a3},
			{Key: "adj_total", Val: adjTot},
			{Key: "adj_total_bps", Val: adjTot * 10_000.0},
			{Key: "a_closed", Val: aCl},
			{Key: "a_rk4", Val: aRK},
			{Key: "a_richardson", Val: aRich},
			{Key: "abs_err_a_rk4", Val: eARK},
			{Key: "abs_err_a_rich", Val: eARich},
		})
	}

	// 2.3 Transisi P: segitiga bawah penuh vs MatExpSeries/VanLoan/GL + MC N=200.000
	fullLowP := afns.DefaultParams()
	fullLowP.KP = [3][3]float64{{0.05, 0, 0}, {0.12, 0.5, 0}, {-0.08, 0.25, 0.4}}
	phiFullCl := fullLowP.TransitionMatrixClosed(0.75)
	negKdt := lag.Zeros(3, 3)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			negKdt[i][j] = -fullLowP.KP[i][j] * 0.75
		}
	}
	phiFullSer := lag.MatExpSeries(negKdt, 35)
	maxPhiMatExpErr := 0.0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if d := math.Abs(phiFullCl[i][j] - phiFullSer[i][j]); d > maxPhiMatExpErr {
				maxPhiMatExpErr = d
			}
		}
	}
	qFullGL := fullLowP.TransitionCovGL(0.75, 32)
	qFullVL := fullLowP.TransitionCovVanLoan(0.75)
	maxQVanLoanGLErr := 0.0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if d := math.Abs(qFullGL[i][j] - qFullVL[i][j]); d > maxQVanLoanGLErr {
				maxQVanLoanGLErr = d
			}
		}
	}
	// Juga cek model default diagonal Q_dt tertutup vs GL vs VanLoan pada dt=1/52
	qDiagCl := afnsP.TransitionCovClosed(kalman.WeeklyDt)
	qDiagVL := afnsP.TransitionCovVanLoan(kalman.WeeklyDt)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if d := math.Abs(qDiagCl[i][j] - qDiagVL[i][j]); d > maxQVanLoanGLErr {
				maxQVanLoanGLErr = d
			}
		}
	}

	x0MC := [3]float64{0.052, -0.015, 0.008}
	transMC := afnsP.VerifyTransitionMC(x0MC, 1.0, 200_000, nChunks, nWorkers, rng.MainSeed)
	transMC1W := afnsP.VerifyTransitionMC(x0MC, 1.0, 20_000, nChunks, 1, rng.MainSeed)
	transMCNW := afnsP.VerifyTransitionMC(x0MC, 1.0, 20_000, nChunks, nWorkers, rng.MainSeed)
	transWorkerEq := math.Float64bits(transMC1W.MaxZScore) == math.Float64bits(transMCNW.MaxZScore)

	// 2.4 Harga risiko & Term Premium 10y
	lam0, lam1 := afnsP.MarketPriceOfRisk()
	tpCl, lrY10, lrR, tpAdjC, tpSlopeC, tpCurvC := afnsP.LongRunTermPremium(10.0)
	tpMC := afnsP.VerifyTermPremiumMC(10.0, 40, 100_000, nChunks, nWorkers, rng.MainSeed)
	elapsedP2 := time.Since(tP2Start).Seconds()

	// =========================================================================
	// BAGIAN 3: Simulasi, Filter Kalman, RTS Smoother, MLE & Replikasi R=50
	// =========================================================================
	tP3Start := time.Now()
	spec := kalman.DefaultSpec()

	// 3.2 Anchor eksak T=5, 3 tenor {1, 5, 10}
	rgAnchor := rng.New(rng.DeriveChunkSeed(rng.MainSeed, 0x3101, 0))
	_, yAnchor := afns.SimulateStationaryPath(spec.AFNS, []float64{1.0, 5.0, 10.0}, 5, kalman.WeeklyDt, spec.SigmaEps, rgAnchor)
	anchorT5, err := kalman.VerifyJointGaussianAnchor(spec, []float64{1.0, 5.0, 10.0}, kalman.WeeklyDt, yAnchor)
	if err != nil {
		panic(err)
	}

	// 3.1 & 3.3 Satu realisasi T=1000 minggu (10 tenor), Filter + Smoother + Single MLE
	rgSingle := rng.New(rng.DeriveChunkSeed(rng.MainSeed, 0x3201, 0))
	trueStates1000, yields1000 := afns.SimulateStationaryPath(spec.AFNS, kalman.StandardTenors, 1000, kalman.WeeklyDt, spec.SigmaEps, rgSingle)
	filtSmooth1000, err := kalman.FilterAndSmooth(spec, kalman.StandardTenors, kalman.WeeklyDt, yields1000)
	if err != nil {
		panic(err)
	}
	initGuess := kalman.PerturbInitialGuess(spec, rgSingle)
	singleMLE := kalman.EstimateMLE(initGuess, spec, kalman.StandardTenors, kalman.WeeklyDt, yields1000)

	// 3.4 Studi replikasi R=50 simulasi (T=1000) + uji determinisme 1 vs nWorkers
	repStudy := kalman.RunReplicationStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 1000, 50, nChunks, nWorkers, rng.MainSeed)
	rep1W := kalman.RunReplicationStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 300, 4, 4, 1, rng.MainSeed)
	repNW := kalman.RunReplicationStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 300, 4, 4, nWorkers, rng.MainSeed)
	kalmanWorkerEq := true
	for i := 0; i < 11; i++ {
		if math.Float64bits(rep1W.Params[i].MeanEstimate) != math.Float64bits(repNW.Params[i].MeanEstimate) {
			kalmanWorkerEq = false
		}
	}

	// 3.5 Kesalahan spesifikasi DNS vs AFNS
	dnsMisspec := kalman.RunDNSMisspecification(spec, kalman.StandardTenors, kalman.WeeklyDt, 1000, rng.MainSeed)
	elapsedP3 := time.Since(tP3Start).Seconds()

	// =========================================================================
	// BAGIAN 4: PCA dan Struktur Kovarians Perubahan Yield (T=5000, R=100)
	// =========================================================================
	tP4Start := time.Now()
	pcaRes, err := pca.RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 5000, 100, nChunks, nWorkers, rng.MainSeed)
	if err != nil {
		panic(err)
	}
	pca1W, err := pca.RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 1000, 8, 4, 1, rng.MainSeed)
	if err != nil {
		panic(err)
	}
	pcaNW, err := pca.RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 1000, 8, 4, nWorkers, rng.MainSeed)
	if err != nil {
		panic(err)
	}
	pcaWorkerEq := math.Float64bits(pca1W.EmpMeanEigenvalues[0]) == math.Float64bits(pcaNW.EmpMeanEigenvalues[0])
	elapsedP4 := time.Since(tP4Start).Seconds()

	workerAllIdentical := interpWorkerEq && transWorkerEq && kalmanWorkerEq && pcaWorkerEq

	// =========================================================================
	// EVALUASI KRITERIA PENERIMAAN (JUJUR, TANPA PENYETELAN)
	// =========================================================================
	crit1 := maxIdentityErr < 1e-12 && limitsOK
	crit2 := maxBootErr < 1e-12
	crit3 := bondRisk.ModDurationYTMRelErr < 1e-6 && bondRisk.ConvexityYTMRelErr < 1e-6 &&
		bondRisk.EffDurationZeroRelErr < 1e-6 && bondRisk.EffConvexityZeroRelErr < 1e-6 &&
		krdReport.SumVsParFDErr < 1e-8 && math.Abs(hedgeRep.ResidualNode1) < 1e-12 && math.Abs(hedgeRep.ResidualNode2) < 1e-12
	crit4 := maxLoadErr < 1e-12 && maxRicBErr < 1e-10 && maxRicAErr < 1e-10 && maxRicAdjErr < 1e-10 && adjNegMono
	crit5 := maxPhiMatExpErr < 1e-12 && maxQVanLoanGLErr < 1e-12 && transMC.MaxZScore <= 3.0 && tpMC.ZScore <= 3.0
	crit6 := anchorT5.LogLikAbsDiff < 1e-9 && anchorT5.MaxSmoothMeanAbsDiff < 1e-9 && anchorT5.MaxSmoothCovAbsDiff < 1e-9
	crit7 := singleMLE.GradRelErrInit < 1e-5 && singleMLE.GradRelErrTrue < 1e-5 && singleMLE.EstLogLik >= singleMLE.TrueLogLik-1e-6

	// Kriteria pra-jalan replikasi R=50 untuk parameter teridentifikasi kuat:
	// lambda (0), sigma1..3 (7, 8, 9), sigma_eps (10): |bias| < 2.0 SE_rep dan cakupan 95% dalam [0.88, 0.99]
	wellIdentifiedIdx := map[int]bool{0: true, 7: true, 8: true, 9: true, 10: true}
	wellIdBiasOK := true
	wellIdCovOK := true
	allParamsBiasOK := true
	allParamsCovOK := true
	for i := 0; i < 11; i++ {
		pr := repStudy.Params[i]
		if pr.BiasInRepSE >= 2.0 {
			allParamsBiasOK = false
			if wellIdentifiedIdx[i] {
				wellIdBiasOK = false
			}
		}
		if pr.Coverage95 < 0.88 || pr.Coverage95 > 0.99 {
			allParamsCovOK = false
			if wellIdentifiedIdx[i] {
				wellIdCovOK = false
			}
		}
	}
	crit8WellId := wellIdBiasOK && wellIdCovOK
	crit8All11 := allParamsBiasOK && allParamsCovOK
	c8Status := "belum tercapai"
	if crit8All11 {
		c8Status = "tercapai"
	} else if crit8WellId {
		c8Status = "sebagian"
	}

	pcaAcc := pcaRes.EvaluateAcceptance()
	allPassed := crit1 && crit2 && crit3 && crit4 && crit5 && crit6 && crit7 && crit8All11 && pcaAcc.AggregatePassed && workerAllIdentical

	// =========================================================================
	// GENERATE 8 SVG FIGURES
	// =========================================================================
	// Fig 1: NS & NSS kurva serta komponen (L, S*g1, C*g2)
	figTaus := make([]float64, 120)
	yNSSeries := make([]float64, 120)
	yNSSSeries := make([]float64, 120)
	compL := make([]float64, 120)
	compSG1 := make([]float64, 120)
	compCG2 := make([]float64, 120)
	for i := 0; i < 120; i++ {
		tau := 0.25 * float64(i+1)
		figTaus[i] = tau
		yNSSeries[i] = nsP.Yield(tau) * 100.0
		yNSSSeries[i] = nssP.Yield(tau) * 100.0
		compL[i] = nsP.L * 100.0
		compSG1[i] = nsP.S * ns.G1(tau, nsP.Lambda) * 100.0
		compCG2[i] = nsP.C * ns.G2(tau, nsP.Lambda) * 100.0
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_ns-curves-and-loadings.svg"), svg.ChartOptions{
		Title:         "Kurva Nelson-Siegel, Svensson & Komponen Faktor (%)",
		XLabel:        "Maturitas tau (Tahun)",
		YLabel:        "Suku Bunga / Kontribusi (%)",
		ZeroLine:      true,
		LegendTopLeft: false,
	}, []svg.Series{
		{Label: "NSS Yield y(tau)", X: figTaus, Y: yNSSSeries, Color: svg.Palette[1]},
		{Label: "NS Yield y(tau)", X: figTaus, Y: yNSSeries, Color: svg.Palette[0]},
		{Label: "Level (L = 4.5%)", X: figTaus, Y: compL, Color: svg.Palette[2], Dashed: true},
		{Label: "Slope S*g1(tau)", X: figTaus, Y: compSG1, Color: svg.Palette[4], Dashed: true},
		{Label: "Curvature C*g2(tau)", X: figTaus, Y: compCG2, Color: svg.Palette[3], Dashed: true},
	}); err != nil {
		panic(err)
	}

	// Fig 2: Yield vs Forward & Stabilitas Interpolasi di Bawah Derau Par (0.5 bp)
	evalTaus := interpRes.LogLinearP.SampleRepEvalTaus
	trueYGrid := make([]float64, len(evalTaus))
	trueFGrid := make([]float64, len(evalTaus))
	llFGrid := make([]float64, len(evalTaus))
	spFGrid := make([]float64, len(evalTaus))
	for i, tau := range evalTaus {
		trueYGrid[i] = nsP.Yield(tau) * 100.0
		trueFGrid[i] = nsP.Forward(tau) * 100.0
		llFGrid[i] = interpRes.LogLinearP.SampleRepLogLinF[i] * 100.0
		spFGrid[i] = interpRes.CubicSpline.SampleRepSplineF[i] * 100.0
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_yield-vs-forward-interp.svg"), svg.ChartOptions{
		Title:        "Forward Seketika vs Yield di Bawah Derau Par (0.5 bp, Contoh 1 Replikasi)",
		XLabel:       "Maturitas tau (Tahun)",
		YLabel:       "Tingkat Bunga (%)",
		LegendBottom: true,
	}, []svg.Series{
		{Label: "Spline Kubik Y -> Forward f(tau)", X: evalTaus, Y: spFGrid, Color: svg.Palette[1]},
		{Label: "Log-Linear P -> Forward f(tau)", X: evalTaus, Y: llFGrid, Color: svg.Palette[4]},
		{Label: "NS Forward Benar f(tau)", X: evalTaus, Y: trueFGrid, Color: svg.Palette[0]},
		{Label: "NS Zero Yield Benar y(tau)", X: evalTaus, Y: trueYGrid, Color: svg.Palette[2], Dashed: true},
	}); err != nil {
		panic(err)
	}

	// Fig 3: Key-Rate Duration & Lindung Nilai Dua Obligasi
	if err := svg.WriteGroupedBarChart(filepath.Join(figDir, "day004a_key-rate-durations.svg"), svg.ChartOptions{
		Title:         "Profil Dollar Key-Rate Duration: Obligasi 10Y 4% vs Terlindung Nilai (5Y+10Y)",
		XLabel:        "Simpul Key-Rate (Tahun)",
		YLabel:        "Dollar KRD (per 1 Nosional)",
		LegendTopLeft: true,
	}, krdCats, []svg.BarSeries{
		{Label: "Target 10Y 4% (Tanpa Hedge)", Values: hedgeRep.TargetDollarKRDs, Color: svg.Palette[0]},
		{Label: "Portofolio Terlindung Nilai (Netral 5Y & 10Y)", Values: hedgeRep.HedgedDollarKRDs, Color: svg.Palette[2]},
	}); err != nil {
		panic(err)
	}

	// Fig 4: Suku Penyesuaian AFNS -A(tau)/tau (bps) vs tau
	adjTaus := make([]float64, 60)
	adjTotBps := make([]float64, 60)
	adjRK4Bps := make([]float64, 60)
	adjL1Bps := make([]float64, 60)
	adjS2Bps := make([]float64, 60)
	adjC3Bps := make([]float64, 60)
	for i := 0; i < 60; i++ {
		tau := 0.5 * float64(i+1)
		adjTaus[i] = tau
		a1, a2, a3, at := afnsP.YieldAdjustmentComponents(tau)
		_, _, _, aRich := afnsP.SolveRiccatiRK4(tau, int(math.Max(400, math.Ceil(tau*300))))
		adjTotBps[i] = at * 10_000.0
		adjRK4Bps[i] = (-aRich / tau) * 10_000.0
		adjL1Bps[i] = a1 * 10_000.0
		adjS2Bps[i] = a2 * 10_000.0
		adjC3Bps[i] = a3 * 10_000.0
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_afns-yield-adjustment.svg"), svg.ChartOptions{
		Title:        "Suku Penyesuaian Bebas-Arbitrase AFNS -A(tau)/tau dan Komponennya",
		XLabel:       "Jatuh Tempo tau (Tahun)",
		YLabel:       "Penyesuaian Yield (bps)",
		ZeroLine:     true,
		LegendBottom: true,
	}, []svg.Series{
		{Label: "Total Tertutup -A(tau)/tau", X: adjTaus, Y: adjTotBps, Color: svg.Palette[0]},
		{Label: "Verifikasi ODE Riccati RK4+Rich", X: adjTaus, Y: adjRK4Bps, Color: svg.Palette[1], Dashed: true},
		{Label: "Komponen Level (-A1/tau)", X: adjTaus, Y: adjL1Bps, Color: svg.Palette[2]},
		{Label: "Komponen Kurvatur (-A3/tau)", X: adjTaus, Y: adjC3Bps, Color: svg.Palette[3]},
		{Label: "Komponen Slope (-A2/tau)", X: adjTaus, Y: adjS2Bps, Color: svg.Palette[4]},
	}); err != nil {
		panic(err)
	}

	// Fig 5: Faktor Keadaan True vs RTS Smoother Kalman (250 minggu pertama)
	nPlotWeeks := 250
	weeks := make([]float64, nPlotWeeks)
	x1True := make([]float64, nPlotWeeks)
	x1Smooth := make([]float64, nPlotWeeks)
	x2True := make([]float64, nPlotWeeks)
	x2Smooth := make([]float64, nPlotWeeks)
	for t := 0; t < nPlotWeeks; t++ {
		weeks[t] = float64(t)
		x1True[t] = trueStates1000[t][0] * 100.0
		x1Smooth[t] = filtSmooth1000.XSmooth[t][0] * 100.0
		x2True[t] = trueStates1000[t][1] * 100.0
		x2Smooth[t] = filtSmooth1000.XSmooth[t][1] * 100.0
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_kalman-smoothed-factors.svg"), svg.ChartOptions{
		Title:         "Faktor Laten AFNS: Simulasi Benar vs RTS Smoother Kalman (Minggu 0-249)",
		XLabel:        "Minggu (t)",
		YLabel:        "Nilai Faktor (%)",
		ZeroLine:      true,
		LegendTopLeft: false,
	}, []svg.Series{
		{Label: "Level X1 Benar", X: weeks, Y: x1True, Color: svg.Palette[0]},
		{Label: "Level X1 RTS Smoothed", X: weeks, Y: x1Smooth, Color: svg.Palette[1], Dashed: true},
		{Label: "Slope X2 Benar", X: weeks, Y: x2True, Color: svg.Palette[2]},
		{Label: "Slope X2 RTS Smoothed", X: weeks, Y: x2Smooth, Color: svg.Palette[4], Dashed: true},
	}); err != nil {
		panic(err)
	}

	// Fig 6: Replikasi MLE R=50 (Rasio Estimasi / Benar +- 1.96 SD)
	pIdxX := make([]float64, 11)
	refOne := make([]float64, 11)
	singleRatio := make([]float64, 11)
	repMeanRatio := make([]float64, 11)
	repLoRatio := make([]float64, 11)
	repHiRatio := make([]float64, 11)
	for i := 0; i < 11; i++ {
		pIdxX[i] = float64(i + 1)
		refOne[i] = 1.0
		tv := repStudy.Params[i].TrueValue
		singleRatio[i] = singleMLE.EstParams[i] / tv
		repMeanRatio[i] = repStudy.Params[i].MeanEstimate / tv
		band := 1.96 * repStudy.Params[i].EmpiricalSD / math.Abs(tv)
		repLoRatio[i] = repMeanRatio[i] - band
		repHiRatio[i] = repMeanRatio[i] + band
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_mle-true-vs-estimated.svg"), svg.ChartOptions{
		Title:         "Studi Replikasi MLE (R=50, T=1000): Rasio Estimasi/Benar (±1.96 SD)",
		XLabel:        "Indeks Parameter (1=lam, 2..4=KP, 5..7=thetaP, 8..10=sigma, 11=sig_eps)",
		YLabel:        "Rasio terhadap Nilai Benar (1.0 = Eksak)",
		LegendTopLeft: true,
	}, []svg.Series{
		{Label: "Mean R=50 (±1.96 SD Empiris)", X: pIdxX, Y: repMeanRatio, Lower: repLoRatio, Upper: repHiRatio, Color: svg.Palette[0]},
		{Label: "Single MLE T=1000", X: pIdxX, Y: singleRatio, Color: svg.Palette[1], Dashed: true},
		{Label: "Referensi Nilai Benar (1.0)", X: pIdxX, Y: refOne, Color: svg.Palette[2], Dashed: true},
	}); err != nil {
		panic(err)
	}

	// Fig 7: PCA Scree Plot (Teoretis vs Empiris R=100, T=5000)
	compX := make([]float64, len(kalman.StandardTenors))
	screeTheo := make([]float64, len(kalman.StandardTenors))
	screeEmp := make([]float64, len(kalman.StandardTenors))
	screeLo := make([]float64, len(kalman.StandardTenors))
	screeHi := make([]float64, len(kalman.StandardTenors))
	for i := range compX {
		compX[i] = float64(i + 1)
		screeTheo[i] = pcaRes.TheoEigenvalues[i] * 1e8 // (bps)^2
		screeEmp[i] = pcaRes.EmpMeanEigenvalues[i] * 1e8
		screeLo[i] = math.Max(1e-3, (pcaRes.EmpMeanEigenvalues[i]-1.96*pcaRes.EmpSDEigenvalues[i])*1e8)
		screeHi[i] = (pcaRes.EmpMeanEigenvalues[i] + 1.96*pcaRes.EmpSDEigenvalues[i]) * 1e8
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_pca-scree-eigenvalues.svg"), svg.ChartOptions{
		Title:         "Scree Plot PCA Perubahan Yield Mingguan (Teoretis vs T=5000, R=100)",
		XLabel:        "Komponen Utama ke-k",
		YLabel:        "Nilai Eigen (bps^2, Skala Log)",
		LogY:          true,
		LegendTopLeft: false,
	}, []svg.Series{
		{Label: "Teoretis Jacobi", X: compX, Y: screeTheo, Color: svg.Palette[0]},
		{Label: "Empiris Mean R=100 (±1.96 SD)", X: compX, Y: screeEmp, Lower: screeLo, Upper: screeHi, Color: svg.Palette[1], Dashed: true},
	}); err != nil {
		panic(err)
	}

	// Fig 8: Tiga Eigenvektor Pertama PCA (Level, Slope, Kurvatur)
	tenorIdxX := make([]float64, len(kalman.StandardTenors))
	for i := range tenorIdxX {
		tenorIdxX[i] = float64(i + 1)
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004a_pca-eigenvectors.svg"), svg.ChartOptions{
		Title:         "Tiga Eigenvektor Pertama PCA Perubahan Yield (Level, Slope, Kurvatur)",
		XLabel:        "Indeks Tenor (1=0.25y, 2=0.5y, 3=1y, 4=2y, 5=3y, 6=5y, 7=7y, 8=10y, 9=20y, 10=30y)",
		YLabel:        "Pemuatan Eigenvektor Satuan",
		ZeroLine:      true,
		LegendTopLeft: true,
	}, []svg.Series{
		{Label: "PC1 Level (Teoretis, 0 flip)", X: tenorIdxX, Y: pcaRes.TheoEigenvectors[0], Color: svg.Palette[0]},
		{Label: "PC1 Level (Empiris R=100)", X: tenorIdxX, Y: pcaRes.EmpMeanTop3Vecs[0], Color: svg.Palette[0], Dashed: true},
		{Label: "PC2 Slope (Teoretis, 1 flip)", X: tenorIdxX, Y: pcaRes.TheoEigenvectors[1], Color: svg.Palette[1]},
		{Label: "PC2 Slope (Empiris R=100)", X: tenorIdxX, Y: pcaRes.EmpMeanTop3Vecs[1], Color: svg.Palette[1], Dashed: true},
		{Label: "PC3 Kurvatur (Teoretis, 2 flip)", X: tenorIdxX, Y: pcaRes.TheoEigenvectors[2], Color: svg.Palette[2]},
		{Label: "PC3 Kurvatur (Empiris R=100)", X: tenorIdxX, Y: pcaRes.EmpMeanTop3Vecs[2], Color: svg.Palette[2], Dashed: true},
	}); err != nil {
		panic(err)
	}

	// =========================================================================
	// SUSUN DAN TULIS OUTPUT JSON (17 DIGIT)
	// =========================================================================
	elapsedTotal := time.Since(tTotalStart).Seconds()

	mleSingleObjs := make([]jsonout.Obj, 11)
	for i, name := range kalman.ParamNames {
		mleSingleObjs[i] = jsonout.Obj{
			{Key: "name", Val: name},
			{Key: "init", Val: singleMLE.InitParams[i]},
			{Key: "true", Val: singleMLE.TrueParams[i]},
			{Key: "est", Val: singleMLE.EstParams[i]},
			{Key: "err", Val: singleMLE.EstParams[i] - singleMLE.TrueParams[i]},
			{Key: "hessian_se", Val: singleMLE.HessianSE[i]},
			{Key: "z_score_vs_true", Val: math.Abs(singleMLE.EstParams[i]-singleMLE.TrueParams[i]) / math.Max(singleMLE.HessianSE[i], 1e-30)},
		}
	}

	repParamObjs := make([]jsonout.Obj, 11)
	for i := 0; i < 11; i++ {
		pr := repStudy.Params[i]
		repParamObjs[i] = jsonout.Obj{
			{Key: "name", Val: pr.Name},
			{Key: "true_value", Val: pr.TrueValue},
			{Key: "mean_estimate", Val: pr.MeanEstimate},
			{Key: "bias", Val: pr.Bias},
			{Key: "rep_se", Val: pr.RepSE},
			{Key: "abs_bias_in_rep_se", Val: pr.BiasInRepSE},
			{Key: "empirical_sd", Val: pr.EmpiricalSD},
			{Key: "rmse", Val: pr.RMSE},
			{Key: "mean_hessian_se", Val: pr.MeanHessianSE},
			{Key: "coverage_95", Val: pr.Coverage95},
			{Key: "covered_count", Val: pr.CoveredCount},
			{Key: "well_identified_group", Val: wellIdentifiedIdx[i]},
		}
	}

	pcaEigObjs := make([]jsonout.Obj, len(kalman.StandardTenors))
	for i := range kalman.StandardTenors {
		pcaEigObjs[i] = jsonout.Obj{
			{Key: "component", Val: i + 1},
			{Key: "theo_eigenvalue", Val: pcaRes.TheoEigenvalues[i]},
			{Key: "approx_innovation_eigenvalue", Val: pcaRes.ApproxEigenvalues[i]},
			{Key: "emp_mean_eigenvalue", Val: pcaRes.EmpMeanEigenvalues[i]},
			{Key: "emp_single_sample_sd", Val: pcaRes.EmpSDEigenvalues[i]},
			{Key: "emp_mean_se", Val: pcaRes.EmpMeanSE[i]},
			{Key: "z_score_mean_vs_theo", Val: pcaRes.ZScoreSampleMean[i]},
			{Key: "single_run_eigenvalue", Val: pcaRes.SingleRunEigenvals[i]},
			{Key: "single_run_z_score", Val: pcaRes.SingleRunZScores[i]},
			{Key: "theo_var_explained", Val: pcaRes.TheoVarExplained[i]},
			{Key: "emp_var_explained", Val: pcaRes.EmpVarExplained[i]},
		}
	}

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day004a"},
			{Key: "date", Val: "2026-11-22"},
			{Key: "language", Val: "go"},
			{Key: "main_seed", Val: rng.MainSeed},
			{Key: "n_workers", Val: nWorkers},
			{Key: "n_chunks", Val: nChunks},
			{Key: "worker_invariance_bit_identical", Val: workerAllIdentical},
			{Key: "elapsed_seconds_part1", Val: elapsedP1},
			{Key: "elapsed_seconds_part2", Val: elapsedP2},
			{Key: "elapsed_seconds_part3", Val: elapsedP3},
			{Key: "elapsed_seconds_part4", Val: elapsedP4},
			{Key: "elapsed_seconds_total", Val: elapsedTotal},
		}},
		{Key: "part1_ns_bootstrap_risk", Val: jsonout.Obj{
			{Key: "identity_checks", Val: idObjs},
			{Key: "max_identity_abs_err", Val: maxIdentityErr},
			{Key: "limits", Val: jsonout.Obj{
				{Key: "ns_y_tau_1e_12", Val: nsLim0},
				{Key: "ns_limit_l_plus_s", Val: nsP.L + nsP.S},
				{Key: "ns_y_tau_1e9", Val: nsLimInf},
				{Key: "ns_limit_l", Val: nsP.L},
				{Key: "nss_y_tau_1e_12", Val: nssLim0},
				{Key: "nss_limit_l_plus_s", Val: nssP.L + nssP.S},
				{Key: "nss_y_tau_1e9", Val: nssLimInf},
				{Key: "nss_limit_l", Val: nssP.L},
			}},
			{Key: "small_tau_stability", Val: smallObjs},
			{Key: "bootstrap", Val: bootSummary},
			{Key: "max_bootstrap_err", Val: maxBootErr},
			{Key: "noisy_interp_r200", Val: jsonout.Obj{
				{Key: "n_rep", Val: interpRes.NumReps},
				{Key: "sigma_par_bps", Val: interpRes.SigmaPar * 10_000.0},
				{Key: "n_knots", Val: interpRes.NumKnots},
				{Key: "n_eval_grid", Val: interpRes.NumEvalGrid},
				{Key: "log_linear_p", Val: jsonout.Obj{
					{Key: "neg_forward_mean_count", Val: interpRes.LogLinearP.NegForwardMeanCount},
					{Key: "neg_forward_max_count", Val: interpRes.LogLinearP.NegForwardMaxCount},
					{Key: "neg_forward_frac_reps", Val: interpRes.LogLinearP.NegForwardFracReps},
					{Key: "forward_tv_mean", Val: interpRes.LogLinearP.ForwardTVMean},
					{Key: "forward_tv_sd", Val: interpRes.LogLinearP.ForwardTVSD},
					{Key: "forward_jump_rms_mean", Val: interpRes.LogLinearP.ForwardJumpRMSMean},
					{Key: "forward_d2_rms_mean", Val: interpRes.LogLinearP.ForwardD2RMSMean},
					{Key: "forward_rmse_mean", Val: interpRes.LogLinearP.ForwardRMSEMean},
					{Key: "yield_rmse_mean", Val: interpRes.LogLinearP.YieldRMSEMean},
					{Key: "yield_rmse_sd", Val: interpRes.LogLinearP.YieldRMSESD},
				}},
				{Key: "natural_cubic_spline_y", Val: jsonout.Obj{
					{Key: "neg_forward_mean_count", Val: interpRes.CubicSpline.NegForwardMeanCount},
					{Key: "neg_forward_max_count", Val: interpRes.CubicSpline.NegForwardMaxCount},
					{Key: "neg_forward_frac_reps", Val: interpRes.CubicSpline.NegForwardFracReps},
					{Key: "forward_tv_mean", Val: interpRes.CubicSpline.ForwardTVMean},
					{Key: "forward_tv_sd", Val: interpRes.CubicSpline.ForwardTVSD},
					{Key: "forward_jump_rms_mean", Val: interpRes.CubicSpline.ForwardJumpRMSMean},
					{Key: "forward_d2_rms_mean", Val: interpRes.CubicSpline.ForwardD2RMSMean},
					{Key: "forward_rmse_mean", Val: interpRes.CubicSpline.ForwardRMSEMean},
					{Key: "yield_rmse_mean", Val: interpRes.CubicSpline.YieldRMSEMean},
					{Key: "yield_rmse_sd", Val: interpRes.CubicSpline.YieldRMSESD},
				}},
			}},
			{Key: "bond_10y_4pct_risk", Val: jsonout.Obj{
				{Key: "price", Val: bondRisk.Price},
				{Key: "ytm", Val: bondRisk.YTM},
				{Key: "mac_duration_ytm", Val: bondRisk.MacDurationYTM},
				{Key: "mod_duration_ytm_closed", Val: bondRisk.ModDurationYTMClosed},
				{Key: "mod_duration_ytm_fd", Val: bondRisk.ModDurationYTMFD},
				{Key: "mod_duration_ytm_rel_err", Val: bondRisk.ModDurationYTMRelErr},
				{Key: "convexity_ytm_closed", Val: bondRisk.ConvexityYTMClosed},
				{Key: "convexity_ytm_fd", Val: bondRisk.ConvexityYTMFD},
				{Key: "convexity_ytm_rel_err", Val: bondRisk.ConvexityYTMRelErr},
				{Key: "eff_duration_zero_closed", Val: bondRisk.EffDurationZeroClosed},
				{Key: "eff_duration_zero_fd", Val: bondRisk.EffDurationZeroFD},
				{Key: "eff_duration_zero_rel_err", Val: bondRisk.EffDurationZeroRelErr},
				{Key: "eff_convexity_zero_closed", Val: bondRisk.EffConvexityZeroClosed},
				{Key: "eff_convexity_zero_fd", Val: bondRisk.EffConvexityZeroFD},
				{Key: "eff_convexity_zero_rel_err", Val: bondRisk.EffConvexityZeroRelErr},
			}},
			{Key: "key_rate_durations", Val: jsonout.Obj{
				{Key: "nodes", Val: krdReport.Nodes},
				{Key: "rows", Val: krdObjs},
				{Key: "sum_krd_fd", Val: krdReport.SumKRDFD},
				{Key: "eff_parallel_fd", Val: krdReport.EffParallelFD},
				{Key: "sum_vs_parallel_fd_err", Val: krdReport.SumVsParFDErr},
				{Key: "sum_vs_parallel_closed_err", Val: krdReport.SumVsParClErr},
			}},
			{Key: "swap_10y_pv01", Val: jsonout.Obj{
				{Key: "par_swap_rate", Val: swapPV01.ParSwapRate},
				{Key: "annuity", Val: swapPV01.Annuity},
				{Key: "annuity_pv01_per_1", Val: swapPV01.AnnuityPV01Per1},
				{Key: "annuity_pv01_per_1m", Val: swapPV01.AnnuityPV01Per1M},
				{Key: "payer_curve_dv01_exact", Val: swapPV01.PayerCurveDV01Ex},
				{Key: "payer_curve_dv01_fd", Val: swapPV01.PayerCurveDV01FD},
				{Key: "payer_curve_dv01_1m", Val: swapPV01.PayerCurveDV011M},
			}},
			{Key: "two_bond_krd_hedge", Val: jsonout.Obj{
				{Key: "target", Val: hedgeRep.TargetName},
				{Key: "hedge_1", Val: hedgeRep.Hedge1Name},
				{Key: "hedge_2", Val: hedgeRep.Hedge2Name},
				{Key: "neutralized_nodes", Val: []float64{hedgeRep.TargetNode1, hedgeRep.TargetNode2}},
				{Key: "weight_hedge_1", Val: hedgeRep.WeightHedge1},
				{Key: "weight_hedge_2", Val: hedgeRep.WeightHedge2},
				{Key: "residual_node_5y", Val: hedgeRep.ResidualNode1},
				{Key: "residual_node_10y", Val: hedgeRep.ResidualNode2},
				{Key: "target_total_par_dollar_dur", Val: hedgeRep.TargetTotalParD},
				{Key: "hedged_residual_total_par_dollar_dur", Val: hedgeRep.ResidualTotalParD},
			}},
		}},
		{Key: "part2_afns_riccati", Val: jsonout.Obj{
			{Key: "max_loadings_vs_ns_err", Val: maxLoadErr},
			{Key: "max_riccati_b_err", Val: maxRicBErr},
			{Key: "max_riccati_a_err", Val: maxRicAErr},
			{Key: "max_riccati_yield_adj_err", Val: maxRicAdjErr},
			{Key: "yield_adj_all_negative_and_monotone", Val: adjNegMono},
			{Key: "adj_tau10", Val: jsonout.Obj{
				{Key: "level", Val: adj10[0]},
				{Key: "slope", Val: adj10[1]},
				{Key: "curvature", Val: adj10[2]},
				{Key: "total", Val: adj10[3]},
				{Key: "total_bps", Val: adj10[3] * 10_000.0},
			}},
			{Key: "adj_tau30", Val: jsonout.Obj{
				{Key: "level", Val: adj30[0]},
				{Key: "slope", Val: adj30[1]},
				{Key: "curvature", Val: adj30[2]},
				{Key: "total", Val: adj30[3]},
				{Key: "total_bps", Val: adj30[3] * 10_000.0},
			}},
			{Key: "riccati_table", Val: ricObjs},
			{Key: "p_dynamics_transition", Val: jsonout.Obj{
				{Key: "lower_tri_phi_vs_matexp_err", Val: maxPhiMatExpErr},
				{Key: "q_vanloan_vs_gl_err", Val: maxQVanLoanGLErr},
				{Key: "phi_weekly", Val: mat3ToSlice(afnsP.TransitionMatrixClosed(kalman.WeeklyDt))},
				{Key: "q_weekly", Val: mat3ToSlice(afnsP.TransitionCovClosed(kalman.WeeklyDt))},
				{Key: "q_stationary", Val: mat3ToSlice(afnsP.StationaryCov())},
				{Key: "mc_transition_1y_200k", Val: jsonout.Obj{
					{Key: "x0", Val: []float64{transMC.X0[0], transMC.X0[1], transMC.X0[2]}},
					{Key: "dt", Val: transMC.Dt},
					{Key: "n_paths", Val: transMC.NumPaths},
					{Key: "theo_mean", Val: []float64{transMC.TheoMean[0], transMC.TheoMean[1], transMC.TheoMean[2]}},
					{Key: "mc_mean", Val: []float64{transMC.MCMean[0], transMC.MCMean[1], transMC.MCMean[2]}},
					{Key: "mean_z_scores", Val: []float64{transMC.MeanZScores[0], transMC.MeanZScores[1], transMC.MeanZScores[2]}},
					{Key: "theo_cov", Val: mat3ToSlice(transMC.TheoCov)},
					{Key: "mc_cov", Val: mat3ToSlice(transMC.MCCov)},
					{Key: "cov_z_scores", Val: mat3ToSlice(transMC.CovZScores)},
					{Key: "max_z_score", Val: transMC.MaxZScore},
				}},
			}},
			{Key: "market_price_of_risk", Val: jsonout.Obj{
				{Key: "lambda0", Val: []float64{lam0[0], lam0[1], lam0[2]}},
				{Key: "lambda1", Val: mat3ToSlice(lam1)},
			}},
			{Key: "term_premium_10y", Val: jsonout.Obj{
				{Key: "closed_form_tp", Val: tpCl},
				{Key: "closed_form_tp_bps", Val: tpCl * 10_000.0},
				{Key: "long_run_yield_10y", Val: lrY10},
				{Key: "long_run_short_rate", Val: lrR},
				{Key: "adj_component", Val: tpAdjC},
				{Key: "slope_component", Val: tpSlopeC},
				{Key: "curvature_component", Val: tpCurvC},
				{Key: "mc_paths", Val: tpMC.NumPaths},
				{Key: "mc_mean_yield_10y", Val: tpMC.MCMeanYield10y},
				{Key: "mc_mean_avg_short_rate", Val: tpMC.MCMeanAvgShort},
				{Key: "mc_tp", Val: tpMC.MCTermPremium},
				{Key: "mc_tp_bps", Val: tpMC.MCTermPremium * 10_000.0},
				{Key: "mc_tp_se_bps", Val: tpMC.MCTermPremiumSE * 10_000.0},
				{Key: "z_score", Val: tpMC.ZScore},
			}},
		}},
		{Key: "part3_kalman_mle", Val: jsonout.Obj{
			{Key: "t5_joint_gaussian_anchor", Val: jsonout.Obj{
				{Key: "t", Val: anchorT5.T},
				{Key: "tenors", Val: []float64{1, 5, 10}},
				{Key: "kalman_loglik", Val: anchorT5.KalmanLogLik},
				{Key: "joint_loglik", Val: anchorT5.JointLogLik},
				{Key: "loglik_abs_diff", Val: anchorT5.LogLikAbsDiff},
				{Key: "max_smooth_mean_abs_diff", Val: anchorT5.MaxSmoothMeanAbsDiff},
				{Key: "max_smooth_cov_abs_diff", Val: anchorT5.MaxSmoothCovAbsDiff},
			}},
			{Key: "single_mle_t1000", Val: jsonout.Obj{
				{Key: "t_steps", Val: 1000},
				{Key: "tenors", Val: kalman.StandardTenors},
				{Key: "init_loglik", Val: singleMLE.InitLogLik},
				{Key: "true_loglik", Val: singleMLE.TrueLogLik},
				{Key: "est_loglik", Val: singleMLE.EstLogLik},
				{Key: "grad_rel_err_init", Val: singleMLE.GradRelErrInit},
				{Key: "grad_rel_err_true", Val: singleMLE.GradRelErrTrue},
				{Key: "grad_abs_diff_opt", Val: singleMLE.GradAbsDiffOpt},
				{Key: "grad_norm_opt", Val: singleMLE.GradNormOpt},
				{Key: "hessian_min_eig", Val: singleMLE.HessianMinEig},
				{Key: "hessian_regularized", Val: singleMLE.HessianRegularized},
				{Key: "iterations", Val: singleMLE.Iterations},
				{Key: "converged", Val: singleMLE.Converged},
				{Key: "hit_bound", Val: singleMLE.HitBound},
				{Key: "parameters", Val: mleSingleObjs},
			}},
			{Key: "replication_study_r50", Val: jsonout.Obj{
				{Key: "n_rep", Val: repStudy.NumReps},
				{Key: "t_steps", Val: repStudy.TSteps},
				{Key: "converged_count", Val: repStudy.ConvergedCount},
				{Key: "ll_ge_true_count", Val: repStudy.LLGeTrueCount},
				{Key: "bound_hit_count", Val: repStudy.BoundHitCount},
				{Key: "pre_specified_well_identified_params", Val: []string{"lambda", "sigma1", "sigma2", "sigma3", "sigma_eps"}},
				{Key: "well_identified_bias_lt_2p0_se", Val: wellIdBiasOK},
				{Key: "well_identified_cov95_in_88_99", Val: wellIdCovOK},
				{Key: "all_11_params_bias_lt_2p0_se", Val: allParamsBiasOK},
				{Key: "all_11_params_cov95_in_88_99", Val: allParamsCovOK},
				{Key: "parameters", Val: repParamObjs},
			}},
			{Key: "dns_misspecification_t1000", Val: jsonout.Obj{
				{Key: "dns_loglik", Val: dnsMisspec.DNSLogLik},
				{Key: "afns_loglik", Val: dnsMisspec.AFNSLogLik},
				{Key: "loglik_drop_afns_vs_dns", Val: dnsMisspec.LogLikDrop},
				{Key: "lambda_est_dns", Val: dnsMisspec.DNSEst[0]},
				{Key: "lambda_est_afns", Val: dnsMisspec.AFNSEst[0]},
				{Key: "lambda_bias_dns", Val: dnsMisspec.LambdaBiasDNS},
				{Key: "lambda_bias_afns", Val: dnsMisspec.LambdaBiasAFNS},
				{Key: "theta1_est_dns", Val: dnsMisspec.DNSEst[4]},
				{Key: "theta1_est_afns", Val: dnsMisspec.AFNSEst[4]},
				{Key: "theta1_bias_dns", Val: dnsMisspec.Theta1BiasDNS},
				{Key: "theta1_bias_afns", Val: dnsMisspec.Theta1BiasAFNS},
			}},
		}},
		{Key: "part4_pca", Val: jsonout.Obj{
			{Key: "tenors", Val: kalman.StandardTenors},
			{Key: "t_steps", Val: pcaRes.TSteps},
			{Key: "n_rep", Val: pcaRes.NumReps},
			{Key: "theo_top3_cum_var", Val: pcaRes.TheoTop3CumVar},
			{Key: "emp_top3_cum_var", Val: pcaRes.EmpTop3CumVar},
			{Key: "theo_sign_changes_top3", Val: []int{pcaRes.TheoSignChanges[0], pcaRes.TheoSignChanges[1], pcaRes.TheoSignChanges[2]}},
			{Key: "emp_sign_changes_top3", Val: []int{pcaRes.EmpSignChanges[0], pcaRes.EmpSignChanges[1], pcaRes.EmpSignChanges[2]}},
			{Key: "emp_mean_top3_cosine", Val: []float64{pcaRes.EmpMeanTop3Cos[0], pcaRes.EmpMeanTop3Cos[1], pcaRes.EmpMeanTop3Cos[2]}},
			{Key: "emp_min_top3_cosine", Val: []float64{pcaRes.EmpMinTop3Cos[0], pcaRes.EmpMinTop3Cos[1], pcaRes.EmpMinTop3Cos[2]}},
			{Key: "single_run_top3_cosine", Val: []float64{pcaRes.SingleRunTop3Cos[0], pcaRes.SingleRunTop3Cos[1], pcaRes.SingleRunTop3Cos[2]}},
			{Key: "theo_top3_eigenvectors", Val: pcaRes.TheoEigenvectors[:3]},
			{Key: "emp_mean_top3_eigenvectors", Val: pcaRes.EmpMeanTop3Vecs},
			{Key: "noise_bulk_pc4_10_theo_mean", Val: pcaRes.NoiseBulkTheoMean},
			{Key: "noise_bulk_pc4_10_emp_mean", Val: pcaRes.NoiseBulkEmpMean},
			{Key: "noise_bulk_pc4_10_emp_mean_se", Val: pcaRes.NoiseBulkMeanSE},
			{Key: "noise_bulk_pc4_10_z_score", Val: pcaRes.NoiseBulkZScore},
			{Key: "eigenvalues", Val: pcaEigObjs},
		}},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_ns_integral_and_limits", Val: crit1},
			{Key: "c2_bootstrap_par_recovery", Val: crit2},
			{Key: "c3_bond_duration_convexity_krd_hedge", Val: crit3},
			{Key: "c4_afns_loadings_and_riccati_rk4", Val: crit4},
			{Key: "c5_p_transition_and_term_premium", Val: crit5},
			{Key: "c6_kalman_joint_t5_and_rts_smoother", Val: crit6},
			{Key: "c7_mle_score_and_loglik", Val: crit7},
			{Key: "c8a_mle_replication_well_identified", Val: crit8WellId},
			{Key: "c8b_mle_replication_all11_params", Val: crit8All11},
			{Key: "c8_status", Val: c8Status},
			{Key: "c9a_pca_top3_eigenvalues_within_3se", Val: pcaAcc.Top3EigenvaluesWithin3SE},
			{Key: "c9b_pca_noise_bulk_mean_within_3se", Val: pcaAcc.NoiseBulkMeanWithin3SE},
			{Key: "c9c_pca_all10_sorted_eigenvalues_within_3se", Val: pcaAcc.DiagnosticAll10SortedWithin3SE},
			{Key: "c9d_pca_top3_cosine_gt_0p99", Val: pcaAcc.Top3CosineGt0p99},
			{Key: "c9e_pca_sign_patterns_012", Val: pcaAcc.SignPatterns012},
			{Key: "c9_pca_eigenvalues_cosine_and_signs", Val: pcaAcc.AggregatePassed},
			{Key: "c9_status", Val: pcaAcc.Status},
			{Key: "worker_invariance_bit_identical", Val: workerAllIdentical},
			{Key: "all_passed", Val: allPassed},
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		panic(err)
	}
	fmt.Printf("day004a selesai dalam %.2fs (P1=%.2fs P2=%.2fs P3=%.2fs P4=%.2fs) | workers=%d | all_passed=%v\n",
		elapsedTotal, elapsedP1, elapsedP2, elapsedP3, elapsedP4, nWorkers, allPassed)
}
