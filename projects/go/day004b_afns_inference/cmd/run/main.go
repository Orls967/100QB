package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/afns"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/jsonout"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/kalman"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/parallel"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/rng"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/svg"
)

func findRepoRoot() string {
	wd, err := os.Getwd()
	if err != nil {
		return "."
	}
	dir := wd
	for i := 0; i < 8; i++ {
		if _, err := os.Stat(filepath.Join(dir, "outputs")); err == nil {
			if _, err2 := os.Stat(filepath.Join(dir, "projects", "go", "day004b_afns_inference")); err2 == nil {
				return dir
			}
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			break
		}
		dir = parent
	}
	return wd
}

func configSpecToJSON(s kalman.ConfigSpec) jsonout.Obj {
	return jsonout.Obj{
		{Key: "name", Val: s.Name},
		{Key: "group", Val: s.Group},
		{Key: "frequency", Val: s.FrequencyLabel},
		{Key: "horizon_years", Val: s.HorizonYears},
		{Key: "dt_years", Val: s.Dt},
		{Key: "t_steps", Val: s.TSteps},
		{Key: "n_rep", Val: s.NumReps},
		{Key: "salt", Val: s.Salt},
		{Key: "shared_prefix_group", Val: s.SharedPrefixGroup},
		{Key: "reused_from", Val: s.ReusedFrom},
	}
}

func paramSummaryToJSON(p kalman.ParamSummary) jsonout.Obj {
	return jsonout.Obj{
		{Key: "name", Val: p.Name},
		{Key: "true_value", Val: p.TrueValue},
		{Key: "mean_estimate", Val: p.MeanEstimate},
		{Key: "bias", Val: p.Bias},
		{Key: "bias_mcse", Val: p.BiasMCSE},
		{Key: "abs_bias_in_mcse", Val: p.BiasInMCSE},
		{Key: "empirical_sd", Val: p.EmpiricalSD},
		{Key: "rmse", Val: p.RMSE},
		{Key: "mean_hessian_se", Val: p.MeanHessianSE},
		{Key: "median_hessian_se", Val: p.MedianHessianSE},
		{Key: "mean_hessian_se_to_emp_sd_ratio", Val: p.MeanHessianSEToEmpSDRatio},
		{Key: "median_hessian_se_to_emp_sd_ratio", Val: p.MedianHessianSEToEmpSDRatio},
		{Key: "coverage_95", Val: p.Coverage95},
		{Key: "coverage_wilson_95_lower", Val: p.CoverageWilsonLower},
		{Key: "coverage_wilson_95_upper", Val: p.CoverageWilsonUpper},
		{Key: "covered_count", Val: p.CoveredCount},
		{Key: "total_reps", Val: p.TotalReps},
		{Key: "bound_hit_count", Val: p.BoundHitCount},
		{Key: "bound_hit_pct", Val: p.BoundHitPct},
	}
}

func mleResultToJSON(repIdx int, r kalman.MLEResult) jsonout.Obj {
	hb := make([]bool, 11)
	for i := 0; i < 11; i++ {
		hb[i] = r.HitBoundByParam[i]
	}
	return jsonout.Obj{
		{Key: "rep_idx", Val: repIdx},
		{Key: "init_params", Val: r.InitParams[:]},
		{Key: "est_params", Val: r.EstParams[:]},
		{Key: "init_loglik", Val: r.InitLogLik},
		{Key: "true_loglik", Val: r.TrueLogLik},
		{Key: "est_loglik", Val: r.EstLogLik},
		{Key: "grad_rel_err_init", Val: r.GradRelErrInit},
		{Key: "grad_rel_err_true", Val: r.GradRelErrTrue},
		{Key: "grad_abs_diff_opt", Val: r.GradAbsDiffOpt},
		{Key: "grad_norm_opt", Val: r.GradNormOpt},
		{Key: "scaled_grad_norm_opt", Val: r.ScaledGradNormOpt},
		{Key: "hessian_se", Val: r.HessianSE[:]},
		{Key: "hessian_min_eig", Val: r.HessianMinEig},
		{Key: "hessian_regularized", Val: r.HessianRegularized},
		{Key: "iterations", Val: r.Iterations},
		{Key: "converged", Val: r.Converged},
		{Key: "hit_bound", Val: r.HitBound},
		{Key: "hit_bound_by_param", Val: hb},
	}
}

func studyResultToJSON(s kalman.StudyConfigResult) jsonout.Obj {
	pObjs := make([]jsonout.Obj, 11)
	for i := 0; i < 11; i++ {
		pObjs[i] = paramSummaryToJSON(s.Params[i])
	}
	repObjs := make([]jsonout.Obj, len(s.Replications))
	for rIdx, rep := range s.Replications {
		repObjs[rIdx] = mleResultToJSON(rIdx, rep)
	}
	return jsonout.Obj{
		{Key: "config", Val: configSpecToJSON(s.Spec)},
		{Key: "workers_used", Val: s.NumWorkers},
		{Key: "elapsed_seconds", Val: s.ElapsedSeconds},
		{Key: "converged_count", Val: s.ConvergedCount},
		{Key: "converged_pct", Val: s.ConvergedPct},
		{Key: "ll_ge_true_count", Val: s.LLGeTrueCount},
		{Key: "ll_ge_true_pct", Val: s.LLGeTruePct},
		{Key: "bound_hit_count", Val: s.BoundHitCount},
		{Key: "bound_hit_pct", Val: s.BoundHitPct},
		{Key: "hessian_regularized_count", Val: s.HessianRegCount},
		{Key: "parameters", Val: pObjs},
		{Key: "replications", Val: repObjs},
	}
}

func trendToJSON(tr kalman.InvHorizonTrendSummary) jsonout.Obj {
	return jsonout.Obj{
		{Key: "param", Val: tr.ParamName},
		{Key: "horizons_years", Val: tr.Horizons},
		{Key: "inv_horizons", Val: tr.InvHorizons},
		{Key: "biases", Val: tr.Biases},
		{Key: "bias_mcse", Val: tr.BiasMCSEs},
		{Key: "rmse", Val: tr.RMSEs},
		{Key: "bias_times_h", Val: tr.BiasTimesH},
		{Key: "bias_times_h_mcse", Val: tr.BiasTimesHMCSE},
		{Key: "ols_intercept", Val: tr.Trend.OLSIntercept},
		{Key: "ols_intercept_prefix_mcse", Val: tr.OLSInterceptPrefixMCSE},
		{Key: "ols_slope_vs_inv_h", Val: tr.Trend.OLSSlope},
		{Key: "ols_slope_prefix_mcse", Val: tr.OLSSlopePrefixMCSE},
		{Key: "ols_slope_indep_mcse", Val: tr.Trend.OLSSlopeMCSE},
		{Key: "ols_slope_residual_se", Val: tr.Trend.OLSSlopeResidualSE},
		{Key: "wls_intercept", Val: tr.Trend.WLSIntercept},
		{Key: "wls_intercept_prefix_mcse", Val: tr.WLSInterceptPrefixMCSE},
		{Key: "wls_slope_vs_inv_h", Val: tr.Trend.WLSSlope},
		{Key: "wls_slope_prefix_mcse", Val: tr.WLSSlopePrefixMCSE},
		{Key: "wls_slope_indep_mcse", Val: tr.Trend.WLSSlopeMCSE},
	}
}

func allFiniteCheck(s kalman.StudyConfigResult, expectedReps int) bool {
	if len(s.Replications) != expectedReps {
		return false
	}
	check := func(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }
	for _, p := range s.Params {
		if p.TotalReps != expectedReps {
			return false
		}
		for _, v := range []float64{
			p.MeanEstimate, p.Bias, p.BiasMCSE, p.EmpiricalSD, p.RMSE,
			p.MeanHessianSE, p.MedianHessianSE, p.Coverage95,
			p.CoverageWilsonLower, p.CoverageWilsonUpper,
		} {
			if !check(v) {
				return false
			}
		}
	}
	return true
}

func buildFigures(
	figDir string,
	kpTrends [3]kalman.InvHorizonTrendSummary,
	infillM, infillW, infillD kalman.StudyConfigResult,
	longSpans []kalman.StudyConfigResult,
	baseline kalman.StudyConfigResult,
	allDistinct []kalman.StudyConfigResult,
) error {
	if err := os.MkdirAll(figDir, 0o755); err != nil {
		return err
	}

	// Grafik 1: Bias K^P_ii vs 1/H beserta ketidakpastian +-1.96 MCSE dan garis pembanding 4/H
	invHs := kpTrends[0].InvHorizons
	heuristic4OverH := make([]float64, len(invHs))
	for i, invH := range invHs {
		heuristic4OverH[i] = 4.0 * invH
	}
	var series1 []svg.Series
	for k := 0; k < 3; k++ {
		tr := kpTrends[k]
		lo := make([]float64, len(tr.Biases))
		hi := make([]float64, len(tr.Biases))
		for i := range tr.Biases {
			lo[i] = tr.Biases[i] - tr.BiasMCSEs[i]
			hi[i] = tr.Biases[i] + tr.BiasMCSEs[i]
		}
		series1 = append(series1, svg.Series{
			Label: tr.ParamName + " bias (±1 MCSE)",
			X:     tr.InvHorizons,
			Y:     tr.Biases,
			Lower: lo,
			Upper: hi,
			Color: svg.Palette[k],
		})
	}
	series1 = append(series1, svg.Series{
		Label:  "Heuristik 4/H (Tang-Chen / Kendall)",
		X:      invHs,
		Y:      heuristic4OverH,
		Color:  "#333333",
		Dashed: true,
	})
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004b_kp-bias-vs-inv-h.svg"), svg.ChartOptions{
		Title:         "Bias Finite-Sample K^P_ii vs 1/H (Weekly Shared-Prefix, R=50)",
		XLabel:        "1 / Horizon Kalender H [tahun^-1]  (H = 80, 40, 20, 10 thn)",
		YLabel:        "Bias MLE E[K^P_est - K^P_true] [tahun^-1]",
		ZeroLine:      true,
		LegendTopLeft: true,
	}, series1); err != nil {
		return err
	}

	// Grafik 2: Relative RMSE (RMSE / |true|) terhadap frekuensi observasi pada H = 20 tahun (log-log)
	freqX := []float64{12.0, 52.0, 252.0}
	keyIdxs := []int{1, 2, 4, 0, 7, 10} // kp11, kp22, theta1, lambda, sigma1, sigma_eps
	var series2 []svg.Series
	for sIdx, pIdx := range keyIdxs {
		pm, pw, pd := infillM.Params[pIdx], infillW.Params[pIdx], infillD.Params[pIdx]
		denom := math.Abs(pm.TrueValue)
		series2 = append(series2, svg.Series{
			Label: pm.Name + " (rel RMSE)",
			X:     freqX,
			Y:     []float64{pm.RMSE / denom, pw.RMSE / denom, pd.RMSE / denom},
			Color: svg.Palette[sIdx%len(svg.Palette)],
		})
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004b_infill-vs-longspan-rmse.svg"), svg.ChartOptions{
		Title:        "Infill Asymptotics pada Horizon Tetap H=20 Tahun (Relative RMSE vs Frekuensi)",
		XLabel:       "Frekuensi Observasi per Tahun [12=Bulanan, 52=Mingguan, 252=Harian]",
		YLabel:       "RMSE Relatif (RMSE / |Nilai Benar|)",
		LogX:         true,
		LogY:         true,
		LegendBottom: true,
	}, series2); err != nil {
		return err
	}

	// Grafik 3: Cakupan Interval Wald 95% + Interval Binomial Wilson 95% untuk ke-11 parameter
	cats := kalman.ParamNames[:]
	covConfigs := []struct {
		label string
		res   kalman.StudyConfigResult
		col   string
	}{
		{"Baseline 1000w (19.2y)", baseline, svg.Palette[0]},
		{"Long-span 10y (520w)", longSpans[0], svg.Palette[1]},
		{"Long-span 80y (4160w)", longSpans[len(longSpans)-1], svg.Palette[2]},
		{"Infill 20y Daily (5040d)", infillD, svg.Palette[3]},
	}
	var series3 []svg.CategoryPointSeries
	for _, cc := range covConfigs {
		vals := make([]float64, 11)
		lo := make([]float64, 11)
		hi := make([]float64, 11)
		for i := 0; i < 11; i++ {
			vals[i] = cc.res.Params[i].Coverage95
			lo[i] = cc.res.Params[i].CoverageWilsonLower
			hi[i] = cc.res.Params[i].CoverageWilsonUpper
		}
		series3 = append(series3, svg.CategoryPointSeries{
			Label:  cc.label,
			Values: vals,
			Lower:  lo,
			Upper:  hi,
			Color:  cc.col,
		})
	}
	y0, y1 := 0.0, 1.04
	if err := svg.WriteCategoryPointChart(filepath.Join(figDir, "day004b_wald-coverage-intervals.svg"), svg.ChartOptions{
		Title:        "Cakupan Empiris Interval Wald 95% & Interval Wilson 95% per Parameter",
		XLabel:       "Parameter Model AFNS (11 Parameter)",
		YLabel:       "Proporsi Cakupan 95% (+/- Wilson 95% CI)",
		YMinOverride: &y0,
		YMaxOverride: &y1,
		LegendBottom: true,
		HRefLines: []svg.HRefLine{
			{Value: 0.95, Label: "Nominal 95%", Color: "#111111"},
		},
	}, cats, series3); err != nil {
		return err
	}

	// Grafik 4: Rasio Median Hessian SE / Empirical SD per parameter
	var series4 []svg.BarSeries
	for _, cc := range covConfigs {
		vals := make([]float64, 11)
		for i := 0; i < 11; i++ {
			vals[i] = cc.res.Params[i].MedianHessianSEToEmpSDRatio
		}
		series4 = append(series4, svg.BarSeries{
			Label:  cc.label,
			Values: vals,
			Color:  cc.col,
		})
	}
	if err := svg.WriteGroupedBarChart(filepath.Join(figDir, "day004b_hessian-se-to-emp-sd-ratio.svg"), svg.ChartOptions{
		Title:         "Rasio Median Hessian SE terhadap SD Empiris Antar-Replikasi",
		XLabel:        "Parameter Model AFNS",
		YLabel:        "Median Hessian SE / Empirical SD",
		LegendTopLeft: true,
		HRefLines: []svg.HRefLine{
			{Value: 1.0, Label: "Kalibrasi Ideal = 1.0", Color: "#111111"},
		},
	}, cats, series4); err != nil {
		return err
	}

	// Grafik 5 (diagnostik pelengkap): Tingkat Konvergensi dan Solusi Menyentuh Batas per Konfigurasi
	cfgLabels := make([]string, len(allDistinct))
	convPcts := make([]float64, len(allDistinct))
	hitPcts := make([]float64, len(allDistinct))
	kp11HitPcts := make([]float64, len(allDistinct))
	for i, st := range allDistinct {
		short := strings.TrimPrefix(st.Spec.Name, "longspan_")
		short = strings.TrimPrefix(short, "infill_")
		cfgLabels[i] = short
		convPcts[i] = st.ConvergedPct
		hitPcts[i] = st.BoundHitPct
		kp11HitPcts[i] = st.Params[1].BoundHitPct
	}
	yB0, yB1 := 0.0, 108.0
	if err := svg.WriteGroupedBarChart(filepath.Join(figDir, "day004b_convergence-and-bound-hits.svg"), svg.ChartOptions{
		Title:        "Diagnostik Konvergensi Optimisasi & Solusi Menyentuh Batas (%) per Studi",
		XLabel:       "Konfigurasi Eksperimen",
		YLabel:       "Persentase dari R Replikasi (%)",
		YMinOverride: &yB0,
		YMaxOverride: &yB1,
		LegendBottom: true,
	}, cfgLabels, []svg.BarSeries{
		{Label: "Konvergen (%)", Values: convPcts, Color: svg.Palette[2]},
		{Label: "Batas Aktif Total (%)", Values: hitPcts, Color: svg.Palette[1]},
		{Label: "Batas Aktif kp11 (%)", Values: kp11HitPcts, Color: svg.Palette[4]},
	}); err != nil {
		return err
	}

	return nil
}

func main() {
	root := findRepoRoot()
	quick := flag.Bool("quick", false, "jalankan mode cepat pengembangan ke direktori sementara tanpa menimpa keluaran penuh")
	workersFlag := flag.Int("workers", 0, "jumlah worker goroutine (0 = min(8, NumCPU))")
	repsFlag := flag.Int("r", 0, "override jumlah replikasi (bawaan: 50 untuk penuh, 4 untuk -quick)")
	outFlag := flag.String("out", "", "override berkas keluaran JSON")
	figDirFlag := flag.String("fig-dir", "", "override folder keluaran grafik SVG")
	flag.Parse()

	nWorkers := *workersFlag
	if nWorkers <= 0 {
		nWorkers = parallel.DefaultWorkers()
	}
	nRep := 50
	if *quick {
		nRep = 4
	}
	if *repsFlag > 0 {
		nRep = *repsFlag
	}

	outPath := filepath.Join(root, "outputs", "day004b", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	if *quick && *outFlag == "" {
		outPath = filepath.Join(os.TempDir(), "day004b_quick_results.json")
	}
	if *quick && *figDirFlag == "" {
		figDir = filepath.Join(os.TempDir(), "day004b_quick_figures")
	}
	if *outFlag != "" {
		outPath = *outFlag
	}
	if *figDirFlag != "" {
		figDir = *figDirFlag
	}

	totalT0 := time.Now()
	trueSpec := kalman.DefaultSpec()

	// Bagian 1a: Anchor Eksak T=5, 3 Tenor (Gaussian gabungan vs Kalman Filter & RTS Smoother)
	tAnc0 := time.Now()
	anchorTenors := []float64{0.5, 2.0, 10.0}
	ancRng := rng.New(rng.DeriveSeed(rng.MainSeed, 0xA001, 0))
	_, ancYields := afns.SimulateStationaryPath(trueSpec.AFNS, anchorTenors, 5, kalman.WeeklyDt, trueSpec.SigmaEps, ancRng)
	anchorRes, err := kalman.VerifyJointGaussianAnchor(trueSpec, anchorTenors, kalman.WeeklyDt, ancYields)
	if err != nil {
		fmt.Fprintf(os.Stderr, "VerifyJointGaussianAnchor error: %v\n", err)
		os.Exit(1)
	}
	fastAncLL, err := kalman.FastLogLikelihood(trueSpec, anchorTenors, kalman.WeeklyDt, ancYields)
	if err != nil {
		fmt.Fprintf(os.Stderr, "FastLogLikelihood anchor error: %v\n", err)
		os.Exit(1)
	}
	fastAncDiff := math.Abs(fastAncLL - anchorRes.JointLogLik)

	// Bagian 1b: Verifikasi invariansi worker 1 vs nWorkers pada fixture cepat
	detCfg := kalman.ConfigSpec{
		Name: "worker_determinism_fixture", Group: "verification", FrequencyLabel: "weekly",
		HorizonYears: 100.0 * kalman.WeeklyDt, Dt: kalman.WeeklyDt, TSteps: 100, NumReps: 4, Salt: 0xDEAD,
	}
	det1 := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, detCfg, 1, rng.MainSeed)
	detN := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, detCfg, nWorkers, rng.MainSeed)
	workerDetPassed := true
	for r := 0; r < detCfg.NumReps; r++ {
		if math.Float64bits(det1.Replications[r].EstLogLik) != math.Float64bits(detN.Replications[r].EstLogLik) {
			workerDetPassed = false
		}
		for i := 0; i < 11; i++ {
			if math.Float64bits(det1.Replications[r].EstParams[i]) != math.Float64bits(detN.Replications[r].EstParams[i]) {
				workerDetPassed = false
			}
		}
	}
	elapsedPart1 := time.Since(tAnc0).Seconds()

	// Bagian 2a: Replikasi baseline 1.000 minggu (~19,2308 tahun)
	baselineCfg := kalman.ConfigSpec{
		Name:           "baseline_1000w",
		Group:          "baseline",
		FrequencyLabel: "weekly",
		HorizonYears:   1000.0 * kalman.WeeklyDt,
		Dt:             kalman.WeeklyDt,
		TSteps:         1000,
		NumReps:        nRep,
		Salt:           0x4B01,
	}
	fmt.Printf("[day004b] Menjalankan baseline 1000w (R=%d, workers=%d)...\n", nRep, nWorkers)
	baselineRes := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, baselineCfg, nWorkers, rng.MainSeed)
	baselineCmp := kalman.CompareBaselineToDay004a(baselineRes)

	// Bagian 2b: Long-span mingguan {10, 20, 40, 80} tahun dengan prefiks 80y bersama
	fmt.Printf("[day004b] Menjalankan long-span mingguan {10,20,40,80} tahun (R=%d, workers=%d)...\n", nRep, nWorkers)
	longSpanRes := kalman.RunLongSpanSharedPrefixStudy(trueSpec, kalman.StandardTenors, []int{10, 20, 40, 80}, nRep, nWorkers, rng.MainSeed, 0x4B02)
	var longSpanElapsed float64
	for _, r := range longSpanRes {
		longSpanElapsed += r.ElapsedSeconds
	}

	// Analisis tren bias K^P_ii vs 1/H + semua 11 parameter
	var kpTrends [3]kalman.InvHorizonTrendSummary
	for k := 0; k < 3; k++ {
		kpTrends[k] = kalman.AnalyzeLongSpanTrend(longSpanRes, k+1)
	}
	var allTrends [11]kalman.InvHorizonTrendSummary
	for i := 0; i < 11; i++ {
		allTrends[i] = kalman.AnalyzeLongSpanTrend(longSpanRes, i)
	}

	// Bagian 2c: Infill pada H=20 tahun: bulanan (T=240), mingguan (REUSED dari longspan_20y_weekly), harian bisnis (T=5040)
	fmt.Printf("[day004b] Menjalankan infill 20y bulanan (T=240) & harian bisnis (T=5040) (R=%d, workers=%d)...\n", nRep, nWorkers)
	infillMonthlyCfg := kalman.ConfigSpec{
		Name:           "infill_20y_monthly",
		Group:          "infill",
		FrequencyLabel: "monthly",
		HorizonYears:   20.0,
		Dt:             kalman.MonthlyDt,
		TSteps:         240,
		NumReps:        nRep,
		Salt:           0x4B03,
	}
	infillMonthlyRes := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, infillMonthlyCfg, nWorkers, rng.MainSeed)
	infillWeeklyRes := kalman.ReuseAsInfillConfig(longSpanRes[1]) // index 1 == 20y
	infillDailyCfg := kalman.ConfigSpec{
		Name:           "infill_20y_daily",
		Group:          "infill",
		FrequencyLabel: "daily_business",
		HorizonYears:   20.0,
		Dt:             kalman.DailyBusinessDt,
		TSteps:         5040,
		NumReps:        nRep,
		Salt:           0x4B04,
	}
	infillDailyRes := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, infillDailyCfg, nWorkers, rng.MainSeed)
	infillScaling := kalman.AnalyzeInfillScaling(infillMonthlyRes, infillWeeklyRes, infillDailyRes)

	allDistinct := []kalman.StudyConfigResult{
		baselineRes,
		longSpanRes[0], longSpanRes[1], longSpanRes[2], longSpanRes[3],
		infillMonthlyRes, infillDailyRes,
	}

	// Evaluasi skor gradien terburuk sepanjang seluruh studi
	maxGradErrInit, maxGradErrTrue := 0.0, 0.0
	for _, st := range allDistinct {
		for _, rep := range st.Replications {
			if rep.GradRelErrInit > maxGradErrInit {
				maxGradErrInit = rep.GradRelErrInit
			}
			if rep.GradRelErrTrue > maxGradErrTrue {
				maxGradErrTrue = rep.GradRelErrTrue
			}
		}
	}

	// Render grafik SVG
	tSvg0 := time.Now()
	if err := buildFigures(figDir, kpTrends, infillMonthlyRes, infillWeeklyRes, infillDailyRes, longSpanRes, baselineRes, allDistinct); err != nil {
		fmt.Fprintf(os.Stderr, "buildFigures error: %v\n", err)
		os.Exit(1)
	}
	elapsedFigures := time.Since(tSvg0).Seconds()
	elapsedTotal := time.Since(totalT0).Seconds()

	// Evaluasi 4 kriteria penerimaan pra-spesifikasi (tidak memaksakan hasil riset)
	c1Anchor := anchorRes.LogLikAbsDiff < 1e-9 &&
		anchorRes.MaxSmoothMeanAbsDiff < 1e-9 &&
		anchorRes.MaxSmoothCovAbsDiff < 1e-9 &&
		fastAncDiff < 1e-9 &&
		maxGradErrInit < 1e-5 &&
		maxGradErrTrue < 1e-5

	c2BaselineRepro := baselineCmp.AllWithin3SE
	c3Completeness := allFiniteCheck(baselineRes, nRep) &&
		allFiniteCheck(infillMonthlyRes, nRep) &&
		allFiniteCheck(infillWeeklyRes, nRep) &&
		allFiniteCheck(infillDailyRes, nRep)
	for _, ls := range longSpanRes {
		c3Completeness = c3Completeness && allFiniteCheck(ls, nRep)
	}
	c4Determinism := workerDetPassed && (elapsedTotal < 600.0 || *quick)
	allMethPassed := c1Anchor && c2BaselineRepro && c3Completeness && c4Determinism

	statusStr := func(ok bool) string {
		if ok {
			return "tercapai"
		}
		return "sebagian"
	}

	// Bangun JSON baseline comparison
	cmpRows := make([]jsonout.Obj, 11)
	for i, c := range baselineCmp.Params {
		cmpRows[i] = jsonout.Obj{
			{Key: "name", Val: c.Name},
			{Key: "bias_day004b", Val: c.Bias4b},
			{Key: "mcse_day004b", Val: c.MCSE4b},
			{Key: "bias_day004a", Val: c.Bias4a},
			{Key: "mcse_day004a", Val: c.MCSE4a},
			{Key: "bias_diff", Val: c.BiasDiff},
			{Key: "combined_mcse", Val: c.CombinedMCSE},
			{Key: "diff_z_score", Val: c.DiffZScore},
			{Key: "within_2_combined_mcse", Val: c.Within2CombinedMCSE},
			{Key: "within_3_combined_mcse", Val: c.Within3CombinedMCSE},
		}
	}

	infillRows := make([]jsonout.Obj, 11)
	for i, sc := range infillScaling {
		infillRows[i] = jsonout.Obj{
			{Key: "name", Val: sc.Name},
			{Key: "monthly_bias", Val: sc.MonthlyBias},
			{Key: "weekly_bias", Val: sc.WeeklyBias},
			{Key: "daily_bias", Val: sc.DailyBias},
			{Key: "monthly_rmse", Val: sc.MonthlyRMSE},
			{Key: "weekly_rmse", Val: sc.WeeklyRMSE},
			{Key: "daily_rmse", Val: sc.DailyRMSE},
			{Key: "monthly_emp_sd", Val: sc.MonthlyEmpSD},
			{Key: "weekly_emp_sd", Val: sc.WeeklyEmpSD},
			{Key: "daily_emp_sd", Val: sc.DailyEmpSD},
			{Key: "daily_over_monthly_rmse", Val: sc.DailyOverMonthlyRMSE},
			{Key: "daily_over_monthly_emp_sd", Val: sc.DailyOverMonthlyEmpSD},
			{Key: "daily_over_monthly_hessian_se", Val: sc.DailyOverMonthlyHessSE},
			{Key: "theoretical_sqrt_n_obs_ratio", Val: sc.SqrtNObsRatio},
		}
	}

	lsTrendObjs := make([]jsonout.Obj, 3)
	for k := 0; k < 3; k++ {
		lsTrendObjs[k] = trendToJSON(kpTrends[k])
	}
	allTrendObjs := make([]jsonout.Obj, 11)
	for i := 0; i < 11; i++ {
		allTrendObjs[i] = trendToJSON(allTrends[i])
	}

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day004b"},
			{Key: "date_wita", Val: "2026-11-22"},
			{Key: "topic", Val: "AFNS — bias sampel terbatas, rentang kalender (long-span), dan frekuensi observasi (infill)"},
			{Key: "quick_mode", Val: *quick},
			{Key: "main_seed", Val: rng.MainSeed},
			{Key: "n_rep", Val: nRep},
			{Key: "workers_used", Val: nWorkers},
			{Key: "units", Val: jsonout.Obj{
				{Key: "yields_and_theta", Val: "desimal per tahun (0.045 = 4.5%)"},
				{Key: "lambda_and_kp", Val: "tahun^-1 (yr^-1)"},
				{Key: "sigma_factors", Val: "tahun^-1/2 (yr^-1/2)"},
				{Key: "sigma_eps", Val: "desimal yield per observasi iid (0.0005 = 5 bp)"},
				{Key: "tenors_and_horizon", Val: "tahun"},
			}},
			{Key: "true_spec", Val: jsonout.Obj{
				{Key: "lambda", Val: trueSpec.AFNS.Lambda},
				{Key: "kp_diag", Val: []float64{trueSpec.AFNS.KP[0][0], trueSpec.AFNS.KP[1][1], trueSpec.AFNS.KP[2][2]}},
				{Key: "theta_p", Val: trueSpec.AFNS.ThetaP[:]},
				{Key: "sigma", Val: trueSpec.AFNS.Sigma[:]},
				{Key: "sigma_eps", Val: trueSpec.SigmaEps},
				{Key: "tenors", Val: kalman.StandardTenors},
			}},
			{Key: "bounds", Val: jsonout.Obj{
				{Key: "param_names", Val: kalman.ParamNames[:]},
				{Key: "lower_original", Val: kalman.LowerBoundsOriginal[:]},
				{Key: "upper_original", Val: kalman.UpperBoundsOriginal[:]},
			}},
			{Key: "elapsed_seconds", Val: jsonout.Obj{
				{Key: "part1_anchor_and_determinism", Val: elapsedPart1},
				{Key: "baseline_1000w", Val: baselineRes.ElapsedSeconds},
				{Key: "longspan_weekly_total", Val: longSpanElapsed},
				{Key: "infill_20y_monthly", Val: infillMonthlyRes.ElapsedSeconds},
				{Key: "infill_20y_weekly_reused", Val: 0.0},
				{Key: "infill_20y_daily", Val: infillDailyRes.ElapsedSeconds},
				{Key: "figures_svg", Val: elapsedFigures},
				{Key: "total", Val: elapsedTotal},
			}},
		}},
		{Key: "part1_exact_anchors_and_verification", Val: jsonout.Obj{
			{Key: "joint_gaussian_t5_3tenors", Val: jsonout.Obj{
				{Key: "t_steps", Val: anchorRes.T},
				{Key: "num_tenors", Val: anchorRes.NumTenors},
				{Key: "tenors", Val: anchorTenors},
				{Key: "kalman_loglik", Val: anchorRes.KalmanLogLik},
				{Key: "joint_gaussian_loglik", Val: anchorRes.JointLogLik},
				{Key: "loglik_abs_diff", Val: anchorRes.LogLikAbsDiff},
				{Key: "fast_ll_vs_joint_abs_diff", Val: fastAncDiff},
				{Key: "rts_smooth_mean_max_abs_diff", Val: anchorRes.MaxSmoothMeanAbsDiff},
				{Key: "rts_smooth_cov_max_abs_diff", Val: anchorRes.MaxSmoothCovAbsDiff},
			}},
			{Key: "score_fd_vs_richardson", Val: jsonout.Obj{
				{Key: "max_rel_err_at_init", Val: maxGradErrInit},
				{Key: "max_rel_err_at_true", Val: maxGradErrTrue},
			}},
			{Key: "worker_determinism_fixture", Val: jsonout.Obj{
				{Key: "fixture_t_steps", Val: detCfg.TSteps},
				{Key: "fixture_n_rep", Val: detCfg.NumReps},
				{Key: "workers_tested", Val: []int{1, nWorkers}},
				{Key: "bit_identical_1_vs_n", Val: workerDetPassed},
			}},
		}},
		{Key: "experiments", Val: jsonout.Obj{
			{Key: "baseline_1000w", Val: studyResultToJSON(baselineRes)},
			{Key: "longspan_10y_weekly", Val: studyResultToJSON(longSpanRes[0])},
			{Key: "longspan_20y_weekly", Val: studyResultToJSON(longSpanRes[1])},
			{Key: "longspan_40y_weekly", Val: studyResultToJSON(longSpanRes[2])},
			{Key: "longspan_80y_weekly", Val: studyResultToJSON(longSpanRes[3])},
			{Key: "infill_20y_monthly", Val: studyResultToJSON(infillMonthlyRes)},
			{Key: "infill_20y_weekly", Val: studyResultToJSON(infillWeeklyRes)},
			{Key: "infill_20y_daily", Val: studyResultToJSON(infillDailyRes)},
		}},
		{Key: "research_analysis", Val: jsonout.Obj{
			{Key: "baseline_vs_day004a", Val: jsonout.Obj{
				{Key: "max_abs_z", Val: baselineCmp.MaxAbsZ},
				{Key: "count_within_2_combined_mcse", Val: baselineCmp.CountWithin2SE},
				{Key: "count_within_3_combined_mcse", Val: baselineCmp.CountWithin3SE},
				{Key: "all_within_3_combined_mcse", Val: baselineCmp.AllWithin3SE},
				{Key: "parameters", Val: cmpRows},
			}},
			{Key: "longspan_kp_bias_vs_inv_h", Val: lsTrendObjs},
			{Key: "longspan_all11_trends_vs_inv_h", Val: allTrendObjs},
			{Key: "infill_20y_scaling", Val: infillRows},
		}},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_exact_anchor_and_score_fd", Val: c1Anchor},
			{Key: "c1_status", Val: statusStr(c1Anchor)},
			{Key: "c2_baseline_reproducible_within_3_mcse", Val: c2BaselineRepro},
			{Key: "c2_status", Val: statusStr(c2BaselineRepro)},
			{Key: "c3_complete_unfiltered_summaries", Val: c3Completeness},
			{Key: "c3_status", Val: statusStr(c3Completeness)},
			{Key: "c4_worker_determinism_and_runtime_lt_600s", Val: c4Determinism},
			{Key: "c4_status", Val: statusStr(c4Determinism)},
			{Key: "all_methodological_passed", Val: allMethPassed},
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		fmt.Fprintf(os.Stderr, "jsonout.WriteFile error: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("day004b selesai dalam %.2fs | workers=%d | c1=%v c2=%v c3=%v c4=%v | all_methodological_passed=%v\n",
		elapsedTotal, nWorkers, c1Anchor, c2BaselineRepro, c3Completeness, c4Determinism, allMethPassed)
	fmt.Printf("Output JSON: %s | Figures dir: %s\n", outPath, figDir)
}
