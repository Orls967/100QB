package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/afns"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/interval"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/jsonout"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/kalman"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/parallel"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/rng"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/stats"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/svg"
)

func findRepoRoot() string {
	wd, err := os.Getwd()
	if err != nil {
		return "."
	}
	dir := wd
	for i := 0; i < 8; i++ {
		if _, err := os.Stat(filepath.Join(dir, "outputs")); err == nil {
			if _, err2 := os.Stat(filepath.Join(dir, "projects", "go", "day004d_afns_interval_calibration")); err2 == nil {
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

func optDiagToJSON(d bootstrap.OptAggregateDiag) jsonout.Obj {
	return jsonout.Obj{
		{Key: "total_fits", Val: d.TotalFits},
		{Key: "finite_count", Val: d.FiniteCount},
		{Key: "non_finite_count", Val: d.NonFiniteCount},
		{Key: "converged_count", Val: d.ConvergedCount},
		{Key: "converged_pct", Val: d.ConvergedPct},
		{Key: "scaled_grad_le_1e6_count", Val: d.ScaledGradLe1e6Count},
		{Key: "scaled_grad_le_1e6_pct", Val: d.ScaledGradLe1e6Pct},
		{Key: "bound_hit_count", Val: d.BoundHitCount},
		{Key: "bound_hit_pct", Val: d.BoundHitPct},
		{Key: "hessian_regularized_count", Val: d.HessianRegularizedCount},
		{Key: "termination_reason_counts", Val: jsonout.Obj{
			{Key: "gradient_tolerance", Val: d.GradTolCount},
			{Key: "boundary_stationary", Val: d.BoundaryStationaryCount},
			{Key: "line_search_failed", Val: d.LineSearchFailedCount},
			{Key: "max_iterations", Val: d.MaxIterCount},
			{Key: "invalid_objective_or_gradient", Val: d.InvalidObjCount},
		}},
		{Key: "line_search_diagnostics", Val: jsonout.Obj{
			{Key: "obj_bit_identical_count", Val: d.LSObjBitIdenticalCount},
			{Key: "armijo_not_met_count", Val: d.LSArmijoNotMetCount},
			{Key: "step_underflow_plateau_count", Val: d.LSStepUnderflowCount},
			{Key: "invalid_objective_count", Val: d.LSInvalidObjCount},
		}},
	}
}

func methodSummaryToJSON(m bootstrap.MethodSummary) jsonout.Obj {
	return jsonout.Obj{
		{Key: "method", Val: m.Method},
		{Key: "applicable", Val: m.Applicable},
		{Key: "total_outer_reps", Val: m.TotalOuterReps},
		{Key: "valid_count", Val: m.ValidCount},
		{Key: "invalid_count", Val: m.InvalidCount},
		{Key: "covered_count", Val: m.CoveredCount},
		{Key: "physically_valid_covered_count", Val: m.PhysicallyValidCoveredCount},
		{Key: "coverage_95", Val: m.Coverage95},
		{Key: "wilson_95_lower", Val: m.WilsonLower},
		{Key: "wilson_95_upper", Val: m.WilsonUpper},
		{Key: "mean_width", Val: m.MeanWidth},
		{Key: "median_width", Val: m.MedianWidth},
		{Key: "q25_width", Val: m.Q25Width},
		{Key: "q75_width", Val: m.Q75Width},
		{Key: "iqr_width", Val: m.IQRWidth},
		{Key: "p95_width", Val: m.P95Width},
		{Key: "max_width", Val: m.MaxWidth},
		{Key: "non_physical_count", Val: m.NonPhysicalCount},
		{Key: "non_physical_pct", Val: m.NonPhysicalPct},
		{Key: "out_of_bounds_count", Val: m.OutOfBoundsCount},
		{Key: "out_of_bounds_pct", Val: m.OutOfBoundsPct},
		{Key: "invalid_reasons", Val: jsonout.Obj{
			{Key: "outer_non_finite", Val: m.InvalidOuterNonFin},
			{Key: "outer_se_non_positive", Val: m.InvalidOuterSENonPos},
			{Key: "outer_non_positive_for_log", Val: m.InvalidOuterNonPosLg},
			{Key: "insufficient_finite_draws", Val: m.InvalidInsuffFinite},
			{Key: "insufficient_valid_se_draws", Val: m.InvalidInsuffValidSE},
			{Key: "insufficient_positive_draws", Val: m.InvalidInsuffPos},
		}},
	}
}

func intervalResToJSON(ir interval.IntervalResult) jsonout.Obj {
	lo, hi, w := ir.Lower, ir.Upper, ir.Width
	if !ir.Valid {
		lo, hi, w = 0.0, 0.0, 0.0
	}
	return jsonout.Obj{
		{Key: "valid", Val: ir.Valid},
		{Key: "invalid_reason", Val: ir.InvalidReason},
		{Key: "lower", Val: lo},
		{Key: "upper", Val: hi},
		{Key: "width", Val: w},
		{Key: "covered", Val: ir.Covered},
		{Key: "non_physical", Val: ir.NonPhysical},
		{Key: "out_of_bounds", Val: ir.OutOfBounds},
		{Key: "valid_draws_count", Val: ir.ValidDrawsCount},
	}
}

func paramSummaryToJSON(p bootstrap.ParamCalibrationSummary) jsonout.Obj {
	return jsonout.Obj{
		{Key: "name", Val: p.Name},
		{Key: "true_value", Val: p.TrueValue},
		{Key: "total_outer_reps", Val: p.TotalOuterReps},
		{Key: "valid_corrected_outer_count", Val: p.ValidCorrectedOuterCount},
		{Key: "raw_mean", Val: p.RawMean},
		{Key: "raw_bias", Val: p.RawBias},
		{Key: "raw_bias_mcse", Val: p.RawBiasMCSE},
		{Key: "raw_empirical_sd", Val: p.RawEmpiricalSD},
		{Key: "raw_rmse", Val: p.RawRMSE},
		{Key: "raw_bound_hit_count", Val: p.RawBoundHitCount},
		{Key: "raw_bound_hit_pct", Val: p.RawBoundHitPct},
		{Key: "boot_mean_conditional_bias", Val: p.BootMeanConditionalBias},
		{Key: "boot_corrected_mean", Val: p.BootCorrectedMean},
		{Key: "boot_corrected_bias", Val: p.BootCorrectedBias},
		{Key: "boot_corrected_bias_mcse", Val: p.BootCorrectedBiasMCSE},
		{Key: "boot_corrected_empirical_sd", Val: p.BootCorrectedEmpiricalSD},
		{Key: "boot_corrected_rmse", Val: p.BootCorrectedRMSE},
		{Key: "rmse_ratio", Val: p.RMSERatio},
		{Key: "intervals", Val: jsonout.Obj{
			{Key: "wald", Val: methodSummaryToJSON(p.Wald)},
			{Key: "percentile", Val: methodSummaryToJSON(p.Percentile)},
			{Key: "basic", Val: methodSummaryToJSON(p.Basic)},
			{Key: "studentized", Val: methodSummaryToJSON(p.Studentized)},
			{Key: "log_basic", Val: methodSummaryToJSON(p.LogBasic)},
		}},
	}
}

func checkComplete(st bootstrap.CalibrationStudyResult, expR, expB int) bool {
	if len(st.OuterReplications) != expR || st.TotalInnerRefits != expR*expB {
		return false
	}
	ok := func(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }
	for _, p := range st.Params {
		if p.TotalOuterReps != expR {
			return false
		}
		for _, m := range []bootstrap.MethodSummary{p.Wald, p.Percentile, p.Basic, p.Studentized, p.LogBasic} {
			if !m.Applicable {
				continue
			}
			if m.ValidCount+m.InvalidCount != expR {
				return false
			}
			if !ok(m.Coverage95) || !ok(m.WilsonLower) || !ok(m.WilsonUpper) {
				return false
			}
		}
	}
	return true
}

func checkAlgebra(st bootstrap.CalibrationStudyResult) bool {
	tol := 1e-11
	for _, p := range st.Params {
		if math.Abs((p.RawMean-p.TrueValue)-p.RawBias) > tol {
			return false
		}
		if math.Abs((p.BootCorrectedMean-p.TrueValue)-p.BootCorrectedBias) > tol {
			return false
		}
		if p.RawRMSE > 0 && math.Abs((p.BootCorrectedRMSE/p.RawRMSE)-p.RMSERatio) > tol {
			return false
		}
		for _, m := range []bootstrap.MethodSummary{p.Wald, p.Percentile, p.Basic, p.Studentized, p.LogBasic} {
			if !m.Applicable {
				continue
			}
			expRate := float64(m.CoveredCount) / float64(m.TotalOuterReps)
			if math.Abs(expRate-m.Coverage95) > tol {
				return false
			}
		}
	}
	return true
}

func buildFigures(figDir string, st bootstrap.CalibrationStudyResult) error {
	if err := os.MkdirAll(figDir, 0o755); err != nil {
		return err
	}

	// Grafik 1: Cakupan 95% per metode untuk kp11 dan theta1 (+ Wilson 95% CI)
	cats1 := []string{"kp11 (true=0.05)", "theta1 (true=0.045)"}
	kp := st.Params[1]
	th := st.Params[4]
	y0, y1 := 0.0, 1.05
	if err := svg.WriteCategoryPointChart(filepath.Join(figDir, "day004d_coverage-by-method-kp11-theta1.svg"), svg.ChartOptions{
		Title:        "Cakupan Interval 95% Antar-Metode Bootstrap pada kp11 & theta1 (+ Wilson 95% CI)",
		XLabel:       "Parameter Sasaran (R=100 Replikasi Luar, B=100 Refit)",
		YLabel:       "Proporsi Cakupan 95% (+/- Wilson 95% CI)",
		YMinOverride: &y0,
		YMaxOverride: &y1,
		LegendBottom: true,
		HRefLines: []svg.HRefLine{
			{Value: 0.95, Label: "Nominal 95%", Color: "#111111"},
		},
	}, cats1, []svg.CategoryPointSeries{
		{
			Label:  "Wald Hessian",
			Values: []float64{kp.Wald.Coverage95, th.Wald.Coverage95},
			Lower:  []float64{kp.Wald.WilsonLower, th.Wald.WilsonLower},
			Upper:  []float64{kp.Wald.WilsonUpper, th.Wald.WilsonUpper},
			Color:  svg.Palette[0],
		},
		{
			Label:  "Percentile (Raw)",
			Values: []float64{kp.Percentile.Coverage95, th.Percentile.Coverage95},
			Lower:  []float64{kp.Percentile.WilsonLower, th.Percentile.WilsonLower},
			Upper:  []float64{kp.Percentile.WilsonUpper, th.Percentile.WilsonUpper},
			Color:  svg.Palette[1],
		},
		{
			Label:  "Basic / Reversed",
			Values: []float64{kp.Basic.Coverage95, th.Basic.Coverage95},
			Lower:  []float64{kp.Basic.WilsonLower, th.Basic.WilsonLower},
			Upper:  []float64{kp.Basic.WilsonUpper, th.Basic.WilsonUpper},
			Color:  svg.Palette[2],
		},
		{
			Label:  "Studentized (Boot-t)",
			Values: []float64{kp.Studentized.Coverage95, th.Studentized.Coverage95},
			Lower:  []float64{kp.Studentized.WilsonLower, th.Studentized.WilsonLower},
			Upper:  []float64{kp.Studentized.WilsonUpper, th.Studentized.WilsonUpper},
			Color:  svg.Palette[3],
		},
		{
			Label:  "Log-Basic (kp11>0)",
			Values: []float64{kp.LogBasic.Coverage95},
			Lower:  []float64{kp.LogBasic.WilsonLower},
			Upper:  []float64{kp.LogBasic.WilsonUpper},
			Color:  svg.Palette[4],
		},
	}); err != nil {
		return err
	}

	// Grafik 2: Distribusi batas bawah interval kp11 skala asli (Basic) vs skala log (Log-Basic)
	var basicLos, logLos []float64
	for _, rep := range st.OuterReplications {
		if rep.BasicIntervals[1].Valid {
			basicLos = append(basicLos, rep.BasicIntervals[1].Lower)
		}
		if rep.LogBasicIntervals[1].Valid {
			logLos = append(logLos, rep.LogBasicIntervals[1].Lower)
		}
	}
	if err := svg.WriteDualHistogramSVG(
		filepath.Join(figDir, "day004d_kp11-endpoints-raw-vs-log.svg"),
		"Distribusi Batas Bawah Interval 95% kp11: Skala Asli (Basic) vs Skala Log (Log-Basic)",
		svg.HistogramPanel{
			Subtitle: fmt.Sprintf("Basic Skala Asli (%d/%d negatif <= 0)", kp.Basic.NonPhysicalCount, st.NumOuterReps),
			XLabel:   "Batas Bawah Basic (Tanpa Clamping)",
			Values:   basicLos,
			ZeroLine: true,
			Bins:     16,
			Markers: []svg.VMarker{
				{Value: 0.0, Label: "Batas Fisik 0", Color: "#d62728", Dashed: true},
				{Value: 0.05, Label: "True kp11=0.05", Color: "#2ca02c"},
			},
		},
		svg.HistogramPanel{
			Subtitle: fmt.Sprintf("Log-Basic (%d/%d negatif <= 0)", kp.LogBasic.NonPhysicalCount, st.NumOuterReps),
			XLabel:   "Batas Bawah Log-Basic exp(2*ln(hat)-q97.5)",
			Values:   logLos,
			ZeroLine: true,
			Bins:     16,
			Markers: []svg.VMarker{
				{Value: 0.0, Label: "Batas Fisik 0", Color: "#d62728", Dashed: true},
				{Value: 0.05, Label: "True kp11=0.05", Color: "#2ca02c"},
			},
		},
	); err != nil {
		return err
	}

	// Grafik 3: Lebar interval antar-metode untuk theta1 (Median dengan rentang IQR [Q25, Q75])
	var wWidths, pWidths, bWidths, sWidths []float64
	for _, rep := range st.OuterReplications {
		if rep.WaldIntervals[4].Valid {
			wWidths = append(wWidths, rep.WaldIntervals[4].Width)
		}
		if rep.PctIntervals[4].Valid {
			pWidths = append(pWidths, rep.PctIntervals[4].Width)
		}
		if rep.BasicIntervals[4].Valid {
			bWidths = append(bWidths, rep.BasicIntervals[4].Width)
		}
		if rep.StudIntervals[4].Valid {
			sWidths = append(sWidths, rep.StudIntervals[4].Width)
		}
	}
	medIQ := func(xs []float64) (med, q25, q75 float64) {
		if len(xs) == 0 {
			return 0, 0, 0
		}
		return stats.Quantile(xs, 0.5), stats.Quantile(xs, 0.25), stats.Quantile(xs, 0.75)
	}
	wM, w25, w75 := medIQ(wWidths)
	pM, p25, p75 := medIQ(pWidths)
	bM, b25, b75 := medIQ(bWidths)
	sM, s25, s75 := medIQ(sWidths)
	theoEmpWidth := 2.0 * stats.Z95 * th.RawEmpiricalSD
	if err := svg.WriteCategoryPointChart(filepath.Join(figDir, "day004d_theta1-interval-widths.svg"), svg.ChartOptions{
		Title:        "Lebar Interval 95% theta1 Antar-Metode (Median & IQR [Q25, Q75])",
		XLabel:       "Metode Interval (H=20y, R=100, B=100)",
		YLabel:       "Lebar Interval theta1",
		LegendBottom: false,
		HRefLines: []svg.HRefLine{
			{Value: theoEmpWidth, Label: fmt.Sprintf("2 * 1.96 * EmpSD = %.4f", theoEmpWidth), Color: "#111111"},
		},
	}, []string{"Wald", "Percentile", "Basic", "Studentized"}, []svg.CategoryPointSeries{
		{
			Label:  "Median & IQR [Q25, Q75]",
			Values: []float64{wM, pM, bM, sM},
			Lower:  []float64{w25, p25, b25, s25},
			Upper:  []float64{w75, p75, b75, s75},
			Color:  svg.Palette[0],
		},
	}); err != nil {
		return err
	}

	// Grafik 4: Hubungan kp11 luar terhadap lebar & cakupan bootstrap Studentized theta1
	var covX, covY, missX, missY []float64
	for _, rep := range st.OuterReplications {
		if !rep.OuterMLE.Finite || !rep.StudIntervals[4].Valid {
			continue
		}
		if rep.StudIntervals[4].Covered {
			covX = append(covX, rep.OuterMLE.EstParams[1])
			covY = append(covY, rep.StudIntervals[4].Width)
		} else {
			missX = append(missX, rep.OuterMLE.EstParams[1])
			missY = append(missY, rep.StudIntervals[4].Width)
		}
	}
	if err := svg.WriteScatterChart(filepath.Join(figDir, "day004d_kp11-vs-theta1-width-coverage.svg"), svg.ChartOptions{
		Title:  "Hubungan Estimasi kp11 Luar dengan Lebar & Cakupan Studentized theta1",
		XLabel: "Estimasi MLE Luar kp11_hat (tahun^-1)",
		YLabel: "Lebar Interval 95% Studentized theta1",
	}, []svg.VMarker{
		{Value: 0.05, Label: "True kp11=0.05", Color: "#2ca02c", Dashed: true},
	}, []svg.ScatterSeries{
		{Label: fmt.Sprintf("Tercakup (%d)", len(covX)), X: covX, Y: covY, Color: svg.Palette[2]},
		{Label: fmt.Sprintf("Tidak Tercakup (%d)", len(missX)), X: missX, Y: missY, Color: svg.Palette[1]},
	}); err != nil {
		return err
	}

	return nil
}

func main() {
	root := findRepoRoot()
	quick := flag.Bool("quick", false, "jalankan mode cepat pengembangan ke direktori sementara tanpa menimpa keluaran penuh")
	workersFlag := flag.Int("workers", 0, "jumlah worker goroutine (0 = min(8, NumCPU))")
	repsFlag := flag.Int("r", 0, "override jumlah replikasi luar R (bawaan: 100 penuh, 4 -quick)")
	bootFlag := flag.Int("b", 0, "override jumlah bootstrap dalam B per data luar (bawaan: 100 penuh, 6 -quick)")
	outFlag := flag.String("out", "", "override berkas keluaran JSON")
	figDirFlag := flag.String("fig-dir", "", "override direktori keluaran SVG")
	flag.Parse()

	nWorkers := *workersFlag
	if nWorkers <= 0 {
		nWorkers = parallel.DefaultWorkers()
	}
	nRep := 100
	nBoot := 100
	if *quick {
		nRep = 4
		nBoot = 6
	}
	if *repsFlag > 0 {
		nRep = *repsFlag
	}
	if *bootFlag > 0 {
		nBoot = *bootFlag
	}

	outPath := filepath.Join(root, "outputs", "day004d", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	if *quick && *outFlag == "" {
		outPath = filepath.Join(os.TempDir(), "day004d_quick_results.json")
	}
	if *quick && *figDirFlag == "" {
		figDir = filepath.Join(os.TempDir(), "day004d_quick_figures")
	}
	if *outFlag != "" {
		outPath = *outFlag
	}
	if *figDirFlag != "" {
		figDir = *figDirFlag
	}

	totalT0 := time.Now()
	trueSpec := kalman.DefaultSpec()

	// Bagian 1a: Anchor eksak Gaussian gabungan T=5 (3 tenor)
	tPart1 := time.Now()
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

	// Bagian 1b: Verifikasi determinisme 1 vs nWorkers pada fixture cepat
	detOuter1 := bootstrap.RunOuterFits(trueSpec, kalman.StandardTenors, kalman.WeeklyDt, 80, 2, 1, rng.MainSeed, 0xD001)
	detOuterN := bootstrap.RunOuterFits(trueSpec, kalman.StandardTenors, kalman.WeeklyDt, 80, 2, nWorkers, rng.MainSeed, 0xD001)
	detB1 := bootstrap.RunCalibrationStudy(trueSpec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, detOuter1, 0, 4, 1, rng.MainSeed, 0xD002)
	detBN := bootstrap.RunCalibrationStudy(trueSpec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, detOuterN, 0, 4, nWorkers, rng.MainSeed, 0xD002)
	workerDetPassed := true
	for r := 0; r < 2; r++ {
		if math.Float64bits(detB1.OuterReplications[r].OuterMLE.EstLogLik) != math.Float64bits(detBN.OuterReplications[r].OuterMLE.EstLogLik) {
			workerDetPassed = false
		}
		for i := 0; i < 11; i++ {
			if math.Float64bits(detB1.OuterReplications[r].StudIntervals[i].Lower) != math.Float64bits(detBN.OuterReplications[r].StudIntervals[i].Lower) ||
				math.Float64bits(detB1.OuterReplications[r].LogBasicIntervals[i].Lower) != math.Float64bits(detBN.OuterReplications[r].LogBasicIntervals[i].Lower) {
				workerDetPassed = false
			}
		}
	}
	elapsedPart1 := time.Since(tPart1).Seconds()

	// Bagian 2: Fit Luar R=100 (H=20y, T=1040 mingguan)
	tOut0 := time.Now()
	fmt.Printf("[day004d] Menjalankan R=%d fit luar MLE 20y mingguan (T=1040, workers=%d)...\n", nRep, nWorkers)
	outerFits := bootstrap.RunOuterFits(trueSpec, kalman.StandardTenors, kalman.WeeklyDt, 1040, nRep, nWorkers, rng.MainSeed, 0x4D01)
	elapsedOuter := time.Since(tOut0).Seconds()

	// Bagian 3: Bootstrap Parametrik + Hessian SE* + 5 Metode Interval (R x B = 10.000 refit)
	fmt.Printf("[day004d] Menjalankan bootstrap parametrik & kalibrasi interval (R=%d, B=%d -> %d refit, workers=%d)...\n", nRep, nBoot, nRep*nBoot, nWorkers)
	study := bootstrap.RunCalibrationStudy(trueSpec, kalman.StandardTenors, 20.0, kalman.WeeklyDt, 1040, outerFits, elapsedOuter, nBoot, nWorkers, rng.MainSeed, 0x4D02)

	maxGradErrInit, maxGradErrTrue := 0.0, 0.0
	for _, rep := range study.OuterReplications {
		if rep.OuterMLE.GradRelErrInit > maxGradErrInit {
			maxGradErrInit = rep.OuterMLE.GradRelErrInit
		}
		if rep.OuterMLE.GradRelErrTrue > maxGradErrTrue {
			maxGradErrTrue = rep.OuterMLE.GradRelErrTrue
		}
	}

	tSvg0 := time.Now()
	if err := buildFigures(figDir, study); err != nil {
		fmt.Fprintf(os.Stderr, "buildFigures error: %v\n", err)
		os.Exit(1)
	}
	elapsedSvg := time.Since(tSvg0).Seconds()
	elapsedTotal := time.Since(totalT0).Seconds()

	c1AnchorScore := anchorRes.LogLikAbsDiff < 1e-9 &&
		anchorRes.MaxSmoothMeanAbsDiff < 1e-9 &&
		anchorRes.MaxSmoothCovAbsDiff < 1e-9 &&
		fastAncDiff < 1e-9 &&
		maxGradErrInit < 1e-5 &&
		maxGradErrTrue < 1e-5
	c2Complete := checkComplete(study, nRep, nBoot)
	c3SummaryAlg := checkAlgebra(study)
	c4Determinism := workerDetPassed
	c5Runtime := elapsedTotal < 600.0 || *quick
	allTechPassed := c1AnchorScore && c2Complete && c3SummaryAlg && c4Determinism && c5Runtime

	statusStr := func(ok bool) string {
		if ok {
			return "tercapai"
		}
		return "belum tercapai"
	}

	pObjs := make([]jsonout.Obj, 11)
	for i := 0; i < 11; i++ {
		pObjs[i] = paramSummaryToJSON(study.Params[i])
	}

	repObjs := make([]jsonout.Obj, len(study.OuterReplications))
	for rIdx, rep := range study.OuterReplications {
		hb := make([]bool, 11)
		for i := 0; i < 11; i++ {
			hb[i] = rep.OuterMLE.HitBoundByParam[i]
		}
		targetIntervals := make([]jsonout.Obj, 4)
		for k, pIdx := range bootstrap.TargetParamIndices {
			targetIntervals[k] = jsonout.Obj{
				{Key: "name", Val: kalman.ParamNames[pIdx]},
				{Key: "outer_est", Val: rep.OuterMLE.EstParams[pIdx]},
				{Key: "outer_hessian_se", Val: rep.OuterMLE.HessianSE[pIdx]},
				{Key: "boot_bias", Val: rep.BootBias[pIdx]},
				{Key: "boot_corrected", Val: rep.BootCorrected[pIdx]},
				{Key: "wald", Val: intervalResToJSON(rep.WaldIntervals[pIdx])},
				{Key: "percentile", Val: intervalResToJSON(rep.PctIntervals[pIdx])},
				{Key: "basic", Val: intervalResToJSON(rep.BasicIntervals[pIdx])},
				{Key: "studentized", Val: intervalResToJSON(rep.StudIntervals[pIdx])},
				{Key: "log_basic", Val: intervalResToJSON(rep.LogBasicIntervals[pIdx])},
			}
		}
		innerObjs := make([]jsonout.Obj, len(rep.InnerDraws))
		for bIdx, d := range rep.InnerDraws {
			ihb := make([]bool, 11)
			for i := 0; i < 11; i++ {
				ihb[i] = d.HitBoundByParam[i]
			}
			innerObjs[bIdx] = jsonout.Obj{
				{Key: "boot_idx", Val: d.BootIdx},
				{Key: "finite", Val: d.Finite},
				{Key: "converged", Val: d.Converged},
				{Key: "termination_reason", Val: d.TerminationReason},
				{Key: "hit_bound", Val: d.HitBound},
				{Key: "hit_bound_by_param", Val: ihb},
				{Key: "hessian_regularized", Val: d.HessianRegularized},
				{Key: "best_stage", Val: d.BestStage},
				{Key: "iterations_best_stage", Val: d.IterationsBest},
				{Key: "iterations_cumulative", Val: d.IterationsCum},
				{Key: "scaled_grad_norm_opt", Val: d.ScaledGradNormOpt},
				{Key: "line_search_failure_cause", Val: d.LineSearchFailureCause},
				{Key: "est_loglik", Val: d.EstLogLik},
				{Key: "est_params", Val: d.EstParams[:]},
				{Key: "hessian_se", Val: d.HessianSE[:]},
			}
		}
		repObjs[rIdx] = jsonout.Obj{
			{Key: "outer_rep_idx", Val: rep.OuterRepIdx},
			{Key: "finite", Val: rep.OuterMLE.Finite},
			{Key: "converged", Val: rep.OuterMLE.Converged},
			{Key: "termination_reason", Val: rep.OuterMLE.TerminationReason},
			{Key: "line_search_failure_cause", Val: rep.OuterMLE.LineSearchFailureCause},
			{Key: "hit_bound", Val: rep.OuterMLE.HitBound},
			{Key: "hit_bound_by_param", Val: hb},
			{Key: "hessian_regularized", Val: rep.OuterMLE.HessianRegularized},
			{Key: "est_loglik", Val: rep.OuterMLE.EstLogLik},
			{Key: "true_loglik", Val: rep.OuterMLE.TrueLogLik},
			{Key: "scaled_grad_norm_opt", Val: rep.OuterMLE.ScaledGradNormOpt},
			{Key: "est_params", Val: rep.OuterMLE.EstParams[:]},
			{Key: "hessian_se", Val: rep.OuterMLE.HessianSE[:]},
			{Key: "finite_draws_count", Val: rep.FiniteDrawsCount},
			{Key: "invalid_draws_count", Val: rep.InvalidDrawsCount},
			{Key: "target_intervals", Val: targetIntervals},
			{Key: "inner_draws", Val: innerObjs},
		}
	}

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day004d"},
			{Key: "date_wita", Val: "2026-11-22"},
			{Key: "topic", Val: "AFNS — interval bootstrap yang peka bias dan batas parameter"},
			{Key: "quick_mode", Val: *quick},
			{Key: "main_seed", Val: rng.MainSeed},
			{Key: "n_rep", Val: nRep},
			{Key: "n_boot", Val: nBoot},
			{Key: "total_inner_refits", Val: study.TotalInnerRefits},
			{Key: "min_required_valid_per_outer", Val: study.MinRequiredValidPerOuter},
			{Key: "min_valid_fraction_threshold", Val: study.MinValidFractionThreshold},
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
				{Key: "dt_years", Val: kalman.WeeklyDt},
				{Key: "t_steps", Val: 1040},
				{Key: "horizon_years", Val: 20.0},
			}},
			{Key: "elapsed_seconds", Val: jsonout.Obj{
				{Key: "part1_anchor_and_determinism", Val: elapsedPart1},
				{Key: "outer_fits", Val: elapsedOuter},
				{Key: "inner_bootstrap_and_intervals", Val: study.InnerElapsedSeconds},
				{Key: "figures_svg", Val: elapsedSvg},
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
				{Key: "fixture_t_steps", Val: 80},
				{Key: "fixture_outer_r", Val: 2},
				{Key: "fixture_boot_b", Val: 4},
				{Key: "workers_tested", Val: []int{1, nWorkers}},
				{Key: "bit_identical_1_vs_n", Val: workerDetPassed},
			}},
		}},
		{Key: "optimizer_diagnostics", Val: jsonout.Obj{
			{Key: "outer_fits", Val: optDiagToJSON(study.OuterOptDiag)},
			{Key: "inner_bootstrap_refits", Val: optDiagToJSON(study.InnerOptDiag)},
		}},
		{Key: "stratification_kp11_vs_theta1", Val: jsonout.Obj{
			{Key: "median_outer_kp11", Val: study.Stratification.MedianOuterKP11},
			{Key: "low_kp11_count", Val: study.Stratification.LowKP11Count},
			{Key: "high_kp11_count", Val: study.Stratification.HighKP11Count},
			{Key: "bound_hit_outer_count", Val: study.Stratification.BoundHitOuterCount},
			{Key: "low_kp11_theta1_wald_cov", Val: study.Stratification.LowKP11Theta1WaldCov},
			{Key: "high_kp11_theta1_wald_cov", Val: study.Stratification.HighKP11Theta1WaldCov},
			{Key: "low_kp11_theta1_stud_cov", Val: study.Stratification.LowKP11Theta1StudCov},
			{Key: "high_kp11_theta1_stud_cov", Val: study.Stratification.HighKP11Theta1StudCov},
			{Key: "low_kp11_theta1_wald_mean_width", Val: study.Stratification.LowKP11Theta1WaldW},
			{Key: "high_kp11_theta1_wald_mean_width", Val: study.Stratification.HighKP11Theta1WaldW},
			{Key: "low_kp11_theta1_stud_mean_width", Val: study.Stratification.LowKP11Theta1StudW},
			{Key: "high_kp11_theta1_stud_mean_width", Val: study.Stratification.HighKP11Theta1StudW},
			{Key: "low_kp11_theta1_wald_median_width", Val: study.Stratification.LowKP11Theta1WaldMedW},
			{Key: "high_kp11_theta1_wald_median_width", Val: study.Stratification.HighKP11Theta1WaldMedW},
			{Key: "low_kp11_theta1_stud_median_width", Val: study.Stratification.LowKP11Theta1StudMedW},
			{Key: "high_kp11_theta1_stud_median_width", Val: study.Stratification.HighKP11Theta1StudMedW},
			{Key: "corr_outer_kp11_with_wald_width", Val: study.Stratification.CorrKP11WithWaldW},
			{Key: "corr_outer_kp11_with_stud_width", Val: study.Stratification.CorrKP11WithStudW},
		}},
		{Key: "parameters", Val: pObjs},
		{Key: "outer_replications", Val: repObjs},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_exact_anchor_and_score", Val: c1AnchorScore},
			{Key: "c1_status", Val: statusStr(c1AnchorScore)},
			{Key: "c2_complete_unfiltered_replications", Val: c2Complete},
			{Key: "c2_status", Val: statusStr(c2Complete)},
			{Key: "c3_summary_math_verified", Val: c3SummaryAlg},
			{Key: "c3_status", Val: statusStr(c3SummaryAlg)},
			{Key: "c4_seed_and_worker_determinism", Val: c4Determinism},
			{Key: "c4_status", Val: statusStr(c4Determinism)},
			{Key: "c5_runtime_within_600s", Val: c5Runtime},
			{Key: "c5_status", Val: statusStr(c5Runtime)},
			{Key: "all_technical_passed", Val: allTechPassed},
			{Key: "technical_scope_note", Val: "all_technical_passed hanya menandakan kriteria teknis C1-C5 terpenuhi dan tidak menyiratkan seluruh metode inferensi terkalibrasi secara statistik"},
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		fmt.Fprintf(os.Stderr, "jsonout.WriteFile error: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("day004d selesai dalam %.2fs | workers=%d | c1=%v c2=%v c3=%v c4=%v c5=%v | all_technical_passed=%v\n",
		elapsedTotal, nWorkers, c1AnchorScore, c2Complete, c3SummaryAlg, c4Determinism, c5Runtime, allTechPassed)
	fmt.Printf("Output JSON: %s | Figures dir: %s\n", outPath, figDir)
}
