package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/afns"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/jackknife"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/jsonout"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/kalman"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/parallel"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/svg"
)

func findRepoRoot() string {
	wd, err := os.Getwd()
	if err != nil {
		return "."
	}
	dir := wd
	for i := 0; i < 8; i++ {
		if _, err := os.Stat(filepath.Join(dir, "outputs")); err == nil {
			if _, err2 := os.Stat(filepath.Join(dir, "projects", "go", "day004c_afns_bias_correction")); err2 == nil {
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

func stageDiagToJSON(s kalman.StageDiagnostic) jsonout.Obj {
	return jsonout.Obj{
		{Key: "stage", Val: s.Stage},
		{Key: "iterations", Val: s.Iterations},
		{Key: "func_evals", Val: s.FuncEvals},
		{Key: "loglik", Val: s.LogLik},
		{Key: "raw_grad_norm", Val: s.RawGradNorm},
		{Key: "grad_norm", Val: s.GradNorm},
		{Key: "scaled_grad_norm", Val: s.ScaledGradNorm},
		{Key: "converged", Val: s.Converged},
		{Key: "hit_bound", Val: s.HitBound},
		{Key: "termination_reason", Val: s.TerminationReason},
	}
}

func mleResultToJSON(r kalman.MLEResult) jsonout.Obj {
	hb := make([]bool, 11)
	for i := 0; i < 11; i++ {
		hb[i] = r.HitBoundByParam[i]
	}
	return jsonout.Obj{
		{Key: "finite", Val: r.Finite},
		{Key: "init_params", Val: r.InitParams[:]},
		{Key: "est_params", Val: r.EstParams[:]},
		{Key: "init_loglik", Val: r.InitLogLik},
		{Key: "true_loglik", Val: r.TrueLogLik},
		{Key: "est_loglik", Val: r.EstLogLik},
		{Key: "grad_rel_err_init", Val: r.GradRelErrInit},
		{Key: "grad_rel_err_true", Val: r.GradRelErrTrue},
		{Key: "grad_abs_diff_opt", Val: r.GradAbsDiffOpt},
		{Key: "raw_grad_norm_opt", Val: r.RawGradNormOpt},
		{Key: "grad_norm_opt", Val: r.GradNormOpt},
		{Key: "scaled_grad_norm_opt", Val: r.ScaledGradNormOpt},
		{Key: "hessian_se", Val: r.HessianSE[:]},
		{Key: "hessian_min_eig", Val: r.HessianMinEig},
		{Key: "hessian_regularized", Val: r.HessianRegularized},
		{Key: "iterations_best_stage", Val: r.BestStageIterations},
		{Key: "iterations_cumulative", Val: r.CumulativeIterations},
		{Key: "best_stage", Val: r.BestStage},
		{Key: "converged", Val: r.Converged},
		{Key: "termination_reason", Val: r.TerminationReason},
		{Key: "hit_bound", Val: r.HitBound},
		{Key: "hit_bound_by_param", Val: hb},
		{Key: "stage1", Val: stageDiagToJSON(r.Stage1)},
		{Key: "stage2", Val: stageDiagToJSON(r.Stage2)},
	}
}

func segmentDiagToJSON(d jackknife.SegmentDiagnostics) jsonout.Obj {
	return jsonout.Obj{
		{Key: "segment", Val: d.Segment},
		{Key: "total_reps", Val: d.TotalReps},
		{Key: "finite_count", Val: d.FiniteCount},
		{Key: "converged_count", Val: d.ConvergedCount},
		{Key: "converged_pct", Val: d.ConvergedPct},
		{Key: "scaled_grad_le_1e6_count", Val: d.ScaledGradLe1e6Count},
		{Key: "scaled_grad_le_1e6_pct", Val: d.ScaledGradLe1e6Pct},
		{Key: "ll_ge_true_count", Val: d.LLGeTrueCount},
		{Key: "bound_hit_count", Val: d.BoundHitCount},
		{Key: "bound_hit_pct", Val: d.BoundHitPct},
		{Key: "hessian_regularized_count", Val: d.HessianRegCount},
		{Key: "termination_reason_counts", Val: jsonout.Obj{
			{Key: "gradient_tolerance", Val: d.GradTolCount},
			{Key: "boundary_stationary", Val: d.BoundaryStationaryCount},
			{Key: "line_search_failed", Val: d.LineSearchFailedCount},
			{Key: "max_iterations", Val: d.MaxIterCount},
			{Key: "invalid_objective_or_gradient", Val: d.InvalidObjCount},
		}},
		{Key: "best_stage_counts", Val: jsonout.Obj{
			{Key: "stage1", Val: d.BestStage1Count},
			{Key: "stage2", Val: d.BestStage2Count},
		}},
	}
}

func pairedParamSummaryToJSON(p jackknife.PairedParamSummary) jsonout.Obj {
	return jsonout.Obj{
		{Key: "name", Val: p.Name},
		{Key: "true_value", Val: p.TrueValue},
		{Key: "total_reps", Val: p.TotalReps},
		{Key: "raw_mean", Val: p.RawMean},
		{Key: "raw_bias", Val: p.RawBias},
		{Key: "raw_bias_mcse", Val: p.RawBiasMCSE},
		{Key: "raw_empirical_sd", Val: p.RawEmpiricalSD},
		{Key: "raw_rmse", Val: p.RawRMSE},
		{Key: "raw_mean_hessian_se", Val: p.RawMeanHessianSE},
		{Key: "raw_median_hessian_se", Val: p.RawMedianHessianSE},
		{Key: "raw_mean_hessian_se_to_emp_sd_ratio", Val: p.RawMeanHessianSEToEmpSDRatio},
		{Key: "raw_median_hessian_se_to_emp_sd_ratio", Val: p.RawMedianHessianSEToEmpSDRatio},
		{Key: "wald_coverage_95", Val: p.WaldCoverage95},
		{Key: "wald_wilson_95_lower", Val: p.WaldWilsonLower},
		{Key: "wald_wilson_95_upper", Val: p.WaldWilsonUpper},
		{Key: "wald_covered_count", Val: p.WaldCoveredCount},
		{Key: "raw_bound_hit_count", Val: p.RawBoundHitCount},
		{Key: "raw_bound_hit_pct", Val: p.RawBoundHitPct},
		{Key: "jackknife_mean", Val: p.JackknifeMean},
		{Key: "jackknife_bias", Val: p.JackknifeBias},
		{Key: "jackknife_bias_mcse", Val: p.JackknifeBiasMCSE},
		{Key: "jackknife_empirical_sd", Val: p.JackknifeEmpiricalSD},
		{Key: "jackknife_rmse", Val: p.JackknifeRMSE},
		{Key: "jackknife_out_of_bounds_count", Val: p.JackknifeOutOfBoundsCount},
		{Key: "jackknife_out_of_bounds_pct", Val: p.JackknifeOutOfBoundsPct},
		{Key: "jackknife_non_physical_count", Val: p.JackknifeNonPhysicalCount},
		{Key: "jackknife_non_physical_pct", Val: p.JackknifeNonPhysicalPct},
		{Key: "paired_bias_diff", Val: p.PairedBiasDiff},
		{Key: "paired_bias_diff_mcse", Val: p.PairedBiasDiffMCSE},
		{Key: "abs_bias_reduction", Val: p.AbsBiasReduction},
		{Key: "rmse_diff", Val: p.RMSEDiff},
		{Key: "rmse_ratio", Val: p.RMSERatio},
	}
}

func jackknifeResultToJSON(jk jackknife.HorizonJackknifeResult) jsonout.Obj {
	pObjs := make([]jsonout.Obj, 11)
	for i := 0; i < 11; i++ {
		pObjs[i] = pairedParamSummaryToJSON(jk.Params[i])
	}
	repObjs := make([]jsonout.Obj, len(jk.Replications))
	for rIdx, rep := range jk.Replications {
		oob := make([]bool, 11)
		nph := make([]bool, 11)
		for i := 0; i < 11; i++ {
			oob[i] = rep.OutOfBoundsByParam[i]
			nph[i] = rep.NonPhysicalByParam[i]
		}
		repObjs[rIdx] = jsonout.Obj{
			{Key: "rep_idx", Val: rep.RepIdx},
			{Key: "full", Val: mleResultToJSON(rep.Full)},
			{Key: "half1", Val: mleResultToJSON(rep.Half1)},
			{Key: "half2", Val: mleResultToJSON(rep.Half2)},
			{Key: "jackknife_est", Val: rep.JackknifeEst[:]},
			{Key: "jackknife_out_of_bounds_by_param", Val: oob},
			{Key: "jackknife_non_physical_by_param", Val: nph},
		}
	}
	return jsonout.Obj{
		{Key: "name", Val: jk.Name},
		{Key: "horizon_years", Val: jk.HorizonYears},
		{Key: "dt_years", Val: jk.Dt},
		{Key: "t_steps", Val: jk.TSteps},
		{Key: "half_t_steps", Val: jk.HalfTSteps},
		{Key: "n_rep", Val: jk.NumReps},
		{Key: "workers_used", Val: jk.NumWorkers},
		{Key: "elapsed_seconds", Val: jk.ElapsedSeconds},
		{Key: "segment_diagnostics", Val: jsonout.Obj{
			{Key: "full", Val: segmentDiagToJSON(jk.FullDiag)},
			{Key: "half1", Val: segmentDiagToJSON(jk.Half1Diag)},
			{Key: "half2", Val: segmentDiagToJSON(jk.Half2Diag)},
		}},
		{Key: "parameters", Val: pObjs},
		{Key: "replications", Val: repObjs},
	}
}

func bootParamSummaryToJSON(p bootstrap.BootParamSummary) jsonout.Obj {
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
		{Key: "boot_mean_conditional_bias", Val: p.BootMeanOfConditionalBias},
		{Key: "boot_corrected_mean", Val: p.BootCorrectedMean},
		{Key: "boot_corrected_bias", Val: p.BootCorrectedBias},
		{Key: "boot_corrected_bias_mcse", Val: p.BootCorrectedBiasMCSE},
		{Key: "boot_corrected_empirical_sd", Val: p.BootCorrectedEmpiricalSD},
		{Key: "boot_corrected_rmse", Val: p.BootCorrectedRMSE},
		{Key: "paired_bias_diff", Val: p.PairedBiasDiff},
		{Key: "paired_bias_diff_mcse", Val: p.PairedBiasDiffMCSE},
		{Key: "abs_bias_reduction", Val: p.AbsBiasReduction},
		{Key: "rmse_diff", Val: p.RMSEDiff},
		{Key: "rmse_ratio", Val: p.RMSERatio},
		{Key: "out_of_bounds_count", Val: p.OutOfBoundsCount},
		{Key: "out_of_bounds_pct", Val: p.OutOfBoundsPct},
		{Key: "non_physical_count", Val: p.NonPhysicalCount},
		{Key: "non_physical_pct", Val: p.NonPhysicalPct},
		{Key: "wald_coverage_95", Val: p.WaldCoverage95},
		{Key: "wald_wilson_95_lower", Val: p.WaldWilsonLower},
		{Key: "wald_wilson_95_upper", Val: p.WaldWilsonUpper},
		{Key: "wald_covered_count", Val: p.WaldCoveredCount},
		{Key: "boot_percentile_coverage_95", Val: p.BootPercentileCoverage95},
		{Key: "boot_percentile_wilson_95_lower", Val: p.BootPercentileWilsonLower},
		{Key: "boot_percentile_wilson_95_upper", Val: p.BootPercentileWilsonUpper},
		{Key: "boot_percentile_covered_count", Val: p.BootPercentileCoveredCount},
		{Key: "boot_percentile_valid_outer", Val: p.BootPercentileValidOuter},
		{Key: "boot_percentile_invalid_outer", Val: p.BootPercentileInvalidOuter},
		{Key: "boot_percentile_mean_width", Val: p.BootPercentileMeanWidth},
		{Key: "boot_percentile_median_width", Val: p.BootPercentileMedianWidth},
	}
}

func bootResultToJSON(b bootstrap.ParametricBootstrapStudyResult) jsonout.Obj {
	pObjs := make([]jsonout.Obj, 11)
	for i := 0; i < 11; i++ {
		pObjs[i] = bootParamSummaryToJSON(b.Params[i])
	}

	outerObjs := make([]jsonout.Obj, len(b.OuterReplications))
	for rIdx, rep := range b.OuterReplications {
		paramBootObjs := make([]jsonout.Obj, 4)
		for k, pIdx := range jackknife.TargetParamIndices {
			bp := rep.Params[pIdx]
			paramBootObjs[k] = jsonout.Obj{
				{Key: "name", Val: kalman.ParamNames[pIdx]},
				{Key: "outer_est", Val: bp.OuterEst},
				{Key: "boot_mean", Val: bp.BootMean},
				{Key: "boot_bias", Val: bp.BootBias},
				{Key: "boot_corrected", Val: bp.BootCorrected},
				{Key: "boot_sd", Val: bp.BootSD},
				{Key: "percentile_95_lower", Val: bp.Percentile95Lower},
				{Key: "percentile_95_upper", Val: bp.Percentile95Upper},
				{Key: "interval_valid", Val: bp.IntervalValid},
				{Key: "percentile_95_covered", Val: bp.Percentile95Covered},
				{Key: "corrected_out_of_bounds", Val: bp.CorrectedOutOfBounds},
				{Key: "corrected_non_physical", Val: bp.CorrectedNonPhysical},
			}
		}

		innerObjs := make([]jsonout.Obj, len(rep.InnerDraws))
		for bIdx, d := range rep.InnerDraws {
			hb := make([]bool, 11)
			for i := 0; i < 11; i++ {
				hb[i] = d.HitBoundByParam[i]
			}
			innerObjs[bIdx] = jsonout.Obj{
				{Key: "boot_idx", Val: d.BootIdx},
				{Key: "finite", Val: d.Finite},
				{Key: "converged", Val: d.Converged},
				{Key: "termination_reason", Val: d.TerminationReason},
				{Key: "hit_bound", Val: d.HitBound},
				{Key: "hit_bound_by_param", Val: hb},
				{Key: "best_stage", Val: d.BestStage},
				{Key: "iterations_best_stage", Val: d.IterationsBest},
				{Key: "iterations_cumulative", Val: d.IterationsCum},
				{Key: "grad_norm_opt", Val: d.GradNormOpt},
				{Key: "scaled_grad_norm_opt", Val: d.ScaledGradNormOpt},
				{Key: "est_loglik", Val: d.EstLogLik},
				{Key: "est_params", Val: d.EstParams[:]},
			}
		}

		outerObjs[rIdx] = jsonout.Obj{
			{Key: "outer_rep_idx", Val: rep.OuterRepIdx},
			{Key: "outer_est", Val: rep.OuterEst[:]},
			{Key: "outer_hessian_se", Val: rep.OuterHessianSE[:]},
			{Key: "total_boot_draws", Val: rep.TotalBootDraws},
			{Key: "finite_draws_count", Val: rep.FiniteDrawsCount},
			{Key: "invalid_draws_count", Val: rep.InvalidDrawsCount},
			{Key: "converged_count", Val: rep.ConvergedCount},
			{Key: "scaled_grad_le_1e6_count", Val: rep.ScaledGradLe1e6Count},
			{Key: "bound_hit_count", Val: rep.BoundHitCount},
			{Key: "termination_reason_counts", Val: jsonout.Obj{
				{Key: "gradient_tolerance", Val: rep.GradTolCount},
				{Key: "boundary_stationary", Val: rep.BoundaryStationaryCount},
				{Key: "line_search_failed", Val: rep.LineSearchFailedCount},
				{Key: "max_iterations", Val: rep.MaxIterCount},
				{Key: "invalid_objective_or_gradient", Val: rep.InvalidObjCount},
			}},
			{Key: "interval_valid", Val: rep.IntervalValid},
			{Key: "target_params_summary", Val: paramBootObjs},
			{Key: "inner_draws", Val: innerObjs},
		}
	}

	return jsonout.Obj{
		{Key: "horizon_years", Val: b.HorizonYears},
		{Key: "dt_years", Val: b.Dt},
		{Key: "t_steps", Val: b.TSteps},
		{Key: "n_outer_reps", Val: b.NumOuterReps},
		{Key: "n_inner_boot", Val: b.NumInnerBoot},
		{Key: "total_refits", Val: b.TotalRefits},
		{Key: "min_required_finite_per_outer", Val: b.MinRequiredFinitePerOuter},
		{Key: "min_valid_fraction_threshold", Val: b.MinValidFractionThreshold},
		{Key: "workers_used", Val: b.NumWorkers},
		{Key: "elapsed_seconds", Val: b.ElapsedSeconds},
		{Key: "aggregate_refit_diagnostics", Val: jsonout.Obj{
			{Key: "finite_refits", Val: b.TotalFiniteRefits},
			{Key: "non_finite_refits", Val: b.TotalNonFiniteRefits},
			{Key: "converged_refits", Val: b.TotalConvergedRefits},
			{Key: "scaled_grad_le_1e6_refits", Val: b.TotalScaledGradLe1e6Refits},
			{Key: "bound_hit_refits", Val: b.TotalBoundHitRefits},
			{Key: "valid_outer_intervals", Val: b.ValidOuterIntervals},
			{Key: "invalid_outer_intervals", Val: b.InvalidOuterIntervals},
			{Key: "termination_reason_counts", Val: jsonout.Obj{
				{Key: "gradient_tolerance", Val: b.TotalGradTolRefits},
				{Key: "boundary_stationary", Val: b.TotalBoundaryStationaryRefit},
				{Key: "line_search_failed", Val: b.TotalLineSearchFailedRefits},
				{Key: "max_iterations", Val: b.TotalMaxIterRefits},
				{Key: "invalid_objective_or_gradient", Val: b.TotalInvalidObjRefits},
			}},
		}},
		{Key: "parameters", Val: pObjs},
		{Key: "outer_replications", Val: outerObjs},
	}
}

func checkAllFiniteJK(jk jackknife.HorizonJackknifeResult, expectedReps int) bool {
	if len(jk.Replications) != expectedReps {
		return false
	}
	ok := func(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }
	for _, p := range jk.Params {
		if p.TotalReps != expectedReps {
			return false
		}
		for _, v := range []float64{
			p.RawMean, p.RawBias, p.RawBiasMCSE, p.RawEmpiricalSD, p.RawRMSE,
			p.JackknifeMean, p.JackknifeBias, p.JackknifeBiasMCSE, p.JackknifeEmpiricalSD, p.JackknifeRMSE,
			p.WaldCoverage95, p.WaldWilsonLower, p.WaldWilsonUpper, p.PairedBiasDiffMCSE,
		} {
			if !ok(v) {
				return false
			}
		}
	}
	return true
}

func checkAllFiniteBoot(b bootstrap.ParametricBootstrapStudyResult, expectedR, expectedB int) bool {
	if len(b.OuterReplications) != expectedR || b.TotalRefits != expectedR*expectedB {
		return false
	}
	if b.ValidOuterIntervals+b.InvalidOuterIntervals != expectedR {
		return false
	}
	ok := func(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }
	for _, p := range b.Params {
		if p.TotalOuterReps != expectedR || p.BootPercentileValidOuter+p.BootPercentileInvalidOuter != expectedR {
			return false
		}
		for _, v := range []float64{
			p.RawMean, p.RawBias, p.RawRMSE,
			p.WaldCoverage95, p.WaldWilsonLower, p.WaldWilsonUpper,
			p.BootPercentileCoverage95, p.BootPercentileWilsonLower, p.BootPercentileWilsonUpper,
		} {
			if !ok(v) {
				return false
			}
		}
	}
	return true
}

func checkSummaryAlgebra(jk20, jk80 jackknife.HorizonJackknifeResult, boot bootstrap.ParametricBootstrapStudyResult) bool {
	tol := 1e-11
	for _, jk := range []jackknife.HorizonJackknifeResult{jk20, jk80} {
		for _, p := range jk.Params {
			if math.Abs((p.RawMean-p.TrueValue)-p.RawBias) > tol {
				return false
			}
			if math.Abs((p.JackknifeMean-p.TrueValue)-p.JackknifeBias) > tol {
				return false
			}
			if math.Abs((math.Abs(p.RawBias)-math.Abs(p.JackknifeBias))-p.AbsBiasReduction) > tol {
				return false
			}
			if math.Abs((p.JackknifeRMSE-p.RawRMSE)-p.RMSEDiff) > tol {
				return false
			}
		}
	}
	for _, p := range boot.Params {
		if math.Abs((math.Abs(p.RawBias)-math.Abs(p.BootCorrectedBias))-p.AbsBiasReduction) > tol {
			return false
		}
		if math.Abs((p.BootCorrectedRMSE-p.RawRMSE)-p.RMSEDiff) > tol {
			return false
		}
		expCov := float64(p.BootPercentileCoveredCount) / float64(p.TotalOuterReps)
		if math.Abs(expCov-p.BootPercentileCoverage95) > tol {
			return false
		}
	}
	return true
}

func buildFigures(
	figDir string,
	jk20, jk80 jackknife.HorizonJackknifeResult,
	boot bootstrap.ParametricBootstrapStudyResult,
	exampleRepIdx int,
) error {
	if err := os.MkdirAll(figDir, 0o755); err != nil {
		return err
	}

	// Grafik 1: Bias & RMSE Relatif (terhadap |nilai benar|) MLE mentah vs Jackknife pada 20y dan 80y
	cats1 := make([]string, 0, 8)
	rawBiasRel := make([]float64, 0, 8)
	jkBiasRel := make([]float64, 0, 8)
	rawRMSERel := make([]float64, 0, 8)
	jkRMSERel := make([]float64, 0, 8)
	for _, hpair := range []struct {
		lbl string
		res jackknife.HorizonJackknifeResult
	}{{"20y", jk20}, {"80y", jk80}} {
		for _, idx := range jackknife.TargetParamIndices {
			p := hpair.res.Params[idx]
			denom := math.Abs(p.TrueValue)
			cats1 = append(cats1, hpair.lbl+":"+p.Name)
			rawBiasRel = append(rawBiasRel, p.RawBias/denom)
			jkBiasRel = append(jkBiasRel, p.JackknifeBias/denom)
			rawRMSERel = append(rawRMSERel, p.RawRMSE/denom)
			jkRMSERel = append(jkRMSERel, p.JackknifeRMSE/denom)
		}
	}
	if err := svg.WriteGroupedBarChart(filepath.Join(figDir, "day004c_jackknife-vs-mle-bias-rmse.svg"), svg.ChartOptions{
		Title:    "MLE Mentah vs Jackknife Temporal: Bias & RMSE Relatif (20y vs 80y)",
		XLabel:   "Horizon & Parameter Sasaran (Dinormalisasi terhadap |Nilai Benar|)",
		YLabel:   "Rasio Bias atau RMSE terhadap |Nilai Benar|",
		ZeroLine: true,
	}, cats1, []svg.BarSeries{
		{Label: "Raw MLE Bias / |true|", Values: rawBiasRel, Color: svg.Palette[0]},
		{Label: "Jackknife Bias / |true|", Values: jkBiasRel, Color: svg.Palette[1]},
		{Label: "Raw MLE RMSE / |true|", Values: rawRMSERel, Color: svg.Palette[2]},
		{Label: "Jackknife RMSE / |true|", Values: jkRMSERel, Color: svg.Palette[3]},
	}); err != nil {
		return err
	}

	// Grafik 2: Perbandingan Bias Relatif (+- 1 MCSE) antara Raw MLE, Jackknife, & Bootstrap pada H=20y
	cats2 := make([]string, 4)
	rVals, rLo, rHi := make([]float64, 4), make([]float64, 4), make([]float64, 4)
	jVals, jLo, jHi := make([]float64, 4), make([]float64, 4), make([]float64, 4)
	bVals, bLo, bHi := make([]float64, 4), make([]float64, 4), make([]float64, 4)
	for k, idx := range jackknife.TargetParamIndices {
		pj := jk20.Params[idx]
		pb := boot.Params[idx]
		denom := math.Abs(pj.TrueValue)
		cats2[k] = pj.Name
		rVals[k] = pj.RawBias / denom
		rLo[k], rHi[k] = (pj.RawBias-pj.RawBiasMCSE)/denom, (pj.RawBias+pj.RawBiasMCSE)/denom
		jVals[k] = pj.JackknifeBias / denom
		jLo[k], jHi[k] = (pj.JackknifeBias-pj.JackknifeBiasMCSE)/denom, (pj.JackknifeBias+pj.JackknifeBiasMCSE)/denom
		bVals[k] = pb.BootCorrectedBias / denom
		bLo[k], bHi[k] = (pb.BootCorrectedBias-pb.BootCorrectedBiasMCSE)/denom, (pb.BootCorrectedBias+pb.BootCorrectedBiasMCSE)/denom
	}
	if err := svg.WriteCategoryPointChart(filepath.Join(figDir, "day004c_bootstrap-vs-jackknife-bias-20y.svg"), svg.ChartOptions{
		Title:    "Perbandingan Bias Relatif pada H=20y (MLE Mentah vs Jackknife vs Bootstrap)",
		XLabel:   "Parameter Sasaran (Berpasangan pada R Replikasi Luar yang Sama)",
		YLabel:   "Bias / |Nilai Benar| (+/- 1 MCSE)",
		ZeroLine: true,
	}, cats2, []svg.CategoryPointSeries{
		{Label: "Raw MLE 20y", Values: rVals, Lower: rLo, Upper: rHi, Color: svg.Palette[0]},
		{Label: "Temporal Jackknife 20y", Values: jVals, Lower: jLo, Upper: jHi, Color: svg.Palette[1]},
		{Label: "Bootstrap Corrected 20y", Values: bVals, Lower: bLo, Upper: bHi, Color: svg.Palette[2]},
	}); err != nil {
		return err
	}

	// Grafik 3: Cakupan Interval 95% Wald vs Bootstrap Persentil + Wilson 95% CI
	w20V, w20L, w20U := make([]float64, 4), make([]float64, 4), make([]float64, 4)
	b20V, b20L, b20U := make([]float64, 4), make([]float64, 4), make([]float64, 4)
	w80V, w80L, w80U := make([]float64, 4), make([]float64, 4), make([]float64, 4)
	for k, idx := range jackknife.TargetParamIndices {
		pb := boot.Params[idx]
		p80 := jk80.Params[idx]
		w20V[k], w20L[k], w20U[k] = pb.WaldCoverage95, pb.WaldWilsonLower, pb.WaldWilsonUpper
		b20V[k], b20L[k], b20U[k] = pb.BootPercentileCoverage95, pb.BootPercentileWilsonLower, pb.BootPercentileWilsonUpper
		w80V[k], w80L[k], w80U[k] = p80.WaldCoverage95, p80.WaldWilsonLower, p80.WaldWilsonUpper
	}
	y0, y1 := 0.0, 1.05
	if err := svg.WriteCategoryPointChart(filepath.Join(figDir, "day004c_wald-vs-bootstrap-coverage.svg"), svg.ChartOptions{
		Title:        "Cakupan Interval 95% Wald vs Bootstrap Persentil (+ Wilson 95% CI)",
		XLabel:       "Parameter Sasaran (Denominator = Seluruh R Replikasi Luar)",
		YLabel:       "Proporsi Cakupan 95% (+/- Wilson 95% CI)",
		YMinOverride: &y0,
		YMaxOverride: &y1,
		LegendBottom: true,
		HRefLines: []svg.HRefLine{
			{Value: 0.95, Label: "Nominal 95%", Color: "#111111"},
		},
	}, cats2, []svg.CategoryPointSeries{
		{Label: "Wald Hessian 20y", Values: w20V, Lower: w20L, Upper: w20U, Color: svg.Palette[0]},
		{Label: "Bootstrap Persentil 20y", Values: b20V, Lower: b20L, Upper: b20U, Color: svg.Palette[2]},
		{Label: "Wald Hessian 80y", Values: w80V, Lower: w80L, Upper: w80U, Color: svg.Palette[4]},
	}); err != nil {
		return err
	}

	// Grafik 4: Distribusi Bootstrap untuk kp11 dan theta1 pada replikasi luar pilihan pra-run (rep_idx = 0)
	exRep := boot.OuterReplications[exampleRepIdx]
	makeHistPanel := func(paramIdx int) svg.HistogramPanel {
		draws := make([]float64, 0, len(exRep.InnerDraws))
		for _, d := range exRep.InnerDraws {
			if d.Finite {
				draws = append(draws, d.EstParams[paramIdx])
			}
		}
		bp := exRep.Params[paramIdx]
		return svg.HistogramPanel{
			Subtitle: fmt.Sprintf("%s (Outer rep %d, %d/%d finite)", kalman.ParamNames[paramIdx], exampleRepIdx, bp.FiniteCount, len(exRep.InnerDraws)),
			XLabel:   fmt.Sprintf("%s | Persentil CI95=[%.4f, %.4f]", kalman.ParamNames[paramIdx], bp.Percentile95Lower, bp.Percentile95Upper),
			Values:   draws,
			CILower:  bp.Percentile95Lower,
			CIUpper:  bp.Percentile95Upper,
			ShowCI:   bp.IntervalValid,
			Bins:     14,
			Markers: []svg.VMarker{
				{Value: bp.TrueValue, Label: fmt.Sprintf("True=%.4f", bp.TrueValue), Color: "#2ca02c", Dashed: true},
				{Value: bp.OuterEst, Label: fmt.Sprintf("OuterMLE=%.4f", bp.OuterEst), Color: "#d62728"},
				{Value: bp.BootCorrected, Label: fmt.Sprintf("BootCorr=%.4f", bp.BootCorrected), Color: "#9467bd", Dashed: true},
			},
		}
	}
	leftPanel := makeHistPanel(1)
	rightPanel := makeHistPanel(4)
	if err := svg.WriteDualHistogramSVG(
		filepath.Join(figDir, "day004c_bootstrap-distributions-rep0.svg"),
		fmt.Sprintf("Distribusi Bootstrap Parametrik 20y (B=%d) pada Replikasi Pra-Run rep_idx=%d", boot.NumInnerBoot, exampleRepIdx),
		leftPanel,
		rightPanel,
	); err != nil {
		return err
	}

	return nil
}

func main() {
	root := findRepoRoot()
	quick := flag.Bool("quick", false, "jalankan mode cepat pengembangan ke direktori sementara tanpa menimpa keluaran penuh")
	workersFlag := flag.Int("workers", 0, "jumlah worker goroutine (0 = min(8, NumCPU))")
	repsFlag := flag.Int("r", 0, "override jumlah replikasi luar R (bawaan: 50 penuh, 3 -quick)")
	bootFlag := flag.Int("b", 0, "override jumlah bootstrap dalam B per data luar (bawaan: 100 penuh, 6 -quick)")
	outFlag := flag.String("out", "", "override berkas keluaran JSON")
	figDirFlag := flag.String("fig-dir", "", "override direktori keluaran SVG")
	flag.Parse()

	nWorkers := *workersFlag
	if nWorkers <= 0 {
		nWorkers = parallel.DefaultWorkers()
	}
	nRep := 50
	nBoot := 100
	if *quick {
		nRep = 3
		nBoot = 6
	}
	if *repsFlag > 0 {
		nRep = *repsFlag
	}
	if *bootFlag > 0 {
		nBoot = *bootFlag
	}

	outPath := filepath.Join(root, "outputs", "day004c", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	if *quick && *outFlag == "" {
		outPath = filepath.Join(os.TempDir(), "day004c_quick_results.json")
	}
	if *quick && *figDirFlag == "" {
		figDir = filepath.Join(os.TempDir(), "day004c_quick_figures")
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

	// Verifikasi aljabar jackknife c/H
	algErr20 := math.Abs(jackknife.CombineSplit(0.05+4.0/20.0, 0.05+4.0/10.0, 0.05+4.0/10.0) - 0.05)
	algErr80 := math.Abs(jackknife.CombineSplit(0.05+4.0/80.0, 0.05+4.0/40.0, 0.05+4.0/40.0) - 0.05)

	// Verifikasi determinisme 1 vs nWorkers untuk jackknife dan bootstrap pada fixture cepat
	detPaths1 := jackknife.GenerateSharedPrefixOuterPaths(trueSpec, kalman.StandardTenors, 80, kalman.WeeklyDt, 2, 1, rng.MainSeed, 0xD001)
	detPathsN := jackknife.GenerateSharedPrefixOuterPaths(trueSpec, kalman.StandardTenors, 80, kalman.WeeklyDt, 2, nWorkers, rng.MainSeed, 0xD001)
	detJK1 := jackknife.RunJackknifeHorizonStudy("det_jk", trueSpec, kalman.StandardTenors, detPaths1, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, 1)
	detJKN := jackknife.RunJackknifeHorizonStudy("det_jk", trueSpec, kalman.StandardTenors, detPathsN, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, nWorkers)
	detOuter1 := []kalman.MLEResult{detJK1.Replications[0].Full, detJK1.Replications[1].Full}
	detOuterN := []kalman.MLEResult{detJKN.Replications[0].Full, detJKN.Replications[1].Full}
	detB1 := bootstrap.RunParametricBootstrapStudy(trueSpec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, detOuter1, 4, 1, rng.MainSeed, 0xD002)
	detBN := bootstrap.RunParametricBootstrapStudy(trueSpec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, detOuterN, 4, nWorkers, rng.MainSeed, 0xD002)

	workerDetPassed := true
	for r := 0; r < 2; r++ {
		if math.Float64bits(detJK1.Replications[r].Full.EstLogLik) != math.Float64bits(detJKN.Replications[r].Full.EstLogLik) {
			workerDetPassed = false
		}
		for i := 0; i < 11; i++ {
			if math.Float64bits(detJK1.Replications[r].JackknifeEst[i]) != math.Float64bits(detJKN.Replications[r].JackknifeEst[i]) {
				workerDetPassed = false
			}
			if math.Float64bits(detB1.Params[i].BootCorrectedMean) != math.Float64bits(detBN.Params[i].BootCorrectedMean) {
				workerDetPassed = false
			}
		}
		for b := 0; b < 4; b++ {
			if math.Float64bits(detB1.OuterReplications[r].InnerDraws[b].EstLogLik) != math.Float64bits(detBN.OuterReplications[r].InnerDraws[b].EstLogLik) {
				workerDetPassed = false
			}
		}
	}
	elapsedPart1 := time.Since(tPart1).Seconds()

	// Bagian 2: Jalur 80 tahun prefiks bersama + Jackknife Temporal pada H=20y (T=1040) dan H=80y (T=4160)
	tPath0 := time.Now()
	fmt.Printf("[day004c] Membangkitkan R=%d lintasan prefiks-bersama 80y (T=4160, workers=%d)...\n", nRep, nWorkers)
	outerPaths := jackknife.GenerateSharedPrefixOuterPaths(trueSpec, kalman.StandardTenors, 4160, kalman.WeeklyDt, nRep, nWorkers, rng.MainSeed, 0x4C01)
	elapsedPaths := time.Since(tPath0).Seconds()

	fmt.Printf("[day004c] Menjalankan jackknife temporal H=20y (T=1040, half=520, R=%d)...\n", nRep)
	jk20 := jackknife.RunJackknifeHorizonStudy("jackknife_20y_weekly", trueSpec, kalman.StandardTenors, outerPaths, 20.0, kalman.WeeklyDt, 1040, nWorkers)

	fmt.Printf("[day004c] Menjalankan jackknife temporal H=80y (T=4160, half=2080, R=%d)...\n", nRep)
	jk80 := jackknife.RunJackknifeHorizonStudy("jackknife_80y_weekly", trueSpec, kalman.StandardTenors, outerPaths, 80.0, kalman.WeeklyDt, 4160, nWorkers)

	// Bagian 3: Bootstrap parametrik 20y (R replikasi luar dari jk20.Full, masing-masing B refit)
	outerFits20 := make([]kalman.MLEResult, nRep)
	for r := 0; r < nRep; r++ {
		outerFits20[r] = jk20.Replications[r].Full
	}
	fmt.Printf("[day004c] Menjalankan bootstrap parametrik 20y (R=%d, B=%d -> %d refit, workers=%d)...\n", nRep, nBoot, nRep*nBoot, nWorkers)
	boot20 := bootstrap.RunParametricBootstrapStudy(trueSpec, kalman.StandardTenors, 20.0, kalman.WeeklyDt, 1040, outerFits20, nBoot, nWorkers, rng.MainSeed, 0x4C02)

	// Hitung galat gradien Richardson terburuk sepanjang seluruh segmen luar
	maxGradErrInit, maxGradErrTrue := 0.0, 0.0
	for _, jk := range []jackknife.HorizonJackknifeResult{jk20, jk80} {
		for _, rep := range jk.Replications {
			for _, seg := range []kalman.MLEResult{rep.Full, rep.Half1, rep.Half2} {
				if seg.GradRelErrInit > maxGradErrInit {
					maxGradErrInit = seg.GradRelErrInit
				}
				if seg.GradRelErrTrue > maxGradErrTrue {
					maxGradErrTrue = seg.GradRelErrTrue
				}
			}
		}
	}

	// Grafik SVG dengan contoh replikasi luar ditetapkan pra-run = rep_idx 0
	const preSpecifiedExampleRep = 0
	tSvg0 := time.Now()
	if err := buildFigures(figDir, jk20, jk80, boot20, preSpecifiedExampleRep); err != nil {
		fmt.Fprintf(os.Stderr, "buildFigures error: %v\n", err)
		os.Exit(1)
	}
	elapsedSvg := time.Since(tSvg0).Seconds()
	elapsedTotal := time.Since(totalT0).Seconds()

	// Evaluasi 5 kriteria teknis pra-spesifikasi (tanpa menjadikan penurunan bias/RMSE atau kenaikan cakupan sebagai syarat lulus)
	c1AnchorScore := anchorRes.LogLikAbsDiff < 1e-9 &&
		anchorRes.MaxSmoothMeanAbsDiff < 1e-9 &&
		anchorRes.MaxSmoothCovAbsDiff < 1e-9 &&
		fastAncDiff < 1e-9 &&
		maxGradErrInit < 1e-5 &&
		maxGradErrTrue < 1e-5 &&
		algErr20 < 1e-13 &&
		algErr80 < 1e-13

	c2Complete := checkAllFiniteJK(jk20, nRep) && checkAllFiniteJK(jk80, nRep) && checkAllFiniteBoot(boot20, nRep, nBoot)
	c3SummaryAlg := checkSummaryAlgebra(jk20, jk80, boot20)
	c4Determinism := workerDetPassed
	c5Runtime := elapsedTotal < 600.0 || *quick
	allTechPassed := c1AnchorScore && c2Complete && c3SummaryAlg && c4Determinism && c5Runtime

	statusStr := func(ok bool) string {
		if ok {
			return "tercapai"
		}
		return "belum tercapai"
	}

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day004c"},
			{Key: "date_wita", Val: "2026-11-22"},
			{Key: "topic", Val: "AFNS — koreksi bias mean-reversion dan inferensi sampel terbatas dengan jackknife serta bootstrap parametrik"},
			{Key: "quick_mode", Val: *quick},
			{Key: "main_seed", Val: rng.MainSeed},
			{Key: "n_rep", Val: nRep},
			{Key: "n_boot", Val: nBoot},
			{Key: "workers_used", Val: nWorkers},
			{Key: "pre_specified_example_rep_idx", Val: preSpecifiedExampleRep},
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
				{Key: "outer_path_generation", Val: elapsedPaths},
				{Key: "jackknife_20y_weekly", Val: jk20.ElapsedSeconds},
				{Key: "jackknife_80y_weekly", Val: jk80.ElapsedSeconds},
				{Key: "bootstrap_20y_parametric", Val: boot20.ElapsedSeconds},
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
			{Key: "jackknife_algebraic_c_over_h_check", Val: jsonout.Obj{
				{Key: "abs_err_h20", Val: algErr20},
				{Key: "abs_err_h80", Val: algErr80},
			}},
			{Key: "worker_determinism_fixture", Val: jsonout.Obj{
				{Key: "fixture_t_steps", Val: 80},
				{Key: "fixture_outer_r", Val: 2},
				{Key: "fixture_boot_b", Val: 4},
				{Key: "workers_tested", Val: []int{1, nWorkers}},
				{Key: "bit_identical_1_vs_n", Val: workerDetPassed},
			}},
		}},
		{Key: "jackknife_studies", Val: jsonout.Obj{
			{Key: "longspan_20y_weekly", Val: jackknifeResultToJSON(jk20)},
			{Key: "longspan_80y_weekly", Val: jackknifeResultToJSON(jk80)},
		}},
		{Key: "bootstrap_study_20y_weekly", Val: bootResultToJSON(boot20)},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_exact_anchor_score_and_algebraic_check", Val: c1AnchorScore},
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
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		fmt.Fprintf(os.Stderr, "jsonout.WriteFile error: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("day004c selesai dalam %.2fs | workers=%d | c1=%v c2=%v c3=%v c4=%v c5=%v | all_technical_passed=%v\n",
		elapsedTotal, nWorkers, c1AnchorScore, c2Complete, c3SummaryAlg, c4Determinism, c5Runtime, allTechPassed)
	fmt.Printf("Output JSON: %s | Figures dir: %s\n", outPath, figDir)
}
