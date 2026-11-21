package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/ak"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/closedform"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/is"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/jsonout"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/lundberg"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/panjer"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/parallel"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/process"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/reinsurance"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/stats"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/svg"
)

const (
	nChunks       = 32
	nPathsMoments = 200_000
	nPathsDeficit = 200_000
	ksNullReps    = 2_000
	nPathsCrudeMC = 1_000_000
	nPathsIS      = 100_000
	nRepsAK       = 200_000
	nPathsHeavyMC = 200_000
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

func main() {
	tTotal := time.Now()
	root := findRepoRoot()
	outPath := filepath.Join(root, "outputs", "day003f", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	nWorkers := parallel.DefaultWorkers()

	lightDists := claims.LightTailedDists()
	heavyDists := claims.HeavyTailedDists()
	allDists := claims.AllDists()

	// =========================================================================
	// BAGIAN 1: Bentuk Tertutup, Eksponen Lundberg, Momen & Defisit Saat Ruin
	// =========================================================================
	t1 := time.Now()
	uClosedGrid := []float64{0, 1, 5, 10, 20, 50, 100}
	closedObjs := make([]jsonout.Obj, 0, len(lightDists))
	lundbergMap := make(map[claims.Kind]lundberg.Result)
	closedMap := make(map[claims.Kind]closedform.Formula)
	part1ClosedOK := true

	for _, d := range lightDists {
		m1, m2, cv := d.Moments()
		resL, err := lundberg.Solve(d)
		if err != nil {
			panic(err)
		}
		cf, err := closedform.Derive(d)
		if err != nil {
			panic(err)
		}
		lundbergMap[d.Kind] = resL
		closedMap[d.Kind] = cf

		rDiff := math.Abs(resL.R - cf.R1)
		cDiff := math.Abs(resL.C - cf.C1)
		if d.Kind == claims.KindExponential {
			rExactA := claims.Theta / ((1.0 + claims.Theta) * claims.Mu)
			if math.Abs(resL.R-rExactA) >= 1e-12 {
				part1ClosedOK = false
			}
		} else if rDiff >= 1e-10 {
			part1ClosedOK = false
		}
		if math.Abs(cf.Eval(0)-claims.Rho) > 1e-14 || cDiff > 1e-12 {
			part1ClosedOK = false
		}

		pts := make([]jsonout.Obj, 0, len(uClosedGrid))
		prev := 1.0
		lundbergOK := true
		for _, u := range uClosedGrid {
			psi := cf.Eval(u)
			bound := lundberg.Bound(resL.R, u)
			clAsym := resL.C * math.Exp(-resL.R*u)
			ratio := psi / clAsym
			if psi > bound+1e-15 || psi >= prev {
				lundbergOK = false
				part1ClosedOK = false
			}
			prev = psi
			pts = append(pts, jsonout.Obj{
				{Key: "u", Val: u},
				{Key: "psi_exact", Val: psi},
				{Key: "lundberg_bound_exp_neg_Ru", Val: bound},
				{Key: "cramer_lundberg_asymp", Val: clAsym},
				{Key: "ratio_exact_over_cl", Val: ratio},
				{Key: "lundberg_bound_satisfied", Val: psi <= bound+1e-15},
			})
		}

		closedObjs = append(closedObjs, jsonout.Obj{
			{Key: "dist_code", Val: d.Code},
			{Key: "dist_name", Val: d.Name},
			{Key: "mean", Val: m1},
			{Key: "second_moment", Val: m2},
			{Key: "cv", Val: cv},
			{Key: "lundberg_R_brent", Val: resL.R},
			{Key: "laplace_R1", Val: cf.R1},
			{Key: "laplace_R2", Val: cf.R2},
			{Key: "abs_err_R", Val: rDiff},
			{Key: "cramer_lundberg_C", Val: resL.C},
			{Key: "residue_C1", Val: cf.C1},
			{Key: "residue_C2", Val: cf.C2},
			{Key: "abs_err_C", Val: cDiff},
			{Key: "psi_zero", Val: cf.Eval(0)},
			{Key: "all_bounds_and_monotone_ok", Val: lundbergOK},
			{Key: "points", Val: pts},
		})
	}

	momRows := process.VerifyMoments(allDists, 10.0, 10.0, nPathsMoments, nChunks, nWorkers, rng.MainSeed)
	momObjs := make([]jsonout.Obj, 0, len(momRows))
	momentsLightOK := true
	momentsMeanAllOK := true
	for _, r := range momRows {
		if r.ZN > 3.5 || r.ZUt > 3.5 {
			momentsMeanAllOK = false
		}
		if r.DistCode != "D" && math.Abs(r.VarStDiff) > 0.05 {
			momentsLightOK = false
		}
		momObjs = append(momObjs, jsonout.Obj{
			{Key: "dist_code", Val: r.DistCode},
			{Key: "dist_name", Val: r.DistName},
			{Key: "theo_E_Nt", Val: r.TheoEN},
			{Key: "sim_E_Nt", Val: r.SimEN},
			{Key: "sim_E_Nt_se", Val: r.SimENSE},
			{Key: "z_Nt", Val: r.ZN},
			{Key: "theo_Var_Nt", Val: r.TheoVarN},
			{Key: "sim_Var_Nt", Val: r.SimVarN},
			{Key: "theo_E_Ut", Val: r.TheoEUt},
			{Key: "sim_E_Ut", Val: r.SimEUt},
			{Key: "sim_E_Ut_se", Val: r.SimEUtSE},
			{Key: "z_Ut", Val: r.ZUt},
			{Key: "theo_Var_St", Val: r.TheoVarSt},
			{Key: "sim_Var_St", Val: r.SimVarSt},
			{Key: "var_St_rel_diff", Val: r.VarStDiff},
		})
	}

	defRows := process.VerifyExponentialDeficit([]float64{2.0, 5.0}, nPathsDeficit, ksNullReps, nChunks, nWorkers, rng.MainSeed)
	defObjs := make([]jsonout.Obj, 0, len(defRows))
	deficitOK := true
	for _, r := range defRows {
		meanZ := math.Abs(r.DeficitMean-1.0) / r.DeficitMeanSE
		if r.PsiZScore > 3.5 || meanZ > 3.5 || !r.KSPass {
			deficitOK = false
		}
		defObjs = append(defObjs, jsonout.Obj{
			{Key: "u", Val: r.U},
			{Key: "psi_exact", Val: r.PsiExact},
			{Key: "psi_mc", Val: r.PsiMC},
			{Key: "psi_mc_se", Val: r.PsiMCSE},
			{Key: "psi_z_score", Val: r.PsiZScore},
			{Key: "n_ruined", Val: r.NRuined},
			{Key: "deficit_mean", Val: r.DeficitMean},
			{Key: "deficit_mean_se", Val: r.DeficitMeanSE},
			{Key: "deficit_mean_z", Val: meanZ},
			{Key: "deficit_sd", Val: r.DeficitSD},
			{Key: "ks_stat", Val: r.KSStat},
			{Key: "ks_crit_95_simulated_null", Val: r.KSCrit95},
			{Key: "ks_p_value", Val: r.KSPValue},
			{Key: "ks_pass", Val: r.KSPass},
		})
	}
	elapsedP1 := time.Since(t1).Seconds()

	// =========================================================================
	// BAGIAN 2: Rekursi Panjer/Beekman, Batas Diskretisasi & Asimtotik
	// =========================================================================
	t2 := time.Now()
	hGrid := []float64{0.1, 0.05, 0.02, 0.01}
	uPanjerLight := []float64{0, 1, 5, 10, 20, 50}
	uPanjerHeavy := []float64{0, 1, 5, 10, 20, 50, 100, 200}

	panjerLightObjs := make([]jsonout.Obj, 0, len(lightDists))
	panjerHeavyObjs := make([]jsonout.Obj, 0, len(heavyDists))
	panjerBracketAllOK := true
	panjerOrderAllOK := true
	panjerHeavyStableOK := true

	// Untuk grafik SVG
	orderSeries := make([]svg.Series, 0, len(allDists))
	panjerHeavyMap := make(map[claims.Kind]panjer.Series)

	for idx, d := range lightDists {
		rows := panjer.EvaluateGrid(d, hGrid, uPanjerLight)
		rowObjs := make([]jsonout.Obj, 0, len(rows))
		logH := make([]float64, 0, len(hGrid))
		logW10 := make([]float64, 0, len(hGrid))
		w10Vals := make([]float64, 0, len(hGrid))
		allBr := true
		for _, r := range rows {
			if !r.Bracketed {
				allBr = false
				panjerBracketAllOK = false
			}
			if r.U == 10.0 {
				logH = append(logH, math.Log(r.H))
				logW10 = append(logW10, math.Log(r.Width))
				w10Vals = append(w10Vals, r.Width)
			}
			rowObjs = append(rowObjs, jsonout.Obj{
				{Key: "h", Val: r.H},
				{Key: "u", Val: r.U},
				{Key: "psi_low", Val: r.PsiLow},
				{Key: "psi_mid", Val: r.PsiMid},
				{Key: "psi_high", Val: r.PsiHigh},
				{Key: "width", Val: r.Width},
				{Key: "psi_exact", Val: r.PsiExact},
				{Key: "bracketed", Val: r.Bracketed},
				{Key: "mid_abs_err", Val: r.MidAbsErr},
			})
		}
		slope := stats.LinearSlope(logH, logW10)
		if math.Abs(slope-1.0) > 0.08 {
			panjerOrderAllOK = false
		}
		orderSeries = append(orderSeries, svg.Series{
			Label: fmt.Sprintf("%s (slope=%.3f)", d.Code, slope),
			X:     hGrid,
			Y:     w10Vals,
			Color: svg.Palette[idx],
		})
		panjerLightObjs = append(panjerLightObjs, jsonout.Obj{
			{Key: "dist_code", Val: d.Code},
			{Key: "dist_name", Val: d.Name},
			{Key: "all_points_bracketed", Val: allBr},
			{Key: "width_order_u10_slope", Val: slope},
			{Key: "grid_rows", Val: rowObjs},
		})
	}

	for idx, d := range heavyDists {
		rows := panjer.EvaluateGrid(d, hGrid, uPanjerHeavy)
		// Simpan seri beresolusi terhalus h = 0.01 untuk dipakai AK di Bagian 4.
		panjerHeavyMap[d.Kind] = panjer.ComputeBounds(claims.Rho, d.LadderSurvival, 200.0, 0.01)

		rowObjs := make([]jsonout.Obj, 0, len(rows))
		logH := make([]float64, 0, len(hGrid))
		logW10 := make([]float64, 0, len(hGrid))
		w10Vals := make([]float64, 0, len(hGrid))
		for _, r := range rows {
			if r.U == 10.0 {
				logH = append(logH, math.Log(r.H))
				logW10 = append(logW10, math.Log(r.Width))
				w10Vals = append(w10Vals, r.Width)
			}
			if r.H == 0.01 && (r.PsiLow <= 0 || r.PsiLow > r.PsiHigh) {
				panjerHeavyStableOK = false
			}
			rowObjs = append(rowObjs, jsonout.Obj{
				{Key: "h", Val: r.H},
				{Key: "u", Val: r.U},
				{Key: "psi_low", Val: r.PsiLow},
				{Key: "psi_mid", Val: r.PsiMid},
				{Key: "psi_high", Val: r.PsiHigh},
				{Key: "width", Val: r.Width},
				{Key: "low_cdf_sum_formula", Val: r.LowCdfSum},
				{Key: "cdf_sum_vs_tail_diff", Val: math.Abs(r.LowCdfSum - r.PsiLow)},
				{Key: "subexp_asymptote", Val: r.Asymptote},
				{Key: "panjer_mid_over_asymp", Val: r.RatioAsym},
			})
		}
		slope := stats.LinearSlope(logH, logW10)
		orderSeries = append(orderSeries, svg.Series{
			Label:  fmt.Sprintf("%s (slope=%.3f)", d.Code, slope),
			X:      hGrid,
			Y:      w10Vals,
			Color:  svg.Palette[3+idx],
			Dashed: true,
		})
		panjerHeavyObjs = append(panjerHeavyObjs, jsonout.Obj{
			{Key: "dist_code", Val: d.Code},
			{Key: "dist_name", Val: d.Name},
			{Key: "width_order_u10_slope", Val: slope},
			{Key: "grid_rows", Val: rowObjs},
		})
	}
	elapsedP2 := time.Since(t2).Seconds()

	// =========================================================================
	// BAGIAN 3: Monte Carlo Biasa & Importance Sampling Eksponensial
	// =========================================================================
	t3 := time.Now()
	uCrude := []float64{0, 5, 10, 20}
	crudeRows := process.SimulateCrudeMC(lightDists, uCrude, nPathsCrudeMC, nChunks, nWorkers, rng.MainSeed)
	crudeObjs := make([]jsonout.Obj, 0, len(crudeRows))
	crudeMap := make(map[string]process.CrudeMCRow)
	crudeOK := true
	for _, r := range crudeRows {
		if r.ZScore > 3.5 || !(r.PsiT50 <= r.PsiT200 && r.PsiT200 <= r.PsiInfMC) {
			crudeOK = false
		}
		crudeMap[fmt.Sprintf("%s_%.0f", r.DistCode, r.U)] = r
		crudeObjs = append(crudeObjs, jsonout.Obj{
			{Key: "dist_code", Val: r.DistCode},
			{Key: "dist_name", Val: r.DistName},
			{Key: "u", Val: r.U},
			{Key: "psi_exact", Val: r.PsiExact},
			{Key: "psi_T50_mc", Val: r.PsiT50},
			{Key: "psi_T200_mc", Val: r.PsiT200},
			{Key: "psi_inf_mc", Val: r.PsiInfMC},
			{Key: "psi_inf_se", Val: r.PsiInfSE},
			{Key: "rel_se", Val: r.RelSE},
			{Key: "z_score", Val: r.ZScore},
			{Key: "elapsed_seconds_dist", Val: r.ElapsedS},
		})
	}

	uISGrid := []float64{0, 5, 10, 20, 50, 100}
	isOptObjs := make([]jsonout.Obj, 0, len(lightDists)*len(uISGrid))
	isOptOK := true
	isMap := make(map[string]is.Row)
	for _, d := range lightDists {
		rOpt := lundbergMap[d.Kind].R
		// Dapatkan waktu rata-rata per titik Crude MC pada distribusi ini untuk VRF waktu-setara.
		tCrudeRef := crudeMap[fmt.Sprintf("%s_10", d.Code)].ElapsedS / 4.0
		for _, u := range uISGrid {
			row := is.Simulate(d, u, rOpt, 1.0, nPathsIS, nChunks, nWorkers, rng.MainSeed)
			if row.ZScore > 3.5 {
				isOptOK = false
			}
			isMap[fmt.Sprintf("%s_%.0f", d.Code, u)] = row

			// VRF teoritis (varians per path): Var_crude / Var_IS
			// VRF waktu-setara: (Var_crude * T_crude_per_path) / (Var_IS * T_IS_per_path)
			vrfTime := 0.0
			if row.ElapsedS > 0 && tCrudeRef > 0 {
				tCrudePerPath := tCrudeRef / float64(nPathsCrudeMC)
				tISPerPath := row.ElapsedS / float64(nPathsIS)
				vrfTime = row.VRFTheory * (tCrudePerPath / tISPerPath)
			}
			isOptObjs = append(isOptObjs, jsonout.Obj{
				{Key: "dist_code", Val: row.DistCode},
				{Key: "dist_name", Val: row.DistName},
				{Key: "u", Val: row.U},
				{Key: "r_tilt", Val: row.RTilt},
				{Key: "psi_exact", Val: row.PsiExact},
				{Key: "psi_is", Val: row.PsiIS},
				{Key: "se_is", Val: row.SEIS},
				{Key: "rel_se_is", Val: row.RelSEIS},
				{Key: "z_score", Val: row.ZScore},
				{Key: "sd_is", Val: row.SDIS},
				{Key: "vrf_per_sample", Val: row.VRFTheory},
				{Key: "vrf_time_adjusted", Val: vrfTime},
				{Key: "elapsed_seconds", Val: row.ElapsedS},
			})
		}
	}

	misspecObjs := make([]jsonout.Obj, 0, len(lightDists)*6)
	misspecOK := true
	for _, d := range lightDists {
		rOpt := lundbergMap[d.Kind].R
		for _, u := range []float64{10, 20} {
			r07 := is.Simulate(d, u, 0.7*rOpt, 0.7, nPathsIS, nChunks, nWorkers, rng.MainSeed)
			r10 := is.Simulate(d, u, 1.0*rOpt, 1.0, nPathsIS, nChunks, nWorkers, rng.MainSeed)
			r13 := is.Simulate(d, u, 1.3*rOpt, 1.3, nPathsIS, nChunks, nWorkers, rng.MainSeed)
			for _, rr := range []is.Row{r07, r10, r13} {
				if rr.ZScore > 3.5 {
					misspecOK = false
				}
				misspecObjs = append(misspecObjs, jsonout.Obj{
					{Key: "dist_code", Val: rr.DistCode},
					{Key: "u", Val: rr.U},
					{Key: "tilt_ratio", Val: rr.TiltRatio},
					{Key: "r_tilt", Val: rr.RTilt},
					{Key: "psi_exact", Val: rr.PsiExact},
					{Key: "psi_is", Val: rr.PsiIS},
					{Key: "se_is", Val: rr.SEIS},
					{Key: "sd_is", Val: rr.SDIS},
					{Key: "z_score", Val: rr.ZScore},
					{Key: "sd_ratio_vs_optimal", Val: rr.SDIS / r10.SDIS},
				})
			}
			if r10.SDIS >= r07.SDIS || r10.SDIS >= r13.SDIS {
				misspecOK = false
			}
		}
	}
	elapsedP3 := time.Since(t3).Seconds()

	// =========================================================================
	// BAGIAN 4: Ekor Berat — Kegagalan IS Eksponensial & Asmussen-Kroese
	// =========================================================================
	t4 := time.Now()
	mgfDivRows := is.VerifyHeavyTailMGFFailure(heavyDists, 0.1, []float64{10, 50, 100, 200, 500})
	mgfDivObjs := make([]jsonout.Obj, 0, len(mgfDivRows))
	for _, r := range mgfDivRows {
		mgfDivObjs = append(mgfDivObjs, jsonout.Obj{
			{Key: "dist_code", Val: r.DistCode},
			{Key: "dist_name", Val: r.DistName},
			{Key: "r_test", Val: r.RTest},
			{Key: "x_max", Val: r.XMax},
			{Key: "partial_mgf", Val: r.PartialMGF},
		})
	}

	uHeavyGrid := []float64{10, 20, 50, 100, 200}
	heavyCrude := process.SimulateCrudeMC(heavyDists, uHeavyGrid, nPathsHeavyMC, nChunks, nWorkers, rng.MainSeed^0x444)
	heavyCrudeMap := make(map[string]process.CrudeMCRow)
	for _, r := range heavyCrude {
		heavyCrudeMap[fmt.Sprintf("%s_%.0f", r.DistCode, r.U)] = r
	}

	akObjs := make([]jsonout.Obj, 0, len(heavyDists)*len(uHeavyGrid))
	akMap := make(map[string]ak.Row)
	akAllOK := true
	for _, d := range heavyDists {
		pSer := panjerHeavyMap[d.Kind]
		for _, u := range uHeavyGrid {
			row := ak.Simulate(d, u, pSer, nRepsAK, nChunks, nWorkers, rng.MainSeed)
			if !row.WithinBand {
				akAllOK = false
			}
			akMap[fmt.Sprintf("%s_%.0f", d.Code, u)] = row
			cmc := heavyCrudeMap[fmt.Sprintf("%s_%.0f", d.Code, u)]
			crudeRelErrVsMid := math.NaN()
			if row.PanjerMid > 0 {
				crudeRelErrVsMid = (cmc.PsiInfMC - row.PanjerMid) / row.PanjerMid
			}
			akRelErrVsMid := (row.PsiAK - row.PanjerMid) / row.PanjerMid
			asymRelErrVsMid := (row.Asymptote - row.PanjerMid) / row.PanjerMid
			vrfTime := 0.0
			if row.ElapsedS > 0 && cmc.ElapsedS > 0 && row.VRFTheory > 0 {
				tCrudePerPath := (cmc.ElapsedS / float64(len(uHeavyGrid))) / float64(nPathsHeavyMC)
				tAKPerPath := row.ElapsedS / float64(nRepsAK)
				vrfTime = row.VRFTheory * (tCrudePerPath / tAKPerPath)
			}
			akObjs = append(akObjs, jsonout.Obj{
				{Key: "dist_code", Val: row.DistCode},
				{Key: "dist_name", Val: row.DistName},
				{Key: "u", Val: row.U},
				{Key: "panjer_low_h001", Val: row.PanjerLow},
				{Key: "panjer_mid_h001", Val: row.PanjerMid},
				{Key: "panjer_high_h001", Val: row.PanjerHigh},
				{Key: "ak_psi", Val: row.PsiAK},
				{Key: "ak_se", Val: row.SEAK},
				{Key: "ak_rel_se", Val: row.RelSEAK},
				{Key: "ak_z_vs_panjer_mid", Val: row.ZVsMid},
				{Key: "ak_within_3se_plus_band", Val: row.WithinBand},
				{Key: "vrf_per_sample", Val: row.VRFTheory},
				{Key: "vrf_time_adjusted", Val: vrfTime},
				{Key: "crude_mc_psi", Val: cmc.PsiInfMC},
				{Key: "crude_mc_se", Val: cmc.PsiInfSE},
				{Key: "subexp_asymptote", Val: row.Asymptote},
				{Key: "ak_over_asymptote", Val: row.AKOverAsym},
				{Key: "rel_err_crude_vs_mid", Val: crudeRelErrVsMid},
				{Key: "rel_err_ak_vs_mid", Val: akRelErrVsMid},
				{Key: "rel_err_asym_vs_mid", Val: asymRelErrVsMid},
				{Key: "elapsed_seconds", Val: row.ElapsedS},
			})
		}
	}
	elapsedP4 := time.Since(t4).Seconds()

	// =========================================================================
	// BAGIAN 5: Reinsurans Excess-of-Loss & Keputusan Modal Awal Minimum
	// =========================================================================
	t5 := time.Now()
	reinsDists := []claims.Dist{claims.Exponential(), claims.Hyperexp2()}
	optResults := make([]reinsurance.OptimalRetentionResult, 0, len(reinsDists))
	capAllRows := make([]reinsurance.CapitalRequirementRow, 0, 6)
	reinsObjs := make([]jsonout.Obj, 0, len(reinsDists))
	reinsOK := true
	epsGrid := []float64{1e-2, 1e-3, 1e-4}

	for _, d := range reinsDists {
		opt, err := reinsurance.FindOptimalRetention(d, 10.0, 0.01)
		if err != nil {
			panic(err)
		}
		if d.Kind == claims.KindExponential && (!opt.MaxClosedVsBand || opt.InfLimitDiff > 1e-10) {
			reinsOK = false
		}
		if opt.MinPsiMid >= opt.PsiNoReins {
			reinsOK = false
		}
		optResults = append(optResults, opt)

		caps, err := reinsurance.SolveMinimumCapital(d, opt.OptimalM, epsGrid, 0.02)
		if err != nil {
			panic(err)
		}
		capAllRows = append(capAllRows, caps...)

		scanObjs := make([]jsonout.Obj, 0, len(opt.ScanM))
		for i, m := range opt.ScanM {
			rp := reinsurance.EvaluateRetention(d, m)
			row := jsonout.Obj{
				{Key: "M", Val: m},
				{Key: "net_premium_c_M", Val: rp.NetPremiumC},
				{Key: "safety_loading_M", Val: rp.SafetyLoading},
				{Key: "psi_panjer_low", Val: opt.ScanPsiLow[i]},
				{Key: "psi_panjer_mid", Val: opt.ScanPsiMid[i]},
				{Key: "psi_panjer_high", Val: opt.ScanPsiHigh[i]},
			}
			if d.Kind == claims.KindExponential {
				row = append(row, jsonout.Field{Key: "psi_closed_form", Val: opt.ScanPsiClosed[i]})
			}
			scanObjs = append(scanObjs, row)
		}
		capObjs := make([]jsonout.Obj, 0, len(caps))
		for _, c := range caps {
			if c.UStarWithReins >= c.UStarNoReins {
				reinsOK = false
			}
			capObjs = append(capObjs, jsonout.Obj{
				{Key: "epsilon", Val: c.Epsilon},
				{Key: "u_star_no_reins", Val: c.UStarNoReins},
				{Key: "u_star_with_reins", Val: c.UStarWithReins},
				{Key: "capital_saving", Val: c.CapitalSaving},
				{Key: "saving_percent", Val: c.SavingPercent},
				{Key: "psi_check_no_reins", Val: c.PsiCheckNo},
				{Key: "psi_check_with_reins", Val: c.PsiCheckWith},
			})
		}

		reinsObjs = append(reinsObjs, jsonout.Obj{
			{Key: "dist_code", Val: opt.DistCode},
			{Key: "dist_name", Val: opt.DistName},
			{Key: "u_eval", Val: opt.U},
			{Key: "m_min_valid", Val: opt.MMin},
			{Key: "optimal_M", Val: opt.OptimalM},
			{Key: "min_psi_panjer_mid", Val: opt.MinPsiMid},
			{Key: "psi_no_reinsurance", Val: opt.PsiNoReins},
			{Key: "closed_form_within_panjer_band", Val: opt.MaxClosedVsBand},
			{Key: "m_inf_limit_abs_diff", Val: opt.InfLimitDiff},
			{Key: "scan_psi_vs_m", Val: scanObjs},
			{Key: "capital_requirements", Val: capObjs},
		})
	}
	elapsedP5 := time.Since(t5).Seconds()

	// Uji determinisme pekerja (1 vs nWorkers)
	w1 := is.Simulate(claims.Exponential(), 20.0, lundbergMap[claims.KindExponential].R, 1.0, 20000, 8, 1, rng.MainSeed)
	wn := is.Simulate(claims.Exponential(), 20.0, lundbergMap[claims.KindExponential].R, 1.0, 20000, 8, nWorkers, rng.MainSeed)
	workerBitIdentical := math.Float64bits(w1.PsiIS) == math.Float64bits(wn.PsiIS)

	// =========================================================================
	// BUAT 7 GRAFIK SVG KE outputs/figures/day003f_*.svg
	// =========================================================================
	// 1. psi(u) eksak vs batas Panjer h=0.01 (A, B, C)
	svg1Series := make([]svg.Series, 0, 6)
	uPlotLight := []float64{0, 1, 5, 10, 20, 50}
	for idx, d := range lightDists {
		ser := panjer.ComputeBounds(claims.Rho, d.LadderSurvival, 50.0, 0.01)
		cf := closedMap[d.Kind]
		exY := make([]float64, len(uPlotLight))
		loY := make([]float64, len(uPlotLight))
		hiY := make([]float64, len(uPlotLight))
		for i, u := range uPlotLight {
			exY[i] = cf.Eval(u)
			loY[i], _, hiY[i], _, _, _ = ser.At(u)
		}
		svg1Series = append(svg1Series,
			svg.Series{Label: d.Code + " Eksak", X: uPlotLight, Y: exY, Color: svg.Palette[idx]},
			svg.Series{Label: d.Code + " Panjer h=0.01 (±)", X: uPlotLight, Y: exY, Lower: loY, Upper: hiY, Color: svg.Palette[idx], Dashed: true},
		)
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_psi-exact-vs-panjer.svg"), svg.ChartOptions{
		Title:  "Probabilitas Ruin psi(u) Eksak vs Batas Panjer (h=0.01)",
		XLabel: "Modal Awal u",
		YLabel: "psi(u) (skala log)",
		LogY:   true,
	}, svg1Series); err != nil {
		panic(err)
	}

	// 2. Orde diskretisasi Panjer log-log (lebar interval vs h di u=10)
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_panjer-discretization-order.svg"), svg.ChartOptions{
		Title:         "Konvergensi Lebar Interval Panjer vs Langkah h (u=10, Log-Log)",
		XLabel:        "Langkah Diskretisasi h",
		YLabel:        "psi_high(10) - psi_low(10)",
		LogX:          true,
		LogY:          true,
		LegendTopLeft: true,
	}, orderSeries); err != nil {
		panic(err)
	}

	// 3. Rasio Cramér-Lundberg psi(u) / (C e^{-Ru}) -> 1
	ratioSeries := make([]svg.Series, 0, len(lightDists))
	for idx, d := range lightDists {
		cf := closedMap[d.Kind]
		rl := lundbergMap[d.Kind]
		ys := make([]float64, len(uClosedGrid))
		for i, u := range uClosedGrid {
			ys[i] = cf.Eval(u) / (rl.C * math.Exp(-rl.R*u))
		}
		ratioSeries = append(ratioSeries, svg.Series{
			Label: fmt.Sprintf("%s (C=%.4f, R=%.4f)", d.Code, rl.C, rl.R),
			X:     uClosedGrid,
			Y:     ys,
			Color: svg.Palette[idx],
		})
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_cramer-lundberg-ratio.svg"), svg.ChartOptions{
		Title:  "Konvergensi Asimtotik Cramer-Lundberg: psi(u) / (C exp(-Ru))",
		XLabel: "Modal Awal u",
		YLabel: "Rasio psi(u) / (C exp(-Ru))",
	}, ratioSeries); err != nil {
		panic(err)
	}

	// 4. SE relatif IS eksponensial vs MC biasa terhadap u (skala log Y)
	relSESeries := make([]svg.Series, 0, 6)
	for idx, d := range lightDists {
		isYs := make([]float64, len(uISGrid))
		crudeTheoYs := make([]float64, len(uISGrid))
		cf := closedMap[d.Kind]
		for i, u := range uISGrid {
			isYs[i] = isMap[fmt.Sprintf("%s_%.0f", d.Code, u)].RelSEIS
			psiEx := cf.Eval(u)
			// SE relatif teoritis MC biasa pada anggaran N_IS yang setara (N=100k) agar perbandingan adil terhadap u
			crudeTheoYs[i] = math.Sqrt((1.0 - psiEx) / (float64(nPathsIS) * psiEx))
		}
		relSESeries = append(relSESeries,
			svg.Series{Label: d.Code + " IS Eksponensial (R)", X: uISGrid, Y: isYs, Color: svg.Palette[idx]},
			svg.Series{Label: d.Code + " MC Biasa (N setara)", X: uISGrid, Y: crudeTheoYs, Color: svg.Palette[idx], Dashed: true},
		)
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_is-vs-mc-rel-se.svg"), svg.ChartOptions{
		Title:         "Efisiensi Logaritmik: SE Relatif IS Eksponensial vs MC Biasa",
		XLabel:        "Modal Awal u",
		YLabel:        "Galat Baku Relatif SE/psi (log)",
		LogY:          true,
		LegendTopLeft: true,
	}, relSESeries); err != nil {
		panic(err)
	}

	// 5. Ekor berat: Asmussen-Kroese vs Panjer h=0.01 vs Asimtotik Subeksponensial
	heavySeries := make([]svg.Series, 0, 6)
	for idx, d := range heavyDists {
		akY := make([]float64, len(uHeavyGrid))
		akLo := make([]float64, len(uHeavyGrid))
		akHi := make([]float64, len(uHeavyGrid))
		asymY := make([]float64, len(uHeavyGrid))
		panY := make([]float64, len(uHeavyGrid))
		for i, u := range uHeavyGrid {
			r := akMap[fmt.Sprintf("%s_%.0f", d.Code, u)]
			akY[i] = r.PsiAK
			akLo[i] = math.Max(r.PsiAK-1.96*r.SEAK, 1e-12)
			akHi[i] = r.PsiAK + 1.96*r.SEAK
			asymY[i] = r.Asymptote
			panY[i] = r.PanjerMid
		}
		heavySeries = append(heavySeries,
			svg.Series{Label: d.Code + " Asmussen-Kroese (±1.96 SE)", X: uHeavyGrid, Y: akY, Lower: akLo, Upper: akHi, Color: svg.Palette[idx]},
			svg.Series{Label: d.Code + " Panjer h=0.01", X: uHeavyGrid, Y: panY, Color: svg.Palette[idx], Dashed: true},
			svg.Series{Label: d.Code + " Asimtotik Subeksponensial", X: uHeavyGrid, Y: asymY, Color: svg.Palette[3+idx], Dashed: true},
		)
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_heavy-tail-ak-vs-asymptotic.svg"), svg.ChartOptions{
		Title:  "Ekor Berat: Asmussen-Kroese vs Panjer vs Asimtotik (u sampai 200)",
		XLabel: "Modal Awal u",
		YLabel: "psi(u) (skala log)",
		LogY:   true,
	}, heavySeries); err != nil {
		panic(err)
	}

	// 6. Reinsurans: psi_M(10) vs retensi M untuk (A) dan (C)
	reinsSeries := make([]svg.Series, 0, 3)
	for idx, opt := range optResults {
		reinsSeries = append(reinsSeries, svg.Series{
			Label: fmt.Sprintf("%s Panjer Mid (M*=%.2f, min=%.4f)", opt.DistCode, opt.OptimalM, opt.MinPsiMid),
			X:     opt.ScanM,
			Y:     opt.ScanPsiMid,
			Lower: opt.ScanPsiLow,
			Upper: opt.ScanPsiHigh,
			Color: svg.Palette[idx],
		})
		if opt.DistCode == "A" {
			reinsSeries = append(reinsSeries, svg.Series{
				Label:  "A Tertutup Eksak",
				X:      opt.ScanM,
				Y:      opt.ScanPsiClosed,
				Color:  svg.Palette[2],
				Dashed: true,
			})
		}
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_reinsurance-psi-vs-m.svg"), svg.ChartOptions{
		Title:  "Optimisasi Retensi Excess-of-Loss (xi=0.4) di Modal u=10",
		XLabel: "Batas Retensi M",
		YLabel: "Probabilitas Ruin psi_M(10)",
	}, reinsSeries); err != nil {
		panic(err)
	}

	// 7. Modal minimum u^*(epsilon) dengan dan tanpa reinsurans
	capSeries := make([]svg.Series, 0, 4)
	for idx, d := range reinsDists {
		noY := make([]float64, 0, len(epsGrid))
		wiY := make([]float64, 0, len(epsGrid))
		for _, r := range capAllRows {
			if r.DistCode == d.Code {
				noY = append(noY, r.UStarNoReins)
				wiY = append(wiY, r.UStarWithReins)
			}
		}
		capSeries = append(capSeries,
			svg.Series{Label: d.Code + " Tanpa Reinsurans", X: epsGrid, Y: noY, Color: svg.Palette[idx], Dashed: true},
			svg.Series{Label: fmt.Sprintf("%s Reinsurans M*=%.2f", d.Code, optResults[idx].OptimalM), X: epsGrid, Y: wiY, Color: svg.Palette[idx]},
		)
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day003f_capital-u-star.svg"), svg.ChartOptions{
		Title:  "Kebutuhan Modal Awal Minimum u*(epsilon) Dengan vs Tanpa Reinsurans",
		XLabel: "Toleransi Ruin epsilon (skala log)",
		YLabel: "Modal Awal Minimum u*",
		LogX:   true,
	}, capSeries); err != nil {
		panic(err)
	}

	elapsedTotal := time.Since(tTotal).Seconds()
	allPassed := part1ClosedOK && momentsMeanAllOK && momentsLightOK && deficitOK &&
		panjerBracketAllOK && panjerOrderAllOK && panjerHeavyStableOK &&
		crudeOK && isOptOK && misspecOK && akAllOK && reinsOK && workerBitIdentical

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day003f"},
			{Key: "date", Val: "2026-11-21"},
			{Key: "language", Val: "go"},
			{Key: "main_seed", Val: rng.MainSeed},
			{Key: "n_workers", Val: nWorkers},
			{Key: "n_chunks", Val: nChunks},
			{Key: "worker_invariance_1_vs_n_bit_identical", Val: workerBitIdentical},
			{Key: "n_paths_moments", Val: nPathsMoments},
			{Key: "n_paths_deficit", Val: nPathsDeficit},
			{Key: "ks_null_replicates", Val: ksNullReps},
			{Key: "n_paths_crude_mc", Val: nPathsCrudeMC},
			{Key: "n_paths_is", Val: nPathsIS},
			{Key: "n_reps_ak", Val: nRepsAK},
			{Key: "n_paths_heavy_crude", Val: nPathsHeavyMC},
			{Key: "timing_seconds", Val: jsonout.Obj{
				{Key: "part1_closed_moments_deficit", Val: elapsedP1},
				{Key: "part2_panjer_cramer_lundberg", Val: elapsedP2},
				{Key: "part3_crude_mc_and_is", Val: elapsedP3},
				{Key: "part4_heavy_tail_mgf_and_ak", Val: elapsedP4},
				{Key: "part5_reinsurance_and_capital", Val: elapsedP5},
				{Key: "total", Val: elapsedTotal},
			}},
		}},
		{Key: "model_parameters", Val: jsonout.Obj{
			{Key: "lambda", Val: claims.Lambda},
			{Key: "mean_claim_mu", Val: claims.Mu},
			{Key: "loading_theta", Val: claims.Theta},
			{Key: "premium_c", Val: claims.Premium},
			{Key: "rho", Val: claims.Rho},
			{Key: "hyperexp2_p1", Val: 0.5},
			{Key: "hyperexp2_beta1", Val: 2.0 / 3.0},
			{Key: "hyperexp2_beta2", Val: 2.0},
			{Key: "pareto_alpha", Val: 2.5},
			{Key: "pareto_scale", Val: 1.5},
			{Key: "weibull_shape", Val: 0.5},
			{Key: "weibull_scale", Val: 0.5},
			{Key: "reinsurance_loading_xi", Val: reinsurance.ReinsurerLoading},
		}},
		{Key: "part1_closed_form_moments_deficit", Val: jsonout.Obj{
			{Key: "closed_forms_and_lundberg", Val: closedObjs},
			{Key: "process_moments_t10", Val: momObjs},
			{Key: "exponential_deficit_ks", Val: defObjs},
		}},
		{Key: "part2_panjer_discretization", Val: jsonout.Obj{
			{Key: "light_tailed_bounds", Val: panjerLightObjs},
			{Key: "heavy_tailed_bounds_to_u200", Val: panjerHeavyObjs},
		}},
		{Key: "part3_crude_mc_and_is", Val: jsonout.Obj{
			{Key: "crude_mc_light_1m", Val: crudeObjs},
			{Key: "exponential_is_optimal", Val: isOptObjs},
			{Key: "exponential_is_misspecified", Val: misspecObjs},
		}},
		{Key: "part4_heavy_tails_ak", Val: jsonout.Obj{
			{Key: "mgf_divergence_proof", Val: mgfDivObjs},
			{Key: "asmussen_kroese_vs_panjer_crude_asymp", Val: akObjs},
		}},
		{Key: "part5_reinsurance", Val: jsonout.Obj{
			{Key: "optimal_retention_and_capital", Val: reinsObjs},
		}},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_closed_form_and_lundberg_R", Val: part1ClosedOK},
			{Key: "c1_process_moments_means_3se", Val: momentsMeanAllOK},
			{Key: "c1_process_var_St_finite_4th_moment", Val: momentsLightOK},
			{Key: "c1_exponential_deficit_ks_pass", Val: deficitOK},
			{Key: "c2_panjer_brackets_exact_all_hu", Val: panjerBracketAllOK},
			{Key: "c2_panjer_order_approx_1", Val: panjerOrderAllOK},
			{Key: "c2_panjer_heavy_tail_stable_u200", Val: panjerHeavyStableOK},
			{Key: "c3_crude_mc_1m_and_horizon_monotone", Val: crudeOK},
			{Key: "c3_exponential_is_unbiased_to_u100", Val: isOptOK},
			{Key: "c3_misspecified_is_unbiased_higher_var", Val: misspecOK},
			{Key: "c4_asmussen_kroese_within_panjer_band", Val: akAllOK},
			{Key: "c5_reinsurance_optimal_retention_and_capital", Val: reinsOK},
			{Key: "worker_bit_identical", Val: workerBitIdentical},
			{Key: "all_passed", Val: allPassed},
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		panic(err)
	}
	fmt.Printf("day003f selesai dalam %.2fs | P1=%.2fs P2=%.2fs P3=%.2fs P4=%.2fs P5=%.2fs | all_passed=%v\n",
		elapsedTotal, elapsedP1, elapsedP2, elapsedP3, elapsedP4, elapsedP5, allPassed)
}
