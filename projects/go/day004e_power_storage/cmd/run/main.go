package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/battery"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/dp"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/jsonout"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/jumps"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/parallel"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/prices"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/stats"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/svg"
)

func findRepoRoot() string {
	wd, err := os.Getwd()
	if err != nil {
		return "."
	}
	dir := wd
	for i := 0; i < 8; i++ {
		if _, err := os.Stat(filepath.Join(dir, "outputs")); err == nil {
			if _, err2 := os.Stat(filepath.Join(dir, "projects", "go", "day004e_power_storage")); err2 == nil {
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

func isFinite(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}

type pathDiag struct {
	TotalHours        int
	NegPriceHours     int
	NegPricePct       float64
	ScarcityHours     int
	OversupplyHours   int
	NormalHours       int
	ScarcityJumps     int
	OversupplyJumps   int
	NormalPoissonCnt  int
	TotalNonZeroJumps int
	HoursWithJump     int
	MultiJumpHours    int
	MinPrice          float64
	Q01Price          float64
	Q05Price          float64
	MedianPrice       float64
	MeanPrice         float64
	Q95Price          float64
	Q99Price          float64
	MaxPrice          float64
	DaysWithNegPrice  int
	DaysWithSpike200  int
	DaysWithBoth      int
}

func analyzePaths(days []prices.DayPath) pathDiag {
	nDays := len(days)
	allP := make([]float64, 0, nDays*prices.HoursPerDay)
	var d pathDiag
	d.MinPrice = math.Inf(1)
	d.MaxPrice = math.Inf(-1)
	for _, day := range days {
		hasNeg, hasSpike := false, false
		for h := 0; h < prices.HoursPerDay; h++ {
			p := day.Prices[h]
			r := day.Regimes[h]
			nj := day.NumJumps[h]
			allP = append(allP, p)
			d.TotalHours++
			if p < 0 {
				d.NegPriceHours++
				hasNeg = true
			}
			if p >= 200.0 {
				hasSpike = true
			}
			if p < d.MinPrice {
				d.MinPrice = p
			}
			if p > d.MaxPrice {
				d.MaxPrice = p
			}
			switch r {
			case regime.Normal:
				d.NormalHours++
				d.NormalPoissonCnt += nj
			case regime.Scarcity:
				d.ScarcityHours++
				d.ScarcityJumps += nj
				if nj > 0 {
					d.HoursWithJump++
					d.TotalNonZeroJumps += nj
				}
				if nj > 1 {
					d.MultiJumpHours++
				}
			case regime.Oversupply:
				d.OversupplyHours++
				d.OversupplyJumps += nj
				if nj > 0 {
					d.HoursWithJump++
					d.TotalNonZeroJumps += nj
				}
				if nj > 1 {
					d.MultiJumpHours++
				}
			}
		}
		if hasNeg {
			d.DaysWithNegPrice++
		}
		if hasSpike {
			d.DaysWithSpike200++
		}
		if hasNeg && hasSpike {
			d.DaysWithBoth++
		}
	}
	d.NegPricePct = 100.0 * float64(d.NegPriceHours) / math.Max(float64(d.TotalHours), 1.0)
	d.MeanPrice = stats.Mean(allP)
	d.Q01Price = stats.Quantile(allP, 0.01)
	d.Q05Price = stats.Quantile(allP, 0.05)
	d.MedianPrice = stats.Quantile(allP, 0.50)
	d.Q95Price = stats.Quantile(allP, 0.95)
	d.Q99Price = stats.Quantile(allP, 0.99)
	return d
}

// selectExampleDay memilih hari contoh menggunakan aturan tetap pra-run:
// hari pertama dengan harga < 0 dan harga >= 200; fallback ke rentang harga terbesar.
func selectExampleDay(days []prices.DayPath) (int, string) {
	for i, day := range days {
		hasNeg, hasHigh := false, false
		for h := 0; h < prices.HoursPerDay; h++ {
			if day.Prices[h] < 0 {
				hasNeg = true
			}
			if day.Prices[h] >= 200.0 {
				hasHigh = true
			}
		}
		if hasNeg && hasHigh {
			return i, "first_day_with_negative_price_and_spike_ge_200"
		}
	}
	bestIdx := 0
	bestSpan := -1.0
	for i, day := range days {
		lo, hi := day.Prices[0], day.Prices[0]
		for h := 1; h < prices.HoursPerDay; h++ {
			if day.Prices[h] < lo {
				lo = day.Prices[h]
			}
			if day.Prices[h] > hi {
				hi = day.Prices[h]
			}
		}
		if hi-lo > bestSpan {
			bestSpan = hi - lo
			bestIdx = i
		}
	}
	return bestIdx, "fallback_max_intraday_price_span"
}

// buildContinuousSampleWindow mensimulasikan jendela kontinu 7 hari (168 jam) setelah burn-in 60 hari
// dan memilih seed indeks pertama yang memuat baik harga negatif maupun jump scarcity untuk visualisasi.
func buildContinuousSampleWindow(spec prices.ModelSpec, pi [regime.NumRegimes]float64) ([]prices.DayPath, int) {
	const hours = 7 * prices.HoursPerDay
	for try := 0; try < 50; try++ {
		rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x5A00, uint64(try)))
		r := regime.SampleStationary(pi, rg)
		p := spec.SeasonalMean(r, 0)
		for step := 0; step < spec.BurnInDays*prices.HoursPerDay; step++ {
			r, p, _, _ = spec.Step(r, step, p, rg)
		}
		win := make([]prices.DayPath, 7)
		hasNeg, hasScarJump := false, false
		for idx := 0; idx < hours; idx++ {
			r, p, nj, jv := spec.Step(r, idx, p, rg)
			d := idx / prices.HoursPerDay
			h := idx % prices.HoursPerDay
			win[d].DayIdx = d
			win[d].Prices[h] = p
			win[d].Regimes[h] = r
			win[d].NumJumps[h] = nj
			win[d].JumpVals[h] = jv
			if p < 0 {
				hasNeg = true
			}
			if r == regime.Scarcity && nj > 0 {
				hasScarJump = true
			}
		}
		if hasNeg && hasScarJump {
			return win, try
		}
	}
	// Fallback ke try 0 jika tidak ditemukan dalam 50 percobaan
	rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x5A00, 0))
	r := regime.SampleStationary(pi, rg)
	p := spec.SeasonalMean(r, 0)
	for step := 0; step < spec.BurnInDays*prices.HoursPerDay; step++ {
		r, p, _, _ = spec.Step(r, step, p, rg)
	}
	win := make([]prices.DayPath, 7)
	for idx := 0; idx < hours; idx++ {
		r, p, nj, jv := spec.Step(r, idx, p, rg)
		d := idx / prices.HoursPerDay
		h := idx % prices.HoursPerDay
		win[d].DayIdx = d
		win[d].Prices[h] = p
		win[d].Regimes[h] = r
		win[d].NumJumps[h] = nj
		win[d].JumpVals[h] = jv
	}
	return win, 0
}

func strategyToJSON(sr dp.StrategyResult, inclDaily bool) jsonout.Obj {
	obj := jsonout.Obj{
		{Key: "strategy", Val: sr.Name},
		{Key: "mean_daily_profit_eur", Val: sr.MeanProfit},
		{Key: "median_daily_profit_eur", Val: sr.MedianProfit},
		{Key: "sd_daily_profit_eur", Val: sr.SDProfit},
		{Key: "mcse_daily_profit_eur", Val: sr.MCSEProfit},
		{Key: "var_5pct_daily_profit_eur", Val: sr.VaR5Profit},
		{Key: "es_5pct_daily_profit_eur", Val: sr.ES5Profit},
		{Key: "mean_gross_value_eur", Val: sr.MeanGrossValue},
		{Key: "mean_cashflow_only_eur", Val: sr.MeanCashflow},
		{Key: "mean_daily_grid_throughput_mwh", Val: sr.MeanThroughput},
		{Key: "mean_daily_efc", Val: sr.MeanEFC},
		{Key: "mean_regret_vs_pf_eur", Val: sr.MeanRegretVsPF},
		{Key: "mcse_regret_vs_pf_eur", Val: sr.MCSERegretVsPF},
		{Key: "max_regret_vs_pf_eur", Val: sr.MaxRegretVsPF},
		{Key: "action_counts", Val: jsonout.Obj{
			{Key: "charge_total", Val: sr.TotalChargeActions},
			{Key: "charge_2mw", Val: sr.Charge2MWCount},
			{Key: "charge_1mw", Val: sr.Charge1MWCount},
			{Key: "idle", Val: sr.TotalIdleActions},
			{Key: "discharge_total", Val: sr.TotalDischargeActions},
			{Key: "discharge_09mw", Val: sr.Discharge09MWCount},
			{Key: "discharge_18mw", Val: sr.Discharge18MWCount},
		}},
		{Key: "total_profit_neg_price_hours_eur", Val: sr.TotalProfitNegPriceHours},
		{Key: "mean_profit_neg_price_per_day_eur", Val: sr.MeanProfitNegPricePerDay},
		{Key: "total_profit_scarcity_or_jump_hours_eur", Val: sr.TotalProfitScarcityOrJumpHrs},
		{Key: "mean_profit_scarcity_or_jump_per_day_eur", Val: sr.MeanProfitScarcityOrJumpPerD},
		{Key: "soc_violation_count", Val: sr.SoCViolationCount},
		{Key: "power_violation_count", Val: sr.PowerViolationCount},
		{Key: "pf_domination_violations", Val: sr.PerfectForesightDomViolations},
		{Key: "max_pf_domination_shortfall_eur", Val: sr.MaxDominationShortfall},
	}
	if inclDaily {
		obj = append(obj, jsonout.Field{Key: "daily_profits_eur", Val: sr.DailyProfits})
		obj = append(obj, jsonout.Field{Key: "daily_regret_vs_pf_eur", Val: sr.DailyRegretVsPF})
		obj = append(obj, jsonout.Field{Key: "daily_efc", Val: sr.DailyEFC})
	}
	return obj
}

func pairedDiffToJSON(label string, pd stats.PairedDiffSummary) jsonout.Obj {
	return jsonout.Obj{
		{Key: "comparison", Val: label},
		{Key: "n_days", Val: pd.N},
		{Key: "mean_diff_eur", Val: pd.MeanDiff},
		{Key: "sd_diff_eur", Val: pd.SDDiff},
		{Key: "mcse_diff_eur", Val: pd.MCSE},
		{Key: "ci_95_level", Val: 0.95},
		{Key: "ci_95_z", Val: stats.Z95},
		{Key: "ci_95_lower_eur", Val: pd.CI95Low},
		{Key: "ci_95_upper_eur", Val: pd.CI95High},
		{Key: "positive_day_fraction", Val: pd.PositiveFrac},
		{Key: "z_stat", Val: pd.ZStat},
	}
}

func kernelToJSON(k dp.TransitionKernel, pol dp.DPPolicy) jsonout.Obj {
	return jsonout.Obj{
		{Key: "grid_name", Val: k.Grid.Name},
		{Key: "num_price_bins", Val: k.Grid.NumBins},
		{Key: "edges", Val: k.Grid.Edges},
		{Key: "representative_prices", Val: k.Grid.RepPrices},
		{Key: "draws_per_state", Val: k.DrawsPerState},
		{Key: "num_states_h_r_b", Val: k.NumStates},
		{Key: "total_conditional_draws", Val: k.TotalDraws},
		{Key: "lower_tail_threshold_eur", Val: k.Grid.LowerTailAt},
		{Key: "upper_tail_threshold_eur", Val: k.Grid.UpperTailAt},
		{Key: "lower_tail_prob_mass", Val: k.LowerTailMass},
		{Key: "upper_tail_prob_mass", Val: k.UpperTailMass},
		{Key: "lower_tail_count", Val: k.LowerTailCount},
		{Key: "upper_tail_count", Val: k.UpperTailCount},
		{Key: "lower_tail_empirical_mean_eur", Val: k.LowerTailMeanP},
		{Key: "upper_tail_empirical_mean_eur", Val: k.UpperTailMeanP},
		{Key: "kernel_build_seconds", Val: k.BuildSeconds},
		{Key: "dp_solve_seconds", Val: pol.SolveSeconds},
		{Key: "max_bellman_residual", Val: pol.MaxBellmanResidual},
		{Key: "terminal_max_error", Val: pol.TerminalMaxError},
	}
}

func buildFigures(
	figDir string,
	sampleWin []prices.DayPath,
	exDay prices.DayPath,
	pol17, pol33 dp.DPPolicy,
	res17, res33, resThresh, resPF dp.StrategyResult,
	bat battery.Spec,
) error {
	if err := os.MkdirAll(figDir, 0o755); err != nil {
		return err
	}

	// Grafik 1: Lintasan harga contoh 7 hari kontinu dengan pita rezim, jump, dan harga negatif
	var xs, ys []float64
	var jumpMask []bool
	var bands []svg.VBand
	for d, day := range sampleWin {
		for h := 0; h < prices.HoursPerDay; h++ {
			t := float64(d*prices.HoursPerDay + h)
			xs = append(xs, t)
			ys = append(ys, day.Prices[h])
			isNonZeroJump := day.Regimes[h] != regime.Normal && day.NumJumps[h] > 0
			jumpMask = append(jumpMask, isNonZeroJump)
			if day.Regimes[h] == regime.Scarcity {
				bands = append(bands, svg.VBand{XStart: t - 0.5, XEnd: t + 0.5, Color: "#fddbc7"})
			} else if day.Regimes[h] == regime.Oversupply {
				bands = append(bands, svg.VBand{XStart: t - 0.5, XEnd: t + 0.5, Color: "#d1e5f0"})
			}
		}
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004e_sample-price-regimes-jumps.svg"), svg.ChartOptions{
		Title:    "Lintasan Harga Listrik Sintetis 7 Hari (Rezim, Harga Negatif & Jump)",
		XLabel:   "Jam Simulasi Kontinu Pasca Burn-In (0..167)",
		YLabel:   "Harga Spot (EUR/MWh)",
		ZeroLine: true,
	}, bands, []svg.LineSeries{
		{
			Label:      "Harga Spot (Merah=Jump Non-Nol, Pita=Rezim)",
			X:          xs,
			Y:          ys,
			Color:      "#1f77b4",
			StrokeW:    1.9,
			PointMask:  jumpMask,
			PointColor: "#d62728",
		},
	}); err != nil {
		return err
	}

	// Grafik 2: Dispatch baterai & SoC pada hari contoh yang dipilih lewat aturan pra-run
	trDP33 := dp.TraceExampleDay(exDay, bat, &pol33, false, false)
	trTh := dp.TraceExampleDay(exDay, bat, nil, true, false)
	trPF := dp.TraceExampleDay(exDay, bat, nil, false, true)
	hrs24 := make([]float64, prices.HoursPerDay)
	hrs25 := make([]float64, prices.HoursPerDay+1)
	for h := 0; h < prices.HoursPerDay; h++ {
		hrs24[h] = float64(h)
		hrs25[h] = float64(h)
	}
	hrs25[24] = 24.0
	if err := svg.WriteDualDayDispatchSVG(
		filepath.Join(figDir, "day004e_example-day-dispatch-soc.svg"),
		fmt.Sprintf("Dispatch Baterai & SoC pada Hari Contoh #%d (Aturan Pra-Run)", exDay.DayIdx),
		hrs24,
		exDay.Prices[:],
		[]svg.LineSeries{
			{Label: fmt.Sprintf("DP-33 (Profit %.1f EUR)", trDP33.NetProfit), X: hrs24, Y: trDP33.PowerMW[:], Color: svg.Palette[0], StrokeW: 2.2},
			{Label: fmt.Sprintf("Threshold (Profit %.1f EUR)", trTh.NetProfit), X: hrs24, Y: trTh.PowerMW[:], Color: svg.Palette[4], StrokeW: 2.0, Dashed: true},
			{Label: fmt.Sprintf("Perfect Foresight (%.1f EUR)", trPF.NetProfit), X: hrs24, Y: trPF.PowerMW[:], Color: svg.Palette[2], StrokeW: 1.8, Dashed: true},
		},
		[]svg.LineSeries{
			{Label: "SoC DP-33 (MWh)", X: hrs25, Y: trDP33.SoCMWh[:], Color: svg.Palette[0], StrokeW: 2.2},
			{Label: "SoC Threshold (MWh)", X: hrs25, Y: trTh.SoCMWh[:], Color: svg.Palette[4], StrokeW: 2.0, Dashed: true},
			{Label: "SoC Perfect Foresight (MWh)", X: hrs25, Y: trPF.SoCMWh[:], Color: svg.Palette[2], StrokeW: 1.8, Dashed: true},
		},
	); err != nil {
		return err
	}

	// Grafik 3: Distribusi profit harian per strategi
	if err := svg.WriteOverlayDistributionSVG(
		filepath.Join(figDir, "day004e_daily-profit-distribution.svg"),
		"Distribusi Profit Arbitrase Harian Antar-Strategi (N=2.000 Hari CRN)",
		"Profit Bersih Harian (EUR, net terhadap nilai terminal awal 81 EUR)",
		[]svg.DistSeries{
			{Label: fmt.Sprintf("DP-33 (Mean=%.1f, VaR5=%.1f)", res33.MeanProfit, res33.VaR5Profit), Values: res33.DailyProfits, Color: svg.Palette[0]},
			{Label: fmt.Sprintf("DP-17 (Mean=%.1f, VaR5=%.1f)", res17.MeanProfit, res17.VaR5Profit), Values: res17.DailyProfits, Color: svg.Palette[3], Dashed: true},
			{Label: fmt.Sprintf("Threshold (Mean=%.1f, VaR5=%.1f)", resThresh.MeanProfit, resThresh.VaR5Profit), Values: resThresh.DailyProfits, Color: svg.Palette[4]},
			{Label: fmt.Sprintf("Perfect Foresight (Mean=%.1f)", resPF.MeanProfit), Values: resPF.DailyProfits, Color: svg.Palette[2], Dashed: true},
		},
		32,
	); err != nil {
		return err
	}

	// Grafik 4: Nilai kebijakan DP terhadap SoC dan rezim pada jam h=12 (bin harga sekitar 40-50 EUR/MWh)
	bMid := pol33.Grid.BinIndex(45.0)
	socVals := make([]float64, bat.NumSoCStates)
	valNorm := make([]float64, bat.NumSoCStates)
	valScar := make([]float64, bat.NumSoCStates)
	valOver := make([]float64, bat.NumSoCStates)
	for s := 0; s < bat.NumSoCStates; s++ {
		socVals[s] = battery.SoCMWh(s)
		valNorm[s] = pol33.ValueAt(12, regime.Normal, bMid, s)
		valScar[s] = pol33.ValueAt(12, regime.Scarcity, bMid, s)
		valOver[s] = pol33.ValueAt(12, regime.Oversupply, bMid, s)
	}
	if err := svg.WriteLineChart(filepath.Join(figDir, "day004e_dp-value-vs-soc-regime.svg"), svg.ChartOptions{
		Title:         fmt.Sprintf("Fungsi Nilai Bellman DP-33 V(h=12, r, bin~%.0f EUR/MWh, SoC)", pol33.Grid.RepPrices[bMid]),
		XLabel:        "State of Charge SoC (MWh)",
		YLabel:        "Nilai Harapan Total V (EUR, termasuk nilai terminal)",
		LegendTopLeft: true,
	}, nil, []svg.LineSeries{
		{Label: "Rezim Normal", X: socVals, Y: valNorm, Color: svg.Palette[0], StrokeW: 2.3},
		{Label: "Rezim Scarcity", X: socVals, Y: valScar, Color: svg.Palette[1], StrokeW: 2.3},
		{Label: "Rezim Oversupply", X: socVals, Y: valOver, Color: svg.Palette[2], StrokeW: 2.3, Dashed: true},
	}); err != nil {
		return err
	}

	// Grafik 5: Sensitivitas hasil stochastic DP terhadap resolusi grid harga (17 vs 33 bin vs Threshold)
	if err := svg.WriteCategoryPointChart(
		filepath.Join(figDir, "day004e_grid-resolution-sensitivity.svg"),
		svg.ChartOptions{
			Title:        "Sensitivitas Kinerja DP terhadap Resolusi Grid Harga (17 vs 33 Bin vs Threshold)",
			XLabel:       "Metrik Profit Harian (EUR, +/- 1.96 MCSE untuk Rata-rata)",
			YLabel:       "EUR per Hari",
			ZeroLine:     true,
			LegendBottom: true,
		},
		[]string{"Rata-rata Profit", "Median Profit", "VaR 5% Profit"},
		[]svg.CategoryPointSeries{
			{
				Label:  "Stochastic DP (17 bin)",
				Values: []float64{res17.MeanProfit, res17.MedianProfit, res17.VaR5Profit},
				Lower:  []float64{res17.MeanProfit - stats.Z95*res17.MCSEProfit, res17.MedianProfit, res17.VaR5Profit},
				Upper:  []float64{res17.MeanProfit + stats.Z95*res17.MCSEProfit, res17.MedianProfit, res17.VaR5Profit},
				Color:  svg.Palette[3],
			},
			{
				Label:  "Stochastic DP (33 bin)",
				Values: []float64{res33.MeanProfit, res33.MedianProfit, res33.VaR5Profit},
				Lower:  []float64{res33.MeanProfit - stats.Z95*res33.MCSEProfit, res33.MedianProfit, res33.VaR5Profit},
				Upper:  []float64{res33.MeanProfit + stats.Z95*res33.MCSEProfit, res33.MedianProfit, res33.VaR5Profit},
				Color:  svg.Palette[0],
			},
			{
				Label:  "Threshold Policy",
				Values: []float64{resThresh.MeanProfit, resThresh.MedianProfit, resThresh.VaR5Profit},
				Lower:  []float64{resThresh.MeanProfit - stats.Z95*resThresh.MCSEProfit, resThresh.MedianProfit, resThresh.VaR5Profit},
				Upper:  []float64{resThresh.MeanProfit + stats.Z95*resThresh.MCSEProfit, resThresh.MedianProfit, resThresh.VaR5Profit},
				Color:  svg.Palette[4],
			},
		},
	); err != nil {
		return err
	}
	return nil
}

func main() {
	root := findRepoRoot()
	quick := flag.Bool("quick", false, "mode cepat pengembangan ke temp dir tanpa menimpa keluaran penuh")
	workersFlag := flag.Int("workers", 0, "jumlah worker goroutine (0 = min(8, NumCPU))")
	daysFlag := flag.Int("days", 0, "override jumlah hari evaluasi N (bawaan: 2000 penuh, 50 -quick)")
	drawsFlag := flag.Int("draws", 0, "override jumlah draw kernel per state (bawaan: 10000 penuh, 500 -quick)")
	outFlag := flag.String("out", "", "override berkas keluaran JSON")
	figDirFlag := flag.String("fig-dir", "", "override direktori keluaran SVG")
	flag.Parse()

	nWorkers := *workersFlag
	if nWorkers <= 0 {
		nWorkers = parallel.DefaultWorkers()
	}
	nDays := 2000
	nDraws := 10000
	if *quick {
		nDays = 50
		nDraws = 500
	}
	if *daysFlag > 0 {
		nDays = *daysFlag
	}
	if *drawsFlag > 0 {
		nDraws = *drawsFlag
	}

	outPath := filepath.Join(root, "outputs", "day004e", "results.json")
	figDir := filepath.Join(root, "outputs", "figures")
	if *quick && *outFlag == "" {
		outPath = filepath.Join(os.TempDir(), "day004e_quick_results.json")
	}
	if *quick && *figDirFlag == "" {
		figDir = filepath.Join(os.TempDir(), "day004e_quick_figures")
	}
	if *outFlag != "" {
		outPath = *outFlag
	}
	if *figDirFlag != "" {
		figDir = *figDirFlag
	}

	tTotal := time.Now()
	spec := prices.DefaultSpec()
	bat := battery.DefaultSpec()

	// 1. Verifikasi model: matriks transisi, distribusi stasioner, momen jump, transisi satu langkah, determinisme worker
	tVerify0 := time.Now()
	pValidErr := regime.ValidateTransitionMatrix(spec.TransitionMatrix)
	pi, piErr := regime.StationaryDistribution(spec.TransitionMatrix)
	piRes := regime.StationaryResidual(spec.TransitionMatrix, pi)
	c1StationaryOK := pValidErr == nil && piErr == nil && piRes < 1e-13

	const mcDraws = 200000
	jumpVerifyObjs := make([]jsonout.Obj, regime.NumRegimes)
	c2JumpsOK := true
	for r := 0; r < regime.NumRegimes; r++ {
		rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x2100+uint64(r), 0))
		var sum, sum2 float64
		for i := 0; i < mcDraws; i++ {
			_, jv := jumps.SampleJump(spec.Jumps, r, rg)
			sum += jv
			sum2 += jv * jv
		}
		empM := sum / mcDraws
		empV := (sum2 - mcDraws*empM*empM) / (mcDraws - 1)
		theoM, theoV := jumps.TheoreticalMoments(spec.Jumps, r)
		var mcseM, mcseV float64
		var ok bool
		if r == regime.Normal {
			ok = empM == 0.0 && empV == 0.0
		} else {
			mcseM = math.Sqrt(theoV / mcDraws)
			rg2 := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x2100+uint64(r), 0))
			var m4 float64
			for i := 0; i < mcDraws; i++ {
				_, jv := jumps.SampleJump(spec.Jumps, r, rg2)
				d := jv - empM
				d2 := d * d
				m4 += d2 * d2
			}
			m4 /= mcDraws
			mcseV = math.Sqrt(math.Max(0, (m4-empV*empV)/mcDraws))
			ok = math.Abs(empM-theoM) <= 4.0*mcseM && math.Abs(empV-theoV) <= 4.0*mcseV
		}
		if !ok {
			c2JumpsOK = false
		}
		jumpVerifyObjs[r] = jsonout.Obj{
			{Key: "regime", Val: regime.Names[r]},
			{Key: "lambda", Val: spec.Jumps[r].Lambda},
			{Key: "mean_jump", Val: spec.Jumps[r].MeanJump},
			{Key: "theoretical_mean", Val: theoM},
			{Key: "empirical_mean", Val: empM},
			{Key: "mcse_mean", Val: mcseM},
			{Key: "theoretical_variance", Val: theoV},
			{Key: "empirical_variance", Val: empV},
			{Key: "mcse_variance", Val: mcseV},
			{Key: "within_4_mcse", Val: ok},
		}
	}

	stepSpecs := []struct {
		r int
		h int
		p float64
	}{
		{regime.Normal, 18, 55.0},
		{regime.Scarcity, 23, 210.0},
		{regime.Oversupply, 4, -85.0},
	}
	stepVerifyObjs := make([]jsonout.Obj, len(stepSpecs))
	c2StepOK := true
	for idx, tc := range stepSpecs {
		rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x3100+uint64(idx), 0))
		var sum, sum2 float64
		for i := 0; i < mcDraws; i++ {
			_, pNext, _, _ := spec.Step(tc.r, tc.h, tc.p, rg)
			sum += pNext
			sum2 += pNext * pNext
		}
		empM := sum / mcDraws
		empV := (sum2 - mcDraws*empM*empM) / (mcDraws - 1)
		theoM, theoV := spec.OneStepConditionalMoments(tc.r, tc.h, tc.p)
		mcseM := math.Sqrt(theoV / mcDraws)
		rg2 := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x3100+uint64(idx), 0))
		var m4 float64
		for i := 0; i < mcDraws; i++ {
			_, pNext, _, _ := spec.Step(tc.r, tc.h, tc.p, rg2)
			d := pNext - empM
			d2 := d * d
			m4 += d2 * d2
		}
		m4 /= mcDraws
		mcseV := math.Sqrt(math.Max(0, (m4-empV*empV)/mcDraws))
		ok := math.Abs(empM-theoM) <= 4.0*mcseM && math.Abs(empV-theoV) <= 4.0*mcseV
		if !ok {
			c2StepOK = false
		}
		stepVerifyObjs[idx] = jsonout.Obj{
			{Key: "r_now", Val: regime.Names[tc.r]},
			{Key: "hour_now", Val: tc.h},
			{Key: "p_now", Val: tc.p},
			{Key: "theoretical_mean", Val: theoM},
			{Key: "empirical_mean", Val: empM},
			{Key: "mcse_mean", Val: mcseM},
			{Key: "theoretical_variance", Val: theoV},
			{Key: "empirical_variance", Val: empV},
			{Key: "mcse_variance", Val: mcseV},
			{Key: "within_4_mcse", Val: ok},
		}
	}
	c2MomentsOK := c2JumpsOK && c2StepOK

	// Determinisme worker 1 vs nWorkers pada fixture cepat
	g17 := dp.NewCoarseGrid17()
	g33 := dp.NewFineGrid33()
	fixK1 := dp.BuildKernel(spec, g17, 200, 1, rng.KernelSeed, 0xD4E1)
	fixKN := dp.BuildKernel(spec, g17, 200, nWorkers, rng.KernelSeed, 0xD4E1)
	fixDays1 := prices.SimulateEvalDays(20, 1, spec, pi, rng.EvalSeed, 0xD4E2)
	fixDaysN := prices.SimulateEvalDays(20, nWorkers, spec, pi, rng.EvalSeed, 0xD4E2)
	workerDetOK := true
	for i := range fixK1.Probs {
		if math.Float64bits(fixK1.Probs[i]) != math.Float64bits(fixKN.Probs[i]) {
			workerDetOK = false
			break
		}
	}
	for d := range fixDays1 {
		for h := 0; h < prices.HoursPerDay; h++ {
			if math.Float64bits(fixDays1[d].Prices[h]) != math.Float64bits(fixDaysN[d].Prices[h]) {
				workerDetOK = false
				break
			}
		}
	}
	verifySeconds := time.Since(tVerify0).Seconds()

	// 2. Uji Bellman vs enumerasi brute-force 4 jam deterministik
	bfPrices := []float64{-40.0, 35.0, 240.0, 10.0}
	dpBFNet, dpBFGross, _, _ := dp.SolveDeterministicDP(bfPrices, bat)
	// Hitung brute force langsung di runner agar teraudit di results.json
	bestBFGross := math.Inf(-1)
	var recBF func(h, s int, cf float64)
	recBF = func(h, s int, cf float64) {
		if h == len(bfPrices) {
			if tot := cf + bat.TerminalValue(s); tot > bestBFGross {
				bestBFGross = tot
			}
			return
		}
		for _, a := range battery.Actions {
			if sn, ok := bat.NextSoCIdx(s, a); ok {
				recBF(h+1, sn, cf+bat.HourlyCashflow(bfPrices[h], a))
			}
		}
	}
	recBF(0, bat.InitialSoCIdx, 0.0)
	bfNet := bestBFGross - bat.TerminalValue(bat.InitialSoCIdx)
	bfDiff := math.Abs(dpBFGross - bestBFGross)

	// 3. Bangun kernel transisi (17 bin & 33 bin) dan selesaikan Bellman DP
	fmt.Printf("[day004e] Membangun kernel DP (17 & 33 bin, draws/state=%d, workers=%d)...\n", nDraws, nWorkers)
	k17 := dp.BuildKernel(spec, g17, nDraws, nWorkers, rng.KernelSeed, 0x4E17)
	pol17 := dp.SolveStochasticDP(k17, bat)

	k33 := dp.BuildKernel(spec, g33, nDraws, nWorkers, rng.KernelSeed, 0x4E33)
	pol33 := dp.SolveStochasticDP(k33, bat)

	// Audit sensitivitas tail diperluas [-1000, 5000] dengan CRN kondisional pada bin interior yang sama
	gExp := dp.NewExpandedTailGrid()
	crnMapper := func(h, r, b int) (uint64, int) {
		if b >= 3 && b <= 33 {
			return 0x4E33, b - 2
		}
		return 0x4E38, b
	}
	kExp := dp.BuildKernelWithSeedMapper(spec, gExp, nDraws, nWorkers, rng.KernelSeed, 0x4E38, crnMapper)
	polExp := dp.SolveStochasticDP(kExp, bat)

	c4BellmanOK := bfDiff < 1e-11 && pol17.MaxBellmanResidual < 1e-11 && pol33.MaxBellmanResidual < 1e-11 &&
		polExp.MaxBellmanResidual < 1e-11 && pol17.TerminalMaxError < 1e-12 && pol33.TerminalMaxError < 1e-12 &&
		polExp.TerminalMaxError < 1e-12

	// 4. Simulasi N hari evaluasi independen (CRN untuk semua strategi)
	tEval0 := time.Now()
	fmt.Printf("[day004e] Mensimulasikan N=%d hari evaluasi (burn-in %d hari) dan mengevaluasi 5 strategi...\n", nDays, spec.BurnInDays)
	evalDays := prices.SimulateEvalDays(nDays, nWorkers, spec, pi, rng.EvalSeed, 0x4E99)
	pDiag := analyzePaths(evalDays)

	res17 := dp.EvalStochasticDP(dp.StrategyDP17, evalDays, pol17, bat)
	res33 := dp.EvalStochasticDP(dp.StrategyDP33, evalDays, pol33, bat)
	resExp := dp.EvalStochasticDP("stochastic_dp_expanded_38bin", evalDays, polExp, bat)
	resThresh := dp.EvalThresholdPolicy(evalDays, bat)
	resPF := dp.EvalPerfectForesightDP(evalDays, bat)
	resNoBat := dp.EvalNoBattery(evalDays)

	const domTol = 1e-9
	dp.ComputeRegretVsPF(&res17, resPF, domTol)
	dp.ComputeRegretVsPF(&res33, resPF, domTol)
	dp.ComputeRegretVsPF(&resExp, resPF, domTol)
	dp.ComputeRegretVsPF(&resThresh, resPF, domTol)
	dp.ComputeRegretVsPF(&resNoBat, resPF, domTol)
	evalSeconds := time.Since(tEval0).Seconds()

	diff17VsThresh := stats.PairedDiffCI(res17.DailyProfits, resThresh.DailyProfits)
	diff33VsThresh := stats.PairedDiffCI(res33.DailyProfits, resThresh.DailyProfits)
	diff33Vs17 := stats.PairedDiffCI(res33.DailyProfits, res17.DailyProfits)
	diffExpVs33 := stats.PairedDiffCI(resExp.DailyProfits, res33.DailyProfits)
	diffExpVs17 := stats.PairedDiffCI(resExp.DailyProfits, res17.DailyProfits)
	evalStepAgree, interiorStateAgree, evalDiffSteps, evalTotSteps, intAgreeCnt, intTotStates := dp.PolicyActionAgreement(evalDays, pol33, polExp, bat)

	// 5. Pilih hari contoh pra-run dan buat 5 grafik SVG
	tSvg0 := time.Now()
	exIdx, exRule := selectExampleDay(evalDays)
	sampleWin, sampleWinTry := buildContinuousSampleWindow(spec, pi)
	if err := buildFigures(figDir, sampleWin, evalDays[exIdx], pol17, pol33, res17, res33, resThresh, resPF, bat); err != nil {
		fmt.Fprintf(os.Stderr, "buildFigures error: %v\n", err)
		os.Exit(1)
	}
	svgSeconds := time.Since(tSvg0).Seconds()
	totalSeconds := time.Since(tTotal).Seconds()

	// 6. Evaluasi 8 kriteria teknis pra-run
	c3BoundsOK := res17.SoCViolationCount == 0 && res17.PowerViolationCount == 0 &&
		res33.SoCViolationCount == 0 && res33.PowerViolationCount == 0 &&
		resExp.SoCViolationCount == 0 && resExp.PowerViolationCount == 0 &&
		resThresh.SoCViolationCount == 0 && resThresh.PowerViolationCount == 0 &&
		resPF.SoCViolationCount == 0 && resPF.PowerViolationCount == 0
	c5PFDomOK := res17.PerfectForesightDomViolations == 0 &&
		res33.PerfectForesightDomViolations == 0 &&
		resExp.PerfectForesightDomViolations == 0 &&
		resThresh.PerfectForesightDomViolations == 0 &&
		resNoBat.PerfectForesightDomViolations == 0
	c6DetOK := workerDetOK
	c7NegAndSpikeOK := pDiag.NegPriceHours > 0 && pDiag.TotalNonZeroJumps > 0 && pDiag.ScarcityHours > 0 && pDiag.OversupplyHours > 0
	c8OutputsFiniteOK := isFinite(res17.MeanProfit) && isFinite(res33.MeanProfit) &&
		isFinite(resThresh.MeanProfit) && isFinite(resPF.MeanProfit) &&
		isFinite(diff33VsThresh.MeanDiff) && isFinite(diff33Vs17.MeanDiff) &&
		(totalSeconds < 600.0 || *quick)
	allTechPassed := c1StationaryOK && c2MomentsOK && c3BoundsOK && c4BellmanOK &&
		c5PFDomOK && c6DetOK && c7NegAndSpikeOK && c8OutputsFiniteOK

	statusStr := func(ok bool) string {
		if ok {
			return "tercapai"
		}
		return "belum tercapai"
	}

	// Susun JSON lengkap berpresisi 17 digit
	regParams := make([]jsonout.Obj, regime.NumRegimes)
	for r := 0; r < regime.NumRegimes; r++ {
		regParams[r] = jsonout.Obj{
			{Key: "name", Val: regime.Names[r]},
			{Key: "mu_eur_mwh", Val: spec.Mu[r]},
			{Key: "phi", Val: spec.Phi[r]},
			{Key: "sigma_eur_mwh_sqrt_h", Val: spec.Sigma[r]},
			{Key: "lambda_per_h", Val: spec.Jumps[r].Lambda},
			{Key: "mean_jump_eur_mwh", Val: spec.Jumps[r].MeanJump},
			{Key: "jump_sign", Val: spec.Jumps[r].Sign},
			{Key: "stationary_prob", Val: pi[r]},
		}
	}

	rootObj := jsonout.Obj{
		{Key: "metadata", Val: jsonout.Obj{
			{Key: "session", Val: "day004e"},
			{Key: "date_wita", Val: "2026-11-22"},
			{Key: "topic", Val: "Harga listrik negatif dan arbitrase baterai (DP stokastik)"},
			{Key: "synthetic_model_disclaimer", Val: "Model harga listrik sintetis berparameter tetap untuk studi eksperimen kontrol stokastik; bukan kalibrasi pasar nyata."},
			{Key: "quick_mode", Val: *quick},
			{Key: "workers_used", Val: nWorkers},
			{Key: "seeds", Val: jsonout.Obj{
				{Key: "main_seed", Val: rng.MainSeed},
				{Key: "kernel_seed", Val: rng.KernelSeed},
				{Key: "eval_seed", Val: rng.EvalSeed},
				{Key: "verify_seed", Val: rng.VerifySeed},
			}},
			{Key: "eval_days_n", Val: nDays},
			{Key: "burn_in_days", Val: spec.BurnInDays},
			{Key: "draws_per_state", Val: nDraws},
			{Key: "elapsed_seconds", Val: jsonout.Obj{
				{Key: "model_verification", Val: verifySeconds},
				{Key: "kernel_17_build", Val: k17.BuildSeconds},
				{Key: "dp_17_solve", Val: pol17.SolveSeconds},
				{Key: "kernel_33_build", Val: k33.BuildSeconds},
				{Key: "dp_33_solve", Val: pol33.SolveSeconds},
				{Key: "simulation_and_eval", Val: evalSeconds},
				{Key: "figures_svg", Val: svgSeconds},
				{Key: "total", Val: totalSeconds},
			}},
		}},
		{Key: "model_parameters", Val: jsonout.Obj{
			{Key: "price_unit", Val: "EUR/MWh"},
			{Key: "time_step", Val: "1 jam (h modulo 24)"},
			{Key: "seasonal_amplitude_eur_mwh", Val: spec.SeasonalAmp},
			{Key: "seasonal_peak_hour", Val: spec.SeasonalPeakHour},
			{Key: "transition_matrix", Val: [][]float64{
				spec.TransitionMatrix[0][:],
				spec.TransitionMatrix[1][:],
				spec.TransitionMatrix[2][:],
			}},
			{Key: "stationary_distribution", Val: pi[:]},
			{Key: "regimes", Val: regParams},
			{Key: "battery", Val: jsonout.Obj{
				{Key: "capacity_mwh", Val: bat.CapacityMWh},
				{Key: "initial_soc_mwh", Val: bat.InitialSoCMWh},
				{Key: "soc_step_mwh", Val: bat.SoCStepMWh},
				{Key: "num_soc_states", Val: bat.NumSoCStates},
				{Key: "eta_charge", Val: bat.EtaC},
				{Key: "eta_discharge", Val: bat.EtaD},
				{Key: "max_charge_mw", Val: bat.MaxChargeMW},
				{Key: "max_discharge_mw", Val: bat.MaxDischargeMW},
				{Key: "degradation_cost_eur_per_mwh_grid", Val: bat.DegradationCost},
				{Key: "terminal_value_formula", Val: "45 * eta_d * soc_end"},
				{Key: "initial_terminal_value_eur", Val: bat.TerminalValue(bat.InitialSoCIdx)},
				{Key: "profit_definition", Val: "net_arbitrage_profit = sum_h cashflow_h + terminalValue(soc_24) - terminalValue(soc_0), sehingga no_battery dan selalu idle bernilai 0 EUR"},
				{Key: "efc_definition", Val: "sum_h |Delta SoC_h (MWh)| / (2 * E), yakni throughput energi tersimpan dibagi 2*4=8 MWh"},
				{Key: "threshold_policy_rule", Val: "charge penuh saat P <= 0, discharge penuh saat P >= 200, idle selain itu (dengan pembatas SoC)"},
			}},
		}},
		{Key: "model_verification", Val: jsonout.Obj{
			{Key: "stationary_residual_max_abs", Val: piRes},
			{Key: "jump_conditional_moments_mc_200k", Val: jumpVerifyObjs},
			{Key: "one_step_price_moments_mc_200k", Val: stepVerifyObjs},
			{Key: "bellman_4h_brute_force", Val: jsonout.Obj{
				{Key: "test_prices", Val: bfPrices},
				{Key: "brute_force_gross_eur", Val: bestBFGross},
				{Key: "dp_gross_eur", Val: dpBFGross},
				{Key: "brute_force_net_eur", Val: bfNet},
				{Key: "dp_net_eur", Val: dpBFNet},
				{Key: "abs_diff_eur", Val: bfDiff},
			}},
			{Key: "worker_determinism_fixture", Val: jsonout.Obj{
				{Key: "workers_tested", Val: []int{1, nWorkers}},
				{Key: "bit_identical", Val: workerDetOK},
			}},
		}},
		{Key: "evaluation_path_diagnostics", Val: jsonout.Obj{
			{Key: "total_eval_hours", Val: pDiag.TotalHours},
			{Key: "negative_price_hours", Val: pDiag.NegPriceHours},
			{Key: "negative_price_pct", Val: pDiag.NegPricePct},
			{Key: "normal_hours", Val: pDiag.NormalHours},
			{Key: "scarcity_hours", Val: pDiag.ScarcityHours},
			{Key: "oversupply_hours", Val: pDiag.OversupplyHours},
			{Key: "scarcity_jumps_total", Val: pDiag.ScarcityJumps},
			{Key: "oversupply_jumps_total", Val: pDiag.OversupplyJumps},
			{Key: "normal_poisson_draws_ignored_in_J", Val: pDiag.NormalPoissonCnt},
			{Key: "total_nonzero_jumps", Val: pDiag.TotalNonZeroJumps},
			{Key: "hours_with_nonzero_jump", Val: pDiag.HoursWithJump},
			{Key: "hours_with_multi_jump_gt_1", Val: pDiag.MultiJumpHours},
			{Key: "days_with_negative_price", Val: pDiag.DaysWithNegPrice},
			{Key: "days_with_spike_ge_200", Val: pDiag.DaysWithSpike200},
			{Key: "days_with_both_neg_and_spike", Val: pDiag.DaysWithBoth},
			{Key: "price_quantiles_eur_mwh", Val: jsonout.Obj{
				{Key: "min", Val: pDiag.MinPrice},
				{Key: "q01", Val: pDiag.Q01Price},
				{Key: "q05", Val: pDiag.Q05Price},
				{Key: "median", Val: pDiag.MedianPrice},
				{Key: "mean", Val: pDiag.MeanPrice},
				{Key: "q95", Val: pDiag.Q95Price},
				{Key: "q99", Val: pDiag.Q99Price},
				{Key: "max", Val: pDiag.MaxPrice},
			}},
			{Key: "example_day_selected_idx", Val: exIdx},
			{Key: "example_day_selection_rule", Val: exRule},
			{Key: "continuous_7day_sample_try", Val: sampleWinTry},
		}},
		{Key: "dp_kernels_and_policies", Val: jsonout.Obj{
			{Key: "coarse_17bin", Val: kernelToJSON(k17, pol17)},
			{Key: "fine_33bin", Val: kernelToJSON(k33, pol33)},
		}},
		{Key: "expanded_tail_sensitivity_audit", Val: jsonout.Obj{
			{Key: "description", Val: "Audit perluasan domain harga menjadi [-1000, 5000] EUR/MWh (38 bin) dengan CRN kondisional pada bin interior fine_33 yang sama"},
			{Key: "kernel_and_policy", Val: kernelToJSON(kExp, polExp)},
			{Key: "strategy_summary", Val: strategyToJSON(resExp, false)},
			{Key: "paired_vs_dp33", Val: pairedDiffToJSON("expanded38_minus_dp33", diffExpVs33)},
			{Key: "paired_vs_dp17", Val: pairedDiffToJSON("expanded38_minus_dp17", diffExpVs17)},
			{Key: "policy_agreement_vs_dp33", Val: jsonout.Obj{
				{Key: "eval_steps_agree_fraction", Val: evalStepAgree},
				{Key: "eval_steps_diff_count", Val: evalDiffSteps},
				{Key: "eval_steps_total", Val: evalTotSteps},
				{Key: "interior_states_agree_fraction", Val: interiorStateAgree},
				{Key: "interior_states_agree_count", Val: intAgreeCnt},
				{Key: "interior_states_total", Val: intTotStates},
			}},
		}},
		{Key: "paired_comparisons", Val: []jsonout.Obj{
			pairedDiffToJSON("dp17_minus_threshold", diff17VsThresh),
			pairedDiffToJSON("dp33_minus_threshold", diff33VsThresh),
			pairedDiffToJSON("dp33_minus_dp17", diff33Vs17),
		}},
		{Key: "strategies", Val: []jsonout.Obj{
			strategyToJSON(res17, true),
			strategyToJSON(res33, true),
			strategyToJSON(resThresh, true),
			strategyToJSON(resPF, true),
			strategyToJSON(resNoBat, false),
		}},
		{Key: "acceptance_criteria", Val: jsonout.Obj{
			{Key: "c1_stationary_distribution_valid", Val: c1StationaryOK},
			{Key: "c1_status", Val: statusStr(c1StationaryOK)},
			{Key: "c2_jump_and_onestep_moments_verified", Val: c2MomentsOK},
			{Key: "c2_status", Val: statusStr(c2MomentsOK)},
			{Key: "c3_zero_soc_and_power_violations", Val: c3BoundsOK},
			{Key: "c3_status", Val: statusStr(c3BoundsOK)},
			{Key: "c4_bellman_matches_bruteforce_and_residual", Val: c4BellmanOK},
			{Key: "c4_status", Val: statusStr(c4BellmanOK)},
			{Key: "c5_perfect_foresight_upper_bound", Val: c5PFDomOK},
			{Key: "c5_status", Val: statusStr(c5PFDomOK)},
			{Key: "c6_seed_and_worker_determinism", Val: c6DetOK},
			{Key: "c6_status", Val: statusStr(c6DetOK)},
			{Key: "c7_negative_prices_and_spikes_realized", Val: c7NegAndSpikeOK},
			{Key: "c7_status", Val: statusStr(c7NegAndSpikeOK)},
			{Key: "c8_outputs_finite_and_within_budget", Val: c8OutputsFiniteOK},
			{Key: "c8_status", Val: statusStr(c8OutputsFiniteOK)},
			{Key: "all_technical_passed", Val: allTechPassed},
			{Key: "technical_scope_note", Val: "all_technical_passed hanya menyatakan verifikasi teknis c1-c8 di model sintetis dan bukan bukti profitabilitas pasar nyata"},
		}},
	}

	if err := jsonout.WriteFile(outPath, rootObj); err != nil {
		fmt.Fprintf(os.Stderr, "jsonout.WriteFile error: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("day004e selesai dalam %.2fs | workers=%d | N=%d | DP17=%.2f | DP33=%.2f | Thresh=%.2f | PF=%.2f | all_technical_passed=%v\n",
		totalSeconds, nWorkers, nDays, res17.MeanProfit, res33.MeanProfit, resThresh.MeanProfit, resPF.MeanProfit, allTechPassed)
}
