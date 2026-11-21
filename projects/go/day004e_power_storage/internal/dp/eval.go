package dp

import (
	"math"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/battery"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/prices"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/stats"
)

const (
	StrategyDP17             = "stochastic_dp_17bin"
	StrategyDP33             = "stochastic_dp_33bin"
	StrategyThreshold        = "threshold_policy"
	StrategyPerfectForesight = "perfect_foresight_dp"
	StrategyNoBattery        = "no_battery"

	ThresholdChargeMaxPrice    = 0.0
	ThresholdDischargeMinPrice = 200.0
)

// DayDispatchTrace mencatat jejak dispatch 24 jam satu hari contoh untuk satu strategi.
type DayDispatchTrace struct {
	DayIdx      int
	Actions     [prices.HoursPerDay]battery.ActionID
	PowerMW     [prices.HoursPerDay]float64
	SoCMWh      [prices.HoursPerDay + 1]float64
	Cashflows   [prices.HoursPerDay]float64
	NetProfit   float64
	GrossProfit float64
}

// StrategyResult menyimpan hasil per hari dan ringkasan ekonomi lengkap satu strategi.
type StrategyResult struct {
	Name                     string
	DailyProfits             []float64 // cashflow_sum + Terminal(soc_24) - Terminal(soc_0), EUR
	DailyGrossValues         []float64 // cashflow_sum + Terminal(soc_24), EUR
	DailyCashflows           []float64 // cashflow_sum saja, EUR
	DailyGridThroughputMWh   []float64
	DailyStoredThroughputMWh []float64
	DailyEFC                 []float64
	DailyRegretVsPF          []float64

	MeanProfit     float64
	MedianProfit   float64
	SDProfit       float64
	MCSEProfit     float64
	VaR5Profit     float64
	ES5Profit      float64
	MeanGrossValue float64
	MeanCashflow   float64
	MeanThroughput float64
	MeanEFC        float64
	MeanRegretVsPF float64
	MaxRegretVsPF  float64
	MCSERegretVsPF float64

	TotalChargeActions    int
	TotalDischargeActions int
	TotalIdleActions      int
	Charge2MWCount        int
	Charge1MWCount        int
	Discharge09MWCount    int
	Discharge18MWCount    int

	TotalProfitNegPriceHours      float64
	MeanProfitNegPricePerDay      float64
	TotalProfitScarcityOrJumpHrs  float64
	MeanProfitScarcityOrJumpPerD  float64
	SoCViolationCount             int
	PowerViolationCount           int
	PerfectForesightDomViolations int
	MaxDominationShortfall        float64
}

func runSingleDayFromDecision(day prices.DayPath, bat battery.Spec, decide func(h, socIdx int) battery.ActionID) (
	netProfit, grossVal, cashSum, gridThrough, storedThrough, negHourCF, scarJumpCF float64,
	socViol, powViol int,
	actionCounts [battery.NumActions]int,
	trace DayDispatchTrace,
) {
	soc := bat.InitialSoCIdx
	trace.DayIdx = day.DayIdx
	trace.SoCMWh[0] = battery.SoCMWh(soc)

	for h := 0; h < prices.HoursPerDay; h++ {
		aid := decide(h, soc)
		act := battery.Actions[aid]
		if !bat.IsFeasible(soc, act) {
			if act.PowerMW > bat.MaxChargeMW+1e-12 || (-act.PowerMW) > bat.MaxDischargeMW+1e-12 {
				powViol++
			} else {
				socViol++
			}
			aid = battery.ActionIdle
			act = battery.Actions[aid]
		}
		p := day.Prices[h]
		cf := bat.HourlyCashflow(p, act)
		cashSum += cf
		gridThrough += act.GridThroughMW
		storedThrough += math.Abs(act.StoredDelta)
		if p < 0.0 {
			negHourCF += cf
		}
		if day.Regimes[h] == regime.Scarcity || day.NumJumps[h] > 0 {
			scarJumpCF += cf
		}
		actionCounts[aid]++
		soc += act.DeltaSoCIdx
		if soc < 0 || soc >= bat.NumSoCStates {
			socViol++
		}
		trace.Actions[h] = aid
		trace.PowerMW[h] = act.PowerMW
		trace.SoCMWh[h+1] = battery.SoCMWh(soc)
		trace.Cashflows[h] = cf
	}
	termEnd := bat.TerminalValue(soc)
	termStart := bat.TerminalValue(bat.InitialSoCIdx)
	grossVal = cashSum + termEnd
	netProfit = grossVal - termStart
	trace.NetProfit = netProfit
	trace.GrossProfit = grossVal
	return
}

func finalizeStrategyResult(name string, days []prices.DayPath, bat battery.Spec, decideDay func(daySliceIdx int, day prices.DayPath, h, socIdx int) battery.ActionID) StrategyResult {
	n := len(days)
	res := StrategyResult{
		Name:                     name,
		DailyProfits:             make([]float64, n),
		DailyGrossValues:         make([]float64, n),
		DailyCashflows:           make([]float64, n),
		DailyGridThroughputMWh:   make([]float64, n),
		DailyStoredThroughputMWh: make([]float64, n),
		DailyEFC:                 make([]float64, n),
		DailyRegretVsPF:          make([]float64, n),
	}
	for i, day := range days {
		netP, grossV, cashS, gTh, sTh, negCF, sjCF, sv, pv, cnts, _ := runSingleDayFromDecision(day, bat, func(h, socIdx int) battery.ActionID {
			return decideDay(i, day, h, socIdx)
		})
		res.DailyProfits[i] = netP
		res.DailyGrossValues[i] = grossV
		res.DailyCashflows[i] = cashS
		res.DailyGridThroughputMWh[i] = gTh
		res.DailyStoredThroughputMWh[i] = sTh
		res.DailyEFC[i] = bat.EquivalentFullCycles(sTh)
		res.TotalProfitNegPriceHours += negCF
		res.TotalProfitScarcityOrJumpHrs += sjCF
		res.SoCViolationCount += sv
		res.PowerViolationCount += pv
		res.Charge2MWCount += cnts[battery.ActionCharge2MW]
		res.Charge1MWCount += cnts[battery.ActionCharge1MW]
		res.TotalIdleActions += cnts[battery.ActionIdle]
		res.Discharge09MWCount += cnts[battery.ActionDischarge09MW]
		res.Discharge18MWCount += cnts[battery.ActionDischarge18MW]
		res.TotalChargeActions += cnts[battery.ActionCharge2MW] + cnts[battery.ActionCharge1MW]
		res.TotalDischargeActions += cnts[battery.ActionDischarge09MW] + cnts[battery.ActionDischarge18MW]
	}

	fn := math.Max(float64(n), 1.0)
	res.MeanProfit = stats.Mean(res.DailyProfits)
	res.MedianProfit = stats.Quantile(res.DailyProfits, 0.5)
	res.SDProfit = stats.SD(res.DailyProfits)
	res.MCSEProfit = stats.SE(res.DailyProfits)
	res.VaR5Profit, res.ES5Profit = stats.VaRAndES(res.DailyProfits, 0.05)
	res.MeanGrossValue = stats.Mean(res.DailyGrossValues)
	res.MeanCashflow = stats.Mean(res.DailyCashflows)
	res.MeanThroughput = stats.Mean(res.DailyGridThroughputMWh)
	res.MeanEFC = stats.Mean(res.DailyEFC)
	res.MeanProfitNegPricePerDay = res.TotalProfitNegPriceHours / fn
	res.MeanProfitScarcityOrJumpPerD = res.TotalProfitScarcityOrJumpHrs / fn
	return res
}

// EvalStochasticDP mengevaluasi kebijakan DP stokastik secara nonanticipative di atas lintasan harga kontinu.
func EvalStochasticDP(name string, days []prices.DayPath, policy DPPolicy, bat battery.Spec) StrategyResult {
	return finalizeStrategyResult(name, days, bat, func(_ int, day prices.DayPath, h, socIdx int) battery.ActionID {
		b := policy.Grid.BinIndex(day.Prices[h])
		return policy.ActionAt(h, day.Regimes[h], b, socIdx)
	})
}

// PolicyActionAgreement mengukur kesepakatan aksi antara dua kebijakan DP pada lintasan evaluasi
// (seluruh langkah jam) serta pada state interior yang memiliki harga representatif identik.
func PolicyActionAgreement(days []prices.DayPath, polRef, polOther DPPolicy, bat battery.Spec) (evalStepAgreeFrac, interiorStateAgreeFrac float64, evalDiffCount, evalTotalSteps, interiorAgreeCount, interiorTotalStates int) {
	for _, day := range days {
		socRef := bat.InitialSoCIdx
		socOth := bat.InitialSoCIdx
		for h := 0; h < prices.HoursPerDay; h++ {
			bRef := polRef.Grid.BinIndex(day.Prices[h])
			bOth := polOther.Grid.BinIndex(day.Prices[h])
			aRef := polRef.ActionAt(h, day.Regimes[h], bRef, socRef)
			aOth := polOther.ActionAt(h, day.Regimes[h], bOth, socOth)
			if aRef == aOth {
				evalAgree := true
				_ = evalAgree
			} else {
				evalDiffCount++
			}
			evalTotalSteps++
			if sn, ok := bat.NextSoCIdx(socRef, battery.Actions[aRef]); ok {
				socRef = sn
			}
			if sn, ok := bat.NextSoCIdx(socOth, battery.Actions[aOth]); ok {
				socOth = sn
			}
		}
	}
	if evalTotalSteps > 0 {
		evalStepAgreeFrac = float64(evalTotalSteps-evalDiffCount) / float64(evalTotalSteps)
	}

	// Cocokkan bin interior polRef [1..NumBins-2] terhadap polOther dengan RepPrices identik.
	for bRef := 1; bRef < polRef.NumBins-1; bRef++ {
		rep := polRef.Grid.RepPrices[bRef]
		bOth := polOther.Grid.BinIndex(rep)
		if math.Abs(polOther.Grid.RepPrices[bOth]-rep) > 1e-12 {
			continue
		}
		for h := 0; h < prices.HoursPerDay; h++ {
			for r := 0; r < regime.NumRegimes; r++ {
				for s := 0; s < bat.NumSoCStates; s++ {
					interiorTotalStates++
					if polRef.ActionAt(h, r, bRef, s) == polOther.ActionAt(h, r, bOth, s) {
						interiorAgreeCount++
					}
				}
			}
		}
	}
	if interiorTotalStates > 0 {
		interiorStateAgreeFrac = float64(interiorAgreeCount) / float64(interiorTotalStates)
	}
	return
}

// ThresholdAction menentukan tindakan kebijakan aturan tetap pra-run dengan pembatas SoC:
// charge penuh saat harga <= 0, discharge penuh saat harga >= 200, idle selain itu.
func ThresholdAction(price float64, socIdx int, bat battery.Spec) battery.ActionID {
	if price <= ThresholdChargeMaxPrice {
		if bat.IsFeasible(socIdx, battery.Actions[battery.ActionCharge2MW]) {
			return battery.ActionCharge2MW
		}
		if bat.IsFeasible(socIdx, battery.Actions[battery.ActionCharge1MW]) {
			return battery.ActionCharge1MW
		}
		return battery.ActionIdle
	}
	if price >= ThresholdDischargeMinPrice {
		if bat.IsFeasible(socIdx, battery.Actions[battery.ActionDischarge18MW]) {
			return battery.ActionDischarge18MW
		}
		if bat.IsFeasible(socIdx, battery.Actions[battery.ActionDischarge09MW]) {
			return battery.ActionDischarge09MW
		}
		return battery.ActionIdle
	}
	return battery.ActionIdle
}

// EvalThresholdPolicy mengevaluasi kebijakan threshold tetap pra-run.
func EvalThresholdPolicy(days []prices.DayPath, bat battery.Spec) StrategyResult {
	return finalizeStrategyResult(StrategyThreshold, days, bat, func(_ int, day prices.DayPath, h, socIdx int) battery.ActionID {
		return ThresholdAction(day.Prices[h], socIdx, bat)
	})
}

// EvalPerfectForesightDP menghitung batas atas bersyarat per lintasan dengan DP deterministik 24 jam penuh.
func EvalPerfectForesightDP(days []prices.DayPath, bat battery.Spec) StrategyResult {
	pfActions := make([][prices.HoursPerDay]battery.ActionID, len(days))
	for i, day := range days {
		_, _, acts, _ := SolveDeterministicDP(day.Prices[:], bat)
		copy(pfActions[i][:], acts)
	}
	return finalizeStrategyResult(StrategyPerfectForesight, days, bat, func(daySliceIdx int, _ prices.DayPath, h, socIdx int) battery.ActionID {
		return pfActions[daySliceIdx][h]
	})
}

// EvalNoBattery mengevaluasi benchmark nol tanpa baterai (seluruh cashflow dan nilai terminal nol).
func EvalNoBattery(days []prices.DayPath) StrategyResult {
	n := len(days)
	return StrategyResult{
		Name:                     StrategyNoBattery,
		DailyProfits:             make([]float64, n),
		DailyGrossValues:         make([]float64, n),
		DailyCashflows:           make([]float64, n),
		DailyGridThroughputMWh:   make([]float64, n),
		DailyStoredThroughputMWh: make([]float64, n),
		DailyEFC:                 make([]float64, n),
		DailyRegretVsPF:          make([]float64, n),
		TotalIdleActions:         n * prices.HoursPerDay,
	}
}

// ComputeRegretVsPF menghitung regret per hari (PF_profit - strategy_profit) dan memeriksa dominasi perfect foresight.
func ComputeRegretVsPF(res *StrategyResult, pf StrategyResult, tol float64) {
	n := len(res.DailyProfits)
	maxReg := math.Inf(-1)
	maxShort := 0.0
	viol := 0
	for i := 0; i < n; i++ {
		reg := pf.DailyProfits[i] - res.DailyProfits[i]
		res.DailyRegretVsPF[i] = reg
		if reg > maxReg {
			maxReg = reg
		}
		if reg < -tol {
			viol++
			if -reg > maxShort {
				maxShort = -reg
			}
		}
	}
	if n == 0 {
		maxReg = 0
	}
	res.MeanRegretVsPF = stats.Mean(res.DailyRegretVsPF)
	res.MaxRegretVsPF = maxReg
	res.MCSERegretVsPF = stats.SE(res.DailyRegretVsPF)
	res.PerfectForesightDomViolations = viol
	res.MaxDominationShortfall = maxShort
}

// TraceExampleDay mengekstrak lintasan aksi dan SoC satu hari contoh untuk strategi tertentu.
func TraceExampleDay(day prices.DayPath, bat battery.Spec, policy *DPPolicy, useThreshold, usePerfectForesight bool) DayDispatchTrace {
	var pfActs []battery.ActionID
	if usePerfectForesight {
		_, _, pfActs, _ = SolveDeterministicDP(day.Prices[:], bat)
	}
	_, _, _, _, _, _, _, _, _, _, trace := runSingleDayFromDecision(day, bat, func(h, socIdx int) battery.ActionID {
		if usePerfectForesight {
			return pfActs[h]
		}
		if useThreshold {
			return ThresholdAction(day.Prices[h], socIdx, bat)
		}
		if policy != nil {
			b := policy.Grid.BinIndex(day.Prices[h])
			return policy.ActionAt(h, day.Regimes[h], b, socIdx)
		}
		return battery.ActionIdle
	})
	return trace
}
