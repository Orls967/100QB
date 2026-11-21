package dp

import (
	"math"
	"time"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/battery"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/prices"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
)

// DPPolicy menyimpan tabel nilai V(h, r, b, s) dan kebijakan optimal A*(h, r, b, s)
// hasil induksi mundur Bellman 24 jam.
type DPPolicy struct {
	Grid               PriceGrid
	Battery            battery.Spec
	NumRegimes         int
	NumBins            int
	NumSoC             int
	Values             []float64          // ukuran 25 * NumRegimes * NumBins * NumSoC
	Actions            []battery.ActionID // ukuran 24 * NumRegimes * NumBins * NumSoC
	MaxBellmanResidual float64
	TerminalMaxError   float64
	SolveSeconds       float64
}

func (p *DPPolicy) valOffset(h, r, b, s int) int {
	return (((h*p.NumRegimes+r)*p.NumBins + b) * p.NumSoC) + s
}

func (p *DPPolicy) actOffset(h, r, b, s int) int {
	return (((h*p.NumRegimes+r)*p.NumBins + b) * p.NumSoC) + s
}

// ValueAt mengembalikan nilai fungsi Bellman V(h, r, b, socIdx).
func (p *DPPolicy) ValueAt(h, r, b, socIdx int) float64 {
	return p.Values[p.valOffset(h, r, b, socIdx)]
}

// ActionAt mengembalikan tindakan optimal A*(h, r, b, socIdx).
func (p *DPPolicy) ActionAt(h, r, b, socIdx int) battery.ActionID {
	return p.Actions[p.actOffset(h, r, b, socIdx)]
}

// SolveStochasticDP menyelesaikan rekursi Bellman 24 jam dengan nilai terminal konservatif yang sama.
func SolveStochasticDP(kernel TransitionKernel, bat battery.Spec) DPPolicy {
	t0 := time.Now()
	nBins := kernel.Grid.NumBins
	nSoC := bat.NumSoCStates
	stride := kernel.TargetStride

	pol := DPPolicy{
		Grid:       kernel.Grid,
		Battery:    bat,
		NumRegimes: regime.NumRegimes,
		NumBins:    nBins,
		NumSoC:     nSoC,
		Values:     make([]float64, (prices.HoursPerDay+1)*regime.NumRegimes*nBins*nSoC),
		Actions:    make([]battery.ActionID, prices.HoursPerDay*regime.NumRegimes*nBins*nSoC),
	}

	// Kondisi batas terminal di h = 24.
	var termVals [battery.NumSoCStates]float64
	for s := 0; s < nSoC; s++ {
		termVals[s] = bat.TerminalValue(s)
	}
	for r := 0; r < regime.NumRegimes; r++ {
		for b := 0; b < nBins; b++ {
			off := pol.valOffset(prices.HoursPerDay, r, b, 0)
			for s := 0; s < nSoC; s++ {
				pol.Values[off+s] = termVals[s]
			}
		}
	}

	// Cache per-action cashflow per bin dan daftar indeks target bernilai > 0.
	evCont := make([]float64, nSoC)

	// Induksi mundur h = 23 .. 0.
	for h := prices.HoursPerDay - 1; h >= 0; h-- {
		nextHourValBase := (h + 1) * stride * nSoC
		curHourValBase := h * stride * nSoC
		curHourActBase := h * stride * nSoC

		for r := 0; r < regime.NumRegimes; r++ {
			for b := 0; b < nBins; b++ {
				stIdx := (h*regime.NumRegimes+r)*nBins + b
				probs := kernel.Probs[stIdx*stride : (stIdx+1)*stride]
				repPrice := kernel.Grid.RepPrices[b]

				// Hitung E[V_{h+1}(r', b', s') | h, r, b] sekali per s' untuk seluruh s dan tindakan.
				for sNext := 0; sNext < nSoC; sNext++ {
					evCont[sNext] = 0.0
				}
				for targetIdx, pr := range probs {
					if pr == 0.0 {
						continue
					}
					vSlice := pol.Values[nextHourValBase+targetIdx*nSoC : nextHourValBase+(targetIdx+1)*nSoC]
					for sNext := 0; sNext < nSoC; sNext++ {
						evCont[sNext] += pr * vSlice[sNext]
					}
				}

				var actCF [battery.NumActions]float64
				for aIdx := 0; aIdx < battery.NumActions; aIdx++ {
					actCF[aIdx] = bat.HourlyCashflow(repPrice, battery.Actions[aIdx])
				}

				tbOffset := (r*nBins + b) * nSoC
				for s := 0; s < nSoC; s++ {
					bestVal := math.Inf(-1)
					bestAct := battery.ActionIdle
					for aIdx := 0; aIdx < battery.NumActions; aIdx++ {
						act := battery.Actions[aIdx]
						sNext, ok := bat.NextSoCIdx(s, act)
						if !ok {
							continue
						}
						cand := actCF[aIdx] + evCont[sNext]
						if cand > bestVal {
							bestVal = cand
							bestAct = act.ID
						}
					}
					pol.Values[curHourValBase+tbOffset+s] = bestVal
					pol.Actions[curHourActBase+tbOffset+s] = bestAct
				}
			}
		}
	}

	// Pemeriksaan konsistensi residual Bellman independen setelah tabel selesai.
	maxTermErr := 0.0
	for r := 0; r < regime.NumRegimes; r++ {
		for b := 0; b < nBins; b++ {
			for s := 0; s < nSoC; s++ {
				err := math.Abs(pol.ValueAt(prices.HoursPerDay, r, b, s) - bat.TerminalValue(s))
				if err > maxTermErr {
					maxTermErr = err
				}
			}
		}
	}
	maxRes := 0.0
	for h := 0; h < prices.HoursPerDay; h++ {
		for r := 0; r < regime.NumRegimes; r++ {
			for b := 0; b < nBins; b++ {
				stIdx := (h*regime.NumRegimes+r)*nBins + b
				probs := kernel.Probs[stIdx*stride : (stIdx+1)*stride]
				repPrice := kernel.Grid.RepPrices[b]
				for s := 0; s < nSoC; s++ {
					bestOpt := math.Inf(-1)
					for aIdx := 0; aIdx < battery.NumActions; aIdx++ {
						act := battery.Actions[aIdx]
						sNext, ok := bat.NextSoCIdx(s, act)
						if !ok {
							continue
						}
						ev := 0.0
						for targetIdx, pr := range probs {
							if pr == 0.0 {
								continue
							}
							rNext := targetIdx / nBins
							bNext := targetIdx % nBins
							ev += pr * pol.ValueAt(h+1, rNext, bNext, sNext)
						}
						cand := bat.HourlyCashflow(repPrice, act) + ev
						if cand > bestOpt {
							bestOpt = cand
						}
					}
					polAct := battery.Actions[pol.ActionAt(h, r, b, s)]
					sNextPol, _ := bat.NextSoCIdx(s, polAct)
					evPol := 0.0
					for targetIdx, pr := range probs {
						if pr == 0.0 {
							continue
						}
						rNext := targetIdx / nBins
						bNext := targetIdx % nBins
						evPol += pr * pol.ValueAt(h+1, rNext, bNext, sNextPol)
					}
					rhsPol := bat.HourlyCashflow(repPrice, polAct) + evPol
					d1 := math.Abs(pol.ValueAt(h, r, b, s) - bestOpt)
					d2 := math.Abs(pol.ValueAt(h, r, b, s) - rhsPol)
					if d1 > maxRes {
						maxRes = d1
					}
					if d2 > maxRes {
						maxRes = d2
					}
				}
			}
		}
	}

	pol.TerminalMaxError = maxTermErr
	pol.MaxBellmanResidual = maxRes
	pol.SolveSeconds = time.Since(t0).Seconds()
	return pol
}

// SolveDeterministicDP menyelesaikan DP mundur eksak untuk satu lintasan harga realisasi
// dengan panjang sembarang H jam (digunakan untuk perfect-foresight DP dan uji brute force).
func SolveDeterministicDP(pricePath []float64, bat battery.Spec) (netProfit, grossTotal float64, actions []battery.ActionID, socPath []int) {
	H := len(pricePath)
	nSoC := bat.NumSoCStates
	V := make([][]float64, H+1)
	actTab := make([][]battery.ActionID, H)
	for h := 0; h <= H; h++ {
		V[h] = make([]float64, nSoC)
	}
	for s := 0; s < nSoC; s++ {
		V[H][s] = bat.TerminalValue(s)
	}
	for h := H - 1; h >= 0; h-- {
		actTab[h] = make([]battery.ActionID, nSoC)
		p := pricePath[h]
		for s := 0; s < nSoC; s++ {
			bestVal := math.Inf(-1)
			bestAct := battery.ActionIdle
			for aIdx := 0; aIdx < battery.NumActions; aIdx++ {
				act := battery.Actions[aIdx]
				sNext, ok := bat.NextSoCIdx(s, act)
				if !ok {
					continue
				}
				cand := bat.HourlyCashflow(p, act) + V[h+1][sNext]
				if cand > bestVal {
					bestVal = cand
					bestAct = act.ID
				}
			}
			V[h][s] = bestVal
			actTab[h][s] = bestAct
		}
	}
	actions = make([]battery.ActionID, H)
	socPath = make([]int, H+1)
	soc := bat.InitialSoCIdx
	socPath[0] = soc
	for h := 0; h < H; h++ {
		aid := actTab[h][soc]
		actions[h] = aid
		soc, _ = bat.NextSoCIdx(soc, battery.Actions[aid])
		socPath[h+1] = soc
	}
	grossTotal = V[0][bat.InitialSoCIdx]
	netProfit = grossTotal - bat.TerminalValue(bat.InitialSoCIdx)
	return netProfit, grossTotal, actions, socPath
}
