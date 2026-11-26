package main

import (
	"fmt"
	"os"
)

func writeBarSVG(outDir, fname, title string, vals []float64) {
	mx := 1.0
	for _, v := range vals {
		if v > mx {
			mx = v
		}
	}
	s := `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>`
	bw := 460.0 / float64(max(len(vals), 1))
	for i, v := range vals {
		h := max((v/mx)*250.0, 1.0)
		s += fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#d62728"/>`, 75.0+float64(i)*bw, 310.0-h, bw*0.72, h)
	}
	s += fmt.Sprintf(`<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">%s</text></svg>`+"\n", title)
	os.WriteFile(outDir+"/figures/"+fname, []byte(s), 0644)
}
func writeS39SVGs(s *Sim, outDir string) {
	writeBarSVG(outDir, "inventory_4sku.svg", "Go Total Sold by SKU (A, B, C, D)", []float64{float64(s.totalSold[0]), float64(s.totalSold[1]), float64(s.totalSold[2]), float64(s.totalSold[3])})
	fr := make([]float64, 10)
	for i := 0; i < 10; i++ {
		fr[i] = 100.0
		if s.dailyFillStore[i][0] > 0 {
			fr[i] = 100.0 * float64(s.dailyFillStore[i][1]) / float64(s.dailyFillStore[i][0])
		}
	}
	writeBarSVG(outDir, "fill_rate_stores.svg", "Go Store Fill Rate (%)", fr)
	ec := make([]float64, 4)
	for d := 0; d < 5; d++ {
		for e := 0; e < 4; e++ {
			ec[e] += float64(s.dailyCostEchelon[d][e])
		}
	}
	writeBarSVG(outDir, "cost_by_echelon.svg", "Go Cost by Echelon", ec)
	var lt []float64
	for i := 0; i < min(len(s.leadtimePairs), 8); i++ {
		lt = append(lt, float64(s.leadtimePairs[i][0]))
	}
	if len(lt) == 0 {
		lt = []float64{1}
	}
	writeBarSVG(outDir, "leadtime_impact.svg", "Go Sample Carrier Lead Times", lt)
	bw := make([]float64, 4)
	for k := 0; k < 4; k++ {
		md := 0.0
		for b := 0; b < 100; b++ {
			md += float64(s.binDemand[k][b])
		}
		md /= 100.0
		m2, m4 := 0.0, 0.0
		for b := 0; b < 100; b++ {
			d := float64(s.binDemand[k][b]) - md
			m2 += d * d
			m4 += d * d * d * d
		}
		if m2 > 0 {
			bw[k] = (m4 / 100.0) / ((m2 / 100.0) * (m2 / 100.0))
		}
	}
	writeBarSVG(outDir, "bullwhip_acf.svg", "Go Demand Kurtosis by SKU", bw)
	var ki []float64
	for i := 0; i < min(len(s.kitImbSeries), 12); i++ {
		ki = append(ki, float64(s.kitImbSeries[i]))
	}
	if len(ki) == 0 {
		ki = []float64{0}
	}
	writeBarSVG(outDir, "kit_imbalance.svg", "Go Kit Component Imbalance", ki)
	sh := make([]float64, 3)
	for i := 0; i < 3; i++ {
		if s.shockCnt[i] > 0 {
			sh[i] = float64(s.shockResp[i]) / float64(s.shockCnt[i])
		}
	}
	writeBarSVG(outDir, "shock_study.svg", "Go Shock Event Study", sh)
	writeBarSVG(outDir, "alloc_vs_cont.svg", "Go Allocation Fill vs Retail Sales", []float64{float64(s.allocTotalFill), float64(s.totalSold[0] + s.totalSold[1] + s.totalSold[2] + s.totalSold[3])})
	writeBarSVG(outDir, "carrier_tiers.svg", "Go Carrier Volume Tiers", []float64{float64(s.carrierTierCnt[0]), float64(s.carrierTierCnt[1]), float64(s.carrierTierCnt[2])})
	frTot := 1.0
	if s.demUnitsTot > 0 {
		frTot = float64(s.demUnitsFilled) / float64(s.demUnitsTot)
	}
	writeBarSVG(outDir, "sweep_safety.svg", "Go Store Fill Rate Summary", []float64{frTot, fr[0] / 100.0, fr[9] / 100.0})
}
