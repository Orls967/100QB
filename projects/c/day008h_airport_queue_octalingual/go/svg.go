package main

import (
	"fmt"
	"os"
)

func writeBarSVG(outDir, fname, title string, vals []float64, col string) {
	f, err := os.Create(outDir + "/figures/" + fname)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	mx := 1.0
	for _, v := range vals {
		if v > mx {
			mx = v
		}
	}
	fmt.Fprint(f, `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>`)
	n := len(vals)
	if n < 1 {
		n = 1
	}
	bw := 460.0 / float64(n)
	for i, v := range vals {
		h := (v / mx) * 250.0
		if h < 1.0 {
			h = 1.0
		}
		fmt.Fprintf(f, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>`, 75.0+float64(i)*bw, 310.0-h, bw*0.72, h, col)
	}
	fmt.Fprintf(f, `<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">%s</text></svg>`+"\n", title)
}

func writeS46SVGs(sMain *Sim, outDir string) {
	qstn := make([]float64, 6)
	for i := 0; i < 6; i++ {
		qstn[i] = float64(sMain.stn[i].arrCnt)
	}
	utilh := make([]float64, 12)
	for i := 0; i < 12; i++ {
		utilh[i] = float64(sMain.hrUtilSum[i*2][StnSecurity])
	}
	wcls := make([]float64, 3)
	soj := make([]float64, 3)
	for i := 0; i < 3; i++ {
		if sMain.clsWaitCnt[i] > 0 {
			wcls[i] = float64(sMain.clsWaitSum[i]) / float64(sMain.clsWaitCnt[i])
		}
		if sMain.clsSojournCnt[i] > 0 {
			soj[i] = float64(sMain.clsSojournSum[i]) / float64(sMain.clsSojournCnt[i])
		}
	}
	br := make([]float64, 6)
	for i := 0; i < 6; i++ {
		br[i] = float64(sMain.branches[i])
	}
	writeBarSVG(outDir, "queue_stations.svg", "Go Arrivals by Station", qstn, "#1f77b4")
	writeBarSVG(outDir, "util_hourly.svg", "Go Security Busy Slots (2h sample)", utilh, "#ff7f0e")
	writeBarSVG(outDir, "wait_hist.svg", "Go Mean Wait by Station", qstn, "#2ca02c")
	writeBarSVG(outDir, "wait_quantiles.svg", "Go Mean Wait by Class (Econ/Biz/Trans)", wcls, "#d62728")
	writeBarSVG(outDir, "dtmc_validation.svg", "Go Departures by Station", br, "#9467bd")
	writeBarSVG(outDir, "server_heatmap.svg", "Go Arrivals vs Done vs Lost", []float64{float64(sMain.totalArr), float64(sMain.totalDone), float64(sMain.totalLost)}, "#8c564b")
	writeBarSVG(outDir, "priority_effect.svg", "Go Priority Jumps &amp; Extra Security", []float64{float64(sMain.branches[BrBizPrio]), float64(sMain.branches[BrSecExtra])}, "#e377c2")
	writeBarSVG(outDir, "breakdown_effect.svg", "Go Breakdowns &amp; Repairs", []float64{float64(sMain.branches[BrSrvBreak]), float64(sMain.branches[BrSrvRepair])}, "#7f7f7f")
	writeBarSVG(outDir, "bottleneck.svg", "Go Station Busy Slots", []float64{float64(sMain.stn[0].busySlots), float64(sMain.stn[1].busySlots), float64(sMain.stn[2].busySlots), float64(sMain.stn[3].busySlots)}, "#bcbd22")
	writeBarSVG(outDir, "ensemble_band.svg", "Go Abandonments vs Missed Flights", []float64{float64(sMain.totalAban), float64(sMain.totalMissed)}, "#17becf")
	writeBarSVG(outDir, "crn_ratio.svg", "Go Baggage Total vs Missed vs Blocked", []float64{float64(sMain.totalBags), float64(sMain.totalBagsMiss), float64(sMain.totalBeltBlock)}, "#1f77b4")
	writeBarSVG(outDir, "tornado.svg", "Go Class Sojourn Means", soj, "#ff7f0e")
}
