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

func writeS38SVGs(sNo, sMain *Sim, outDir string) {
	dCnt := int(sMain.rp.days)
	step := dCnt / 12
	if step < 1 {
		step = 1
	}
	epi := make([]float64, 12)
	inc := make([]float64, 12)
	rt := make([]float64, 12)
	hosp := make([]float64, 12)
	ens := make([]float64, 12)
	for i := 0; i < 12; i++ {
		d := i * step
		if d >= dCnt {
			d = dCnt - 1
		}
		epi[i] = float64(sMain.dailyCnt[d][StIA] + sMain.dailyCnt[d][StIS])
		inc[i] = float64(sMain.dailyIncAge[d][0] + sMain.dailyIncAge[d][1] + sMain.dailyIncAge[d][2])
		rt[i] = float64(sNo.dailyCnt[d][StIA] + sNo.dailyCnt[d][StIS])
		hosp[i] = float64(sMain.dailyHosp[d][0])
		ens[i] = float64(sMain.dailyVar[d][1])
	}
	writeBarSVG(outDir, "epi_curves.svg", "Go Active Infectious (Interv)", epi, "#1f77b4")
	writeBarSVG(outDir, "incidence_age.svg", "Go Daily Incidence Sample", inc, "#ff7f0e")
	writeBarSVG(outDir, "rt_series.svg", "Go Active Infectious (No Interv)", rt, "#d62728")
	writeBarSVG(outDir, "hospital_cap.svg", "Go Hospital Occupancy Sample", hosp, "#9467bd")
	writeBarSVG(outDir, "layer_attack.svg", "Go Infections by Layer (HH/SW/Comm)", []float64{float64(sMain.layerInf[0]), float64(sMain.layerInf[1]), float64(sMain.layerInf[2])}, "#2ca02c")
	writeBarSVG(outDir, "offspring_dist.svg", "Go Infections by Age Group", []float64{float64(sMain.ageInf[0]), float64(sMain.ageInf[1]), float64(sMain.ageInf[2])}, "#8c564b")
	writeBarSVG(outDir, "sweep_heatmap.svg", "Go No-Interv vs Interv Total Infections", []float64{float64(sNo.totalTrans), float64(sMain.totalTrans)}, "#e377c2")
	writeBarSVG(outDir, "vaccine_effect.svg", "Go Vaccine Doses &amp; Blocked Transmissions", []float64{float64(sMain.branches[18]), float64(sMain.branches[3])}, "#17becf")
	writeBarSVG(outDir, "variant_share.svg", "Go Variant 1 Daily Incidence", ens, "#bcbd22")
	writeBarSVG(outDir, "ensemble_band.svg", "Go Final Recovered vs Deaths", []float64{float64(sMain.dailyCnt[dCnt-1][StR]), float64(sMain.dailyCnt[dCnt-1][StD])}, "#7f7f7f")
}
