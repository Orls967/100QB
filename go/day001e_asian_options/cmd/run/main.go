package main

import (
	"encoding/json"
	"fmt"
	"image/color"
	"log"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/Orls967/100QB/go/day001e_asian_options/options"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

type sanityCheck struct {
	BSCall float64                 `json:"bs_call"`
	MCCall options.EstimatorResult `json:"mc_call"`
	Diff   float64                 `json:"diff"`
	ZScore float64                 `json:"z_score"`
}

type estimatorEntry struct {
	Name   string   `json:"name"`
	Price  float64  `json:"price"`
	SE     float64  `json:"se"`
	CILow  float64  `json:"ci_low"`
	CIHigh float64  `json:"ci_high"`
	VRF    float64  `json:"vrf"`
	Beta   *float64 `json:"beta,omitempty"`
	YExact *float64 `json:"y_exact,omitempty"`
}

type deltaComparison struct {
	Pathwise options.DeltaResult `json:"pathwise"`
	Bump1Pct options.DeltaResult `json:"bump_1pct"`
	Diff     float64             `json:"diff"`
	ZScore   float64             `json:"z_score"`
}

type convergencePoint struct {
	N                int     `json:"n"`
	SENaive          float64 `json:"se_naive"`
	SEAntithetic     float64 `json:"se_antithetic"`
	SEControlVariate float64 `json:"se_control_variate"`
}

type slopes struct {
	Naive          float64 `json:"naive"`
	Antithetic     float64 `json:"antithetic"`
	ControlVariate float64 `json:"control_variate"`
}

type timingsMs struct {
	SanityCheck    float64 `json:"sanity_check"`
	Naive          float64 `json:"naive"`
	Antithetic     float64 `json:"antithetic"`
	ControlVariate float64 `json:"control_variate"`
	DeltaPathwise  float64 `json:"delta_pathwise"`
	DeltaBump      float64 `json:"delta_bump"`
	Convergence    float64 `json:"convergence"`
	Total200kCore  float64 `json:"total_200k_core"`
}

type fullResults struct {
	Params      options.Params     `json:"params"`
	NPaths      int                `json:"n_paths"`
	Seed        uint64             `json:"seed"`
	SanityCheck sanityCheck        `json:"sanity_check"`
	Estimators  []estimatorEntry   `json:"estimators"`
	Delta       deltaComparison    `json:"delta"`
	Convergence []convergencePoint `json:"convergence"`
	Slopes      slopes             `json:"slopes"`
	TimingsMs   timingsMs          `json:"timings_ms"`
}

func findRepoRoot() (string, error) {
	dir, err := os.Getwd()
	if err != nil {
		return "", err
	}
	for {
		if _, err := os.Stat(filepath.Join(dir, "src", "qb100")); err == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", fmt.Errorf("repository root not found")
		}
		dir = parent
	}
}

func plotConvergence(outPath string, pts []convergencePoint, sl slopes) error {
	p := plot.New()
	p.Title.Text = "Konvergensi Standard Error vs Jumlah Path (Log-Log, Go)"
	p.X.Label.Text = "Jumlah Path N"
	p.Y.Label.Text = "Standard Error (SE)"
	p.X.Scale = plot.LogScale{}
	p.Y.Scale = plot.LogScale{}
	p.X.Tick.Marker = plot.LogTicks{}
	p.Y.Tick.Marker = plot.LogTicks{}
	p.Add(plotter.NewGrid())

	makeXY := func(getY func(convergencePoint) float64) plotter.XYs {
		xys := make(plotter.XYs, len(pts))
		for i, pt := range pts {
			xys[i].X = float64(pt.N)
			xys[i].Y = getY(pt)
		}
		return xys
	}

	series := []struct {
		label string
		col   color.RGBA
		xys   plotter.XYs
	}{
		{fmt.Sprintf("Naive MC (slope=%.3f)", sl.Naive), color.RGBA{R: 31, G: 111, B: 235, A: 255}, makeXY(func(pt convergencePoint) float64 { return pt.SENaive })},
		{fmt.Sprintf("Antithetic (slope=%.3f)", sl.Antithetic), color.RGBA{R: 209, G: 73, B: 91, A: 255}, makeXY(func(pt convergencePoint) float64 { return pt.SEAntithetic })},
		{fmt.Sprintf("Control Variate (slope=%.3f)", sl.ControlVariate), color.RGBA{R: 42, G: 157, B: 143, A: 255}, makeXY(func(pt convergencePoint) float64 { return pt.SEControlVariate })},
	}

	for _, s := range series {
		line, points, err := plotter.NewLinePoints(s.xys)
		if err != nil {
			return err
		}
		line.Color = s.col
		line.Width = vg.Points(2)
		points.Color = s.col
		points.Radius = vg.Points(3)
		p.Add(line, points)
		p.Legend.Add(s.label, line)
	}
	p.Legend.Top = true

	return p.Save(9*vg.Inch, 5.2*vg.Inch, outPath)
}

func plotEstimatorsSE(outPath string, ests []estimatorEntry) error {
	p := plot.New()
	p.Title.Text = "Perbandingan Standard Error Estimator (200k Path, Go)"
	p.Y.Label.Text = "Standard Error (SE)"
	p.Add(plotter.NewGrid())

	vals := make(plotter.Values, len(ests))
	names := make([]string, len(ests))
	for i, e := range ests {
		vals[i] = e.SE
		names[i] = fmt.Sprintf("%s\nSE=%.5f (VRF %.1fx)", e.Name, e.SE, e.VRF)
	}
	bars, err := plotter.NewBarChart(vals, vg.Points(45))
	if err != nil {
		return err
	}
	bars.Color = color.RGBA{R: 31, G: 111, B: 235, A: 220}
	p.Add(bars)
	p.NominalX(names...)

	return p.Save(8.5*vg.Inch, 4.8*vg.Inch, outPath)
}

func main() {
	params := options.Params{
		S0:     100.0,
		K:      100.0,
		R:      0.05,
		Sigma:  0.20,
		T:      1.0,
		NSteps: 252,
	}
	const nPaths = 200_000
	const seed uint64 = 42

	tTotal := time.Now()

	t0 := time.Now()
	bsEuro := options.BSCall(params)
	mcEuro, err := options.PriceEuropeanMC(params, nPaths, seed)
	if err != nil {
		log.Fatal(err)
	}
	tSanity := float64(time.Since(t0).Microseconds()) / 1000.0
	diffEuro := mcEuro.Price - bsEuro
	zEuro := math.Abs(diffEuro) / mcEuro.SE

	t0 = time.Now()
	resNaive, err := options.PriceAsianNaive(params, nPaths, seed)
	if err != nil {
		log.Fatal(err)
	}
	tNaive := float64(time.Since(t0).Microseconds()) / 1000.0

	t0 = time.Now()
	resAnti, err := options.PriceAsianAntithetic(params, nPaths, seed)
	if err != nil {
		log.Fatal(err)
	}
	tAnti := float64(time.Since(t0).Microseconds()) / 1000.0

	t0 = time.Now()
	resCV, err := options.PriceAsianControlVariate(params, nPaths, seed)
	if err != nil {
		log.Fatal(err)
	}
	tCV := float64(time.Since(t0).Microseconds()) / 1000.0

	t0 = time.Now()
	deltaPW, err := options.DeltaAsianPathwise(params, nPaths, seed)
	if err != nil {
		log.Fatal(err)
	}
	tPW := float64(time.Since(t0).Microseconds()) / 1000.0

	t0 = time.Now()
	deltaBmp, err := options.DeltaAsianBump(params, nPaths, 0.01, seed)
	if err != nil {
		log.Fatal(err)
	}
	tBmp := float64(time.Since(t0).Microseconds()) / 1000.0

	coreMs := float64(time.Since(tTotal).Microseconds()) / 1000.0

	deltaDiff := deltaPW.Delta - deltaBmp.Delta
	deltaZ := math.Abs(deltaDiff) / math.Hypot(deltaPW.SE, deltaBmp.SE)

	varNaive := resNaive.SE * resNaive.SE
	betaVal := resCV.Beta
	yExactVal := resCV.YExact
	ests := []estimatorEntry{
		{
			Name:   "Naive MC (200k)",
			Price:  resNaive.Price,
			SE:     resNaive.SE,
			CILow:  resNaive.CILow,
			CIHigh: resNaive.CIHigh,
			VRF:    1.0,
		},
		{
			Name:   "Antithetic (100k pasang)",
			Price:  resAnti.Price,
			SE:     resAnti.SE,
			CILow:  resAnti.CILow,
			CIHigh: resAnti.CIHigh,
			VRF:    varNaive / (resAnti.SE * resAnti.SE),
		},
		{
			Name:   "Control Variate (Geom. Asia)",
			Price:  resCV.Price,
			SE:     resCV.SE,
			CILow:  resCV.CILow,
			CIHigh: resCV.CIHigh,
			VRF:    varNaive / (resCV.SE * resCV.SE),
			Beta:   &betaVal,
			YExact: &yExactVal,
		},
	}

	t0 = time.Now()
	nGrid := []int{1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000}
	convPts := make([]convergencePoint, len(nGrid))
	seN := make([]float64, len(nGrid))
	seA := make([]float64, len(nGrid))
	seC := make([]float64, len(nGrid))
	for i, n := range nGrid {
		rn, err := options.PriceAsianNaive(params, n, seed)
		if err != nil {
			log.Fatal(err)
		}
		ra, err := options.PriceAsianAntithetic(params, n, seed)
		if err != nil {
			log.Fatal(err)
		}
		rc, err := options.PriceAsianControlVariate(params, n, seed)
		if err != nil {
			log.Fatal(err)
		}
		seN[i] = rn.SE
		seA[i] = ra.SE
		seC[i] = rc.SE
		convPts[i] = convergencePoint{
			N:                n,
			SENaive:          rn.SE,
			SEAntithetic:     ra.SE,
			SEControlVariate: rc.SE,
		}
	}
	tConv := float64(time.Since(t0).Microseconds()) / 1000.0

	slN, err := options.FitLogLogSlope(nGrid, seN)
	if err != nil {
		log.Fatal(err)
	}
	slA, err := options.FitLogLogSlope(nGrid, seA)
	if err != nil {
		log.Fatal(err)
	}
	slC, err := options.FitLogLogSlope(nGrid, seC)
	if err != nil {
		log.Fatal(err)
	}
	sl := slopes{Naive: slN, Antithetic: slA, ControlVariate: slC}

	res := fullResults{
		Params: params,
		NPaths: nPaths,
		Seed:   seed,
		SanityCheck: sanityCheck{
			BSCall: bsEuro,
			MCCall: mcEuro,
			Diff:   diffEuro,
			ZScore: zEuro,
		},
		Estimators: ests,
		Delta: deltaComparison{
			Pathwise: deltaPW,
			Bump1Pct: deltaBmp,
			Diff:     deltaDiff,
			ZScore:   deltaZ,
		},
		Convergence: convPts,
		Slopes:      sl,
		TimingsMs: timingsMs{
			SanityCheck:    tSanity,
			Naive:          tNaive,
			Antithetic:     tAnti,
			ControlVariate: tCV,
			DeltaPathwise:  tPW,
			DeltaBump:      tBmp,
			Convergence:    tConv,
			Total200kCore:  coreMs,
		},
	}

	root, err := findRepoRoot()
	if err != nil {
		log.Fatal(err)
	}
	outDir := filepath.Join(root, "outputs", "day001e")
	figDir := filepath.Join(root, "outputs", "figures")
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		log.Fatal(err)
	}
	if err := os.MkdirAll(figDir, 0o755); err != nil {
		log.Fatal(err)
	}

	payload, err := json.MarshalIndent(res, "", "  ")
	if err != nil {
		log.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(outDir, "results.json"), payload, 0o644); err != nil {
		log.Fatal(err)
	}

	if err := plotConvergence(filepath.Join(figDir, "day001e_asian-convergence-loglog.png"), convPts, sl); err != nil {
		log.Fatal(err)
	}
	if err := plotEstimatorsSE(filepath.Join(figDir, "day001e_asian-estimators-ci.png"), ests); err != nil {
		log.Fatal(err)
	}

	fmt.Printf("BS Eropa: %.6f | MC Eropa: %.6f ± %.6f (%.2f SE)\n", bsEuro, mcEuro.Price, mcEuro.SE, zEuro)
	fmt.Printf("Geom. Asian exact: %.6f | beta: %.6f\n", resCV.YExact, resCV.Beta)
	for _, e := range ests {
		fmt.Printf("%-28s: %.6f ± %.6f (CI [%.6f, %.6f]) | VRF = %.2fx\n", e.Name, e.Price, e.SE, e.CILow, e.CIHigh, e.VRF)
	}
	fmt.Printf("Delta PW: %.6f ± %.6f | Bump: %.6f ± %.6f | diff = %+.6f (%.2f SE)\n", deltaPW.Delta, deltaPW.SE, deltaBmp.Delta, deltaBmp.SE, deltaDiff, deltaZ)
	fmt.Printf("Slopes -> Naive: %.4f, Antithetic: %.4f, CV: %.4f\n", sl.Naive, sl.Antithetic, sl.ControlVariate)
	fmt.Printf("Total 200k core time: %.1f ms\n", coreMs)
}
