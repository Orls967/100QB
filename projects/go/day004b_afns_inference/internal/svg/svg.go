package svg

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// Palette warna bersih ramah buta-warna (Okabe-Ito / Tableau).
var Palette = []string{
	"#1f77b4",
	"#d62728",
	"#2ca02c",
	"#9467bd",
	"#ff7f0e",
	"#17becf",
}

// Series merepresentasikan satu seri kurva atau titik berbatang-galat pada grafik SVG.
type Series struct {
	Label  string
	X      []float64
	Y      []float64
	Lower  []float64
	Upper  []float64
	Color  string
	Dashed bool
}

// HRefLine adalah garis referensi horizontal konstan (misal 0.95 atau 1.0).
type HRefLine struct {
	Value float64
	Label string
	Color string
}

// ChartOptions mengonfigurasi judul, label sumbu, skala, garis referensi, dan legenda.
type ChartOptions struct {
	Title         string
	XLabel        string
	YLabel        string
	LogX          bool
	LogY          bool
	ZeroLine      bool
	LegendTopLeft bool
	LegendBottom  bool
	YMinOverride  *float64
	YMaxOverride  *float64
	HRefLines     []HRefLine
}

func esc(s string) string {
	s = strings.ReplaceAll(s, "&", "&amp;")
	s = strings.ReplaceAll(s, "<", "&lt;")
	s = strings.ReplaceAll(s, ">", "&gt;")
	return s
}

func tx(v float64, logScale bool) float64 {
	if !logScale {
		return v
	}
	return math.Log10(math.Max(v, 1e-16))
}

func formatTick(v float64, logScale bool) string {
	if logScale {
		return fmt.Sprintf("%.2g", math.Pow(10, v))
	}
	if math.Abs(v) < 1e-12 {
		return "0"
	}
	return fmt.Sprintf("%.4g", v)
}

// WriteLineChart menggambar grafik garis multi-seri dengan pita galat dan garis referensi opsional.
func WriteLineChart(path string, opts ChartOptions, series []Series) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 860.0, 480.0
	ml, mr, mt, mb := 88.0, 32.0, 54.0, 64.0
	pw, ph := w-ml-mr, h-mt-mb

	xmin, xmax := math.Inf(1), math.Inf(-1)
	ymin, ymax := math.Inf(1), math.Inf(-1)
	for _, s := range series {
		for i := range s.X {
			xv := tx(s.X[i], opts.LogX)
			if xv < xmin {
				xmin = xv
			}
			if xv > xmax {
				xmax = xv
			}
		}
		for i, rawY := range s.Y {
			yv := tx(rawY, opts.LogY)
			if yv < ymin {
				ymin = yv
			}
			if yv > ymax {
				ymax = yv
			}
			if len(s.Lower) == len(s.Y) {
				if lv := tx(s.Lower[i], opts.LogY); lv < ymin {
					ymin = lv
				}
			}
			if len(s.Upper) == len(s.Y) {
				if uv := tx(s.Upper[i], opts.LogY); uv > ymax {
					ymax = uv
				}
			}
		}
	}
	for _, hl := range opts.HRefLines {
		hv := tx(hl.Value, opts.LogY)
		if hv < ymin {
			ymin = hv
		}
		if hv > ymax {
			ymax = hv
		}
	}
	if !opts.LogY && ymin > 0 && ymax > 0 {
		ymin = 0
	}
	if opts.YMinOverride != nil {
		ymin = tx(*opts.YMinOverride, opts.LogY)
	}
	if opts.YMaxOverride != nil {
		ymax = tx(*opts.YMaxOverride, opts.LogY)
	}
	if math.IsInf(xmin, 1) || xmax <= xmin {
		xmin, xmax = 0, 1
	}
	if math.IsInf(ymin, 1) || ymax <= ymin {
		ymin--
		ymax++
	}
	ypad := 0.08 * (ymax - ymin)
	ymax += ypad
	if (ymin != 0 || opts.LogY) && opts.YMinOverride == nil {
		ymin -= 0.05 * (ymax - ymin)
	}
	xpad := 0.04 * (xmax - xmin)
	xmin -= xpad
	xmax += xpad

	px := func(rawX float64) float64 {
		return ml + ((tx(rawX, opts.LogX)-xmin)/(xmax-xmin))*pw
	}
	py := func(rawY float64) float64 {
		return mt + ph - ((tx(rawY, opts.LogY)-ymin)/(ymax-ymin))*ph
	}

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="30" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(opts.Title))
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)

	for k := 0; k <= 5; k++ {
		frac := float64(k) / 5.0
		gx := ml + frac*pw
		xv := xmin + frac*(xmax-xmin)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", gx, mt, gx, mt+ph)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#444">%s</text>`+"\n", gx, mt+ph+18, esc(formatTick(xv, opts.LogX)))

		gy := mt + ph - frac*ph
		yv := ymin + frac*(ymax-ymin)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", ml, gy, ml+pw, gy)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" fill="#444">%s</text>`+"\n", ml-8, gy+4, esc(formatTick(yv, opts.LogY)))
	}

	if opts.ZeroLine && !opts.LogY && ymin < 0 && ymax > 0 {
		y0 := py(0)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888" stroke-width="1" stroke-dasharray="3,3"/>`+"\n", ml, y0, ml+pw, y0)
	}
	for _, hl := range opts.HRefLines {
		yh := py(hl.Value)
		col := hl.Color
		if col == "" {
			col = "#555"
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.3" stroke-dasharray="5,4"/>`+"\n", ml, yh, ml+pw, yh, col)
		if hl.Label != "" {
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>`+"\n", ml+6, yh-5, col, esc(hl.Label))
		}
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))

	for idx, s := range series {
		col := s.Color
		if col == "" {
			col = Palette[idx%len(Palette)]
		}
		if len(s.Lower) == len(s.Y) && len(s.Upper) == len(s.Y) {
			for i := range s.X {
				cx := px(s.X[i])
				yLo := py(s.Lower[i])
				yHi := py(s.Upper[i])
				fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.4" stroke-opacity="0.8"/>`+"\n", cx, yLo, cx, yHi, col)
				fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.4"/>`+"\n", cx-4, yLo, cx+4, yLo, col)
				fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.4"/>`+"\n", cx-4, yHi, cx+4, yHi, col)
			}
		}
		var pts strings.Builder
		for i := range s.X {
			fmt.Fprintf(&pts, "%.2f,%.2f ", px(s.X[i]), py(s.Y[i]))
		}
		dash := ""
		if s.Dashed {
			dash = ` stroke-dasharray="6,4"`
		}
		fmt.Fprintf(&sb, `<polyline fill="none" stroke="%s" stroke-width="2.2"%s points="%s"/>`+"\n", col, dash, pts.String())
		if len(s.X) <= 16 {
			for i := range s.X {
				fmt.Fprintf(&sb, `<circle cx="%.2f" cy="%.2f" r="3.4" fill="%s"/>`+"\n", px(s.X[i]), py(s.Y[i]), col)
			}
		}
	}

	legW := 284.0
	legH := float64(len(series))*18.0 + 12.0
	legX := ml + pw - legW - 10.0
	if opts.LegendTopLeft {
		legX = ml + 12.0
	}
	legY := mt + 10.0
	if opts.LegendBottom {
		legY = mt + ph - legH - 10.0
	}
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffffff" fill-opacity="0.92" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legW, legH)
	for i, s := range series {
		col := s.Color
		if col == "" {
			col = Palette[i%len(Palette)]
		}
		lx := legX + 10.0
		ly := legY + 15.0 + float64(i)*18.0
		dash := ""
		if s.Dashed {
			dash = ` stroke-dasharray="5,3"`
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.4"%s/>`+"\n", lx, ly, lx+22, ly, col, dash)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="11" fill="#222">%s</text>`+"\n", lx+28, ly+4, esc(s.Label))
	}

	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}

// CategoryPointSeries merepresentasikan titik estimasi dan interval ketidakpastian pada kategori diskret.
type CategoryPointSeries struct {
	Label  string
	Values []float64
	Lower  []float64
	Upper  []float64
	Color  string
}

// WriteCategoryPointChart menggambar grafik titik berkelompok + error bar pada kategori diskret
// (misal cakupan Wald + Wilson 95% CI atau rasio SE Hessian/SD empiris per parameter).
func WriteCategoryPointChart(path string, opts ChartOptions, categories []string, series []CategoryPointSeries) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 900.0, 480.0
	ml, mr, mt, mb := 82.0, 30.0, 54.0, 78.0
	pw, ph := w-ml-mr, h-mt-mb

	ymin, ymax := math.Inf(1), math.Inf(-1)
	for _, s := range series {
		for i, v := range s.Values {
			if v < ymin {
				ymin = v
			}
			if v > ymax {
				ymax = v
			}
			if len(s.Lower) == len(s.Values) && s.Lower[i] < ymin {
				ymin = s.Lower[i]
			}
			if len(s.Upper) == len(s.Values) && s.Upper[i] > ymax {
				ymax = s.Upper[i]
			}
		}
	}
	for _, hl := range opts.HRefLines {
		if hl.Value < ymin {
			ymin = hl.Value
		}
		if hl.Value > ymax {
			ymax = hl.Value
		}
	}
	if opts.YMinOverride != nil {
		ymin = *opts.YMinOverride
	}
	if opts.YMaxOverride != nil {
		ymax = *opts.YMaxOverride
	}
	if math.IsInf(ymin, 1) || ymax <= ymin {
		ymin, ymax = 0, 1
	}
	ypad := 0.06 * (ymax - ymin)
	if opts.YMaxOverride == nil {
		ymax += ypad
	}
	if opts.YMinOverride == nil {
		ymin -= ypad
	}

	py := func(v float64) float64 {
		return mt + ph - ((v-ymin)/(ymax-ymin))*ph
	}

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="30" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(opts.Title))
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)

	for k := 0; k <= 5; k++ {
		frac := float64(k) / 5.0
		gy := mt + ph - frac*ph
		yv := ymin + frac*(ymax-ymin)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", ml, gy, ml+pw, gy)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" fill="#444">%s</text>`+"\n", ml-8, gy+4, esc(formatTick(yv, false)))
	}

	for _, hl := range opts.HRefLines {
		yh := py(hl.Value)
		col := hl.Color
		if col == "" {
			col = "#555"
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4" stroke-dasharray="5,4"/>`+"\n", ml, yh, ml+pw, yh, col)
		if hl.Label != "" {
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>`+"\n", ml+6, yh-5, col, esc(hl.Label))
		}
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))

	nCat := len(categories)
	nSer := len(series)
	if nCat > 0 && nSer > 0 {
		groupW := pw / float64(nCat)
		stepX := (groupW * 0.70) / math.Max(float64(nSer), 1.0)
		for cIdx, cat := range categories {
			gx := ml + float64(cIdx)*groupW
			if cIdx > 0 {
				fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#f0f0f0" stroke-width="1"/>`+"\n", gx, mt, gx, mt+ph)
			}
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" fill="#333">%s</text>`+"\n", gx+groupW/2, mt+ph+20, esc(cat))
			for sIdx, s := range series {
				if cIdx >= len(s.Values) {
					continue
				}
				col := s.Color
				if col == "" {
					col = Palette[sIdx%len(Palette)]
				}
				cx := gx + groupW*0.15 + (float64(sIdx)+0.5)*stepX
				yv := py(s.Values[cIdx])
				if len(s.Lower) == len(s.Values) && len(s.Upper) == len(s.Values) {
					yLo := py(s.Lower[cIdx])
					yHi := py(s.Upper[cIdx])
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.6"/>`+"\n", cx, yLo, cx, yHi, col)
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx-3.5, yLo, cx+3.5, yLo, col)
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx-3.5, yHi, cx+3.5, yHi, col)
				}
				fmt.Fprintf(&sb, `<circle cx="%.2f" cy="%.2f" r="3.6" fill="%s"/>`+"\n", cx, yv, col)
			}
		}
	}

	legW := 260.0
	legH := float64(len(series))*18.0 + 12.0
	legX := ml + pw - legW - 10.0
	if opts.LegendTopLeft {
		legX = ml + 12.0
	}
	legY := mt + 10.0
	if opts.LegendBottom {
		legY = mt + ph - legH - 10.0
	}
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffffff" fill-opacity="0.92" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legW, legH)
	for i, s := range series {
		col := s.Color
		if col == "" {
			col = Palette[i%len(Palette)]
		}
		lx := legX + 10.0
		ly := legY + 15.0 + float64(i)*18.0
		fmt.Fprintf(&sb, `<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>`+"\n", lx+8, ly, col)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="11" fill="#222">%s</text>`+"\n", lx+20, ly+4, esc(s.Label))
	}

	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}

// BarSeries merepresentasikan satu kelompok batang pada grafik batang berkelompok.
type BarSeries struct {
	Label  string
	Values []float64
	Color  string
}

// WriteGroupedBarChart menggambar grafik batang berkelompok untuk kategori diskret.
func WriteGroupedBarChart(path string, opts ChartOptions, categories []string, series []BarSeries) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 880.0, 470.0
	ml, mr, mt, mb := 88.0, 32.0, 54.0, 74.0
	pw, ph := w-ml-mr, h-mt-mb

	ymin, ymax := 0.0, 0.0
	hasVal := false
	for _, s := range series {
		for _, v := range s.Values {
			if !hasVal {
				ymin, ymax = math.Min(0, v), math.Max(0, v)
				hasVal = true
			} else {
				if v < ymin {
					ymin = v
				}
				if v > ymax {
					ymax = v
				}
			}
		}
	}
	for _, hl := range opts.HRefLines {
		if hl.Value < ymin {
			ymin = hl.Value
		}
		if hl.Value > ymax {
			ymax = hl.Value
		}
	}
	if opts.YMinOverride != nil {
		ymin = *opts.YMinOverride
	}
	if opts.YMaxOverride != nil {
		ymax = *opts.YMaxOverride
	}
	if ymax <= ymin {
		ymin, ymax = -1, 1
	}
	ypad := 0.10 * (ymax - ymin)
	ymax += ypad
	if ymin < 0 {
		ymin -= ypad
	}

	py := func(v float64) float64 {
		return mt + ph - ((v-ymin)/(ymax-ymin))*ph
	}

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="30" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(opts.Title))
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)

	for k := 0; k <= 5; k++ {
		frac := float64(k) / 5.0
		gy := mt + ph - frac*ph
		yv := ymin + frac*(ymax-ymin)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", ml, gy, ml+pw, gy)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" fill="#444">%s</text>`+"\n", ml-8, gy+4, esc(formatTick(yv, false)))
	}

	yZero := py(0)
	if ymin <= 0 && ymax >= 0 {
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#666" stroke-width="1"/>`+"\n", ml, yZero, ml+pw, yZero)
	}
	for _, hl := range opts.HRefLines {
		yh := py(hl.Value)
		col := hl.Color
		if col == "" {
			col = "#555"
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4" stroke-dasharray="5,4"/>`+"\n", ml, yh, ml+pw, yh, col)
		if hl.Label != "" {
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>`+"\n", ml+6, yh-5, col, esc(hl.Label))
		}
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))

	nCat := len(categories)
	nSer := len(series)
	if nCat > 0 && nSer > 0 {
		groupW := pw / float64(nCat)
		barW := (groupW * 0.72) / float64(nSer)
		offset := groupW * 0.14
		for cIdx, cat := range categories {
			gx := ml + float64(cIdx)*groupW
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" fill="#333">%s</text>`+"\n", gx+groupW/2, mt+ph+20, esc(cat))
			for sIdx, s := range series {
				if cIdx >= len(s.Values) {
					continue
				}
				col := s.Color
				if col == "" {
					col = Palette[sIdx%len(Palette)]
				}
				v := s.Values[cIdx]
				bx := gx + offset + float64(sIdx)*barW
				yv := py(v)
				top := math.Min(yv, yZero)
				bh := math.Abs(yZero - yv)
				fmt.Fprintf(&sb, `<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s" fill-opacity="0.88"/>`+"\n", bx, top, math.Max(barW-1.5, 2.0), bh, col)
			}
		}
	}

	legW := 268.0
	legH := float64(len(series))*18.0 + 12.0
	legX := ml + pw - legW - 10.0
	if opts.LegendTopLeft {
		legX = ml + 12.0
	}
	legY := mt + 10.0
	if opts.LegendBottom {
		legY = mt + ph - legH - 10.0
	}
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffffff" fill-opacity="0.92" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legW, legH)
	for i, s := range series {
		col := s.Color
		if col == "" {
			col = Palette[i%len(Palette)]
		}
		lx := legX + 10.0
		ly := legY + 12.0 + float64(i)*18.0
		fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="14" height="10" fill="%s"/>`+"\n", lx, ly-4, col)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="11" fill="#222">%s</text>`+"\n", lx+20, ly+5, esc(s.Label))
	}

	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}
