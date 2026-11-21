package svg

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// Palette warna bersih ramah buta-warna.
var Palette = []string{
	"#1f77b4",
	"#d62728",
	"#2ca02c",
	"#9467bd",
	"#ff7f0e",
	"#17becf",
}

// HRefLine menggambar garis referensi horizontal konstan.
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

// CategoryPointSeries merepresentasikan titik estimasi beserta pita ketidakpastian pada kategori diskret.
type CategoryPointSeries struct {
	Label  string
	Values []float64
	Lower  []float64
	Upper  []float64
	Color  string
}

// WriteCategoryPointChart menggambar grafik titik berkelompok dengan error bar atau selang kepercayaan.
func WriteCategoryPointChart(path string, opts ChartOptions, categories []string, series []CategoryPointSeries) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 900.0, 480.0
	ml, mr, mt, mb := 84.0, 32.0, 54.0, 78.0
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
	if opts.ZeroLine {
		if 0 < ymin {
			ymin = 0
		}
		if 0 > ymax {
			ymax = 0
		}
	}
	if opts.YMinOverride != nil {
		ymin = *opts.YMinOverride
	}
	if opts.YMaxOverride != nil {
		ymax = *opts.YMaxOverride
	}
	if math.IsInf(ymin, 1) || ymax <= ymin {
		ymin, ymax = -0.1, 0.1
	}
	ypad := 0.08 * (ymax - ymin)
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

	if opts.ZeroLine && ymin < 0 && ymax > 0 {
		y0 := py(0)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#777" stroke-width="1.2" stroke-dasharray="3,3"/>`+"\n", ml, y0, ml+pw, y0)
	}
	for _, hl := range opts.HRefLines {
		yh := py(hl.Value)
		col := hl.Color
		if col == "" {
			col = "#444"
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4" stroke-dasharray="5,4"/>`+"\n", ml, yh, ml+pw, yh, col)
		if hl.Label != "" {
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>`+"\n", ml+8, yh-5, col, esc(hl.Label))
		}
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))

	nCat := len(categories)
	nSer := len(series)
	if nCat > 0 && nSer > 0 {
		groupW := pw / float64(nCat)
		stepX := (groupW * 0.72) / math.Max(float64(nSer), 1.0)
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
				cx := gx + groupW*0.14 + (float64(sIdx)+0.5)*stepX
				yv := py(s.Values[cIdx])
				if len(s.Lower) == len(s.Values) && len(s.Upper) == len(s.Values) {
					yLo := py(s.Lower[cIdx])
					yHi := py(s.Upper[cIdx])
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.6"/>`+"\n", cx, yLo, cx, yHi, col)
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx-3.5, yLo, cx+3.5, yLo, col)
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx-3.5, yHi, cx+3.5, yHi, col)
				}
				fmt.Fprintf(&sb, `<circle cx="%.2f" cy="%.2f" r="3.8" fill="%s"/>`+"\n", cx, yv, col)
			}
		}
	}

	legW := 282.0
	legH := float64(len(series))*18.0 + 12.0
	legX := ml + pw - legW - 10.0
	if opts.LegendTopLeft {
		legX = ml + 12.0
	}
	legY := mt + 10.0
	if opts.LegendBottom {
		legY = mt + ph - legH - 10.0
	}
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffffff" fill-opacity="0.93" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legW, legH)
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

// BarSeries merepresentasikan kelompok batang pada grafik batang berkelompok.
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
	w, h := 900.0, 480.0
	ml, mr, mt, mb := 86.0, 32.0, 54.0, 76.0
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
	if opts.YMaxOverride == nil {
		ymax += ypad
	}
	if ymin < 0 && opts.YMinOverride == nil {
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
			col = "#444"
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4" stroke-dasharray="5,4"/>`+"\n", ml, yh, ml+pw, yh, col)
		if hl.Label != "" {
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>`+"\n", ml+8, yh-5, col, esc(hl.Label))
		}
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))

	nCat := len(categories)
	nSer := len(series)
	if nCat > 0 && nSer > 0 {
		groupW := pw / float64(nCat)
		barW := (groupW * 0.74) / float64(nSer)
		offset := groupW * 0.13
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

	legW := 280.0
	legH := float64(len(series))*18.0 + 12.0
	legX := ml + pw - legW - 10.0
	if opts.LegendTopLeft {
		legX = ml + 12.0
	}
	legY := mt + 10.0
	if opts.LegendBottom {
		legY = mt + ph - legH - 10.0
	}
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffffff" fill-opacity="0.93" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legW, legH)
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

// VMarker menggambar penanda garis vertikal pada panel histogram.
type VMarker struct {
	Value  float64
	Label  string
	Color  string
	Dashed bool
}

// HistogramPanel menampung data dan penanda satu panel histogram.
type HistogramPanel struct {
	Subtitle string
	XLabel   string
	Values   []float64
	CILower  float64
	CIUpper  float64
	ShowCI   bool
	Markers  []VMarker
	Bins     int
}

// WriteDualHistogramSVG menggambar dua panel histogram berdampingan dalam satu SVG untuk distribusi bootstrap.
func WriteDualHistogramSVG(path, mainTitle string, leftPanel, rightPanel HistogramPanel) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 960.0, 470.0
	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(mainTitle))

	drawPanel := func(p HistogramPanel, x0, panelW float64) {
		ml, mr, mt, mb := x0+52.0, 18.0, 72.0, 74.0
		pw, ph := panelW-68.0-mr+16.0, h-mt-mb

		xmin, xmax := math.Inf(1), math.Inf(-1)
		for _, v := range p.Values {
			if v < xmin {
				xmin = v
			}
			if v > xmax {
				xmax = v
			}
		}
		for _, mk := range p.Markers {
			if mk.Value < xmin {
				xmin = mk.Value
			}
			if mk.Value > xmax {
				xmax = mk.Value
			}
		}
		if math.IsInf(xmin, 1) || xmax <= xmin {
			xmin, xmax = 0, 1
		}
		span := xmax - xmin
		xmin -= 0.06 * span
		xmax += 0.06 * span

		nBins := p.Bins
		if nBins <= 0 {
			nBins = 14
		}
		counts := make([]int, nBins)
		binW := (xmax - xmin) / float64(nBins)
		maxC := 1
		for _, v := range p.Values {
			b := int((v - xmin) / binW)
			if b >= nBins {
				b = nBins - 1
			}
			if b < 0 {
				b = 0
			}
			counts[b]++
			if counts[b] > maxC {
				maxC = counts[b]
			}
		}
		ymax := float64(maxC) * 1.18

		px := func(v float64) float64 {
			return ml + ((v-xmin)/(xmax-xmin))*pw
		}
		py := func(c float64) float64 {
			return mt + ph - (c/ymax)*ph
		}

		fmt.Fprintf(&sb, `<text x="%.1f" y="54" text-anchor="middle" font-size="13" font-weight="bold" fill="#222">%s</text>`+"\n", ml+pw/2, esc(p.Subtitle))

		if p.ShowCI {
			cx0 := math.Max(ml, px(p.CILower))
			cx1 := math.Min(ml+pw, px(p.CIUpper))
			if cx1 > cx0 {
				fmt.Fprintf(&sb, `<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#1f77b4" fill-opacity="0.11"/>`+"\n", cx0, mt, cx1-cx0, ph)
			}
		}

		fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)
		for k := 0; k <= 4; k++ {
			frac := float64(k) / 4.0
			gx := ml + frac*pw
			xv := xmin + frac*(xmax-xmin)
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" fill="#444">%s</text>`+"\n", gx, mt+ph+17, esc(formatTick(xv, false)))
			gy := mt + ph - frac*ph
			yv := frac * ymax
			fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", ml, gy, ml+pw, gy)
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="10.5" fill="#444">%.0f</text>`+"\n", ml-6, gy+4, yv)
		}

		barPxW := pw / float64(nBins)
		for b := 0; b < nBins; b++ {
			if counts[b] == 0 {
				continue
			}
			bx := ml + float64(b)*barPxW
			by := py(float64(counts[b]))
			bh := mt + ph - by
			fmt.Fprintf(&sb, `<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#1f77b4" fill-opacity="0.68" stroke="#ffffff" stroke-width="0.8"/>`+"\n", bx, by, barPxW, bh)
		}

		for mIdx, mk := range p.Markers {
			mx := px(mk.Value)
			dash := ""
			if mk.Dashed {
				dash = ` stroke-dasharray="5,3"`
			}
			fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="2.0"%s/>`+"\n", mx, mt, mx, mt+ph, mk.Color, dash)
			lx := ml + 8.0
			ly := mt + 14.0 + float64(mIdx)*15.0
			fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.2"%s/>`+"\n", lx, ly-3, lx+16, ly-3, mk.Color, dash)
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="10" fill="#222">%s</text>`+"\n", lx+20, ly, esc(mk.Label))
		}

		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11.5" fill="#222">%s</text>`+"\n", ml+pw/2, h-20, esc(p.XLabel))
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#222" transform="rotate(-90 %.1f %.1f)">Frekuensi Bootstrap</text>`+"\n", x0+14, mt+ph/2, x0+14, mt+ph/2)
	}

	drawPanel(leftPanel, 12.0, (w-24.0)/2.0)
	drawPanel(rightPanel, 12.0+(w-24.0)/2.0, (w-24.0)/2.0)

	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}
