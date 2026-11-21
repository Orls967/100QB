package svg

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// Palette warna ramah buta-warna.
var Palette = []string{
	"#1f77b4",
	"#d62728",
	"#2ca02c",
	"#9467bd",
	"#ff7f0e",
	"#17becf",
}

// HRefLine menggambar garis referensi horizontal.
type HRefLine struct {
	Value float64
	Label string
	Color string
}

// VBand memberi arsiran latar belakang vertikal pada rentang indeks/nilai [XStart, XEnd].
type VBand struct {
	XStart float64
	XEnd   float64
	Color  string
}

// LineSeries menyimpan satu seri garis X-Y dengan penanda opsional.
type LineSeries struct {
	Label      string
	X          []float64
	Y          []float64
	Color      string
	StrokeW    float64
	Dashed     bool
	PointMask  []bool // bila non-nil, menggambar titik di indeks true
	PointColor string
}

// ChartOptions mengonfigurasi judul, label, batas sumbu, dan legenda grafik.
type ChartOptions struct {
	Title         string
	XLabel        string
	YLabel        string
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

func formatTick(v float64) string {
	if math.Abs(v) < 1e-12 {
		return "0"
	}
	return fmt.Sprintf("%.4g", v)
}

// WriteLineChart menggambar grafik garis multi-seri dengan pita latar belakang dan penanda titik opsional.
func WriteLineChart(path string, opts ChartOptions, bands []VBand, series []LineSeries) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 940.0, 480.0
	ml, mr, mt, mb := 84.0, 30.0, 56.0, 74.0
	pw, ph := w-ml-mr, h-mt-mb

	xmin, xmax := math.Inf(1), math.Inf(-1)
	ymin, ymax := math.Inf(1), math.Inf(-1)
	for _, s := range series {
		for i := range s.X {
			if i >= len(s.Y) {
				break
			}
			xv, yv := s.X[i], s.Y[i]
			if xv < xmin {
				xmin = xv
			}
			if xv > xmax {
				xmax = xv
			}
			if yv < ymin {
				ymin = yv
			}
			if yv > ymax {
				ymax = yv
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
	if math.IsInf(xmin, 1) || xmax <= xmin {
		xmin, xmax = 0, 1
	}
	if math.IsInf(ymin, 1) || ymax <= ymin {
		ymin, ymax = -1, 1
	}
	ypad := 0.08 * (ymax - ymin)
	if opts.YMinOverride == nil {
		ymin -= ypad
	}
	if opts.YMaxOverride == nil {
		ymax += ypad
	}

	px := func(v float64) float64 {
		return ml + ((v-xmin)/(xmax-xmin))*pw
	}
	py := func(v float64) float64 {
		return mt + ph - ((v-ymin)/(ymax-ymin))*ph
	}

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="30" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(opts.Title))

	for _, b := range bands {
		x1 := math.Max(ml, px(b.XStart))
		x2 := math.Min(ml+pw, px(b.XEnd))
		if x2 > x1 {
			fmt.Fprintf(&sb, `<rect x="%.2f" y="%.1f" width="%.2f" height="%.1f" fill="%s" fill-opacity="0.45"/>`+"\n", x1, mt, x2-x1, ph, b.Color)
		}
	}

	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)
	for k := 0; k <= 5; k++ {
		frac := float64(k) / 5.0
		gx := ml + frac*pw
		xv := xmin + frac*(xmax-xmin)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", gx, mt, gx, mt+ph)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#444">%s</text>`+"\n", gx, mt+ph+18, esc(formatTick(xv)))

		gy := mt + ph - frac*ph
		yv := ymin + frac*(ymax-ymin)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", ml, gy, ml+pw, gy)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" fill="#444">%s</text>`+"\n", ml-8, gy+4, esc(formatTick(yv)))
	}

	if opts.ZeroLine && ymin < 0 && ymax > 0 {
		y0 := py(0)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#666" stroke-width="1.2" stroke-dasharray="4,3"/>`+"\n", ml, y0, ml+pw, y0)
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

	for sIdx, s := range series {
		col := s.Color
		if col == "" {
			col = Palette[sIdx%len(Palette)]
		}
		sw := s.StrokeW
		if sw <= 0 {
			sw = 2.0
		}
		dash := ""
		if s.Dashed {
			dash = ` stroke-dasharray="6,4"`
		}
		var pts strings.Builder
		for i := range s.X {
			if i >= len(s.Y) {
				break
			}
			if i > 0 {
				pts.WriteByte(' ')
			}
			fmt.Fprintf(&pts, "%.2f,%.2f", px(s.X[i]), py(s.Y[i]))
		}
		fmt.Fprintf(&sb, `<polyline fill="none" stroke="%s" stroke-width="%.2f"%s points="%s"/>`+"\n", col, sw, dash, pts.String())
		if len(s.PointMask) == len(s.X) {
			pCol := s.PointColor
			if pCol == "" {
				pCol = "#d62728"
			}
			for i, ok := range s.PointMask {
				if ok && i < len(s.Y) {
					fmt.Fprintf(&sb, `<circle cx="%.2f" cy="%.2f" r="4.2" fill="%s" stroke="#fff" stroke-width="0.9"/>`+"\n", px(s.X[i]), py(s.Y[i]), pCol)
				}
			}
		}
	}

	legW := 285.0
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
			dash = ` stroke-dasharray="4,2"`
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.2"%s/>`+"\n", lx, ly-3, lx+18, ly-3, col, dash)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="11" fill="#222">%s</text>`+"\n", lx+24, ly+1, esc(s.Label))
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))
	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}

// WriteDualDayDispatchSVG menggambar panel atas harga & aksi daya dan panel bawah lintasan SoC 24 jam.
func WriteDualDayDispatchSVG(path, title string, hours []float64, prices []float64, topSeries []LineSeries, socSeries []LineSeries) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 940.0, 620.0
	ml, mr := 82.0, 30.0
	pw := w - ml - mr

	drawLinePanel := func(sb *strings.Builder, yTop, ph float64, yLabel string, zeroLine bool, yMinOvr, yMaxOvr *float64, sers []LineSeries, drawXAxis bool) {
		xmin, xmax := 0.0, 24.0
		ymin, ymax := math.Inf(1), math.Inf(-1)
		for _, s := range sers {
			for i := range s.X {
				if i >= len(s.Y) {
					break
				}
				if s.Y[i] < ymin {
					ymin = s.Y[i]
				}
				if s.Y[i] > ymax {
					ymax = s.Y[i]
				}
			}
		}
		if zeroLine {
			if 0 < ymin {
				ymin = 0
			}
			if 0 > ymax {
				ymax = 0
			}
		}
		if yMinOvr != nil {
			ymin = *yMinOvr
		}
		if yMaxOvr != nil {
			ymax = *yMaxOvr
		}
		if math.IsInf(ymin, 1) || ymax <= ymin {
			ymin, ymax = -1, 1
		}
		ypad := 0.08 * (ymax - ymin)
		if yMinOvr == nil {
			ymin -= ypad
		}
		if yMaxOvr == nil {
			ymax += ypad
		}
		px := func(v float64) float64 { return ml + ((v-xmin)/(xmax-xmin))*pw }
		py := func(v float64) float64 { return yTop + ph - ((v-ymin)/(ymax-ymin))*ph }

		fmt.Fprintf(sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, yTop, pw, ph)
		for k := 0; k <= 4; k++ {
			frac := float64(k) / 4.0
			gy := yTop + ph - frac*ph
			yv := ymin + frac*(ymax-ymin)
			fmt.Fprintf(sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee" stroke-width="1"/>`+"\n", ml, gy, ml+pw, gy)
			fmt.Fprintf(sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="10.5" fill="#444">%s</text>`+"\n", ml-7, gy+4, esc(formatTick(yv)))
		}
		for hr := 0; hr <= 24; hr += 4 {
			gx := px(float64(hr))
			fmt.Fprintf(sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#f0f0f0" stroke-width="1"/>`+"\n", gx, yTop, gx, yTop+ph)
			if drawXAxis {
				fmt.Fprintf(sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" fill="#444">%d</text>`+"\n", gx, yTop+ph+16, hr)
			}
		}
		if zeroLine && ymin < 0 && ymax > 0 {
			y0 := py(0)
			fmt.Fprintf(sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#777" stroke-width="1.1" stroke-dasharray="4,3"/>`+"\n", ml, y0, ml+pw, y0)
		}
		for sIdx, s := range sers {
			col := s.Color
			if col == "" {
				col = Palette[sIdx%len(Palette)]
			}
			sw := s.StrokeW
			if sw <= 0 {
				sw = 2.0
			}
			dash := ""
			if s.Dashed {
				dash = ` stroke-dasharray="5,3"`
			}
			var pts strings.Builder
			for i := range s.X {
				if i >= len(s.Y) {
					break
				}
				if i > 0 {
					pts.WriteByte(' ')
				}
				fmt.Fprintf(&pts, "%.2f,%.2f", px(s.X[i]), py(s.Y[i]))
			}
			fmt.Fprintf(sb, `<polyline fill="none" stroke="%s" stroke-width="%.2f"%s points="%s"/>`+"\n", col, sw, dash, pts.String())
		}
		legW := 240.0
		legH := float64(len(sers))*17.0 + 10.0
		legX := ml + 10.0
		legY := yTop + 8.0
		fmt.Fprintf(sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffffff" fill-opacity="0.90" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legW, legH)
		for i, s := range sers {
			col := s.Color
			if col == "" {
				col = Palette[i%len(Palette)]
			}
			lx := legX + 8.0
			ly := legY + 14.0 + float64(i)*17.0
			fmt.Fprintf(sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.2"/>`+"\n", lx, ly-3, lx+16, ly-3, col)
			fmt.Fprintf(sb, `<text x="%.1f" y="%.1f" font-size="10.5" fill="#222">%s</text>`+"\n", lx+21, ly+1, esc(s.Label))
		}
		fmt.Fprintf(sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="11.5" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", yTop+ph/2, yTop+ph/2, esc(yLabel))
	}

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(title))

	// Panel 1: Harga per jam (kiri)
	priceSer := []LineSeries{{Label: "Harga Spot P_h (EUR/MWh)", X: hours, Y: prices, Color: "#111111", StrokeW: 2.2}}
	drawLinePanel(&sb, 48.0, 150.0, "Harga (EUR/MWh)", true, nil, nil, priceSer, false)

	// Panel 2: Daya Dispatch MW (+ charge, - discharge)
	pMin, pMax := -2.0, 2.2
	drawLinePanel(&sb, 220.0, 150.0, "Daya Grid (MW)", true, &pMin, &pMax, topSeries, false)

	// Panel 3: State of Charge (0..4 MWh)
	sMin, sMax := 0.0, 4.1
	drawLinePanel(&sb, 392.0, 165.0, "SoC Baterai (MWh)", false, &sMin, &sMax, socSeries, true)

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">Jam dalam Hari Evaluasi (h = 0..24)</text>`+"\n", ml+pw/2, h-16)
	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}

// DistSeries merepresentasikan sampel profit harian satu strategi untuk plot distribusi.
type DistSeries struct {
	Label  string
	Values []float64
	Color  string
	Dashed bool
}

// WriteOverlayDistributionSVG menggambar poligon frekuensi distribusi profit harian antar-strategi.
func WriteOverlayDistributionSVG(path, title, xLabel string, series []DistSeries, nBins int) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	var all []float64
	for _, s := range series {
		all = append(all, s.Values...)
	}
	sort.Float64s(all)
	if len(all) == 0 {
		all = []float64{0, 1}
	}
	// Potong tampilan sumbu-X ke [P1, P99] agar puncak distribusi terbaca jelas sambil mencakup massa utama.
	xmin := all[int(0.01*float64(len(all)-1))]
	xmax := all[int(0.99*float64(len(all)-1))]
	if xmax <= xmin {
		xmin, xmax = all[0]-1, all[len(all)-1]+1
	}
	if nBins <= 0 {
		nBins = 32
	}
	binW := (xmax - xmin) / float64(nBins)

	xs := make([]float64, nBins)
	for b := 0; b < nBins; b++ {
		xs[b] = xmin + (float64(b)+0.5)*binW
	}
	lineSeries := make([]LineSeries, len(series))
	for i, s := range series {
		ys := make([]float64, nBins)
		invN := 1.0 / math.Max(float64(len(s.Values)), 1.0)
		for _, v := range s.Values {
			b := int((v - xmin) / binW)
			if b < 0 {
				b = 0
			}
			if b >= nBins {
				b = nBins - 1
			}
			ys[b] += 100.0 * invN // persen hari per bin
		}
		lineSeries[i] = LineSeries{
			Label:   s.Label,
			X:       xs,
			Y:       ys,
			Color:   s.Color,
			StrokeW: 2.2,
			Dashed:  s.Dashed,
		}
	}
	y0 := 0.0
	return WriteLineChart(path, ChartOptions{
		Title:        title,
		XLabel:       xLabel,
		YLabel:       "Persentase Hari per Bin (%)",
		ZeroLine:     false,
		YMinOverride: &y0,
	}, nil, lineSeries)
}

// CategoryPointSeries merepresentasikan seri estimasi titik beserta selang ketidakpastian.
type CategoryPointSeries struct {
	Label  string
	Values []float64
	Lower  []float64
	Upper  []float64
	Color  string
}

// WriteCategoryPointChart menggambar grafik titik berkelompok dengan error bar.
func WriteCategoryPointChart(path string, opts ChartOptions, categories []string, series []CategoryPointSeries) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 900.0, 480.0
	ml, mr, mt, mb := 84.0, 32.0, 56.0, 78.0
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
	if opts.ZeroLine {
		if 0 < ymin {
			ymin = 0
		}
		if 0 > ymax {
			ymax = 0
		}
	}
	if math.IsInf(ymin, 1) || ymax <= ymin {
		ymin, ymax = -0.1, 0.1
	}
	ypad := 0.10 * (ymax - ymin)
	ymin -= ypad
	ymax += ypad

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
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" fill="#444">%s</text>`+"\n", ml-8, gy+4, esc(formatTick(yv)))
	}
	if opts.ZeroLine && ymin < 0 && ymax > 0 {
		y0 := py(0)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#777" stroke-width="1.2" stroke-dasharray="3,3"/>`+"\n", ml, y0, ml+pw, y0)
	}

	nCat := len(categories)
	nSer := len(series)
	if nCat > 0 && nSer > 0 {
		groupW := pw / float64(nCat)
		stepX := (groupW * 0.70) / float64(nSer)
		for cIdx, cat := range categories {
			gx := ml + float64(cIdx)*groupW
			if cIdx > 0 {
				fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#f0f0f0" stroke-width="1"/>`+"\n", gx, mt, gx, mt+ph)
			}
			fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#333">%s</text>`+"\n", gx+groupW/2, mt+ph+20, esc(cat))
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
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.8"/>`+"\n", cx, yLo, cx, yHi, col)
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.6"/>`+"\n", cx-4, yLo, cx+4, yLo, col)
					fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.6"/>`+"\n", cx-4, yHi, cx+4, yHi, col)
				}
				fmt.Fprintf(&sb, `<circle cx="%.2f" cy="%.2f" r="4.2" fill="%s"/>`+"\n", cx, yv, col)
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

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-16, esc(opts.XLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(opts.YLabel))
	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}
