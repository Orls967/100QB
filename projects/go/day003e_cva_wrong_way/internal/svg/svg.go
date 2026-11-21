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

// Series merepresentasikan satu seri garis pada grafik SVG, dengan pita galat opsional (Upper/Lower).
type Series struct {
	Label  string
	X      []float64
	Y      []float64
	Lower  []float64
	Upper  []float64
	Color  string
	Dashed bool
}

// ChartOptions mengonfigurasi judul, label sumbu, skala log, dan posisi legenda.
type ChartOptions struct {
	Title         string
	XLabel        string
	YLabel        string
	LogX          bool
	LogY          bool
	ZeroLine      bool
	LegendTopLeft bool
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
	return fmt.Sprintf("%.4g", v)
}

// WriteLineChart menggambar grafik garis multi-seri ke berkas SVG, mendukung skala log dan galat bar.
func WriteLineChart(path string, opts ChartOptions, series []Series) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 820.0, 460.0
	ml, mr, mt, mb := 86.0, 30.0, 52.0, 62.0
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
				lv := tx(s.Lower[i], opts.LogY)
				if lv < ymin {
					ymin = lv
				}
			}
			if len(s.Upper) == len(s.Y) {
				uv := tx(s.Upper[i], opts.LogY)
				if uv > ymax {
					ymax = uv
				}
			}
		}
	}
	if !opts.LogY && ymin > 0 && ymax > 0 {
		ymin = 0
	}
	if math.IsInf(xmin, 1) || xmax <= xmin {
		xmin, xmax = 0, 1
	}
	if math.IsInf(ymin, 1) || ymax <= ymin {
		ymin = ymin - 1
		ymax = ymax + 1
	}
	ypad := 0.08 * (ymax - ymin)
	ymax += ypad
	if ymin != 0 || opts.LogY {
		ymin -= 0.04 * (ymax - ymin)
	}
	xpad := 0.03 * (xmax - xmin)
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
	fmt.Fprintf(&sb, `<text x="%.1f" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(opts.Title))
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
		y0 := mt + ph - ((0-ymin)/(ymax-ymin))*ph
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999" stroke-width="1" stroke-dasharray="3,3"/>`+"\n", ml, y0, ml+pw, y0)
	}

	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-14, esc(opts.XLabel))
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
				fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx, yLo, cx, yHi, col)
				fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx-4, yLo, cx+4, yLo, col)
				fmt.Fprintf(&sb, `<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1.5"/>`+"\n", cx-4, yHi, cx+4, yHi, col)
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
		if len(s.X) <= 12 {
			for i := range s.X {
				fmt.Fprintf(&sb, `<circle cx="%.2f" cy="%.2f" r="3.4" fill="%s"/>`+"\n", px(s.X[i]), py(s.Y[i]), col)
			}
		}
	}

	// Kotak legenda semi-transparan di kanan atas atau kiri atas area grafik.
	legX := ml + pw - 255.0
	if opts.LegendTopLeft {
		legX = ml + 12.0
	}
	legY := mt + 10.0
	legH := float64(len(series))*18.0 + 12.0
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="245" height="%.1f" fill="#ffffff" fill-opacity="0.9" stroke="#ccc" rx="4"/>`+"\n", legX, legY, legH)
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
