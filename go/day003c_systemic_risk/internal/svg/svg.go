package svg

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// Palette warna standar 100QB.
var Palette = []string{"#1f6feb", "#d1495b", "#2a9d8f", "#e9c46a", "#6a4c93", "#f4a261"}

// Series merepresentasikan satu kurva (X, Y) pada grafik garis.
type Series struct {
	Label  string
	X      []float64
	Y      []float64
	Color  string
	Dashed bool
}

func esc(s string) string {
	r := strings.NewReplacer("&", "&amp;", "<", "&lt;", ">", "&gt;", `"`, "&quot;")
	return r.Replace(s)
}

// WriteLineChart menggambar grafik garis 2D (linear atau log-log) ke berkas SVG.
func WriteLineChart(path, title, xLabel, yLabel string, series []Series, logX, logY bool) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 820.0, 460.0
	ml, mr, mt, mb := 84.0, 30.0, 48.0, 64.0
	pw, ph := w-ml-mr, h-mt-mb

	tx := func(v float64) float64 {
		if logX {
			if v <= 1e-30 {
				v = 1e-30
			}
			return math.Log10(v)
		}
		return v
	}
	ty := func(v float64) float64 {
		if logY {
			if v <= 1e-30 {
				v = 1e-30
			}
			return math.Log10(v)
		}
		return v
	}

	xmin, xmax := math.Inf(1), math.Inf(-1)
	ymin, ymax := math.Inf(1), math.Inf(-1)
	for _, s := range series {
		for i := range s.X {
			vx, vy := tx(s.X[i]), ty(s.Y[i])
			if !math.IsNaN(vx) && !math.IsInf(vx, 0) && !math.IsNaN(vy) && !math.IsInf(vy, 0) {
				if vx < xmin {
					xmin = vx
				}
				if vx > xmax {
					xmax = vx
				}
				if vy < ymin {
					ymin = vy
				}
				if vy > ymax {
					ymax = vy
				}
			}
		}
	}
	if math.Abs(xmax-xmin) < 1e-15 {
		xmin -= 1
		xmax += 1
	}
	if math.Abs(ymax-ymin) < 1e-15 {
		ymin -= 1
		ymax += 1
	}
	pad := 0.06 * (ymax - ymin)
	ymin -= pad
	ymax += pad

	px := func(v float64) float64 { return ml + (tx(v)-xmin)/(xmax-xmin)*pw }
	py := func(v float64) float64 { return mt + (1.0-(ty(v)-ymin)/(ymax-ymin))*ph }

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(title))

	for k := 0; k <= 4; k++ {
		fr := float64(k) / 4.0
		gx := ml + fr*pw
		gy := mt + (1.0-fr)*ph
		xv := xmin + fr*(xmax-xmin)
		yv := ymin + fr*(ymax-ymin)
		if logX {
			xv = math.Pow(10, xv)
		}
		if logY {
			yv = math.Pow(10, yv)
		}
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e5e7eb" stroke-dasharray="3,3"/>`+"\n", gx, mt, gx, mt+ph)
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e5e7eb" stroke-dasharray="3,3"/>`+"\n", ml, gy, ml+pw, gy)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#444">%.3g</text>`+"\n", gx, mt+ph+18, xv)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" fill="#444">%.3g</text>`+"\n", ml-8, gy+4, yv)
	}

	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)
	fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" fill="#222">%s</text>`+"\n", ml+pw/2, h-14, esc(xLabel))
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(yLabel))

	for _, s := range series {
		var pts strings.Builder
		for i := range s.X {
			fmt.Fprintf(&pts, "%.2f,%.2f ", px(s.X[i]), py(s.Y[i]))
		}
		dash := ""
		if s.Dashed {
			dash = ` stroke-dasharray="6,4"`
		}
		fmt.Fprintf(&sb, `<polyline fill="none" stroke="%s" stroke-width="2.2"%s points="%s"/>`+"\n", s.Color, dash, pts.String())
	}

	for i, s := range series {
		lx := ml + 14.0
		ly := mt + 18.0 + float64(i)*18.0
		fmt.Fprintf(&sb, `<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5"/>`+"\n", lx, ly, lx+20, ly, s.Color)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" font-size="11" fill="#222">%s</text>`+"\n", lx+26, ly+4, esc(s.Label))
	}

	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}

// WriteBarChart menggambar grafik batang terkelompok/tunggal ke berkas SVG.
func WriteBarChart(path, title, yLabel string, labels []string, values []float64) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	w, h := 780.0, 430.0
	ml, mr, mt, mb := 82.0, 30.0, 48.0, 72.0
	pw, ph := w-ml-mr, h-mt-mb

	ymax := 1e-12
	for _, v := range values {
		if v > ymax {
			ymax = v
		}
	}
	ymax *= 1.15
	n := float64(len(labels))
	if n < 1 {
		n = 1
	}
	slot := pw / n
	bw := slot * 0.56

	var sb strings.Builder
	fmt.Fprintf(&sb, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f" style="background:#ffffff;font-family:system-ui,sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&sb, `<text x="%.1f" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">%s</text>`+"\n", w/2, esc(title))
	fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="#333" stroke-width="1"/>`+"\n", ml, mt, pw, ph)

	for i, lab := range labels {
		val := values[i]
		if val < 0 {
			val = 0
		}
		bx := ml + (float64(i)+0.5)*slot - 0.5*bw
		bh := (val / ymax) * ph
		by := mt + ph - bh
		col := Palette[i%len(Palette)]
		fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>`+"\n", bx, by, bw, bh, col)
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#222">%.4g</text>`+"\n", bx+0.5*bw, by-5, values[i])
		fmt.Fprintf(&sb, `<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#333">%s</text>`+"\n", bx+0.5*bw, mt+ph+18, esc(lab))
	}
	fmt.Fprintf(&sb, `<text x="18" y="%.1f" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 %.1f)">%s</text>`+"\n", mt+ph/2, mt+ph/2, esc(yLabel))
	sb.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}
