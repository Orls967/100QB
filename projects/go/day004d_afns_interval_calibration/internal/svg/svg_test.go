package svg

import (
	"encoding/xml"
	"os"
	"path/filepath"
	"testing"
)

func TestSVGChartsValidXML(t *testing.T) {
	dir := t.TempDir()
	p1 := filepath.Join(dir, "pt.svg")
	if err := WriteCategoryPointChart(p1, ChartOptions{Title: "Point", XLabel: "X", YLabel: "Y"}, []string{"a", "b"}, []CategoryPointSeries{
		{Label: "S1", Values: []float64{0.5, 0.9}, Lower: []float64{0.4, 0.8}, Upper: []float64{0.6, 0.95}},
	}); err != nil {
		t.Fatal(err)
	}
	p2 := filepath.Join(dir, "hist.svg")
	if err := WriteDualHistogramSVG(p2, "Dual", HistogramPanel{Subtitle: "L", XLabel: "x", Values: []float64{0.1, 0.2, 0.3}}, HistogramPanel{Subtitle: "R", XLabel: "x", Values: []float64{-0.1, 0.0, 0.1}}); err != nil {
		t.Fatal(err)
	}
	p3 := filepath.Join(dir, "scatter.svg")
	if err := WriteScatterChart(p3, ChartOptions{Title: "Scatter", XLabel: "X", YLabel: "Y"}, []VMarker{{Value: 0.05, Label: "True", Color: "#2ca02c"}}, []ScatterSeries{
		{Label: "Cov", X: []float64{0.05, 0.2}, Y: []float64{0.03, 0.02}},
	}); err != nil {
		t.Fatal(err)
	}
	for _, p := range []string{p1, p2, p3} {
		b, err := os.ReadFile(p)
		if err != nil {
			t.Fatal(err)
		}
		var v any
		if err := xml.Unmarshal(b, &v); err != nil {
			t.Fatalf("invalid XML in %s: %v", p, err)
		}
	}
}
