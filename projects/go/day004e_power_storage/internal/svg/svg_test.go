package svg

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSVGWritersProduceValidXML(t *testing.T) {
	dir := t.TempDir()
	p1 := filepath.Join(dir, "line.svg")
	if err := WriteLineChart(p1, ChartOptions{Title: "Line", XLabel: "X", YLabel: "Y", ZeroLine: true},
		[]VBand{{XStart: 1, XEnd: 2, Color: "#fde0dd"}},
		[]LineSeries{{Label: "S1", X: []float64{0, 1, 2, 3}, Y: []float64{-10, 50, 200, 30}, PointMask: []bool{false, false, true, false}}},
	); err != nil {
		t.Fatalf("WriteLineChart error: %v", err)
	}
	b, err := os.ReadFile(p1)
	if err != nil || !strings.Contains(string(b), "<svg") || !strings.Contains(string(b), "</svg>") {
		t.Fatalf("invalid svg output")
	}

	p2 := filepath.Join(dir, "dist.svg")
	if err := WriteOverlayDistributionSVG(p2, "Dist", "Profit", []DistSeries{
		{Label: "A", Values: []float64{10, 20, 30, 40, 50}},
	}, 10); err != nil {
		t.Fatalf("WriteOverlayDistributionSVG error: %v", err)
	}

	p3 := filepath.Join(dir, "cat.svg")
	if err := WriteCategoryPointChart(p3, ChartOptions{Title: "Cat", XLabel: "C", YLabel: "V"}, []string{"M1"}, []CategoryPointSeries{
		{Label: "S", Values: []float64{12.5}, Lower: []float64{10.0}, Upper: []float64{15.0}},
	}); err != nil {
		t.Fatalf("WriteCategoryPointChart error: %v", err)
	}
}
