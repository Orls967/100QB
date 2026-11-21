package svg

import (
	"encoding/xml"
	"os"
	"path/filepath"
	"testing"
)

func TestSVGChartsValidXML(t *testing.T) {
	dir := t.TempDir()
	ptPath := filepath.Join(dir, "points.svg")
	if err := WriteCategoryPointChart(ptPath, ChartOptions{
		Title:     "Bias Point Chart",
		XLabel:    "Parameter",
		YLabel:    "Bias",
		ZeroLine:  true,
		HRefLines: []HRefLine{{Value: 0.95, Label: "95%"}},
	}, []string{"kp11", "theta1"}, []CategoryPointSeries{{
		Label:  "Raw MLE",
		Values: []float64{0.22, -0.001},
		Lower:  []float64{0.19, -0.005},
		Upper:  []float64{0.25, 0.003},
	}}); err != nil {
		t.Fatalf("WriteCategoryPointChart error: %v", err)
	}

	barPath := filepath.Join(dir, "bars.svg")
	if err := WriteGroupedBarChart(barPath, ChartOptions{
		Title:  "Grouped Bars",
		XLabel: "Config",
		YLabel: "RMSE",
	}, []string{"20y", "80y"}, []BarSeries{{
		Label:  "MLE",
		Values: []float64{0.30, 0.07},
	}}); err != nil {
		t.Fatalf("WriteGroupedBarChart error: %v", err)
	}

	histPath := filepath.Join(dir, "hist.svg")
	pL := HistogramPanel{
		Subtitle: "kp11 draws",
		XLabel:   "kp11",
		Values:   []float64{0.1, 0.15, 0.2, 0.22, 0.25, 0.3, 0.35},
		ShowCI:   true,
		CILower:  0.11,
		CIUpper:  0.34,
		Markers:  []VMarker{{Value: 0.05, Label: "True", Color: "#2ca02c", Dashed: true}},
		Bins:     6,
	}
	pR := pL
	pR.Subtitle = "theta1 draws"
	if err := WriteDualHistogramSVG(histPath, "Bootstrap Draws", pL, pR); err != nil {
		t.Fatalf("WriteDualHistogramSVG error: %v", err)
	}

	for _, p := range []string{ptPath, barPath, histPath} {
		data, err := os.ReadFile(p)
		if err != nil {
			t.Fatalf("read %s: %v", p, err)
		}
		var doc any
		if err := xml.Unmarshal(data, &doc); err != nil {
			t.Fatalf("invalid XML in %s: %v", p, err)
		}
	}
}
