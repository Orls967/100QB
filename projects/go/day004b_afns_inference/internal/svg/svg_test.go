package svg

import (
	"encoding/xml"
	"os"
	"path/filepath"
	"testing"
)

func TestSVGChartsValidXML(t *testing.T) {
	dir := t.TempDir()
	linePath := filepath.Join(dir, "line.svg")
	if err := WriteLineChart(linePath, ChartOptions{
		Title:     "Test Line",
		XLabel:    "X",
		YLabel:    "Y",
		ZeroLine:  true,
		HRefLines: []HRefLine{{Value: 0.5, Label: "Ref", Color: "#555"}},
	}, []Series{{
		Label: "S1",
		X:     []float64{1, 2, 3},
		Y:     []float64{0.2, 0.5, 0.8},
		Lower: []float64{0.1, 0.4, 0.7},
		Upper: []float64{0.3, 0.6, 0.9},
	}}); err != nil {
		t.Fatalf("WriteLineChart failed: %v", err)
	}

	ptPath := filepath.Join(dir, "points.svg")
	if err := WriteCategoryPointChart(ptPath, ChartOptions{
		Title:     "Test Points",
		XLabel:    "Cat",
		YLabel:    "Val",
		HRefLines: []HRefLine{{Value: 0.95, Label: "95%"}},
	}, []string{"A", "B"}, []CategoryPointSeries{{
		Label:  "Cov",
		Values: []float64{0.92, 0.80},
		Lower:  []float64{0.81, 0.67},
		Upper:  []float64{0.97, 0.89},
	}}); err != nil {
		t.Fatalf("WriteCategoryPointChart failed: %v", err)
	}

	barPath := filepath.Join(dir, "bars.svg")
	if err := WriteGroupedBarChart(barPath, ChartOptions{
		Title:  "Test Bars",
		XLabel: "Cat",
		YLabel: "Val",
	}, []string{"A", "B"}, []BarSeries{{
		Label:  "B1",
		Values: []float64{1.2, 0.4},
	}}); err != nil {
		t.Fatalf("WriteGroupedBarChart failed: %v", err)
	}

	for _, p := range []string{linePath, ptPath, barPath} {
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
