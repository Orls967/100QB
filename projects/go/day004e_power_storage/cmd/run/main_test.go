package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func TestQuickRunnerEndToEndAndSVGs(t *testing.T) {
	tmp := t.TempDir()
	outJSON := filepath.Join(tmp, "results.json")
	figDir := filepath.Join(tmp, "figures")

	cmd := exec.Command("go", "run", ".", "-quick", "-workers", "2", "-out", outJSON, "-fig-dir", figDir)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("quick runner failed: %v\n%s", err, string(out))
	}
	if !strings.Contains(string(out), "all_technical_passed=true") {
		t.Fatalf("expected all_technical_passed=true in output:\n%s", string(out))
	}
	jb, err := os.ReadFile(outJSON)
	if err != nil || len(jb) < 500 {
		t.Fatalf("missing or empty results.json: %v", err)
	}

	expectedSVGs := []string{
		"day004e_sample-price-regimes-jumps.svg",
		"day004e_example-day-dispatch-soc.svg",
		"day004e_daily-profit-distribution.svg",
		"day004e_dp-value-vs-soc-regime.svg",
		"day004e_grid-resolution-sensitivity.svg",
	}
	for _, name := range expectedSVGs {
		b, err := os.ReadFile(filepath.Join(figDir, name))
		if err != nil || !strings.Contains(string(b), "<svg") || !strings.Contains(string(b), "</svg>") {
			t.Fatalf("invalid or missing SVG %s: %v", name, err)
		}
	}
}
