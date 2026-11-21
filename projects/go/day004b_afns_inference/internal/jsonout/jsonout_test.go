package jsonout

import (
	"encoding/json"
	"math"
	"path/filepath"
	"strings"
	"testing"
)

func TestMarshalPrecisionAndValidation(t *testing.T) {
	obj := Obj{
		{Key: "pi", Val: math.Pi},
		{Key: "flags", Val: []bool{true, false}},
		{Key: "nested", Val: []Obj{{{Key: "k", Val: 42}}}},
	}
	b, err := Marshal(obj)
	if err != nil {
		t.Fatalf("Marshal failed: %v", err)
	}
	if !strings.Contains(string(b), "3.141592653589793") {
		t.Fatalf("expected 17-digit float formatting, got:\n%s", string(b))
	}
	var parsed map[string]any
	if err := json.Unmarshal(b, &parsed); err != nil {
		t.Fatalf("invalid JSON produced: %v", err)
	}

	if _, err := Marshal(Obj{{Key: "bad", Val: math.NaN()}}); err == nil {
		t.Fatalf("expected error on NaN")
	}

	outPath := filepath.Join(t.TempDir(), "test.json")
	if err := WriteFile(outPath, obj); err != nil {
		t.Fatalf("WriteFile failed: %v", err)
	}
}
