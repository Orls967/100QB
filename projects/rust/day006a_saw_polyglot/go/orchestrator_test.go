package main

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
	"time"
)

func ensureRustBinary(t *testing.T) string {
	t.Helper()
	if envBin := os.Getenv("SAW_BIN"); envBin != "" {
		if _, err := os.Stat(envBin); err == nil {
			return envBin
		}
	}
	cand := filepath.Clean("../rust/target/release/saw")
	if _, err := os.Stat(cand); err == nil {
		abs, _ := filepath.Abs(cand)
		return abs
	}
	cmd := exec.Command("cargo", "build", "--release", "--manifest-path", "../rust/Cargo.toml")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("cargo build failed: %v\n%s", err, string(out))
	}
	abs, err := filepath.Abs(cand)
	if err != nil {
		t.Fatalf("abs path error: %v", err)
	}
	return abs
}

func TestPrefixGenerationCountK8(t *testing.T) {
	shortCounts, tasks := GeneratePrefixesAndShortCounts(12, 8)
	// c_8 = 5916 => with first step forced East (weight 4), there are 5916 / 4 = 1479 prefixes.
	if len(tasks) != 1479 {
		t.Fatalf("expected 1479 prefixes at K=8, got %d", len(tasks))
	}
	if shortCounts[1] != 4 || shortCounts[4] != 100 || shortCounts[7] != 2172 {
		t.Fatalf("unexpected shortCounts: %v", shortCounts[:8])
	}
}

func TestShardedMatchesSingleProcessSmallN(t *testing.T) {
	sawBin := ensureRustBinary(t)
	ctx, cancel := context.WithTimeout(context.Background(), 15*time.Second)
	defer cancel()

	single, _, err := RunRustSingleCount(ctx, sawBin, 14, 0)
	if err != nil {
		t.Fatalf("single count failed: %v", err)
	}
	sharded, nTasks, err := RunShardedEnumeration(ctx, sawBin, 14, 6, 3)
	if err != nil {
		t.Fatalf("sharded count failed: %v", err)
	}
	if nTasks != 780/4 {
		t.Fatalf("expected 195 tasks at K=6, got %d", nTasks)
	}
	if !SlicesEqualU64(single.Counts, sharded) {
		t.Fatalf("mismatch:\nsingle=%v\nsharded=%v", single.Counts, sharded)
	}
}

func TestWorkerCountInvariance(t *testing.T) {
	sawBin := ensureRustBinary(t)
	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Second)
	defer cancel()

	var base []uint64
	for _, w := range []int{1, 2, 4, 8} {
		got, _, err := RunShardedEnumeration(ctx, sawBin, 15, 7, w)
		if err != nil {
			t.Fatalf("workers=%d failed: %v", w, err)
		}
		if base == nil {
			base = got
		} else if !SlicesEqualU64(base, got) {
			t.Fatalf("worker invariance violated for w=%d: %v vs %v", w, got, base)
		}
	}
}

func TestDeterministicTaskOrderAndSeeds(t *testing.T) {
	_, t1 := GeneratePrefixesAndShortCounts(10, 5)
	_, t2 := GeneratePrefixesAndShortCounts(10, 5)
	if len(t1) != len(t2) {
		t.Fatalf("task count mismatch")
	}
	for i := range t1 {
		if t1[i] != t2[i] {
			t.Fatalf("non-deterministic prefix order at %d: %+v vs %+v", i, t1[i], t2[i])
		}
	}
	if DeriveSeed(20261124, 400) != DeriveSeed(20261124, 400) {
		t.Fatalf("DeriveSeed must be deterministic")
	}
	if DeriveSeed(20261124, 400) == DeriveSeed(20261124, 800) {
		t.Fatalf("DeriveSeed collision for distinct N")
	}
}

func TestContextTimeoutCancelsRun(t *testing.T) {
	sawBin := ensureRustBinary(t)
	ctx, cancel := context.WithTimeout(context.Background(), 1*time.Nanosecond)
	defer cancel()
	time.Sleep(1 * time.Millisecond)

	_, _, err := RunShardedEnumeration(ctx, sawBin, 20, 8, 2)
	if err == nil {
		t.Fatalf("expected error on expired context, got nil")
	}
}

func TestParseWorkerResponseJSON(t *testing.T) {
	okLine := []byte(`{"id":7,"ok":true,"counts":[1,4,12,36]}`)
	resp, err := ParseWorkerResponseLine(okLine)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if resp.ID != 7 || len(resp.Counts) != 4 || resp.Counts[3] != 36 {
		t.Fatalf("unexpected parsed response: %+v", resp)
	}

	errLine := []byte(`{"id":0,"ok":false,"error":"prefix is not self-avoiding"}`)
	if _, err := ParseWorkerResponseLine(errLine); err == nil {
		t.Fatalf("expected error for ok=false response")
	}
}
