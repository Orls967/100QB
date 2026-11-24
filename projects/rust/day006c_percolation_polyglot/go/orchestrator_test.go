package main

import (
	"context"
	"encoding/json"
	"os"
	"path/filepath"
	"testing"
	"time"
)

func TestFixedChunksIndependentOfWorkers(t *testing.T) {
	j1 := BuildJobs(20261124, []int{16, 32}, map[int]int{16: 1000, 32: 1000}, 10)
	j2 := BuildJobs(20261124, []int{16, 32}, map[int]int{16: 1000, 32: 1000}, 10)
	for i := range j1 {
		if j1[i] != j2[i] {
			t.Fatalf("mismatch at %d", i)
		}
	}
}

func TestDeterministicReductionOrder(t *testing.T) {
	bin := filepath.Join("..", "rust", "target", "release", "percolate")
	if _, err := os.Stat(bin); err != nil {
		t.Skip("no binary")
	}
	ctx, cancel := context.WithTimeout(context.Background(), 15*time.Second)
	defer cancel()
	jobs := BuildJobs(20261124, []int{8, 16}, map[int]int{8: 800, 16: 800}, 8)
	r1, e1 := RunPool(ctx, bin, jobs, 1)
	r4, e4 := RunPool(ctx, bin, jobs, 4)
	if e1 != nil || e4 != nil {
		t.Fatalf("%v / %v", e1, e4)
	}
	b1, _ := json.Marshal(r1)
	b4, _ := json.Marshal(r4)
	if string(b1) != string(b4) {
		t.Fatal("1 != 4 workers")
	}
}

func TestContextTimeout(t *testing.T) {
	slow := filepath.Join(t.TempDir(), "slow.sh")
	_ = os.WriteFile(slow, []byte("#!/bin/sh\nsleep 5\n"), 0o755)
	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Millisecond)
	defer cancel()
	if _, err := RunPool(ctx, slow, []Job{{Idx: 0, L: 4, Seed: 1, Samples: 10}}, 2); err == nil {
		t.Fatal("expected timeout")
	}
}
