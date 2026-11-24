package main

import (
	"context"
	"os"
	"path/filepath"
	"reflect"
	"testing"
	"time"
)

func sawBinPath(t *testing.T) string {
	t.Helper()
	if p := os.Getenv("SAW_BIN"); p != "" {
		return p
	}
	abs, _ := filepath.Abs("../rust/target/release/saw")
	return abs
}

func TestShardedEqualsSingleAndWorkerInvariance(t *testing.T) {
	ctx, cancel := context.WithTimeout(context.Background(), 15*time.Second)
	defer cancel()
	bin := sawBinPath(t)
	c1, nTasks, err := RunSharded(ctx, bin, 12, 7, 1)
	if err != nil {
		t.Fatal(err)
	}
	if nTasks != 2172/4 {
		t.Fatalf("want 543 tasks for K=7, got %d", nTasks)
	}
	c4, _, err := RunSharded(ctx, bin, 12, 7, 4)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(c1, c4) || c1[12] != 324932 {
		t.Fatalf("mismatch: c1=%v c4=%v", c1, c4)
	}
}

func TestDeterministicPrefixOrder(t *testing.T) {
	_, a := GenPrefixes(10, 6)
	_, b := GenPrefixes(10, 6)
	if !reflect.DeepEqual(a, b) {
		t.Fatal("prefix order not deterministic")
	}
}

func TestContextTimeout(t *testing.T) {
	ctx, cancel := context.WithTimeout(context.Background(), 1*time.Nanosecond)
	defer cancel()
	time.Sleep(time.Millisecond)
	if _, _, err := RunSharded(ctx, sawBinPath(t), 18, 7, 2); err == nil {
		t.Fatal("expected timeout error")
	}
}
