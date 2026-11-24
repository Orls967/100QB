package main

import (
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"time"
)

func main() {
	sawBin := flag.String("saw-bin", "", "Path to compiled Rust saw binary")
	nmax := flag.Int("nmax", 22, "Maximum SAW length for exact enumeration")
	prefixK := flag.Int("k", 8, "Prefix length K for Go task sharding")
	heatmapN := flag.Int("heatmap-n", 16, "Length n for endpoint heatmap")
	seed := flag.Uint64("seed", 20261124, "Master seed")
	timeoutSec := flag.Int("timeout-sec", 120, "Context timeout in seconds")
	outPath := flag.String("out", "", "Output JSON file path")
	flag.Parse()

	if *sawBin == "" {
		fmt.Fprintln(os.Stderr, "error: --saw-bin is required")
		os.Exit(1)
	}

	ctx, cancel := context.WithTimeout(context.Background(), time.Duration(*timeoutSec)*time.Second)
	defer cancel()

	singleOut, singleSec, err := RunRustSingleCount(ctx, *sawBin, *nmax, *heatmapN)
	if err != nil {
		fmt.Fprintf(os.Stderr, "single count error: %v\n", err)
		os.Exit(1)
	}

	maxW := runtime.NumCPU()
	if maxW <= 4 {
		maxW = 8
	}
	workerGrid := []int{1, 2, 4, maxW}
	timings := make([]WorkerTiming, 0, len(workerGrid))
	var refCounts []uint64
	numTasks := 0
	allInvariant := true

	for idx, w := range workerGrid {
		t0 := time.Now()
		counts, nTasks, err := RunShardedEnumeration(ctx, *sawBin, *nmax, *prefixK, w)
		elapsed := time.Since(t0).Seconds()
		if err != nil {
			fmt.Fprintf(os.Stderr, "sharded count (workers=%d) error: %v\n", w, err)
			os.Exit(1)
		}
		numTasks = nTasks
		if idx == 0 {
			refCounts = counts
		}
		match := SlicesEqualU64(counts, refCounts) && SlicesEqualU64(counts, singleOut.Counts)
		if !match {
			allInvariant = false
		}
		speedup := 1.0
		if idx > 0 && elapsed > 0 {
			speedup = timings[0].ElapsedSec / elapsed
		}
		timings = append(timings, WorkerTiming{
			Workers:     w,
			ElapsedSec:  elapsed,
			Speedup:     speedup,
			CountsMatch: match,
			Counts:      counts,
		})
	}

	pivots, err := RunPivotSuite(ctx, *sawBin, *seed, []int{50, 100, 200, 400, 800})
	if err != nil {
		fmt.Fprintf(os.Stderr, "pivot suite error: %v\n", err)
		os.Exit(1)
	}

	bundle := OrchestratedBundle{
		Nmax:            *nmax,
		PrefixK:         *prefixK,
		NumPrefixTasks:  numTasks,
		BaseSeed:        *seed,
		RustSingle:      singleOut,
		RustSingleSec:   singleSec,
		GoShardedCounts: refCounts,
		WorkerTimings:   timings,
		WorkerInvariant: allInvariant,
		Pivots:          pivots,
	}

	encoded, err := json.MarshalIndent(bundle, "", "  ")
	if err != nil {
		fmt.Fprintf(os.Stderr, "json marshal error: %v\n", err)
		os.Exit(1)
	}

	if *outPath != "" {
		if err := os.MkdirAll(filepath.Dir(*outPath), 0o755); err != nil {
			fmt.Fprintf(os.Stderr, "mkdir error: %v\n", err)
			os.Exit(1)
		}
		if err := os.WriteFile(*outPath, encoded, 0o644); err != nil {
			fmt.Fprintf(os.Stderr, "write error: %v\n", err)
			os.Exit(1)
		}
	} else {
		fmt.Println(string(encoded))
	}
}
