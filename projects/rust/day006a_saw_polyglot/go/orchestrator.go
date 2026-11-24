package main

import (
	"bufio"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os/exec"
	"sync"
	"time"
)

type PrefixTask struct {
	ID     uint64 `json:"id"`
	Nmax   int    `json:"nmax"`
	Weight uint64 `json:"weight"`
	Steps  string `json:"steps"`
}

type WorkerResponse struct {
	ID     uint64   `json:"id"`
	OK     bool     `json:"ok"`
	Counts []uint64 `json:"counts"`
	Error  string   `json:"error,omitempty"`
}

type SingleCountOutput struct {
	Nmax     int        `json:"nmax"`
	Counts   []uint64   `json:"counts"`
	HeatmapN int        `json:"heatmap_n,omitempty"`
	Heatmap  [][]uint64 `json:"heatmap,omitempty"`
}

type PivotOutput struct {
	N              int        `json:"n"`
	Seed           uint64     `json:"seed"`
	R2Mean         float64    `json:"r2_mean"`
	R2SE           float64    `json:"r2_se"`
	AcceptanceRate float64    `json:"acceptance_rate"`
	TotalSteps     int        `json:"total_steps"`
	Samples        [][][2]int `json:"samples"`
}

type WorkerTiming struct {
	Workers     int      `json:"workers"`
	ElapsedSec  float64  `json:"elapsed_sec"`
	Speedup     float64  `json:"speedup"`
	CountsMatch bool     `json:"counts_match"`
	Counts      []uint64 `json:"counts"`
}

type OrchestratedBundle struct {
	Nmax            int               `json:"nmax"`
	PrefixK         int               `json:"prefix_k"`
	NumPrefixTasks  int               `json:"num_prefix_tasks"`
	BaseSeed        uint64            `json:"base_seed"`
	RustSingle      SingleCountOutput `json:"rust_single"`
	RustSingleSec   float64           `json:"rust_single_sec"`
	GoShardedCounts []uint64          `json:"go_sharded_counts"`
	WorkerTimings   []WorkerTiming    `json:"worker_timings"`
	WorkerInvariant bool              `json:"worker_invariant"`
	Pivots          []PivotOutput     `json:"pivots"`
}

var dirs = [4]struct {
	dx, dy int
	code   byte
}{
	{1, 0, 'E'},
	{0, 1, 'N'},
	{-1, 0, 'W'},
	{0, -1, 'S'},
}

// GeneratePrefixesAndShortCounts enumerates all SAWs up to kEff with first step forced East (weight 4).
// It returns exact c_n for n < kEff in shortCounts, and all prefix strings of exact length kEff.
func GeneratePrefixesAndShortCounts(nmax, k int) ([]uint64, []PrefixTask) {
	shortCounts := make([]uint64, nmax+1)
	shortCounts[0] = 1
	if nmax <= 0 {
		return shortCounts, nil
	}
	kEff := k
	if kEff > nmax {
		kEff = nmax
	}
	if kEff < 1 {
		kEff = 1
	}

	const offset = 32
	var visited [64]uint64
	setBit := func(x, y int) { visited[y+offset] |= uint64(1) << uint(x+offset) }
	clrBit := func(x, y int) { visited[y+offset] &^= uint64(1) << uint(x+offset) }
	hasBit := func(x, y int) bool { return (visited[y+offset] & (uint64(1) << uint(x+offset))) != 0 }

	setBit(0, 0)
	setBit(1, 0)
	buf := make([]byte, 1, kEff)
	buf[0] = 'E'
	var tasks []PrefixTask

	var dfs func(x, y, step int)
	dfs = func(x, y, step int) {
		if step < kEff {
			shortCounts[step] += 4
		}
		if step == kEff {
			tasks = append(tasks, PrefixTask{
				ID:     uint64(len(tasks)),
				Nmax:   nmax,
				Weight: 4,
				Steps:  string(buf),
			})
			return
		}
		for _, d := range dirs {
			nx, ny := x+d.dx, y+d.dy
			if !hasBit(nx, ny) {
				setBit(nx, ny)
				buf = append(buf, d.code)
				dfs(nx, ny, step+1)
				buf = buf[:len(buf)-1]
				clrBit(nx, ny)
			}
		}
	}

	dfs(1, 0, 1)
	return shortCounts, tasks
}

func ParseWorkerResponseLine(line []byte) (WorkerResponse, error) {
	var resp WorkerResponse
	if err := json.Unmarshal(line, &resp); err != nil {
		return WorkerResponse{}, err
	}
	if !resp.OK {
		return resp, fmt.Errorf("worker error: %s", resp.Error)
	}
	return resp, nil
}

type persistentWorker struct {
	cmd    *exec.Cmd
	stdin  io.WriteCloser
	reader *bufio.Reader
}

func startPersistentWorker(ctx context.Context, sawBin string) (*persistentWorker, error) {
	cmd := exec.CommandContext(ctx, sawBin, "count", "--stdin")
	stdin, err := cmd.StdinPipe()
	if err != nil {
		return nil, err
	}
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return nil, err
	}
	if err := cmd.Start(); err != nil {
		return nil, err
	}
	return &persistentWorker{
		cmd:    cmd,
		stdin:  stdin,
		reader: bufio.NewReaderSize(stdout, 64*1024),
	}, nil
}

func (w *persistentWorker) close() error {
	_ = w.stdin.Close()
	return w.cmd.Wait()
}

// RunShardedEnumeration runs persistent Rust worker processes over JSON stdin/stdout pipes
// and sums integer counts in deterministic fixed task-ID order.
func RunShardedEnumeration(ctx context.Context, sawBin string, nmax, k, numWorkers int) ([]uint64, int, error) {
	if numWorkers < 1 {
		return nil, 0, errors.New("numWorkers must be >= 1")
	}
	shortCounts, tasks := GeneratePrefixesAndShortCounts(nmax, k)
	if len(tasks) == 0 {
		return shortCounts, 0, nil
	}

	wCount := numWorkers
	if wCount > len(tasks) {
		wCount = len(tasks)
	}

	taskCh := make(chan PrefixTask, len(tasks))
	for _, t := range tasks {
		taskCh <- t
	}
	close(taskCh)

	results := make([][]uint64, len(tasks))
	errCh := make(chan error, wCount)
	var wg sync.WaitGroup

	for wi := 0; wi < wCount; wi++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			pw, err := startPersistentWorker(ctx, sawBin)
			if err != nil {
				errCh <- err
				return
			}
			defer func() { _ = pw.close() }()

			for task := range taskCh {
				if err := ctx.Err(); err != nil {
					errCh <- err
					return
				}
				payload, err := json.Marshal(task)
				if err != nil {
					errCh <- err
					return
				}
				payload = append(payload, '\n')
				if _, err := pw.stdin.Write(payload); err != nil {
					errCh <- err
					return
				}
				line, err := pw.reader.ReadBytes('\n')
				if err != nil {
					errCh <- err
					return
				}
				resp, err := ParseWorkerResponseLine(line)
				if err != nil {
					errCh <- err
					return
				}
				if resp.ID != task.ID {
					errCh <- fmt.Errorf("mismatched response ID: got %d want %d", resp.ID, task.ID)
					return
				}
				results[task.ID] = resp.Counts
			}
		}()
	}

	wg.Wait()
	close(errCh)
	for err := range errCh {
		if err != nil {
			return nil, len(tasks), err
		}
	}

	kEff := k
	if kEff > nmax {
		kEff = nmax
	}
	total := make([]uint64, nmax+1)
	for n := 0; n < kEff; n++ {
		total[n] = shortCounts[n]
	}
	// Deterministic integer summation in fixed task ID order (0..len(tasks)-1), zero floats.
	for id := 0; id < len(tasks); id++ {
		row := results[id]
		if len(row) != nmax+1 {
			return nil, len(tasks), fmt.Errorf("missing or incomplete result at task %d", id)
		}
		for n := kEff; n <= nmax; n++ {
			total[n] += row[n]
		}
	}
	return total, len(tasks), nil
}

func RunRustSingleCount(ctx context.Context, sawBin string, nmax, heatmapN int) (SingleCountOutput, float64, error) {
	t0 := time.Now()
	args := []string{"count", "--nmax", fmt.Sprintf("%d", nmax)}
	if heatmapN > 0 {
		args = append(args, "--heatmap-n", fmt.Sprintf("%d", heatmapN))
	}
	cmd := exec.CommandContext(ctx, sawBin, args...)
	out, err := cmd.Output()
	elapsed := time.Since(t0).Seconds()
	if err != nil {
		return SingleCountOutput{}, elapsed, err
	}
	var parsed SingleCountOutput
	if err := json.Unmarshal(out, &parsed); err != nil {
		return SingleCountOutput{}, elapsed, err
	}
	return parsed, elapsed, nil
}

// DeriveSeed deterministically derives a per-N seed from baseSeed using SplitMix64 mixing.
func DeriveSeed(baseSeed uint64, n int) uint64 {
	z := baseSeed + uint64(n)*0x9E3779B97F4A7C15
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
	z = (z ^ (z >> 27)) * 0x94D049BB133111EB
	return z ^ (z >> 31)
}

func RunPivotSuite(ctx context.Context, sawBin string, baseSeed uint64, ns []int) ([]PivotOutput, error) {
	outputs := make([]PivotOutput, len(ns))
	errCh := make(chan error, len(ns))
	var wg sync.WaitGroup

	for idx, nVal := range ns {
		wg.Add(1)
		go func(i, n int) {
			defer wg.Done()
			seed := DeriveSeed(baseSeed, n)
			samples := 8
			stride := 800
			if n == 400 {
				samples = 120
				stride = 1200
			}
			burnIn := 30 * n
			if burnIn < 3000 {
				burnIn = 3000
			}
			args := []string{
				"pivot",
				"--n", fmt.Sprintf("%d", n),
				"--seed", fmt.Sprintf("%d", seed),
				"--burn-in", fmt.Sprintf("%d", burnIn),
				"--batches", "64",
				"--batch-size", "4000",
				"--samples", fmt.Sprintf("%d", samples),
				"--sample-stride", fmt.Sprintf("%d", stride),
			}
			cmd := exec.CommandContext(ctx, sawBin, args...)
			raw, err := cmd.Output()
			if err != nil {
				errCh <- fmt.Errorf("pivot N=%d failed: %w", n, err)
				return
			}
			var po PivotOutput
			if err := json.Unmarshal(raw, &po); err != nil {
				errCh <- fmt.Errorf("pivot N=%d json error: %w", n, err)
				return
			}
			outputs[i] = po
		}(idx, nVal)
	}

	wg.Wait()
	close(errCh)
	for err := range errCh {
		if err != nil {
			return nil, err
		}
	}
	return outputs, nil
}

func SlicesEqualU64(a, b []uint64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
