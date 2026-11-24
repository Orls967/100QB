package main

import (
	"bufio"
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"time"
)

type PrefixTask struct {
	ID     uint64 `json:"id"`
	Nmax   int    `json:"nmax"`
	Weight uint64 `json:"weight"`
	Steps  string `json:"steps"`
}

type WorkerResp struct {
	ID     uint64   `json:"id"`
	OK     bool     `json:"ok"`
	Counts []uint64 `json:"counts"`
	Error  string   `json:"error,omitempty"`
}

type SingleOut struct {
	Nmax     int        `json:"nmax"`
	Counts   []uint64   `json:"counts"`
	HeatmapN int        `json:"heatmap_n,omitempty"`
	Heatmap  [][]uint64 `json:"heatmap,omitempty"`
}

type PivotOut struct {
	N       int        `json:"n"`
	Seed    uint64     `json:"seed"`
	R2Mean  float64    `json:"r2_mean"`
	R2SE    float64    `json:"r2_se"`
	Acc     float64    `json:"acc"`
	Samples [][][2]int `json:"samples"`
}

func GenPrefixes(nmax, k int) ([]uint64, []PrefixTask) {
	short := make([]uint64, nmax+1)
	short[0] = 1
	if nmax <= 0 {
		return short, nil
	}
	ke := min(max(k, 1), nmax)
	var vis [64]uint64
	set := func(x, y int) { vis[y+32] |= 1 << uint(x+32) }
	clr := func(x, y int) { vis[y+32] &^= 1 << uint(x+32) }
	has := func(x, y int) bool { return (vis[y+32] & (1 << uint(x+32))) != 0 }
	set(0, 0)
	set(1, 0)
	buf := []byte{'E'}
	var tasks []PrefixTask
	dirs := [4]struct {
		dx, dy int
		c      byte
	}{{1, 0, 'E'}, {0, 1, 'N'}, {-1, 0, 'W'}, {0, -1, 'S'}}
	var dfs func(x, y, step int)
	dfs = func(x, y, step int) {
		if step < ke {
			short[step] += 4
		}
		if step == ke {
			tasks = append(tasks, PrefixTask{ID: uint64(len(tasks)), Nmax: nmax, Weight: 4, Steps: string(buf)})
			return
		}
		for _, d := range dirs {
			nx, ny := x+d.dx, y+d.dy
			if !has(nx, ny) {
				set(nx, ny)
				buf = append(buf, d.c)
				dfs(nx, ny, step+1)
				buf = buf[:len(buf)-1]
				clr(nx, ny)
			}
		}
	}
	dfs(1, 0, 1)
	return short, tasks
}

func RunSharded(ctx context.Context, bin string, nmax, k, workers int) ([]uint64, int, error) {
	short, tasks := GenPrefixes(nmax, k)
	if len(tasks) == 0 {
		return short, 0, nil
	}
	w := min(max(workers, 1), len(tasks))
	ch := make(chan PrefixTask, len(tasks))
	for _, t := range tasks {
		ch <- t
	}
	close(ch)
	res := make([][]uint64, len(tasks))
	errCh := make(chan error, w)
	var wg sync.WaitGroup
	for i := 0; i < w; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			cmd := exec.CommandContext(ctx, bin, "count", "--stdin")
			inPipe, err := cmd.StdinPipe()
			if err != nil {
				errCh <- err
				return
			}
			outPipe, err := cmd.StdoutPipe()
			if err != nil {
				errCh <- err
				return
			}
			if err := cmd.Start(); err != nil {
				errCh <- err
				return
			}
			defer func() { _ = inPipe.Close(); _ = cmd.Wait() }()
			rd := bufio.NewReaderSize(outPipe, 32*1024)
			for t := range ch {
				if err := ctx.Err(); err != nil {
					errCh <- err
					return
				}
				b, _ := json.Marshal(t)
				if _, err := inPipe.Write(append(b, '\n')); err != nil {
					errCh <- err
					return
				}
				line, err := rd.ReadBytes('\n')
				if err != nil {
					errCh <- err
					return
				}
				var rp WorkerResp
				if err := json.Unmarshal(line, &rp); err != nil || !rp.OK || rp.ID != t.ID {
					errCh <- fmt.Errorf("bad worker resp: %v %s", err, rp.Error)
					return
				}
				res[t.ID] = rp.Counts
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
	ke := min(max(k, 1), nmax)
	tot := make([]uint64, nmax+1)
	copy(tot[:ke], short[:ke])
	for id := 0; id < len(tasks); id++ {
		for n := ke; n <= nmax; n++ {
			tot[n] += res[id][n]
		}
	}
	return tot, len(tasks), nil
}

func DeriveSeed(base uint64, n int) uint64 {
	z := base + uint64(n)*0x9E3779B97F4A7C15
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
	z = (z ^ (z >> 27)) * 0x94D049BB133111EB
	return z ^ (z >> 31)
}

func runPivotOne(ctx context.Context, bin string, base uint64, n, samples int) (PivotOut, error) {
	cmd := exec.CommandContext(ctx, bin, "pivot", "--n", fmt.Sprint(n),
		"--seed", fmt.Sprint(DeriveSeed(base, n)), "--batches", "50", "--batch-size", "4000",
		"--samples", fmt.Sprint(samples), "--stride", "1000")
	raw, err := cmd.Output()
	if err != nil {
		return PivotOut{}, err
	}
	var po PivotOut
	return po, json.Unmarshal(raw, &po)
}

func main() {
	bin := flag.String("saw-bin", "", "saw binary")
	nmax := flag.Int("nmax", 18, "Nmax")
	k := flag.Int("k", 7, "prefix K")
	hmN := flag.Int("heatmap-n", 14, "heatmap n")
	seed := flag.Uint64("seed", 20261124, "seed")
	out := flag.String("out", "", "output json")
	flag.Parse()

	ctx, cancel := context.WithTimeout(context.Background(), 90*time.Second)
	defer cancel()

	rawSingle, err := exec.CommandContext(ctx, *bin, "count", "--nmax", fmt.Sprint(*nmax), "--heatmap-n", fmt.Sprint(*hmN)).Output()
	if err != nil {
		panic(err)
	}
	var single SingleOut
	if err := json.Unmarshal(rawSingle, &single); err != nil {
		panic(err)
	}

	c1, nTasks, err := RunSharded(ctx, *bin, *nmax, *k, 1)
	if err != nil {
		panic(err)
	}
	nw := max(runtime.NumCPU(), 4)
	cn, _, err := RunSharded(ctx, *bin, *nmax, *k, nw)
	if err != nil {
		panic(err)
	}
	inv := len(c1) == len(cn)
	for i := range c1 {
		if c1[i] != cn[i] {
			inv = false
		}
	}

	ns := []int{50, 100, 200, 400}
	pivs := make([]PivotOut, len(ns))
	var wg sync.WaitGroup
	for i, n := range ns {
		wg.Add(1)
		go func(idx, nVal int) {
			defer wg.Done()
			po, e := runPivotOne(ctx, *bin, *seed, nVal, 0)
			if e != nil {
				panic(e)
			}
			pivs[idx] = po
		}(i, n)
	}
	wg.Wait()
	samp300, err := runPivotOne(ctx, *bin, *seed, 300, 100)
	if err != nil {
		panic(err)
	}

	payload := map[string]any{
		"nmax": *nmax, "prefix_k": *k, "num_tasks": nTasks, "seed": *seed,
		"rust_single": single, "go_counts": cn, "worker_invariant": inv, "workers_n": nw,
		"pivots": pivs, "samples_n300": samp300,
	}
	b, _ := json.MarshalIndent(payload, "", "  ")
	_ = os.MkdirAll(filepath.Dir(*out), 0o755)
	if err := os.WriteFile(*out, b, 0o644); err != nil {
		panic(err)
	}
}
