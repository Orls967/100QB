package main

import (
	"bufio"
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"sync"
	"time"
)

type Job struct {
	Idx     int    `json:"-"`
	L       int    `json:"L"`
	Seed    uint64 `json:"seed"`
	Samples int    `json:"samples"`
}

type JobResp struct {
	L            int      `json:"L"`
	Samples      int      `json:"samples"`
	SumLargest   uint64   `json:"sum_largest"`
	SumLargestSq uint64   `json:"sum_largest_sq"`
	Hist         []uint64 `json:"hist"`
}

type LatticeSummary struct {
	L            int      `json:"L"`
	NBonds       int      `json:"n_bonds"`
	Samples      int      `json:"samples"`
	Chunks       int      `json:"chunks"`
	SumLargest   uint64   `json:"sum_largest"`
	SumLargestSq uint64   `json:"sum_largest_sq"`
	Hist         []uint64 `json:"hist"`
}

func splitMix64(x uint64) uint64 {
	x += 0x9e3779b97f4a7c15
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9
	x = (x ^ (x >> 27)) * 0x94d049bb133111eb
	return x ^ (x >> 31)
}

func ChunkSeed(master uint64, l, chunk int) uint64 {
	h := splitMix64(master ^ (uint64(l) * 0x517cc1b727220a95))
	return splitMix64(h ^ (uint64(chunk) * 0x94d049bb133111eb))
}

func BuildJobs(master uint64, sizes []int, smap map[int]int, chunks int) []Job {
	var jobs []Job
	for _, l := range sizes {
		per := smap[l] / chunks
		for c := 0; c < chunks; c++ {
			jobs = append(jobs, Job{Idx: len(jobs), L: l, Seed: ChunkSeed(master, l, c), Samples: per})
		}
	}
	return jobs
}

func RunPool(ctx context.Context, binPath string, jobs []Job, workers int) ([]LatticeSummary, error) {
	if workers < 1 {
		workers = 1
	}
	resps := make([]JobResp, len(jobs))
	ch := make(chan Job, len(jobs))
	for _, jb := range jobs {
		ch <- jb
	}
	close(ch)
	errCh := make(chan error, workers)
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			cmd := exec.CommandContext(ctx, binPath, "job")
			stdin, err := cmd.StdinPipe()
			if err != nil {
				errCh <- err
				return
			}
			stdout, err := cmd.StdoutPipe()
			if err != nil {
				errCh <- err
				return
			}
			if err := cmd.Start(); err != nil {
				errCh <- err
				return
			}
			defer func() { _ = stdin.Close(); _ = cmd.Wait() }()
			sc := bufio.NewScanner(stdout)
			sc.Buffer(make([]byte, 64*1024), 8*1024*1024)
			for jb := range ch {
				b, _ := json.Marshal(jb)
				if _, err := stdin.Write(append(b, '\n')); err != nil {
					errCh <- err
					return
				}
				if !sc.Scan() {
					if ctx.Err() != nil {
						errCh <- ctx.Err()
					} else {
						errCh <- fmt.Errorf("worker eof: %v", sc.Err())
					}
					return
				}
				if err := json.Unmarshal(sc.Bytes(), &resps[jb.Idx]); err != nil {
					errCh <- err
					return
				}
			}
		}()
	}
	wg.Wait()
	close(errCh)
	if err := <-errCh; err != nil {
		return nil, err
	}
	byL, order := map[int]*LatticeSummary{}, []int{}
	for i, jb := range jobs {
		r := resps[i]
		s, ok := byL[jb.L]
		if !ok {
			nb := 2*jb.L*jb.L - 1
			s = &LatticeSummary{L: jb.L, NBonds: nb, Hist: make([]uint64, nb+1)}
			byL[jb.L], order = s, append(order, jb.L)
		}
		s.Samples += r.Samples
		s.Chunks++
		s.SumLargest += r.SumLargest
		s.SumLargestSq += r.SumLargestSq
		for k, v := range r.Hist {
			s.Hist[k] += v
		}
	}
	out := make([]LatticeSummary, len(order))
	for i, l := range order {
		out[i] = *byL[l]
	}
	return out, nil
}

func main() {
	bin := flag.String("bin", "", "percolate binary")
	outPath := flag.String("out", "", "output json")
	seed := flag.Uint64("seed", 20261124, "master seed")
	tsec := flag.Int("timeout", 55, "timeout s")
	flag.Parse()
	ctx, cancel := context.WithTimeout(context.Background(), time.Duration(*tsec)*time.Second)
	defer cancel()
	nw := min(runtime.NumCPU(), 8)
	smap := map[int]int{16: 20000, 32: 20000, 64: 20000, 128: 6000}
	sum, err := RunPool(ctx, *bin, BuildJobs(*seed, []int{16, 32, 64, 128}, smap, 20), nw)
	if err != nil {
		os.Exit(1)
	}
	ij := BuildJobs(*seed^0xabcdef, []int{16, 32}, map[int]int{16: 2000, 32: 2000}, 10)
	s1, e1 := RunPool(ctx, *bin, ij, 1)
	sn, en := RunPool(ctx, *bin, ij, nw)
	if e1 != nil || en != nil {
		os.Exit(1)
	}
	b1, _ := json.Marshal(s1)
	bn, _ := json.Marshal(sn)
	raw, _ := json.MarshalIndent(map[string]any{"seed": *seed, "workers": nw, "worker_invariant": string(b1) == string(bn), "lattices": sum}, "", "  ")
	if os.WriteFile(*outPath, raw, 0o644) != nil {
		os.Exit(1)
	}
}
