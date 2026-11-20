package parallel

import (
	"runtime"
	"sync"
	"sync/atomic"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
)

// DefaultWorkers mengembalikan jumlah CPU logis.
func DefaultWorkers() int {
	w := runtime.NumCPU()
	if w < 1 {
		return 1
	}
	return w
}

// RunChunks mengeksekusi total item yang dibagi ke nChunks chunk tetap menggunakan nWorkers goroutine.
// Karena batas chunk dan seed hanya bergantung pada indeks chunk (bukan nWorkers),
// hasil reduksi identik bit-demi-bit untuk 1 maupun banyak worker.
func RunChunks[T any](
	total int,
	nChunks int,
	nWorkers int,
	baseSeed uint64,
	workerFn func(chunkIdx, start, end int, seed uint64) T,
) []T {
	if total <= 0 {
		return nil
	}
	if nChunks < 1 {
		nChunks = 1
	}
	if nChunks > total {
		nChunks = total
	}
	if nWorkers < 1 {
		nWorkers = 1
	}
	if nWorkers > nChunks {
		nWorkers = nChunks
	}

	chunkSize := (total + nChunks - 1) / nChunks
	results := make([]T, nChunks)
	var nextIdx int64 = -1
	var wg sync.WaitGroup

	for w := 0; w < nWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for {
				idx := int(atomic.AddInt64(&nextIdx, 1))
				if idx >= nChunks {
					return
				}
				start := idx * chunkSize
				end := start + chunkSize
				if end > total {
					end = total
				}
				if start >= end {
					continue
				}
				cSeed := rng.SplitMix64(baseSeed ^ uint64(idx+1))
				results[idx] = workerFn(idx, start, end, cSeed)
			}
		}()
	}
	wg.Wait()
	return results
}
