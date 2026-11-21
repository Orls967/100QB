package parallel

import (
	"runtime"
	"sync"
)

// DefaultWorkers mengembalikan jumlah goroutine pekerja bawaan (min(8, NumCPU)).
func DefaultWorkers() int {
	w := runtime.NumCPU()
	if w < 1 {
		return 1
	}
	if w > 8 {
		return 8
	}
	return w
}

// RunChunks menjalankan nChunks unit kerja secara paralel memakai nWorkers goroutine
// dan mengembalikan irisan hasil berurut sesuai chunkIdx agar reduksi bit-identik.
func RunChunks[T any](nChunks, nWorkers int, fn func(chunkIdx int) T) []T {
	if nChunks <= 0 {
		return nil
	}
	out := make([]T, nChunks)
	if nWorkers <= 1 || nChunks == 1 {
		for i := 0; i < nChunks; i++ {
			out[i] = fn(i)
		}
		return out
	}
	if nWorkers > nChunks {
		nWorkers = nChunks
	}
	jobs := make(chan int, nChunks)
	for i := 0; i < nChunks; i++ {
		jobs <- i
	}
	close(jobs)

	var wg sync.WaitGroup
	wg.Add(nWorkers)
	for w := 0; w < nWorkers; w++ {
		go func() {
			defer wg.Done()
			for idx := range jobs {
				out[idx] = fn(idx)
			}
		}()
	}
	wg.Wait()
	return out
}
