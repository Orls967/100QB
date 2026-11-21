package parallel

import "testing"

func TestRunChunksWorkerInvariance(t *testing.T) {
	fn := func(i int) int { return i*i + 3*i + 1 }
	a := RunChunks(50, 1, fn)
	b := RunChunks(50, 8, fn)
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("mismatch at %d: %d vs %d", i, a[i], b[i])
		}
	}
}
