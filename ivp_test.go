package ivp

import (
	"testing"
	"math"
)

func benchmarkLeapfrog(b *testing.B, nSteps int) {
	// This benchmark just runs something similar to the commented-out
	// main function
	field := &PointMass{ []float64{ 0, 0, 0 }, 1, 1 }
	x0, v0 := []float64{ 1, 0, 0 }, []float64{ 0, 1, 0 }
	x, v := MakeVectors(2, 3), MakeVectors(2, 3)
	stepsPerOutput := nSteps
	t0, dt := 0.0, 2*math.Pi/1e3

	// All the actual benchmarking code is below here

	// Start the timer here
	b.ResetTimer()
	// b.N will be adaptively increased by Go's benchmarking framework
	for i := 0; i < b.N; i++ {
		// If I wanted to pause the timer to do some setup here, I would have done
		// b.StartTimer()
		// SetupFunction()
		// b.StopTimer()

		// Run the code you're benchmarking here.
		Leapfrog(field, x0, v0, t0, dt, stepsPerOutput, x, v)
	}
}


func BenchmarkLeapfrog1e1Steps(b *testing.B) { benchmarkLeapfrog(b, 10) }
func BenchmarkLeapfrog1e2Steps(b *testing.B) { benchmarkLeapfrog(b, 100) }
func BenchmarkLeapfrog1e3Steps(b *testing.B) { benchmarkLeapfrog(b, 1000) }
func BenchmarkLeapfrog1e4Steps(b *testing.B) { benchmarkLeapfrog(b, 10000) }
