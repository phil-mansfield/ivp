package ivp

import (
	"fmt"
	"math"
)

// AccelerationField is an interface representing objects which can compute
// accelerations. (Note: I'm not sure how big of a performance impact the use
// of an interface here has. Worth testing!)
type AccelerationField interface {
	// Acc computes the acceleration at point x at time t and writes the
	// acceleration to out.
	Acc(t float64, x, out []float64)
}

// PlummerMass is a the acceleration field of a single point mass with a position
// X, mass M, and gravitational constant G.
type PointMass struct {
	X []float64
	M, G float64
}

// Acc computtes the acceleration from a point mass at point x at time t and
// writes the acceleration to out. Note that this function takes a performance
// hit because we made the interface take []flaot64s as inputs and not
// [3]float64s, since it leads to more pointer hopping and stops loops form being
// unrolled.
func (field *PointMass) Acc(t float64, x, out []float64) {
	r2 := 0.0 // r^2
	for k := range field.X {
		dx := field.X[k] - x[k]
		r2 += dx*dx
	}

	r := math.Sqrt(r2)
	aNorm := field.M*field.G/(r*r*r)

	for k := range field.X {
		out[k] = (field.X[k] - x[k]) * aNorm
	}
}

// Leapfrog performs a leapfrog (2nd order velocity Verlet) integration. acc
// is the acceleration field, x0 and v0 is the initial coordinates and velocity,
// t0 is the initial time, dt is the timestep, and xOut/vOut are the outputs.
// Not every timestep needs to be output, and you can control how many timesteps
// occur betwen every output with stepsPerOutput.
//
// Note that for leapfrog to work, all you need is for something that looks
// like an acceleration field, not an actual acceleration field, i.e. you just
// need the differential equation to be d^2x/dt^2 = a(t, x), where x is allowed
// to be a vector. So, for example, you could construct your AccelerationField
// so it packs a bunch of particles into he x and v fields simultaneously if
// you're getting killed by interface and function pointer overhead.
func Leapfrog(
	field AccelerationField, x0, v0 []float64, // System description
	t0, dt float64, stepsPerOutput int, // Integration specifics
	xOut, vOut [][]float64, // Output
) {
	// The first output x and v will always be the input x and v
	xOut[0], vOut[0] = x0, v0

	// Allocate various internal vectors. "prev" means (...)_i and "next"
	// means (...)_i+1 in the standard mathematical notation.
	dim := len(x0)
	xPrev, xNext := make([]float64, dim), make([]float64, dim)
	vPrev, vHalf := make([]float64, dim), make([]float64, dim)
	vNext        := make([]float64, dim)
	aPrev, aNext := make([]float64, dim), make([]float64, dim)

	// Intiailize all "prev" variables
	copy(xPrev, x0)
	copy(vPrev, v0)
	field.Acc(t0, x0, aPrev)
	j := 1 // Index into the output arrays
	nSteps := stepsPerOutput*(len(xOut) - 1) + 1
	for i := 1; i < nSteps; i++ {
		for k := 0; k < dim; k++ {
			vHalf[k] = vPrev[k] + aPrev[k]*dt/2
			xNext[k] = xPrev[k] + vHalf[k]*dt
		}
		tNext := t0 + dt*float64(i+1)
		field.Acc(tNext, xNext, aNext)

		for k := 0; k < 3; k++ {
			vNext[k] = vHalf[k] + aNext[k]*dt/2
		}

		copy(aPrev, aNext)
		copy(xPrev, xNext)
		copy(vPrev, vNext)

		if i == j*stepsPerOutput {
			// Originally I had a bug here where I wrote
			// xOut[j], vOut[j] = xPrev, vPrev
			// Quiz: why was that wrong?
			copy(xOut[j], xPrev)
			copy(vOut[j], vPrev)
			j++
		}
	}
}

// MakeVectors allocates a slice of uniform-size vectors. These are all slices
// into the same underlying array to keep it cache-local. 
func MakeVectors(n, dim int) [][]float64 {
	// This will be less efficient than just allocating an array of, say,
	// [][3]float64 vectors, since you need an extra level of indirection.
	buf := make([]float64, n*dim)
	out := make([][]float64, n)
	for i := range out {
		out[i] = buf[i*dim: (i+1)*dim]
	}
	return out
}


func PrintPhaseSpaceTable(x, v [][]float64, verb string) {
	for i := range x {
		for k := range x[i] {
			fmt.Printf(verb, x[i][k])
			fmt.Printf(" ")
		}
		for k := range v[i] {
			fmt.Printf(verb, v[i][k])
			fmt.Printf(" ")
		}
		fmt.Println()
	}
}

// I didn't have time to write real tests, so I just wrote a main() function
// that tests it.

/*
func main() {
	field := &PointMass{ []float64{ 0, 0, 0 }, 1, 1 }
	x0, v0 := []float64{ 1, 0, 0 }, []float64{ 0, 1, 0 }

	x, v := MakeVectors(5, 3), MakeVectors(5, 3)
	stepsPerOutput := 250
	t0, dt := 0.0, 2*math.Pi/1e3

	Leapfrog(field, x0, v0, t0, dt, stepsPerOutput, x, v)
	PrintPhaseSpaceTable(x, v, "%10.5f")
}
*/	