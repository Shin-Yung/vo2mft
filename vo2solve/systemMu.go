package vo2solve

import (
	"github.com/tflovorn/cmatrix"
	"github.com/tflovorn/scExplorer/bzone"
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Return the absolute error and gradient of the Mu equation w.r.t. the given
// variables (which should be fixed to ["M", "W", "Mu"] for this case).
func AbsErrorMu(env *Environment, variables []string) solve.Diffable {
	F := func(v vec.Vector) (float64, error) {
		env.Set(v, variables)
		L := env.BZPointsPerDim
		innerClosure := func(k vec.Vector) float64 {
			return innerMu(env, k)
		}
		lhs := 1.0
		rhs := 0.5 * bzone.Avg(L, 3, innerClosure)
		return lhs - rhs, nil
	}
	h := 1e-6
	epsabs := 1e-4
	return solve.SimpleDiffable(F, len(variables), h, epsabs)
}

func innerMu(env *Environment, k vec.Vector) float64 {
	H := ElHamiltonian(env, k)
	dim, _ := H.Dims()
	evals, _ := cmatrix.Eigensystem(H)
	sum := 0.0
	for alpha := 0; alpha < dim; alpha++ {
		// Mu is included in H, so not included here.
		sum += env.Fermi(evals[alpha])
	}
	// Multiply by 2 for spin degeneracy.
	return 2.0 * sum
}
