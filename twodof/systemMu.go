package twodof

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
		H := cmatrix.NewCMatrixGSL(4, 4)
		work, evals, evecs := cmatrix.HermEigensystemSetup(H)
		innerClosure := func(k vec.Vector) float64 {
			return innerMu(env, k, H, work, evals, evecs)
		}
		lhs := 1.0
		rhs := bzone.Avg(L, 3, innerClosure)

		H.Destroy()
		cmatrix.HermEigensystemCleanup(work, evals, evecs)
		return lhs - rhs, nil
	}
	h := 1e-6
	epsabs := 1e-4
	return solve.SimpleDiffable(F, len(variables), h, epsabs)
}

func innerMu(env *Environment, k vec.Vector, H *cmatrix.CMatrixGSL, work *cmatrix.HermWorkGSL, evals *cmatrix.VectorGSL, evecs *cmatrix.CMatrixGSL) float64 {
	ElHamiltonian(env, k, H)
	dim, _ := H.Dims()
	cmatrix.HermEigensystem(H, work, evals, evecs)
	sum := 0.0
	for alpha := 0; alpha < dim; alpha++ {
		// Mu is included in H, so not included here.
		sum += env.Fermi(evals.At(alpha))
	}
	// Multiply by 2 for spin degeneracy.
	return 2.0 * sum
}
