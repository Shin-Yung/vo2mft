package twodof

import (
	"encoding/json"
	"math"
	"math/cmplx"
)
import (
	"github.com/tflovorn/cmatrix"
	"github.com/tflovorn/scExplorer/bzone"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Expect some e.v.'s to be pure real or imaginary - panic if imag/real part
// is greater than zero_threshold.
const zero_threshold = 1e-12

type HoppingEV struct {
	// Value of M01 for which the contained hopping e.v.'s have been calculated.
	m01_cached map[string]float64
	// Value of M12 for which the contained hopping e.v.'s have been calculated.
	m12_cached map[string]float64
	// Value of Mu for which the contained hopping e.v.'s have been calculated.
	mu_cached map[string]float64
	// If hopping e.v.'s have not been calculated yet, init = false.
	init map[string]bool
	// Hopping e.v.'s for odd symmetry (pre-calculated).
	dco float64
}

func NewHoppingEV() *HoppingEV {
	names := []string{"dco"}

	Ds := new(HoppingEV)
	Ds.m01_cached = make(map[string]float64)
	Ds.m12_cached = make(map[string]float64)
	Ds.mu_cached = make(map[string]float64)
	Ds.init = make(map[string]bool)

	for _, name := range names {
		Ds.init[name] = false
	}

	return Ds
}

func (Ds *HoppingEV) Dco(env *Environment) float64 {
	if Ds.cacheOk(env, "dco") {
		return Ds.dco
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	H := cmatrix.NewCMatrixGSL(4, 4)
	work, evals, evecs := cmatrix.HermEigensystemSetup(H)
	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_KQ0(env, k, H, work, evals, evecs)
		// make sure that ev is pure imaginary
		//if math.Abs(real(ev)) > zero_threshold {
		//	panic("Expected pure imaginary value for <c_{k+Q,0}^{\\dagger} c_{k,0}>, got finite real part")
		//}
		// -2i * ev = 2 * imag(ev)
		return 2.0 * math.Sin(k[2]) * imag(ev)
	}
	dco := bzone.Avg(env.BZPointsPerDim, 3, inner)

	// Uncomment to verify that Dco is real.
	/*
		inner_re := func(k vec.Vector) float64 {
			ev := GetEV_K0_KQ0(env, k, H, work, evals, evecs)
			return -2.0 * math.Sin(k[2]) * real(ev)
		}
		dco_re := bzone.Avg(env.BZPointsPerDim, 3, inner_re)
		if math.Abs(dco_re) > zero_threshold {
			panic("Expected real value for Dco, got finite imaginary part.")
		}
	*/

	Ds.init["dco"] = true
	Ds.m01_cached["dco"] = env.M01
	Ds.m12_cached["dco"] = env.M12
	Ds.mu_cached["dco"] = env.Mu
	Ds.dco = dco

	H.Destroy()
	cmatrix.HermEigensystemCleanup(work, evals, evecs)
	return dco
}

// Return true iff the cached evaluation of the given D value is still
// OK to use.
func (Ds *HoppingEV) cacheOk(env *Environment, dname string) bool {
	if !Ds.init[dname] {
		return false
	}
	M01_ok := env.M01 == Ds.m01_cached[dname]
	M12_ok := env.M12 == Ds.m12_cached[dname]
	Mu_ok := env.Mu == Ds.mu_cached[dname]

	return M01_ok && M12_ok && Mu_ok
}

// Convert to string by marshalling to JSON.
// Leave out internal cache data.
func (Ds *HoppingEV) StringEnv(env *Environment) string {
	marshalled := Ds.MarshalEnv(env)
	return marshalled
}

func (Ds *HoppingEV) MarshalEnv(env *Environment) string {
	repr := make(map[string]float64)
	repr["Dco"] = Ds.Dco(env)
	marshalled, err := json.Marshal(repr)
	if err != nil {
		panic(err)
	}
	return string(marshalled)
}

// Evaluate <c^{\dagger}_{k,0} c_{k+Q,0}>.
func GetEV_K0_KQ0(env *Environment, k vec.Vector, H *cmatrix.CMatrixGSL, work *cmatrix.HermWorkGSL, evals *cmatrix.VectorGSL, evecs *cmatrix.CMatrixGSL) complex128 {
	return evalEV(env, k, 1, 2, H, work, evals, evecs)
}

// Evaluate <c^{\dagger}_{indexL} c_{indexR}> where the index values have
// the following correspondence:
// 	1 <--> k, 0 ; 2 <--> k+Q, 0 ; 3 <--> k, 1 ; 4 <--> k+Q, 1
func evalEV(env *Environment, k vec.Vector, indexL, indexR int, H *cmatrix.CMatrixGSL, work *cmatrix.HermWorkGSL, evals *cmatrix.VectorGSL, evecs *cmatrix.CMatrixGSL) complex128 {
	ElHamiltonian(env, k, H)
	dim, _ := H.Dims()
	cmatrix.HermEigensystem(H, work, evals, evecs)
	sum := complex(0.0, 0.0)
	for alpha := 0; alpha < dim; alpha++ {
		// Coefficients psi^*_{alpha} psi_{alpha}.
		// Indices reversed relative to matrix since Eigensystem returns
		// a slice of eigenvectors (i.e. eigenvectors in rows instead of columns).
		// Shifted by 1 since the slice is zero-indexed.
		//left := cmplx.Conj(evecs[alpha][indexL-1])
		//right := evecs[alpha][indexR-1]
		left := cmplx.Conj(evecs.At(indexL-1, alpha))
		right := evecs.At(indexR-1, alpha)
		// Fermi-Dirac occupation.
		// Mu is included in H, so not included here.
		occ := env.Fermi(evals.At(alpha))
		// alpha'th eigenvector contribution to EV.
		sum += left * right * complex(occ, 0.0)
	}
	return sum
}
