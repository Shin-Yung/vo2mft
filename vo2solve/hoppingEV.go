package vo2solve

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
	// Value of M for which the contained hopping e.v.'s have been calculated.
	m_dae, m_dce, m_dbe, m_dao, m_dco, m_dbo float64
	// Value of Mu for which the contained hopping e.v.'s have been calculated.
	mu_dae, mu_dce, mu_dbe, mu_dao, mu_dco, mu_dbo float64
	// If hopping e.v.'s have not been calculated yet, init = false.
	init_dae, init_dce, init_dbe, init_dao, init_dco, init_dbo bool
	// Hopping e.v.'s for even symmetry (pre-calculated).
	dae, dce, dbe float64
	// Hopping e.v.'s for odd symmetry (pre-calculated).
	dao, dco, dbo float64
}

func (Ds *HoppingEV) Dae(env *Environment) float64 {
	if Ds.init_dae && (env.M == Ds.m_dae) && (env.Mu == Ds.mu_dae) {
		return Ds.dae
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_K0(env, k)
		// make sure that ev is pure real
		//if math.Abs(imag(ev)) > zero_threshold {
		//	panic("Expected pure real value for <c_{k,0}^{\\dagger} c_{k,0}>, got finite imaginary part")
		//}
		return 4.0 * math.Cos(k[0]) * real(ev)
	}
	dae := 0.5 * bzone.Avg(env.BZPointsPerDim, 3, inner)

	inner_im := func(k vec.Vector) float64 {
		ev := GetEV_K0_K0(env, k)
		return -4.0 * math.Cos(k[0]) * imag(ev)
	}
	dae_im := bzone.Avg(env.BZPointsPerDim, 3, inner_im)
	if math.Abs(dae_im) > zero_threshold {
		panic("Expected real value for Dae, got finite imaginary part.")
	}

	Ds.init_dae = true
	Ds.m_dae = env.M
	Ds.mu_dae = env.Mu
	Ds.dae = dae
	return dae
}

func (Ds *HoppingEV) Dce(env *Environment) float64 {
	if Ds.init_dce && (env.M == Ds.m_dce) && (env.Mu == Ds.mu_dce) {
		return Ds.dce
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_K0(env, k)
		// make sure that ev is pure real
		//if math.Abs(imag(ev)) > zero_threshold {
		//	panic("Expected pure real value for <c_{k,0}^{\\dagger} c_{k,0}>, got finite imaginary part")
		//}
		return 4.0 * math.Cos(k[2]) * real(ev)
	}
	dce := 0.5 * bzone.Avg(env.BZPointsPerDim, 3, inner)

	inner_im := func(k vec.Vector) float64 {
		ev := GetEV_K0_K0(env, k)
		return -4.0 * math.Cos(k[2]) * imag(ev)
	}
	dce_im := bzone.Avg(env.BZPointsPerDim, 3, inner_im)
	if math.Abs(dce_im) > zero_threshold {
		panic("Expected real value for Dce, got finite imaginary part.")
	}

	Ds.init_dce = true
	Ds.m_dce = env.M
	Ds.mu_dce = env.Mu
	Ds.dce = dce
	return dce
}

func (Ds *HoppingEV) Dbe(env *Environment) float64 {
	if Ds.init_dbe && (env.M == Ds.m_dbe) && (env.Mu == Ds.mu_dbe) {
		return Ds.dbe
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_K1(env, k)
		return 2.0 * real(ev+cmplx.Conj(ev))
	}
	dbe := 0.5 * bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dbe = true
	Ds.m_dbe = env.M
	Ds.mu_dbe = env.Mu
	Ds.dbe = dbe
	return dbe
}

func (Ds *HoppingEV) Dao(env *Environment) float64 {
	if Ds.init_dao && (env.M == Ds.m_dao) && (env.Mu == Ds.mu_dao) {
		return Ds.dao
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K0(env, k)
		//println("Dao", ev)
		// make sure that ev is pure imaginary
		//if math.Abs(real(ev)) > zero_threshold {
		//	panic("Expected pure imaginary value for <c_{k+Q,0}^{\\dagger} c_{k,0}>, got finite real part")
		//}
		// 2i * ev = -2 * imag(ev)
		return -2.0 * math.Sin(k[0]) * imag(ev)
	}
	dao := 0.5 * bzone.Avg(env.BZPointsPerDim, 3, inner)

	inner_re := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K0(env, k)
		return 2.0 * math.Sin(k[0]) * real(ev)
	}
	dao_re := bzone.Avg(env.BZPointsPerDim, 3, inner_re)
	if math.Abs(dao_re) > zero_threshold {
		panic("Expected real value for Dao, got finite imaginary part.")
	}

	Ds.init_dao = true
	Ds.m_dao = env.M
	Ds.mu_dao = env.Mu
	Ds.dao = dao
	return dao
}

func (Ds *HoppingEV) Dco(env *Environment) float64 {
	if Ds.init_dco && (env.M == Ds.m_dco) && (env.Mu == Ds.mu_dco) {
		return Ds.dco
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K0(env, k)
		// make sure that ev is pure imaginary
		//if math.Abs(real(ev)) > zero_threshold {
		//	panic("Expected pure imaginary value for <c_{k+Q,0}^{\\dagger} c_{k,0}>, got finite real part")
		//}
		// 2i * ev = -2 * imag(ev)
		return -2.0 * math.Sin(k[2]) * imag(ev)
	}
	dco := 0.5 * bzone.Avg(env.BZPointsPerDim, 3, inner)

	inner_re := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K0(env, k)
		return 2.0 * math.Sin(k[2]) * real(ev)
	}
	dco_re := bzone.Avg(env.BZPointsPerDim, 3, inner_re)
	if math.Abs(dco_re) > zero_threshold {
		panic("Expected real value for Dco, got finite imaginary part.")
	}

	Ds.init_dco = true
	Ds.m_dco = env.M
	Ds.mu_dco = env.Mu
	Ds.dco = dco
	return dco
}

func (Ds *HoppingEV) Dbo(env *Environment) float64 {
	if Ds.init_dbo && (env.M == Ds.m_dbo) && (env.Mu == Ds.mu_dbo) {
		return Ds.dbo
	}
	if !env.FiniteHoppings() {
		return 0.0
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K1(env, k)
		//println("Dbo", ev)
		return real(ev + cmplx.Conj(ev))
	}
	dbo := 0.5 * bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dbo = true
	Ds.m_dbo = env.M
	Ds.mu_dbo = env.Mu
	Ds.dbo = dbo
	return dbo
}

// Convert to string by marshalling to JSON.
// Leave out internal cache data.
func (Ds *HoppingEV) StringEnv(env *Environment) string {
	marshalled := Ds.MarshalEnv(env)
	return marshalled
}

func (Ds *HoppingEV) MarshalEnv(env *Environment) string {
	repr := make(map[string]float64)
	repr["Dae"] = Ds.Dae(env)
	repr["Dce"] = Ds.Dce(env)
	repr["Dbe"] = Ds.Dbe(env)
	repr["Dao"] = Ds.Dao(env)
	repr["Dco"] = Ds.Dco(env)
	repr["Dbo"] = Ds.Dbo(env)
	marshalled, err := json.Marshal(repr)
	if err != nil {
		panic(err)
	}
	return string(marshalled)
}

// Evaluate <c^{\dagger}_{k,0} c_{k,0}>.
func GetEV_K0_K0(env *Environment, k vec.Vector) complex128 {
	return evalEV(env, k, 1, 1)
}

// Evaluate <c^{\dagger}_{k+Q,0} c_{k,0}>.
func GetEV_KQ0_K0(env *Environment, k vec.Vector) complex128 {
	return evalEV(env, k, 2, 1)
}

// Evaluate <c^{\dagger}_{k,0} c_{k,1}>.
func GetEV_K0_K1(env *Environment, k vec.Vector) complex128 {
	return evalEV(env, k, 1, 3)
}

// Evaluate <c^{\dagger}_{k+Q,0} c_{k,1}>.
func GetEV_KQ0_K1(env *Environment, k vec.Vector) complex128 {
	return evalEV(env, k, 2, 3)
}

// Evaluate <c^{\dagger}_{indexL} c_{indexR}> where the index values have
// the following correspondence:
// 	1 <--> k, 0 ; 2 <--> k+Q, 0 ; 3 <--> k, 1 ; 4 <--> k+Q, 1
func evalEV(env *Environment, k vec.Vector, indexL, indexR int) complex128 {
	H := ElHamiltonian(env, k)
	dim, _ := H.Dims()
	evals, evecs := cmatrix.Eigensystem(H)
	sum := complex(0.0, 0.0)
	for alpha := 0; alpha < dim; alpha++ {
		// Coefficients psi^*_{alpha} psi_{alpha}.
		// Indices reversed relative to matrix since Eigensystem returns
		// a slice of eigenvectors (i.e. eigenvectors in rows instead of columns).
		// Shifted by 1 since the slice is zero-indexed.
		left := cmplx.Conj(evecs[alpha][indexL-1])
		right := evecs[alpha][indexR-1]
		// Fermi-Dirac occupation.
		// Mu is included in H, so not included here.
		occ := env.Fermi(evals[alpha])
		// alpha'th eigenvector contribution to EV.
		sum += left * right * complex(occ, 0.0)
	}
	return sum
}
