package vo2mft

import (
	"math"
	"math/cmplx"
)
import (
	"github.com/tflovorn/scExplorer/bzone"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Expect some e.v.'s to be pure real or imaginary - panic if imag/real part
// is greater than zero_threshold.
const zero_threshold = 1e-12

type HoppingEV struct {
	// Value of M for which the contained hopping e.v.'s have been calculated.
	m_dae, m_dce, m_dbe, m_dao, m_dco, m_dbo float64
	// If hopping e.v.'s have not been calculated yet, init = false.
	init_dae, init_dce, init_dbe, init_dao, init_dco, init_dbo bool
	// Hopping e.v.'s for even symmetry (pre-calculated).
	dae, dce, dbe float64
	// Hopping e.v.'s for odd symmetry (pre-calculated).
	dao, dco, dbo float64
}

func (Ds *HoppingEV) Dae(env *Environment) float64 {
	if Ds.init_dae && (env.M == Ds.m_dae) {
		return Ds.dae
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_K0(env, k)
		// make sure that ev is pure real
		if math.Abs(imag(ev)) > zero_threshold {
			panic("Expected pure real value for <c_{k,0}^{\\dagger} c_{k,0}>, got finite imaginary part")
		}
		return 4.0*math.Cos(k[0])*real(ev)
	}
	dae := bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dae = true
	Ds.m_dae = env.M
	Ds.dae = dae
	return dae
}

func (Ds *HoppingEV) Dce(env *Environment) float64 {
	if Ds.init_dce && (env.M == Ds.m_dce) {
		return Ds.dce
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_K0(env, k)
		// make sure that ev is pure real
		if math.Abs(imag(ev)) > zero_threshold {
			panic("Expected pure real value for <c_{k,0}^{\\dagger} c_{k,0}>, got finite imaginary part")
		}
		return 4.0*math.Cos(k[2])*real(ev)
	}
	dce := bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dce = true
	Ds.m_dce = env.M
	Ds.dce = dce
	return dce
}

func (Ds *HoppingEV) Dbe(env *Environment) float64 {
	if Ds.init_dbe && (env.M == Ds.m_dbe) {
		return Ds.dbe
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_K0_K1(env, k)
		return 2.0*real(ev + cmplx.Conj(ev))
	}
	dbe := bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dbe = true
	Ds.m_dbe = env.M
	Ds.dbe = dbe
	return dbe
}

func (Ds *HoppingEV) Dao(env *Environment) float64 {
	if Ds.init_dao && (env.M == Ds.m_dao) {
		return Ds.dao
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K0(env, k)
		// make sure that ev is pure imaginary
		if math.Abs(real(ev)) > zero_threshold {
			panic("Expected pure imaginary value for <c_{k+Q,0}^{\\dagger} c_{k,0}>, got finite real part")
		}
		// 2i * ev = -2 * imag(ev)
		return -2.0*math.Sin(k[0])*imag(ev)
	}
	dao := bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dao = true
	Ds.m_dao = env.M
	Ds.dao = dao
	return dao
}

func (Ds *HoppingEV) Dco(env *Environment) float64 {
	if Ds.init_dco && (env.M == Ds.m_dco) {
		return Ds.dco
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K0(env, k)
		// make sure that ev is pure imaginary
		if math.Abs(real(ev)) > zero_threshold {
			panic("Expected pure imaginary value for <c_{k+Q,0}^{\\dagger} c_{k,0}>, got finite real part")
		}
		// 2i * ev = -2 * imag(ev)
		return -2.0*math.Sin(k[2])*imag(ev)
	}
	dco := bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dco = true
	Ds.m_dco = env.M
	Ds.dco = dco
	return dco
}

func (Ds *HoppingEV) Dbo(env *Environment) float64 {
	if Ds.init_dbo && (env.M == Ds.m_dbo) {
		return Ds.dbo
	}

	inner := func(k vec.Vector) float64 {
		ev := GetEV_KQ0_K1(env, k)
		return real(ev + cmplx.Conj(ev))
	}
	dbo := bzone.Avg(env.BZPointsPerDim, 3, inner)

	Ds.init_dbo = true
	Ds.m_dbo = env.M
	Ds.dbo = dbo
	return dbo
}

func GetEV_K0_K0(emv *Environment, k vec.Vector) complex128 {
	return complex(0.0, 0.0)
}

func GetEV_KQ0_K0(emv *Environment, k vec.Vector) complex128 {
	return complex(0.0, 0.0)
}

func GetEV_K0_K1(emv *Environment, k vec.Vector) complex128 {
	return complex(0.0, 0.0)
}

func GetEV_KQ0_K1(emv *Environment, k vec.Vector) complex128 {
	return complex(0.0, 0.0)
}
