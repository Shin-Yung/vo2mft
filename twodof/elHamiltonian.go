package twodof

import (
	"math"
	"math/cmplx"
)
import (
	"github.com/tflovorn/cmatrix"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Calculate 4x4 electronic Hamiltonian.
// k is in the Cartesian basis, with each component scaled by the corresponding
// lattice constant; i.e. k = (a kx, a ky, c kz) and a kx, a ky, c kz range
// over [-pi, pi) and periodic copies of this interval.
func ElHamiltonian(env *Environment, k vec.Vector) cmatrix.CMatrix {
	KQ := vec.Vector{0.0, math.Pi, math.Pi}
	k.Add(&KQ) // now KQ = k + Q

	EpsAE := EpsilonAE(env, k)
	EpsBE := EpsilonBE(env, k)
	EpsBE_KQ := EpsilonBE(env, KQ)
	EpsAO := EpsilonAO(env, k)
	mu := complex(env.Mu, 0.0)
	ident_part := -0.5 * mu
	m01 := complex(env.M01, 0.0)
	m12 := complex(env.M12, 0.0)

	H := cmatrix.InitSliceCMatrix(4, 4)

	H[0][0] = EpsAE + ident_part
	H[1][0] = -m01 * EpsAO
	H[2][0] = cmplx.Conj(EpsBE)
	H[3][0] = 0.0

	H[0][1] = m01 * EpsAO
	H[1][1] = -EpsAE + ident_part
	H[2][1] = 0.0
	H[3][1] = cmplx.Conj(EpsBE_KQ)

	H[0][2] = EpsBE
	H[1][2] = 0.0
	H[2][2] = EpsAE + ident_part
	H[3][2] = -m12 * EpsAO

	H[0][3] = 0.0
	H[1][3] = EpsBE_KQ
	H[2][3] = m12 * EpsAO
	H[3][3] = -EpsAE + ident_part

	return H
}

// Cubic axes, even symmetry (k, p; k, p)
func EpsilonAE(env *Environment, k vec.Vector) complex128 {
	rp := -env.Tce * math.Cos(k[2])
	return complex(rp, 0.0)
}

// Body diagonal, even symmetry (k, p; k, pbar)
func EpsilonBE(env *Environment, k vec.Vector) complex128 {
	e1 := cmplx.Exp(complex(0.0, 0.0))
	e2 := cmplx.Exp(complex(0.0, -k[0]))
	e3 := cmplx.Exp(complex(0.0, -k[1]))
	e4 := cmplx.Exp(complex(0.0, -k[2]))
	return -complex(env.Tbe, 0.0) * (e1 + e2 + e3 + e4)
}

// Cubic axes, odd symmetry (k, p; k+Q, p)
func EpsilonAO(env *Environment, k vec.Vector) complex128 {
	ip := -2.0 * env.Tco * math.Sin(k[2])
	return complex(0.0, ip)
}
