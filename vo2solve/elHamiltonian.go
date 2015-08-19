package vo2solve

import (
	"math"
	"math/cmplx"
)
import (
	"github.com/tflovorn/cmatrix"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Calculate 4x4 electronic Hamiltonian.
// Assumes that k is scaled such that kx*a, ky*a, kz*c --> kx, ky, kz.
func ElHamiltonian(env *Environment, k vec.Vector) cmatrix.CMatrix {
	KQ := vec.Vector{math.Pi, math.Pi, math.Pi}
	k.Add(&KQ) // now KQ = k + Q

	EpsAE := EpsilonAE(env, k)
	EpsBE := EpsilonBE(env, k)
	EpsBE_KQ := EpsilonBE(env, KQ)
	EpsAO := EpsilonAO(env, k)
	EpsBO := EpsilonBO(env, k)
	ikd := complex(0.0, k[0]/2.0+k[1]/2.0+k[2]/2.0)
	mu := complex(env.Mu, 0.0)

	H := cmatrix.InitSliceCMatrix(4, 4)
	/*
		//H[0][0] = 0.5*EpsAE - 0.5*mu
		H[0][0] = 0.5*EpsAE
		H[1][0] = EpsAO
		H[2][0] = 0.5*EpsBE*cmplx.Exp(-ikd)
		H[3][0] = -0.5*cmplx.Conj(EpsBO)*cmplx.Exp(-ikd)

		H[0][1] = -EpsAO
		//H[1][1] = -0.5*EpsAE - 0.5*mu
		H[1][1] = -0.5*EpsAE
		H[2][1] = 0.5*cmplx.Conj(EpsBO)*cmplx.Exp(-ikd)
		H[3][1] = complex(0.0, 0.5)*EpsBE_KQ*cmplx.Exp(-ikd)

		H[0][2] = 0.5*EpsBE*cmplx.Exp(ikd)
		H[1][2] = 0.5*EpsBO*cmplx.Exp(ikd)
		//H[2][2] = 0.5*EpsAE - 0.5*mu
		H[2][2] = 0.5*EpsAE
		H[3][2] = EpsAO

		H[0][3] = -0.5*EpsBO*cmplx.Exp(ikd)
		H[1][3] = complex(0.0, -0.5)*EpsBE_KQ*cmplx.Exp(ikd)
		H[2][3] = -EpsAO
		//H[3][3] = -0.5*EpsAE - 0.5*mu
		H[3][3] = -0.5*EpsAE
	*/
	H[0][0] = EpsAE - mu
	H[1][0] = 2.0 * EpsAO
	H[2][0] = EpsBE * cmplx.Exp(-ikd)
	H[3][0] = -cmplx.Conj(EpsBO) * cmplx.Exp(-ikd)

	H[0][1] = -2.0 * EpsAO
	H[1][1] = -EpsAE - mu
	H[2][1] = cmplx.Conj(EpsBO) * cmplx.Exp(-ikd)
	H[3][1] = complex(0.0, 1.0) * EpsBE_KQ * cmplx.Exp(-ikd)

	H[0][2] = EpsBE * cmplx.Exp(ikd)
	H[1][2] = EpsBO * cmplx.Exp(ikd)
	H[2][2] = EpsAE - mu
	H[3][2] = 2.0 * EpsAO

	H[0][3] = -EpsBO * cmplx.Exp(ikd)
	H[1][3] = complex(0.0, -1.0) * EpsBE_KQ * cmplx.Exp(ikd)
	H[2][3] = -2.0 * EpsAO
	H[3][3] = -EpsAE - mu

	return H
}

// Cubic axes, even symmetry (k, p; k, p)
func EpsilonAE(env *Environment, k vec.Vector) complex128 {
	rp := -2.0 * (env.Tae*(math.Cos(k[0])+math.Cos(k[1])) + env.Tce*math.Cos(k[2]))
	return complex(rp, 0.0)
}

// Body diagonal, even symmetry (k, p; k, pbar)
func EpsilonBE(env *Environment, k vec.Vector) complex128 {
	rp := -8.0 * env.Tbe * math.Cos(k[0]/2.0) * math.Cos(k[1]/2.0) * math.Cos(k[2]/2.0)
	return complex(rp, 0.0)
}

// Cubic axes, odd symmetry (k, p; k+Q, p)
func EpsilonAO(env *Environment, k vec.Vector) complex128 {
	ip := -2.0 * env.M * (env.Tao*(math.Sin(k[0])+math.Sin(k[1])) + env.Tco*math.Sin(k[2]))
	return complex(0.0, ip)
}

// Body diagonal, odd symmetry (k, p; k+Q, pbar)
func EpsilonBO(env *Environment, k vec.Vector) complex128 {
	rp := -8.0 * env.M * env.Tbo * math.Cos(k[0]/2.0) * math.Cos(k[1]/2.0) * math.Cos(k[2]/2.0)
	ip := 8.0 * env.M * env.Tbo * math.Sin(k[0]/2.0) * math.Sin(k[1]/2.0) * math.Sin(k[2]/2.0)
	return complex(rp, ip)
}
