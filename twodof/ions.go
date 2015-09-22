package twodof

import (
	"math"
)

var cached_all_S = [][]int{}

func (env *Environment) Z1(Ds *HoppingEV) float64 {
	all_S := all_S_configs()
	val := 0.0
	for _, S := range all_S {
		val += math.Exp(-env.Beta * env.H_Ion(S, Ds))
	}
	return val
}

func (env *Environment) Mpa(p, alpha int, Ds *HoppingEV) float64 {
	S_index := p + 2*(alpha-1)
	all_S := all_S_configs()
	val := 0.0
	for _, S := range all_S {
		S_part := float64(S[S_index])
		val += S_part * math.Exp(-env.Beta*env.H_Ion(S, Ds))
	}
	return val / env.Z1(Ds)
}

func (env *Environment) Wpa(p, alpha int, Ds *HoppingEV) float64 {
	S_index := p + 2*(alpha-1)
	all_S := all_S_configs()
	val := 0.0
	for _, S := range all_S {
		S_part := float64(S[S_index] * S[S_index])
		val += S_part * math.Exp(-env.Beta*env.H_Ion(S, Ds))
	}
	return val / env.Z1(Ds)
}

// Single-site ionic Hamiltonian (local and ion-ion parts).
// S = [S01, S11, S02, S12].
func (env *Environment) H_Ion(S []int, Ds *HoppingEV) float64 {
	S01, S11, S02, S12 := float64(S[0]), float64(S[1]), float64(S[2]), float64(S[3])
	Bxy, Bzz, Bxz, Jb, Jc := env.Bxy(), env.Bzz(), env.Bxz(), env.Jb(), env.Jc()
	Dco := Ds.Dco(env)

	S01_part := Bzz*S01*S01 + Bxz*S01*S01*env.W02 - (4.0*Jb*env.M11+2.0*Jc*env.M01+2.0*Dco)*S01
	S11_part := Bxy*S11*S11 + Bxz*S11*S11*env.W12 - 4.0*Jb*env.M01*S11
	S02_part := Bxy*S02*S02 + Bxz*S02*S02*env.W01 - 4.0*Jb*env.M12*S02
	S12_part := Bzz*S12*S12 + Bxz*S12*S12*env.W11 - (4.0*Jb*env.M02+2.0*Jc*env.M12+2.0*Dco)*S12
	return S01_part + S11_part + S02_part + S12_part
}

// Constant part of the ionic Hamiltonian (no S dependence).
func (env *Environment) EConst_Ion() float64 {
	dimer := env.Jc() * (math.Pow(env.M01, 2.0) + math.Pow(env.M12, 2.0))
	cb := 4.0 * env.Jb() * (env.M01*env.M11 + env.M02*env.M12)
	onsite := -env.Bxz() * (env.W01*env.W02 + env.W11*env.W12)
	return dimer + cb + onsite
}

// Constant part of the electron-ion Hamiltonian.
func (env *Environment) EConst_IonEl(Ds *HoppingEV) float64 {
	return 2.0 * env.Tco * (env.M01 + env.M12) * Ds.Dco(env)
}

func all_S_configs() [][]int {
	if len(cached_all_S) > 0 {
		return cached_all_S
	}

	all_S := [][]int{}
	for i := 0; i < 4; i++ {
		add_S_elems(&all_S)
	}
	cached_all_S = all_S
	return all_S
}

func add_S_elems(all_S *[][]int) {
	if len(*all_S) == 0 {
		*all_S = [][]int{[]int{-1}, []int{0}, []int{1}}
	} else {
		new_all_S := make([][]int, 3*len(*all_S))
		for i, v := range *all_S {
			base := []int{-1, 0, 1}
			for j, S := range base {
				new_v := make([]int, len(v))
				copy(new_v, v)
				new_v = append(new_v, S)
				new_all_S[3*i+j] = new_v
			}
		}
		*all_S = new_all_S
	}
}
