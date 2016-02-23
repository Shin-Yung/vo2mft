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
		//fmt.Println(env.Beta*env.H_Ion(S, Ds), math.Exp(-env.Beta*env.H_Ion(S, Ds)), val)
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
	Kbe := 4.0 * env.Kb()
	Kcxxe, Kczze, Kcxz := 2.0*env.Kcxx(), 2.0*env.Kczz(), env.Kcxz()
	Dco := Ds.Dco(env)

	S01_part := (Bzz+Kbe*env.W11+Kczze*env.W01+Kcxz*env.W02)*S01*S01 - (4.0*Jb*env.M11+2.0*Jc*env.M01+2.0*Dco)*S01
	S11_part := (Bxy+Kbe*env.W01+Kcxxe*env.W11+Kcxz*env.W12)*S11*S11 - 4.0*Jb*env.M01*S11
	S02_part := (Bxy+Kbe*env.W12+Kcxxe*env.W02+Kcxz*env.W01)*S02*S02 - 4.0*Jb*env.M12*S02
	S12_part := (Bzz+Kbe*env.W02+Kczze*env.W12+Kcxz*env.W11)*S12*S12 - (4.0*Jb*env.M02+2.0*Jc*env.M12+2.0*Dco)*S12
	S02_S01_part := Bxz * S02 * S02 * S01 * S01
	S11_S12_part := Bxz * S11 * S11 * S12 * S12
	return S01_part + S11_part + S02_part + S12_part + S02_S01_part + S11_S12_part
}

// Constant part of the ionic Hamiltonian (no S dependence).
func (env *Environment) EConst_Ion() float64 {
	dimer_quad := env.Jc() * (math.Pow(env.M01, 2.0) + math.Pow(env.M12, 2.0))
	dimer_quart := -env.Kcxx()*(env.W02*env.W02+env.W11*env.W11) + env.Kczz()*(env.W01*env.W01+env.W12*env.W12) + env.Kcxz()*(env.W02*env.W01+env.W11*env.W12)
	cb := 4.0 * env.Jb() * (env.M01*env.M11 + env.M02*env.M12)
	ccbb := -4.0 * env.Kb() * (env.W01*env.W11 + env.W02*env.W12)
	return dimer_quad + dimer_quart + cb + ccbb
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
