package vo2mft

// Contains parameters necessary to characterize electronic and ionic systems.
// The order parameters M and W must be determined self-consistently.
type Environment struct {
	// Hopping parameters, even symmetry (a, c, diagonal axes).
	Tae, Tce, Tbe float64
	// Hopping parameters, odd symmetry (a, c, diagonal axes).
	Tao, Tco, Tbo float64
	// Order parameter <S>.
	M float64
	// Order parameter <S^2>.
	W float64
	// Chemical potential.
	Mu float64
	// Inverse temperature, 1 / (k_B * T).
	Beta float64
	// One-spin term for BEG model: coefficient for (S_i)^2.
	B float64
	// Exchange parameters for BEG model: coefficients to S_i dot S_j.
	// Jb is excluded since it does not contribute to results.
	Ja, Jc float64
	// Biquadratic exchange parameters for BEG model: coefficients to (S_i)^2 * (S_j)^2.
	Ka, Kc, Kb float64
	// On-site energies in M and R phases.
	EpsilonM, EpsilonR float64
}

func (env *Environment) DeltaS() float64 {
	return env.B + env.EpsilonM - env.EpsilonR
}

func (env *Environment) QK() float64 {
	return 4.0*env.Ka + 2.0*env.Kc + 8.0*env.Kb
}

func (env *Environment) QJ(Ds *HoppingEV) float64 {
	Dao, Dco := Ds.Dao(env), Ds.Dco(env)
	return 4.0*(env.Ja + env.Tao*Dao) + 2.0*(env.Jc + env.Tco*Dco)
}

func (env *Environment) Qele(Ds *HoppingEV) float64 {
	Dae, Dce, Dbe := Ds.Dae(env), Ds.Dce(env), Ds.Dbe(env)
	return 4.0*env.Tae*Dae + 2.0*env.Tce*Dce + 8.0*env.Tbe*Dbe
}
