package vo2mft

type HoppingEV struct {
	// Value of M for which the contained hopping e.v.'s have been calculated.
	M float64
	// If hopping e.v.'s have not been calculated yet, intialized = false.
	initialized bool
	// Hopping e.v.'s for even symmetry (pre-calculated).
	dae, dce, dbe float64
	// Hopping e.v.'s for odd symmetry (pre-calculated).
	dao, dco, dbo float64
}

func (Ds *HoppingEV) Dae(env *Environment) float64 {
	return 0.0
}

func (Ds *HoppingEV) Dce(env *Environment) float64 {
	return 0.0
}

func (Ds *HoppingEV) Dbe(env *Environment) float64 {
	return 0.0
}

func (Ds *HoppingEV) Dao(env *Environment) float64 {
	return 0.0
}

func (Ds *HoppingEV) Dco(env *Environment) float64 {
	return 0.0
}

func (Ds *HoppingEV) Dbo(env *Environment) float64 {
	return 0.0
}
