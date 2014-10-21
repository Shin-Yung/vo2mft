package vo2mft

import (
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

func MWMuSystem(env *Environment, Ds *HoppingEV) (solve.DiffSystem, []float64) {
	variables := []string{"M", "W", "Mu"}
	diffM := AbsErrorM(env, Ds, variables)
	diffW := AbsErrorW(env, Ds, variables)
	diffMu := AbsErrorMu(env, variables)
	system := solve.Combine([]solve.Diffable{diffM, diffW, diffMu})
	start := []float64{env.M, env.W, env.Mu}
	return system, start
}

func MWMuSolve(env *Environment, Ds *HoppingEV, epsAbs, epsRel float64) (vec.Vector, error) {
	system, start := MWMuSystem(env, Ds)
	solution, err := solve.MultiDim(system, start, epsAbs, epsRel)
	if err != nil {
		return nil, err
	}
	return solution, nil
}
