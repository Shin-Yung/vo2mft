package vo2mft

import (
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

func MWSystem(env *Environment) (solve.DiffSystem, []float64) {
	variables := []string{"M", "W"}
	diffM := AbsErrorM(env, variables)
	diffW := AbsErrorW(env, variables)
	system := solve.Combine([]solve.Diffable{diffM, diffW})
	start := []float64{env.M, env.W}
	return system, start
}

func MWSolve(env *Environment, epsAbs, epsRel float64) (vec.Vector, error) {
	system, start := MWSystem(env)
	solution, err := solve.MultiDim(system, start, epsAbs, epsRel)
	if err != nil {
		return nil, err
	}
	return solution, nil
}
