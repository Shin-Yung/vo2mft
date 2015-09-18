package twodof

import (
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

func MSystem(env *Environment) (solve.DiffSystem, []float64) {
	variables := []string{"M01", "M11", "M02", "M12"}
	diffM01 := AbsErrorM(env, variables, 0, 1)
	diffM11 := AbsErrorM(env, variables, 1, 1)
	diffM02 := AbsErrorM(env, variables, 0, 2)
	diffM12 := AbsErrorM(env, variables, 1, 2)
	system := solve.Combine([]solve.Diffable{diffM01, diffM11, diffM02, diffM12})
	start := []float64{env.M01, env.M11, env.M02, env.M12}
	return system, start
}

func MSolve(env *Environment, epsAbs, epsRel float64) (vec.Vector, error) {
	system, start := MSystem(env)
	solution, err := solve.MultiDim(system, start, epsAbs, epsRel)
	if err != nil {
		return nil, err
	}
	return solution, nil
}
