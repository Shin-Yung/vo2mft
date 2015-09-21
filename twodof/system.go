package twodof

import (
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

func MSystem(env *Environment, Ds *HoppingEV, m01_0, m11_0, m02_0, m12_0 bool) (solve.DiffSystem, []float64) {
	variables, start := []string{}, []float64{}
	diff_list := []solve.Diffable{}

	if !m01_0 {
		variables = append(variables, "M01")
		start = append(start, env.M01)
	}
	if !m11_0 {
		variables = append(variables, "M11")
		start = append(start, env.M11)
	}
	if !m02_0 {
		variables = append(variables, "M02")
		start = append(start, env.M02)
	}
	if !m12_0 {
		variables = append(variables, "M12")
		start = append(start, env.M12)
	}

	if !m01_0 {
		diffM01 := AbsErrorM(env, Ds, variables, 0, 1)
		diff_list = append(diff_list, diffM01)
	}
	if !m11_0 {
		diffM11 := AbsErrorM(env, Ds, variables, 1, 1)
		diff_list = append(diff_list, diffM11)
	}
	if !m02_0 {
		diffM02 := AbsErrorM(env, Ds, variables, 0, 2)
		diff_list = append(diff_list, diffM02)
	}
	if !m12_0 {
		diffM12 := AbsErrorM(env, Ds, variables, 1, 2)
		diff_list = append(diff_list, diffM12)
	}

	system := solve.Combine(diff_list)
	return system, start
}

func MMuSystem(env *Environment, Ds *HoppingEV, m01_0, m11_0, m02_0, m12_0 bool) (solve.DiffSystem, []float64) {
	variables, start := []string{}, []float64{}
	diff_list := []solve.Diffable{}

	if !m01_0 {
		variables = append(variables, "M01")
		start = append(start, env.M01)
	}
	if !m11_0 {
		variables = append(variables, "M11")
		start = append(start, env.M11)
	}
	if !m02_0 {
		variables = append(variables, "M02")
		start = append(start, env.M02)
	}
	if !m12_0 {
		variables = append(variables, "M12")
		start = append(start, env.M12)
	}
	variables = append(variables, "Mu")
	start = append(start, env.Mu)

	if !m01_0 {
		diffM01 := AbsErrorM(env, Ds, variables, 0, 1)
		diff_list = append(diff_list, diffM01)
	}
	if !m11_0 {
		diffM11 := AbsErrorM(env, Ds, variables, 1, 1)
		diff_list = append(diff_list, diffM11)
	}
	if !m02_0 {
		diffM02 := AbsErrorM(env, Ds, variables, 0, 2)
		diff_list = append(diff_list, diffM02)
	}
	if !m12_0 {
		diffM12 := AbsErrorM(env, Ds, variables, 1, 2)
		diff_list = append(diff_list, diffM12)
	}
	diffMu := AbsErrorMu(env, variables)
	diff_list = append(diff_list, diffMu)

	system := solve.Combine(diff_list)
	return system, start
}

func MSolve(env *Environment, Ds *HoppingEV, epsAbs, epsRel float64, m01_0, m11_0, m02_0, m12_0 bool) (vec.Vector, error) {
	if m01_0 && m11_0 && m02_0 && m12_0 {
		return []float64{}, nil
	}
	system, start := MSystem(env, Ds, m01_0, m11_0, m02_0, m12_0)
	solution, err := solve.MultiDim(system, start, epsAbs, epsRel)
	if err != nil {
		return nil, err
	}
	return solution, nil
}

func MMuSolve(env *Environment, Ds *HoppingEV, epsAbs, epsRel float64, m01_0, m11_0, m02_0, m12_0 bool) (vec.Vector, error) {
	system, start := MMuSystem(env, Ds, m01_0, m11_0, m02_0, m12_0)
	solution, err := solve.MultiDim(system, start, epsAbs, epsRel)
	if err != nil {
		return nil, err
	}
	return solution, nil
}
