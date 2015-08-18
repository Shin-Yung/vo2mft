package vo2solve

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

func MWSystem(env *Environment, Ds *HoppingEV) (solve.DiffSystem, []float64) {
	variables := []string{"M", "W"}
	diffM := AbsErrorM(env, Ds, variables)
	diffW := AbsErrorW(env, Ds, variables)
	system := solve.Combine([]solve.Diffable{diffM, diffW})
	start := []float64{env.M, env.W}
	return system, start
}

func MWSolve(env *Environment, Ds *HoppingEV, epsAbs, epsRel float64) (vec.Vector, error) {
	system, start := MWSystem(env, Ds)
	solution, err := solve.MultiDim(system, start, epsAbs, epsRel)
	if err != nil {
		return nil, err
	}
	return solution, nil
}

func MuSystem(env *Environment) (solve.DiffSystem, []float64) {
	variables := []string{"Mu"}
	diffMu := AbsErrorMu(env, variables)
	system := solve.Combine([]solve.Diffable{diffMu})
	start := []float64{env.Mu}
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

func MWMuSolve_Iterative(env *Environment, Ds *HoppingEV, epsAbs, epsRel float64) (vec.Vector, error) {
	Mu_system, Mu_start := MuSystem(env)
	MW_system, MW_start := MWSystem(env, Ds)
	stages := []solve.DiffSystem{Mu_system, MW_system}
	start := []vec.Vector{Mu_start, MW_start}
	accept := func(x []vec.Vector) {
		env.Mu = x[0][0]
		env.M = x[1][0]
		env.W = x[1][1]
	}
	epsAbsV := []float64{epsAbs, epsAbs}
	epsRelV := []float64{epsRel, epsRel}
	solutions, err := solve.Iterative(stages, start, epsAbsV, epsRelV, accept)
	if err != nil {
		return nil, err
	}
	result := vec.Vector{solutions[0][0], solutions[1][0], solutions[1][1]}
	return result, nil
}
