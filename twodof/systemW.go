package twodof

import (
//"fmt"
)
import (
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Return the absolute error and gradient of the W_{p, alpha} equation w.r.t.
// the given variables.
func AbsErrorW(env *Environment, Ds *HoppingEV, variables []string, p, alpha int) solve.Diffable {
	W_name := []string{"W01", "W11", "W02", "W12"}
	name_index := p + 2*(alpha-1)

	F := func(v vec.Vector) (float64, error) {
		env.Set(v, variables)
		W_env := env.GetFloat(W_name[name_index])
		W_eq := env.Wpa(p, alpha, Ds)
		//fmt.Printf("p = %d, alpha = %d, M_env = %f, M_eq = %f\n", p, alpha, M_env, M_eq)
		return W_env - W_eq, nil
	}
	h := 1e-6
	epsabs := 1e-4
	return solve.SimpleDiffable(F, len(variables), h, epsabs)
}
