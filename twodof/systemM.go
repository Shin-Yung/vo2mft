package twodof

import (
//"fmt"
)
import (
	"github.com/tflovorn/scExplorer/solve"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Return the absolute error and gradient of the M_{p, alpha} equation w.r.t.
// the given variables.
func AbsErrorM(env *Environment, variables []string, p, alpha int) solve.Diffable {
	M_name := []string{"M01", "M11", "M02", "M12"}
	name_index := p + 2*(alpha-1)

	F := func(v vec.Vector) (float64, error) {
		env.Set(v, variables)
		M_env := env.GetFloat(M_name[name_index])
		M_eq := env.Mpa(p, alpha)
		//fmt.Printf("p = %d, alpha = %d, M_env = %f, M_eq = %f\n", p, alpha, M_env, M_eq)
		return M_env - M_eq, nil
	}
	h := 1e-6
	epsabs := 1e-4
	return solve.SimpleDiffable(F, len(variables), h, epsabs)
}
