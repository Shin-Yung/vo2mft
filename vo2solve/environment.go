package vo2solve

import (
	"fmt"
	"math"
	"reflect"
)
import (
	"github.com/tflovorn/scExplorer/serialize"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Contains parameters necessary to characterize electronic and ionic systems.
// The order parameters M and W must be determined self-consistently.
type Environment struct {
	// Size of BZ on one edge (total number of BZ points is this cubed).
	BZPointsPerDim int
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

// Combined biquadratic coefficient (S_i^2 S_j^2).
func (env *Environment) QK() float64 {
	return 4.0*env.Ka + 2.0*env.Kc + 8.0*env.Kb
}

// Combined renormalized 'exchange' coefficient (S_i S_j) favoring dimers.
func (env *Environment) QJ(Ds *HoppingEV) float64 {
	Dao, Dco := Ds.Dao(env), Ds.Dco(env)
	return 4.0*(env.Ja+env.Tao*Dao) + 2.0*(env.Jc+env.Tco*Dco)
}

func (env *Environment) Qele(Ds *HoppingEV) float64 {
	Dae, Dce, Dbe := Ds.Dae(env), Ds.Dce(env), Ds.Dbe(env)
	return 4.0*env.Tae*Dae + 2.0*env.Tce*Dce + 8.0*env.Tbe*Dbe
}

func (env *Environment) Z1(Ds *HoppingEV) float64 {
	exp := math.Exp(-env.Beta * (env.DeltaS() - env.W*env.QK()))
	return 1.0 + 2.0*exp*math.Cosh(env.Beta*env.M*env.QJ(Ds))
}

// Are electronic hopping finite?
// If not, don't need to calculate D's.
func (env *Environment) FiniteHoppings() bool {
	eps := 1e-9
	even := (math.Abs(env.Tae) > eps) || (math.Abs(env.Tce) > eps) || (math.Abs(env.Tbe) > eps)
	odd := (math.Abs(env.Tao) > eps) || (math.Abs(env.Tco) > eps) || (math.Abs(env.Tbo) > eps)
	return even || odd
}

// Fermi distribution function.
func (env *Environment) Fermi(energy float64) float64 {
	// Need to make this check to be sure we're dividing by a nonzero energy in the next step.
	if energy == 0.0 {
		return 0.5
	}
	// Temperature is 0 or e^(Beta*energy) is too big to calculate
	if env.Beta == math.Inf(1) || env.Beta >= math.Abs(math.MaxFloat64/energy) || math.Abs(env.Beta*energy) >= math.Log(math.MaxFloat64) {
		if energy <= 0 {
			return 1.0
		}
		return 0.0
	}
	// nonzero temperature
	return 1.0 / (math.Exp(energy*env.Beta) + 1.0)
}

// Create an Environment from the given serialized data.
func NewEnvironment(jsonData string) (*Environment, error) {
	// initialize env with input data
	env := new(Environment)
	err := serialize.CopyFromJSON(jsonData, env)
	if err != nil {
		return nil, err
	}

	return env, nil
}

// Convert to string by marshalling to JSON
func (env *Environment) String() string {
	if env.Beta == math.Inf(1) {
		// hack to get around JSON's choice to not allow Inf
		env.Beta = math.MaxFloat64
	}
	marshalled, err := serialize.MakeJSON(env)
	if err != nil {
		panic(err)
	}
	if env.Beta == math.MaxFloat64 {
		env.Beta = math.Inf(1)
	}
	return marshalled
}

// Iterate through v and vars simultaneously. vars specifies the names of
// fields to change in env (they are set to the values given in v).
// Panics if vars specifies a field not contained in env (or a field of
// non-float type).
func (env *Environment) Set(v vec.Vector, vars []string) {
	ev := reflect.ValueOf(env).Elem()
	for i := 0; i < len(vars); i++ {
		field := ev.FieldByName(vars[i])
		if field == reflect.Zero(reflect.TypeOf(env)) {
			panic(fmt.Sprintf("Field %v not present in Environment", vars[i]))
		}
		if field.Type().Kind() != reflect.Float64 {
			panic(fmt.Sprintf("Field %v is non-float", vars[i]))
		}
		field.SetFloat(v[i])
	}
}