package twodof

import (
	"fmt"
	"io/ioutil"
	"math"
	"reflect"
)
import (
	//	"github.com/tflovorn/cmatrix"
	//	"github.com/tflovorn/scExplorer/bzone"
	"github.com/tflovorn/scExplorer/serialize"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Contains parameters necessary to characterize electronic and ionic systems.
// The ionic order parameters M and W and the electronic chemical potential Mu
// must be determined self-consistently.
type Environment struct {
	// Order parameter <S_{p,alpha}>.
	M01, M11, M02, M12 float64
	// Inverse temperature, 1 / (k_B * T).
	Beta float64
	// One-spin term for BEG model: coefficient for (S_i)^2.
	Bxy0, Bzz0 float64
	// Exchange parameters for BEG model: coefficients to S_i dot S_j.
	// Jb is excluded since it does not contribute to results.
	Jb0, Jc0 float64
}

func (env *Environment) Bxy() float64 {
	// TODO - add strain dependence
	return env.Bxy0
}

func (env *Environment) Bzz() float64 {
	return env.Bzz0
}

func (env *Environment) Jb() float64 {
	return env.Jb0
}

func (env *Environment) Jc() float64 {
	return env.Jc0
}

// Environment with all self-consistent values converged.
// Includes additional data for exporting to outside programs.
type FinalEnvironment struct {
	Environment
	FreeEnergy float64
}

// Free energy per cell value (Ncell = 2Nsite).
// Points on the phase diagram include the state with minimum free energy
// (may not reach this state, depending on initial conditions - need to
// consider a set of initial conditions and look for minimum).
func (env *Environment) FreeEnergy() float64 {
	ion_part := env.FreeEnergyIons()
	// avg_avg_part includes <S><S> terms.
	avg_avg_part := env.EConst_Ion()

	return ion_part + avg_avg_part
}

func (env *Environment) FreeEnergyIons() float64 {
	T := 1.0 / env.Beta
	return -T * math.Log(env.Z1())
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

// Create a FinalEnvironment from the given solved Environment and associated
// HoppingEV.
func NewFinalEnvironment(env *Environment) *FinalEnvironment {
	FreeEnergy := env.FreeEnergy()
	fenv := FinalEnvironment{*env, FreeEnergy}
	return &fenv
}

// Load an Environment from the JSON file at envFilePath.
func LoadEnv(envFilePath string) (*Environment, error) {
	data, err := ioutil.ReadFile(envFilePath)
	if err != nil {
		return nil, err
	}
	env, err := NewEnvironment(string(data))
	if err != nil {
		return nil, err
	}
	return env, nil
}

// Convert to string by marshalling to JSON
func (env *Environment) String() string {
	marshalled := env.Marshal()
	return marshalled
}

func (env *Environment) Marshal() string {
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

func (env *FinalEnvironment) String() string {
	marshalled := env.Marshal()
	return marshalled
}

func (env *FinalEnvironment) Marshal() string {
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

// Return the value of the env variable with type float64 with the given name.
func (env *Environment) GetFloat(var_name string) float64 {
	ev := reflect.ValueOf(env).Elem()
	return ev.FieldByName(var_name).Float()
}
