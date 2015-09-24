package twodof

import (
	"fmt"
	"io/ioutil"
	"math"
	"reflect"
)
import (
	"github.com/tflovorn/cmatrix"
	"github.com/tflovorn/scExplorer/bzone"
	"github.com/tflovorn/scExplorer/serialize"
	vec "github.com/tflovorn/scExplorer/vector"
)

// Contains parameters necessary to characterize electronic and ionic systems.
// The ionic order parameters M and W and the electronic chemical potential Mu
// must be determined self-consistently.
type Environment struct {
	// Size of k mesh.
	BZPointsPerDim int
	// Order parameter <S_{p,alpha}>.
	M01, M11, M02, M12 float64
	// Order parameter <S^2_{p,alpha}>.
	W01, W11, W02, W12 float64
	// Inverse temperature, 1 / (k_B * T).
	Beta float64
	// One-spin term for BEG model: coefficient for (S_i)^2.
	Bxy0, Bzz0, Bxz0 float64
	// Exchange parameters for BEG model: coefficients to S_i dot S_j.
	// Jb is excluded since it does not contribute to results.
	Jb0, Jc0 float64
	// Quartic contribution corresponding to Jb term.
	Kb0 float64
	// Quartic contributions along dimers.
	Kcxx0, Kczz0, Kcxz0 float64
	// Hopping along c axis. TODO - strain dependence.
	Tce, Tco float64
	// Hopping along body diagonal.
	Tbe float64
	// Electron chemical potential.
	Mu float64
	// Only do ionic part of calculation (all electronic quantities --> 0)
	IonsOnly bool
}

func (env *Environment) Bxy() float64 {
	// TODO - add strain dependence
	return env.Bxy0
}

func (env *Environment) Bzz() float64 {
	return env.Bzz0
}

func (env *Environment) Bxz() float64 {
	return env.Bxz0
}

func (env *Environment) Jb() float64 {
	return env.Jb0
}

func (env *Environment) Jc() float64 {
	return env.Jc0
}

func (env *Environment) Kb() float64 {
	return env.Kb0
}

func (env *Environment) Kcxx() float64 {
	return env.Kcxx0
}

func (env *Environment) Kczz() float64 {
	return env.Kczz0
}

func (env *Environment) Kcxz() float64 {
	return env.Kcxz0
}

// Are electronic hopping finite?
// If not, don't need to calculate D's.
func (env *Environment) FiniteHoppings() bool {
	eps := 1e-9
	even := (math.Abs(env.Tce) > eps) || (math.Abs(env.Tbe) > eps)
	odd := math.Abs(env.Tco) > eps
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

// Environment with all self-consistent values converged.
// Includes additional data for exporting to outside programs.
type FinalEnvironment struct {
	Environment
	Dco        float64
	FreeEnergy float64
}

// Free energy per cell value (Ncell = 2Nsite).
// Points on the phase diagram include the state with minimum free energy
// (may not reach this state, depending on initial conditions - need to
// consider a set of initial conditions and look for minimum).
func (env *Environment) FreeEnergy(Ds *HoppingEV) float64 {
	ion_part := env.FreeEnergyIons(Ds)
	// avg_avg_part includes <S><S> terms.
	avg_avg_part := env.EConst_Ion() + env.EConst_IonEl(Ds)

	if env.IonsOnly {
		return ion_part + avg_avg_part
	} else {
		electron_part := env.FreeEnergyElectrons()
		return ion_part + electron_part + avg_avg_part
	}
}

func (env *Environment) FreeEnergyIons(Ds *HoppingEV) float64 {
	T := 1.0 / env.Beta
	return -T * math.Log(env.Z1(Ds))
}

func (env *Environment) FreeEnergyElectrons() float64 {
	inner := func(k vec.Vector) float64 {
		H := ElHamiltonian(env, k)
		dim, _ := H.Dims()
		evals, _ := cmatrix.Eigensystem(H)
		sum := 0.0
		for alpha := 0; alpha < dim; alpha++ {
			eps_ka := evals[alpha]
			// Mu excluded from exp argument here since it is
			// included in H.
			val := 1.0 + math.Exp(-env.Beta*eps_ka)
			// Factor of 2 for spins.
			sum += 2.0 * math.Log(val)
		}
		return sum
	}
	L := env.BZPointsPerDim
	T := 1.0 / env.Beta
	band_part := -T * bzone.Avg(L, 3, inner)

	return band_part
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
func NewFinalEnvironment(env *Environment, Ds *HoppingEV) *FinalEnvironment {
	Dco := Ds.Dco(env)
	FreeEnergy := env.FreeEnergy(Ds)
	fenv := FinalEnvironment{*env, Dco, FreeEnergy}
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

// Load an Environment from the JSON file at envFilePath.
// Set all electronic parameters to 0 to restrict to ionic system.
func LoadIonEnv(envFilePath string) (*Environment, error) {
	env, err := LoadEnv(envFilePath)
	if err != nil {
		return nil, err
	}
	env.Tce = 0.0
	env.Tbe = 0.0
	env.Tco = 0.0
	env.Mu = 0.0
	env.IonsOnly = true
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
