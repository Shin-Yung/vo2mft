package vo2mft

import (
	"testing"
	"fmt"
	"io/ioutil"
)
import (
        "github.com/tflovorn/scExplorer/solve"
)

func TestSolveSystem(t *testing.T) {
	solve.DebugReport(true)

	env, err := loadEnv("system_test_env.json")
	if err != nil {
		t.Fatal(err)
	}
	Ds := new(HoppingEV)
	eps := 1e-6
	result, err := MWMuSolve(env, Ds, eps, eps)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Println(result)
}

func TestSolveSystemIons(t *testing.T) {
	solve.DebugReport(true)

	env, err := loadEnv("system_test_env_ions.json")
	if err != nil {
		t.Fatal(err)
	}
	Ds := new(HoppingEV)
	eps := 1e-6
	result, err := MWSolve(env, Ds, eps, eps)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Println(result)
}

func loadEnv(envFilePath string) (*Environment, error) {
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
