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

	env, err := defaultEnv()
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

func defaultEnv() (*Environment, error) {
        data, err := ioutil.ReadFile("system_test_env.json")
        if err != nil {
                return nil, err
        }
        env, err := NewEnvironment(string(data))
        if err != nil {
                return nil, err
        }
        return env, nil
}
