package twodof

import (
	"fmt"
	"testing"
)
import (
	"github.com/tflovorn/scExplorer/solve"
)

func TestSolveSystem(t *testing.T) {
	solve.DebugReport(true)

	env, err := LoadEnv("system_test_env.json")
	if err != nil {
		t.Fatal(err)
	}
	Ds := NewHoppingEV()

	eps := 1e-9
	result, err := MWMuSolve(env, Ds, eps, eps, false, false, false, false)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Println(result)
}

func TestSolveSystemIons(t *testing.T) {
	solve.DebugReport(true)

	env, err := LoadIonEnv("system_test_env.json")
	if err != nil {
		t.Fatal(err)
	}
	Ds := NewHoppingEV()

	eps := 1e-9
	result, err := MWSolve(env, Ds, eps, eps, false, false, false, false)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Println(result)
}
