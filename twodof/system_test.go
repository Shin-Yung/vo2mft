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
	eps := 1e-9
	result, err := MSolve(env, eps, eps)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Println(result)
}
