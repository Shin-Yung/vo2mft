package vo2solve

import (
	"flag"
	"fmt"
	"math"
	"testing"
)
import (
	"github.com/tflovorn/scExplorer/solve"
)

var regression_vals = flag.Bool("regression_vals", false, "Run all regression tests, printing output without checking for errors")

func TestRegressionSolveSystem(t *testing.T) {
	flag.Parse()

	// Initial conditions: M = 1.0; W = 1.0.
	// Beta = 100.0, 10.0, and 0.1.
	// expected_result = [m, w, mu]; expected_Ds = [Dao, Dco].
	expected_result_largeM_b100 := []float64{1.0, 1.0, -2.5674589838322244}
	expected_Ds_largeM_b100 := []float64{0.07204174183966662, 0.04728135302972381}
	expected_result_largeM_b10 := []float64{0.9999999999971557, 0.9999999999971557, -2.538739633495934}
	expected_Ds_largeM_b10 := []float64{0.072591514675989, 0.04677945708493083}
	expected_result_largeM_b0_1 := []float64{-5.895244303210709e-10, 0.6664444023084701, -11.224011853014702}
	expected_Ds_largeM_b0_1 := []float64{-1.0835533804262066e-11, -5.415817658614926e-12}
	// Initial conditions: M = 0.1; W = 0.01.
	// Beta = 100.0, 10.0, and 0.1.
	expected_result_smallM_b100 := []float64{1.0, 1.0, -2.56745699439314}
	expected_Ds_smallM_b100 := []float64{0.07204176342035805, 0.047281379960874975}
	expected_result_smallM_b10 := []float64{-2.562304229335293e-10, 0.6440869234501366, -1.483771702424886}
	expected_Ds_smallM_b10 := []float64{-2.724957953642525e-11, -1.2360228333600988e-11}
	expected_result_smallM_b0_1 := []float64{-8.463415267568632e-14, 0.666444407419614, -11.224023907935099}
	expected_Ds_smallM_b0_1 := []float64{-1.5555838004397371e-15, -7.775630601417968e-16}

	all_setup := []map[string]float64{map[string]float64{"M": 1.0, "W": 1.0, "Beta": 100.0},
		map[string]float64{"M": 1.0, "W": 1.0, "Beta": 10.0},
		map[string]float64{"M": 1.0, "W": 1.0, "Beta": 0.1},
		map[string]float64{"M": 0.1, "W": 0.01, "Beta": 100.0},
		map[string]float64{"M": 0.1, "W": 0.01, "Beta": 10.0},
		map[string]float64{"M": 0.1, "W": 0.01, "Beta": 0.1}}

	all_expected := [][][]float64{[][]float64{expected_result_largeM_b100, expected_Ds_largeM_b100},
		[][]float64{expected_result_largeM_b10, expected_Ds_largeM_b10},
		[][]float64{expected_result_largeM_b0_1, expected_Ds_largeM_b0_1},
		[][]float64{expected_result_smallM_b100, expected_Ds_smallM_b100},
		[][]float64{expected_result_smallM_b10, expected_Ds_smallM_b10},
		[][]float64{expected_result_smallM_b0_1, expected_Ds_smallM_b0_1}}

	for i := 0; i < len(all_setup); i++ {
		// Start with fresh env and Ds.
		// (Ds assumes that only M/W/Mu in env are changing - need to reinitialize it
		// when we change Beta).
		env, err := LoadEnv("system_test_regression_env.json")
		if err != nil {
			t.Fatal(err)
		}
		Ds := new(HoppingEV)

		// Set variables for this test.
		setup := all_setup[i]
		for k, v := range setup {
			env.Set([]float64{v}, []string{k})
		}
		//fmt.Println(env.String())

		// Solve for (M, W, Mu).
		eps := 1e-6
		result, err := MWMuSolve(env, Ds, eps, eps)
		if err != nil {
			t.Fatal(err)
		}
		this_Dao := Ds.Dao(env)
		this_Dco := Ds.Dco(env)

		if *regression_vals {
			fmt.Println("result = ", result)
			fmt.Println("Dao = ", this_Dao)
			fmt.Println("Dco = ", this_Dco)
		} else {
			// Check results.
			expected := all_expected[i]
			if diff_float(result[0], expected[0][0]) {
				t.Fatalf("Incorrect M = %f; expected %f.", result[0], expected[0][0])
			} else if diff_float(result[1], expected[0][1]) {
				t.Fatalf("Incorrect W = %f; expected %f.", result[1], expected[0][1])
			} else if diff_float(result[2], expected[0][2]) {
				t.Fatalf("Incorrect Mu = %f; expected %f.", result[2], expected[0][2])
			} else if diff_float(this_Dao, expected[1][0]) {
				t.Fatalf("Incorrect Dao = %f; expected %f.", this_Dao, expected[1][0])
			} else if diff_float(this_Dco, expected[1][1]) {
				t.Fatalf("Incorrect Dco = %f; expected %f.", this_Dco, expected[1][1])
			}
		}
	}
}

func diff_float(x, y float64) bool {
	eps := 1e-9
	return math.Abs(x-y) > eps
}

func TestSolveSystem(t *testing.T) {
	solve.DebugReport(true)

	env, err := LoadEnv("system_test_env.json")
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

	env, err := LoadEnv("system_test_env_ions.json")
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
