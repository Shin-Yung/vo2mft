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
	expected_result_largeM_b100 := []float64{1.0, 1.0, -2.484757449884048}
	expected_Ds_largeM_b100 := []float64{0.072234636374735, 0.047713408035940756}
	expected_result_largeM_b10 := []float64{0.9999999999971739, 0.9999999999971739, -2.459610225368738}
	expected_Ds_largeM_b10 := []float64{0.07262525926112415, 0.04792141776134653}
	expected_result_largeM_b0_1 := []float64{-5.895211371059444e-10, 0.6664444023085105, -11.224011853221501}
	expected_Ds_largeM_b0_1 := []float64{-1.0835473471881098e-11, -5.415787167984652e-12}
	// Initial conditions: M = 0.1; W = 0.01.
	// Beta = 100.0, 10.0, and 0.1.
	expected_result_smallM_b100 := []float64{1.0, 1.0, -2.484757677675485}
	expected_Ds_smallM_b100 := []float64{0.07223462146639055, 0.04771339475446586}
	expected_result_smallM_b10 := []float64{-1.8376863666856174e-09, 0.6440868923032077, -1.3826195271154291}
	expected_Ds_smallM_b10 := []float64{-1.9388440095595153e-10, -7.70737994809756e-11}
	expected_result_smallM_b0_1 := []float64{-8.463415247009317e-14, 0.666444407419614, -11.224023908144346}
	expected_Ds_smallM_b0_1 := []float64{-1.5557672979419076e-15, -7.776962415169131e-16}

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
