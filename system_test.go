package vo2mft

import (
	"fmt"
	"io/ioutil"
	"math"
	"testing"
)
import (
	"github.com/tflovorn/scExplorer/solve"
)

func TestRegressionSolveSystem(t *testing.T) {
	// Initial conditions: M = 1.0; W = 1.0.
	// Beta = 100.0, 10.0, and 0.1.
	// expected_result = [m, w, mu]; expected_Ds = [Dao, Dco].
	expected_result_largeM_b100 := []float64{1.0, 1.0, -2.6174549060199843}
	expected_Ds_largeM_b100 := []float64{0.07204178607541058, 0.047281408231201555}
	expected_result_largeM_b10 := []float64{0.9999999999971557, 0.9999999999971557, -2.588738149985736}
	expected_Ds_largeM_b10 := []float64{0.07259153285442703, 0.04677947795619789}
	expected_result_largeM_b0_1 := []float64{-5.779155911014687e-10, 0.6664444023617359, -11.274011799110193}
	expected_Ds_largeM_b0_1 := []float64{-1.0622161855274949e-11, -5.309170076421088e-12}
	// Initial conditions: M = 0.1; W = 0.01.
	// Beta = 100.0, 10.0, and 0.1.
	expected_result_smallM_b100 := []float64{1.0, 1.0, -2.61745487843183}
	expected_Ds_smallM_b100 := []float64{0.07204178637469986, 0.047281408604661956}
	expected_result_smallM_b10 := []float64{-3.549391415937292e-10, 0.6440869201456757, -1.5337717118249259}
	expected_Ds_smallM_b10 := []float64{-3.7747049265576935e-11, -1.712181115574212e-11}
	expected_result_smallM_b0_1 := []float64{-8.494291904455557e-14, 0.6664444074196134, -11.274023907908134}
	expected_Ds_smallM_b0_1 := []float64{-1.5612399322102982e-15, -7.803446757606022e-16}

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
		env, err := loadEnv("system_test_regression_env.json")
		if err != nil {
			t.Fatal(err)
		}
		Ds := new(HoppingEV)

		// Set variables for this test.
		setup, expected := all_setup[i], all_expected[i]
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

		//fmt.Println("result = ", result)
		//fmt.Println("Dao = ", this_Dao)
		//fmt.Println("Dco = ", this_Dco)

		// Check results.
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

func diff_float(x, y float64) bool {
	eps := 1e-9
	return math.Abs(x-y) > eps
}

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
