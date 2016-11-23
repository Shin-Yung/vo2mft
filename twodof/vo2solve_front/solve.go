package main

import (
	"github.com/tflovorn/vo2mft/twodof"
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
)

var eps = flag.Float64("eps", 1e-6, "Converged when error below eps")
var ions = flag.Bool("ions", false, "Solve only ionic system")
var m01_0 = flag.Bool("m01_0", false, "Fix m_01 = 0")
var m11_0 = flag.Bool("m11_0", false, "Fix m_11 = 0")
var m02_0 = flag.Bool("m02_0", false, "Fix m_02 = 0")
var m12_0 = flag.Bool("m12_0", false, "Fix m_12 = 0")

//var ions = flag.Bool("ions", false, "Solve only ionic system")

func main() {
	flag.Parse()
	args := flag.Args()
	if len(args) < 2 {
		fmt.Println("Usage: vo2solve_front [--eps EPS] in_path out_path")
		fmt.Println("For flag descriptions, use: vo2solve_front --help")
		os.Exit(2)
	}
	in_path := args[0]
	out_path := args[1]

	Ds := twodof.NewHoppingEV()

	// Load Environment from in_path, then solve system
	// (throw away result from Solve - env is modified in-place).
	var solved_env *twodof.Environment
	if !*ions {
		env, err := twodof.LoadEnv(in_path)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		if *m01_0 {
			env.M01 = 0.0
		}
		if *m11_0 {
			env.M11 = 0.0
		}
		if *m02_0 {
			env.M02 = 0.0
		}
		if *m12_0 {
			env.M12 = 0.0
		}

		_, err = twodof.MWMuSolve(env, Ds, *eps, *eps, *m01_0, *m11_0, *m02_0, *m12_0)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		solved_env = env
	} else {
		env, err := twodof.LoadIonEnv(in_path)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		if *m01_0 {
			env.M01 = 0.0
		}
		if *m11_0 {
			env.M11 = 0.0
		}
		if *m02_0 {
			env.M02 = 0.0
		}
		if *m12_0 {
			env.M12 = 0.0
		}

		_, err = twodof.MWSolve(env, Ds, *eps, *eps, *m01_0, *m11_0, *m02_0, *m12_0)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		solved_env = env
	}

	// Calculate additional data for export from solved Environment.
	fenv := twodof.NewFinalEnvironment(solved_env, Ds)

	// Write output system.
	fenv_out_buf := bytes.NewBufferString(fenv.Marshal())
	fenv_out_data := fenv_out_buf.Bytes()
	ioutil.WriteFile(out_path+"_fenv.json", fenv_out_data, 0644) // u=rw;go=r
}
