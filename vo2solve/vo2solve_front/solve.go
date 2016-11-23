package main

import (
	"github.com/tflovorn/vo2mft/vo2solve"
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
)

var eps = flag.Float64("eps", 1e-6, "Converged when error below eps")
var ions = flag.Bool("ions", false, "Solve only ionic system")

func main() {
	flag.Parse()
	args := flag.Args()
	if len(args) < 2 {
		fmt.Println("Usage: vo2solve_front [--eps EPS] [--ions] in_path out_path")
		fmt.Println("For flag descriptions, use: vo2solve_front --help")
		os.Exit(2)
	}
	in_path := args[0]
	out_path := args[1]

	// Initialize Ds cache.
	Ds := vo2solve.NewHoppingEV()

	// Load Environment from in_path, then solve system
	// (throw away result from Solve - env is modified in-place).
	var solved_env *vo2solve.Environment
	if !*ions {
		env, err := vo2solve.LoadEnv(in_path)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		_, err = vo2solve.MWMuSolve(env, Ds, *eps, *eps)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		solved_env = env
	} else {
		env, err := vo2solve.LoadIonEnv(in_path)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		_, err = vo2solve.MWSolve(env, Ds, *eps, *eps)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		solved_env = env
	}

	// Calculate additional data for export from solved Environment.
	fenv := vo2solve.NewFinalEnvironment(solved_env, Ds)

	// Write output system.
	fenv_out_buf := bytes.NewBufferString(fenv.Marshal())
	fenv_out_data := fenv_out_buf.Bytes()
	ioutil.WriteFile(out_path+"_fenv.json", fenv_out_data, 0644) // u=rw;go=r
}
