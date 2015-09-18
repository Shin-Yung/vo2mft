package main

import (
	"bitbucket.org/tflovorn/vo2mft/twodof"
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
)

var eps = flag.Float64("eps", 1e-6, "Converged when error below eps")

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

	// Load Environment from in_path, then solve system
	// (throw away result from Solve - env is modified in-place).
	var solved_env *twodof.Environment
	env, err := twodof.LoadEnv(in_path)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	_, err = twodof.MSolve(env, *eps, *eps)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	solved_env = env

	// Calculate additional data for export from solved Environment.
	fenv := twodof.NewFinalEnvironment(solved_env)

	// Write output system.
	fenv_out_buf := bytes.NewBufferString(fenv.Marshal())
	fenv_out_data := fenv_out_buf.Bytes()
	ioutil.WriteFile(out_path+"_fenv.json", fenv_out_data, 0644) // u=rw;go=r
}
