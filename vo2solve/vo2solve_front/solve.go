package main

import (
	"bitbucket.org/tflovorn/vo2mft/vo2solve"
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
)

var eps = flag.Float64("eps", 1e-6, "Converged when error below eps")

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

	// Load Environment from in_path.
	env, err := vo2solve.LoadEnv(in_path)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	// Initialize Ds cache.
	Ds := new(vo2solve.HoppingEV)

	// Solve system (throw away result -- env is modified in-place).
	_, err = vo2solve.MWMuSolve(env, Ds, *eps, *eps)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	// Calculate additional data for export from solved Environment.
	fenv := vo2solve.NewFinalEnvironment(env, Ds)

	// Write output system.
	fenv_out_buf := bytes.NewBufferString(fenv.Marshal())
	fenv_out_data := fenv_out_buf.Bytes()
	ioutil.WriteFile(out_path+"_fenv.json", fenv_out_data, 0644) // u=rw;go=r
}
