# Dependencies

Requires Go, scipy, matplotlib, and GSL. On recent Debian-based distributions, obtain with:

    sudo apt-get install golang python3-setuptools python3-numpy python3-matplotlib python3-tk python3-scipy libgsl-dev

Note - for older versions of Debian-based distributions, instead of installing libgsl-dev, obtain GSL with:

    sudo apt-get install gsl-bin libgsl0ldbl libgsl0-dev

If you don't have a GOPATH:

    mkdir ~/gopath

Then add the following to ~/.bashrc:

    export GOPATH=$HOME/gopath

Requires scExplorer and cmatrix:

    go get github.com/tflovorn/scExplorer
    go get github.com/tflovorn/cmatrix

The root-finding implementation in scExplorer uses cgo to interface with GSL.
In Go 1.6, a new rule was implemented which forbids passing a Go object
through C which contains a pointer to another Go object.
For now, we need to turn off this rule.
Add the following to ~/.bashrc:

    export GODEBUG=cgocheck=0

Requires tetra (Python implementation of tetrahedron method):

    cd ~
    git clone https://github.com/tflovorn/tetra.git
    cd tetra
    python3 setup.py develop --user

# Installation

Get this repository:

    go get github.com/tflovorn/vo2mft
    cd $GOPATH/github.com/tflovorn/vo2mft

Build the solver:

    cd vo2solve/vo2solve_front/
    go build
    cd ../../twodof/vo2solve_front
    go build
    cd ../..

Get libraries ctetra (C implementation of tetrahedron method) and bstrlib, included as submodules:

    git submodule init
    git submodule update

Build the DOS calculator:

    cd tetra_dos/ctetra
    make
    cd ..
    make
    cd ..

Get setuptools and install using setup.py:

    python3 setup.py develop --user

(NOTE - the mechanism used to traverse the directory structure of the repo is broken
if setup.py install is used instead of develop. TODO - fix this?)

# Usage

To build Fig. 2 of "Complex quasi two-dimensional crystalline order embedded in VO2 and other crystals":

    cd vo2mft
    python3 twodof_phase_diagram.py --multi_b_cutoff
