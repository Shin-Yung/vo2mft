# Dependencies

Requires scipy, matplotlib, and GSL:

    sudo apt-get install python3-numpy python3-matplotlib python3-tk python3-scipy gsl-bin libgsl0ldbl libgsl0-dev

Requires scExplorer and cmatrix:

    go get github.com/tflovorn/scExplorer
    go get github.com/tflovorn/cmatrix

Requires tetra (Python implementation of tetrahedron method):

    cd ~
    git clone https://bitbucket.org/tflovorn/tetra.git

Requires ctetra (C implementation of tetrahedron method) and bstrlib:

    git submodule init
    git submodule update

# Installation

Get this repository:

    mkdir -p ~/gopath/src/bitbucket.org/tflovorn
    cd ~/gopath/src/bitbucket.org/tflovorn/vo2mft
    git clone https://bitbucket.org/tflovorn/vo2mft.git

Add GOPATH to ~/.bashrc:

    export GOPATH=$HOME/gopath

Get setuptools and install using setup.py:

    sudo apt-get install python3-setuptools
    sudo python3 setup.py install
    cd ~/tetra
    sudo python3 setup.py install

To have changes to the source reflected immediately instead:

    sudo python3 setup.py develop

To run setup.py without root, create ~/local and add $HOME/local to $PYTHONPATH.
Then run setup.py with:

    python3 setup.py install --prefix "/home/username/local"
