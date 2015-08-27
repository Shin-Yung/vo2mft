from uuid import uuid4
import subprocess
import os
import numpy as np
from tetra.dos import DosValues_AllE
from vo2mft.elHamiltonian import ElHamiltonian_Recip
from vo2mft.lattice import _cubic_R
from vo2mft.util import _run_dos_path

def Dos(env, num_dos, n0, use_ctetra=True):
    '''Return two lists, dos_vals and E_vals. dos_vals contains the density
    of states D(E) at num_dos energies E between the minimum and maximum energy
    eigenvalues. E_vals contains the E values at which the corresponding element
    of dos_vals was evaluated.

    n0 gives the number of k-points to use to obtain D(E).
    '''
    if not use_ctetra:
        def Hk(k):
            return ElHamiltonian_Recip(env, k)
        return PytetraDos(Hk, num_dos, n0)

    rundos_path = _run_dos_path()
    out_name = str(uuid4())

    rundos_call = [rundos_path, out_name, str(n0), str(num_dos), str(env["Tae"]), str(env["Tce"]),
            str(env["Tbe"]), str(env["Tao"]), str(env["Tco"]), str(env["Tbo"]), str(env["EpsilonR"]),
            str(env["EpsilonM"]), str(env["M"]), str(env["W"]), str(env["Mu"])]
    subprocess.call(rundos_call)

    dos_vals, E_vals = _get_dos_vals(out_name)
    os.remove(out_name)

    return dos_vals, E_vals

def _get_dos_vals(dos_path):
    dos_vals, E_vals = [], []
    lines = None
    with open(dos_path, 'r') as fp:
        lines = fp.readlines()

    for i, line in enumerate(lines):
        # Skip header.
        if i == 0:
            continue

        split = line.strip().split('\t')
        E_vals.append(float(split[0]))
        dos_vals.append(float(split[1]))

    return dos_vals, E_vals

def PytetraDos(Hk, num_dos, n0, R=None):
    '''Return two lists, dos_vals and E_vals. dos_vals contains the density
    of states D(E) at num_dos energies E between the minimum and maximum energy
    eigenvalues. E_vals contains the E values at which the corresponding element
    of dos_vals was evaluated.

    Hk gives the Hamiltonian as a function of k, where k is in the reciprocal
    lattice basis (i.e. k = (k_1, k_2, k_3) with corresponding Cartesian
    representation k_Cart = k_1 b_1 + k_2 b_2 + k_3 b_3).

    n0 gives the number of k-points to use to obtain D(E).
    R is a numpy matrix with rows equal to the reciprocal lattice vectors.
    If R is not specified, the lattice is assumed to be simple cubic.
    '''
    if R == None:
        R = _cubic_R(1.0)

    def Efn(k):
        evals = sorted(np.linalg.eigvalsh(Hk(k)))
        return evals

    dos_vals, E_vals, tetras, Eks = DosValues_AllE(num_dos, n0, Efn, R)

    return dos_vals, E_vals

def FindGaps(dos_values, Es):
    '''Search for gaps in the density of states given by the list dos_values.
    The elements of dos_values correspond to the energies given by Es.

    Return a list with elements (gap_start, gap_stop) for each gap detected.
    If no gaps are detected, return a 0-element list.
    '''
    gaps = []
    last_dos_nonzero = False
    in_gap = False
    eps = 1e-12
    gap_start = None
    for i, dos in enumerate(dos_values):
        E = Es[i]
        zero_dos = abs(dos) < eps

        # Detect the beginning of a gap.
        # Need to come from a nonzero DOS value to a zero value.
        if not in_gap and last_dos_nonzero and zero_dos:
            in_gap = True
            gap_start = E
        # Detect the end of a gap.
        # Need to come from a zero dos value to a nonzero value
        elif in_gap and not zero_dos:
            in_gap = False
            gap_stop = E
            gaps.append((gap_start, gap_stop))
            gap_start = None
        # Detect the beginning of the bands (i.e. lowest-energy nonzero dos value).
        elif not zero_dos:
            last_dos_nonzero = True

    return gaps
