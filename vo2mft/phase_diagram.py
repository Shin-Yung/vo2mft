from argparse import ArgumentParser
from copy import deepcopy
from multiprocessing import Pool
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from vo2mft.environment import QJ_ion
from vo2mft.min_free_energy import minimize_free_energy
from vo2mft.solve import read_env_file
from vo2mft.dos import Dos, FindGaps
from vo2mft.plot_spectrum import plot_spectrum

def phase_sample(base_env, num_Bs, num_Ts):
    '''Generate a set of envs over which the phase diagram will be sampled,
    with constant properties taken from base_env.
    '''
    # Generate set of (B, T) values.
    # TODO - generalize away from EpsilonM - EpsilonR = 0.
    base_QJ_ion = QJ_ion(base_env)
    #Bratio_start, Bratio_stop = 0.01, 0.7  # ions
    Bratio_start, Bratio_stop = 0.01, 1.2   # electrons
    #Tratio_start, Tratio_stop = 0.1, 1.0
    Tratio_start, Tratio_stop = 0.01, 0.8

    Bratios = np.linspace(Bratio_start, Bratio_stop, num_Bs)
    Tratios = np.linspace(Tratio_start, Tratio_stop, num_Ts)

    B_T_vals = []
    for Br in Bratios:
        for Tr in Tratios:
            B_val = Br*base_QJ_ion
            T_val = Tr*base_QJ_ion
            B_T_vals.append([B_val, T_val])

    sample_envs = []
    for B, T in B_T_vals:
        this_env = deepcopy(base_env)
        this_env["B"] = B
        this_env["Beta"] = 1.0/T
        sample_envs.append(this_env)

    return sample_envs

def min_envs_from_sample(sample_envs, ions, npar):
    eps = 1e-6
    sample_env_eps = []
    for env in sample_envs:
        sample_env_eps.append((env, eps, ions))

    # Find minimal free energy solutions of sample_envs in parallel.
    # Pool() uses number of workers = os.cpu_count() by default.
    pool_size = os.cpu_count()
    if npar != None:
        pool_size = npar

    minimize_result = None
    with Pool(pool_size) as pool:
        minimize_result = pool.starmap(minimize_free_energy, sample_env_eps)

    all_min_envs, all_final_envs = [], []
    for min_env, final_envs in minimize_result:
        all_min_envs.append(min_env)
        all_final_envs.append(final_envs)

    return all_min_envs, all_final_envs

def min_envs_from_base(base_path, ions, num_Bs, num_Ts, npar):
    base_env = read_env_file(base_path)
    sample = phase_sample(base_env, num_Bs, num_Ts)
    min_envs, all_fenvs = min_envs_from_sample(sample, ions, npar)
    return min_envs, all_fenvs

def _min_env_filename(prefix):
    return prefix + "_min_data"

def _save_min_envs(envs, prefix):
    min_envs_strl = []
    for env in envs:
        min_envs_strl.append(json.dumps(env))
    with open(_min_env_filename(prefix), 'w') as fp:
        fp.write('\n'.join(min_envs_strl))

def _read_min_envs(prefix):
    min_envs_lines = None
    with open(_min_env_filename(prefix), 'r') as fp:
        min_envs_lines = fp.readlines()

    min_envs = []
    for line in min_envs_lines:
        min_envs.append(json.loads(line.strip()))

    return min_envs

def _save_all_envs(all_envs, prefix):
    all_envs_strl = []
    for env_list in all_envs:
        all_envs_strl.append(json.dumps(env_list))
    with open(prefix + "_all_data", 'w') as fp:
        fp.write('\n'.join(all_envs_strl))

def _collect_BT_env_val(min_envs, val_name):
    xs, ys, vals = [], [], []
    for this_env in min_envs:
        # May not have found a solution for all envs.
        if this_env == None:
            continue
        # This env was solved -- add it to plot.
        this_QJ_ion = QJ_ion(this_env)
        Bratio = this_env["B"] / this_QJ_ion
        Tratio = (1.0 / this_env["Beta"]) / this_QJ_ion
        xs.append(Bratio)
        ys.append(Tratio)
        vals.append(this_env[val_name])

    return xs, ys, vals

def _collect_BTgaps(min_envs, npar, num_dos=3000, n0=8):
    BTgap_args = []
    for this_env in min_envs:
        BTgap_args.append((this_env, num_dos, n0))

    pool_size = os.cpu_count()
    if npar != None:
        pool_size = npar

    BTgap_result = None
    with Pool() as pool:
        BTgap_result = pool.starmap(_one_BTgap, BTgap_args)

    xs, ys, gaps, dos_vals, E_vals = [], [], [], [], []
    for this_result in BTgap_result:
        if this_result == None:
            continue
        Bratio, Tratio, gap_norm, this_dos_vals, this_E_vals = this_result
        xs.append(Bratio)
        ys.append(Tratio)
        gaps.append(gap_norm)
        dos_vals.append(this_dos_vals)
        E_vals.append(this_E_vals)

    return xs, ys, gaps, dos_vals, E_vals

def _one_BTgap(this_env, num_dos, n0):
    # May not have found a solution for all envs.
    if this_env == None:
        return None
    # This env was solved -- add it to plot.
    this_QJ_ion = QJ_ion(this_env)
    Bratio = this_env["B"] / this_QJ_ion
    Tratio = (1.0 / this_env["Beta"]) / this_QJ_ion

    this_dos_vals, this_E_vals = Dos(this_env, num_dos, n0)
    this_gaps = FindGaps(this_dos_vals, this_E_vals)

    if len(this_gaps) > 1:
        print("WARNING: found more than one gap;")
        print("assuming first gap is around E_F -- this may be incorrect.")

    gap_size = None
    if len(this_gaps) == 0:
        gap_size = 0.0
    else:
        gap_size = this_gaps[0][1] - this_gaps[0][0]

    gap_norm = gap_size / this_QJ_ion

    return (Bratio, Tratio, gap_norm, this_dos_vals, this_E_vals)


def _dos_data_path(out_prefix):
    return out_prefix + "_dos_data"

def _save_dos(Bs, Ts, gaps, dos_vals, E_vals, out_prefix):
    dos_data = []
    for B, T, gap, this_dos_vals, this_E_vals in zip(Bs, Ts, gaps, dos_vals, E_vals):
        dos_dict = {"B": B, "T": T, "dos_vals": this_dos_vals, "E_vals": list(this_E_vals)}
        dos_data.append(json.dumps(dos_dict))

    with open(_dos_data_path(out_prefix), 'w') as fp:
        fp.write('\n'.join(dos_data))

def _make_val_diagram(Bs, Ts, vals, val_label, out_prefix):
    plt.xlabel("$b/q_J^{ion}$", fontsize='x-large')
    plt.ylabel("$T/q_J^{ion}$", fontsize='x-large')
    plt.title(val_label, fontsize='x-large')

    plt.xlim(0.0, max(Bs))
    plt.ylim(0.0, max(Ts))

    #plt.scatter(Bs, Ts, c=vals, cmap='gnuplot', s=100, edgecolors="none") # 10x10
    plt.scatter(Bs, Ts, c=vals, cmap='gnuplot', s=15, edgecolors="none") # 100x100
    plt.colorbar()

    if out_prefix == None:
        plt.show()
    else:
        plt.savefig(out_prefix + '.png', bbox_inches='tight', dpi=500)

    plt.clf()

def _make_BT_plot(min_envs, out_prefix, env_val, env_val_label):
    Bs, Ts, Ms = _collect_BT_env_val(min_envs, env_val)
    _make_val_diagram(Bs, Ts, Ms, env_val_label, "{}_{}".format(out_prefix, env_val))

def _plot_dos(Bs, Ts, dos_vals, E_vals, out_prefix):
    for B, T, dos, Es in zip(Bs, Ts, dos_vals, E_vals):
        out_name = out_prefix + "_B_{}_T_{}_dos.png".format(str(B), str(T))

        plt.xlabel("$E$")
        plt.ylabel("DOS")
        plt.xlim(min(Es), max(Es))
        plt.ylim(min(dos), max(dos))

        plt.plot(Es, dos, 'k-')

        plt.savefig(out_name, bbox_inches='tight', dpi=500)
        plt.clf()

def _plot_spectra(min_envs, out_prefix):
    for env in min_envs:
        B, Beta = env["B"], env["Beta"]
        T = 1.0/Beta
        this_QJ_ion = QJ_ion(env)
        Bratio, Tratio = B / this_QJ_ion, T / this_QJ_ion
        out_name = out_prefix + "_B_{}_T_{}_spectrum".format(str(Bratio), str(Tratio))
        plot_spectrum(env, out_name)

def _main():
    parser = ArgumentParser(description="Construct phase diagram")
    parser.add_argument('--base_env_path', type=str, help="Base environment file path",
            default="phase_diagram_env.json")
    parser.add_argument('--read_prefix', type=str, default=None,
            help="If specified, read solved envs from here instead of solving new systems")
    parser.add_argument('--out_prefix', type=str, help="Output file path prefix",
            default="out_phase_diagram")
    parser.add_argument('--ions', action='store_true', help="Consider only ionic part")
    parser.add_argument('--num_Bs', type=int, help="Number of B points",
            default=20)
    parser.add_argument('--num_Ts', type=int, help="Number of temperature points",
            default=20)
    parser.add_argument('--npar', type=int, help="Number of parallel processes",
            default=None)
    parser.add_argument('--plot_spectrum', action='store_true',
            help="If specified, plot spectrum for each (B, T) point")
    parser.add_argument('--plot_dos', action='store_true',
            help="If specified, plot DOS for each (B, T) point")
    args = parser.parse_args()

    # TODO -- add only_B, only_T options to fix (B, T) point for generating
    # a few spectra/dos

    # TODO - don't assume run in own directory

    min_envs, all_fenvs = None, None
    if args.read_prefix == None:
        min_envs, all_fenvs = min_envs_from_base(args.base_env_path, args.ions, args.num_Bs, args.num_Ts, args.npar)
        _save_min_envs(min_envs, args.out_prefix)
        _save_all_envs(all_fenvs, args.out_prefix)
    else:
        min_envs = _read_min_envs(args.read_prefix)

    _make_BT_plot(min_envs, args.out_prefix, "M", "$m$")
    _make_BT_plot(min_envs, args.out_prefix, "W", "$w$")

    if not args.ions:
        _make_BT_plot(min_envs, args.out_prefix, "Mu", "$\mu$")
        _make_BT_plot(min_envs, args.out_prefix, "Dae", "$D_{ae}$")
        _make_BT_plot(min_envs, args.out_prefix, "Dce", "$D_{ce}$")
        _make_BT_plot(min_envs, args.out_prefix, "Dbe", "$D_{be}$")
        _make_BT_plot(min_envs, args.out_prefix, "Dao", "$D_{ao}$")
        _make_BT_plot(min_envs, args.out_prefix, "Dco", "$D_{co}$")
        _make_BT_plot(min_envs, args.out_prefix, "Dbo", "$D_{bo}$")

        if args.plot_spectrum:
            _plot_spectra(min_envs, args.out_prefix)

        Bs, Ts, gaps, dos_vals, E_vals = _collect_BTgaps(min_envs, args.npar)
        _make_val_diagram(Bs, Ts, gaps, "Gap / $q_J^{ion}$", args.out_prefix + "_gap")

        _save_dos(Bs, Ts, gaps, dos_vals, E_vals, args.out_prefix)
        if args.plot_dos:
            _plot_dos(Bs, Ts, dos_vals, E_vals, args.out_prefix)

if __name__ == "__main__":
    _main()
