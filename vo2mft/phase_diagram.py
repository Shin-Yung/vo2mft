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

def _save_min_envs(envs, prefix):
    min_envs_strl = []
    for env in envs:
        min_envs_strl.append(json.dumps(env))
    with open(prefix + "_min_data", 'w') as fp:
        fp.write('\n'.join(min_envs_strl))

def _save_all_envs(all_envs, prefix):
    all_envs_strl = []
    for env_list in all_envs:
        all_envs_strl.append(json.dumps(env_list))
    with open(prefix + "_all_data", 'w') as fp:
        fp.write('\n'.join(all_envs_strl))

def _collect_BTM(min_envs):
    xs, ys, Ms = [], [], []
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
        Ms.append(this_env["M"])

    return xs, ys, Ms

def _make_M_diagram(Bs, Ts, Ms, out_prefix):
    plt.xlabel("$b/q_J^{ion}$")
    plt.ylabel("$T/q_J^{ion}$")
    plt.title("$m$")

    plt.xlim(0.0, max(Bs))
    plt.ylim(0.0, max(Ts))

    #plt.scatter(Bs, Ts, c=Ms, cmap='gnuplot', s=100, edgecolors="none") # 10x10
    plt.scatter(Bs, Ts, c=Ms, cmap='gnuplot', s=15, edgecolors="none") # 100x100
    plt.colorbar()

    if out_prefix == None:
        plt.show()
    else:
        plt.savefig(out_prefix + '.png', bbox_inches='tight', dpi=500)

    plt.clf()

def _main():
    parser = ArgumentParser(description="Construct phase diagram")
    parser.add_argument('--base_env_path', type=str, help="Base environment file path",
            default="phase_diagram_env.json")
    parser.add_argument('--out_prefix', type=str, help="Output file path prefix",
            default="out_phase_diagram")
    parser.add_argument('--ions', action='store_true', help="Consider only ionic part")
    parser.add_argument('--num_Bs', type=int, help="Number of B points",
            default=20)
    parser.add_argument('--num_Ts', type=int, help="Number of temperature points",
            default=20)
    parser.add_argument('--npar', type=int, help="Number of parallel processes",
            default=None)
    args = parser.parse_args()

    # TODO - don't assume run in own directory
    min_envs, all_fenvs = min_envs_from_base(args.base_env_path, args.ions, args.num_Bs, args.num_Ts, args.npar)
    _save_min_envs(min_envs, args.out_prefix)
    _save_all_envs(all_fenvs, args.out_prefix)

    Bs, Ts, Ms = _collect_BTM(min_envs)
    _make_M_diagram(Bs, Ts, Ms, args.out_prefix)

if __name__ == "__main__":
    _main()
