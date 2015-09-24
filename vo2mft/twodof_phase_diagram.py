from argparse import ArgumentParser
from copy import deepcopy
from multiprocessing import Pool
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from vo2mft.twodof_environment import Jbe
from vo2mft.min_free_energy import minimize_free_energy
from vo2mft.solve import read_env_file
from vo2mft.dos import Dos, FindGaps
from vo2mft.plot_spectrum import plot_spectrum

def phase_sample(base_env, num_Bs, num_Ts):
    '''Generate a set of envs over which the phase diagram will be sampled,
    with constant properties taken from base_env.
    '''
    # Generate set of (B, T) values.
    base_Jbe = Jbe(base_env)

    #Bratio_start, Bratio_stop = 0.01, 0.7
    #Tratio_start, Tratio_stop = 0.01, 0.8
    Bratio_start, Bratio_stop = 0.01, 1.0
    Tratio_start, Tratio_stop = 0.01, 1.0

    Bratios = np.linspace(Bratio_start, Bratio_stop, num_Bs)
    Tratios = np.linspace(Tratio_start, Tratio_stop, num_Ts)

    B_T_vals = []
    for Br in Bratios:
        for Tr in Tratios:
            B_val = Br*base_Jbe
            T_val = Tr*base_Jbe
            B_T_vals.append([B_val, T_val])

    xy_zz_ratio = base_env["Bxy0"] / base_env["Bzz0"]
    sample_envs = []
    for B, T in B_T_vals:
        this_env = deepcopy(base_env)
        this_env["Bxy0"] = xy_zz_ratio * B
        this_env["Bzz0"] = B
        this_env["Beta"] = 1.0/T
        sample_envs.append(this_env)

    return sample_envs

def min_envs_from_sample(sample_envs, ions, npar):
    eps = 1e-8
    sample_env_eps = []
    for env in sample_envs:
        twodof, twodof_body_indep = True, True
        sample_env_eps.append((env, eps, ions, twodof, twodof_body_indep))

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

def _collect_BT_env_val(min_envs, val_name, func=None):
    xs, ys, vals = [], [], []
    for this_env in min_envs:
        # May not have found a solution for all envs.
        if this_env == None:
            continue
        # This env was solved -- add it to plot.
        this_Jbe = Jbe(this_env)
        Bratio = this_env["Bxy0"] / this_Jbe
        Tratio = (1.0 / this_env["Beta"]) / this_Jbe
        xs.append(Bratio)
        ys.append(Tratio)

        if func == None:
            vals.append(this_env[val_name])
        else:
            if func == "abs":
                vals.append(abs(this_env[val_name]))
            elif func == "phase":
                def phase_func(x):
                    eps = 1e-6
                    if abs(x) > eps:
                        return 1.0
                    else:
                        return 0.0
                vals.append(phase_func(this_env[val_name]))

    return xs, ys, vals

def _make_val_diagram(Bs, Ts, vals, val_label, out_prefix, cbar_format=None):
    plt.xlabel("$b_{xx}/4J_{b}$", fontsize='x-large')
    plt.ylabel("$T/4J_{b}$", fontsize='x-large')
    plt.title(val_label, fontsize='x-large')

    plt.xlim(0.0, max(Bs))
    plt.ylim(0.0, max(Ts))

    #plt.scatter(Bs, Ts, c=vals, cmap='gnuplot', s=100, edgecolors="none") # 10x10
    #plt.scatter(Bs, Ts, c=vals, cmap='gnuplot', s=15, edgecolors="none") # 100x100
    #plt.scatter(Bs, Ts, c=vals, cmap='Set1', s=15, edgecolors="none") # 100x100
    plt.scatter(Bs, Ts, c=vals, cmap='Paired', s=15, edgecolors="none") # 100x100

    if cbar_format != None:
        plt.colorbar(format=cbar_format)
    else:
        plt.colorbar()

    if out_prefix == None:
        plt.show()
    else:
        plt.savefig(out_prefix + '.png', bbox_inches='tight', dpi=500)

    plt.clf()

def _make_BT_plot(min_envs, out_prefix, env_val, env_val_label, func=None, cbar_format=None):
    Bs, Ts, Ms = _collect_BT_env_val(min_envs, env_val, func)
    _make_val_diagram(Bs, Ts, Ms, env_val_label, "{}_{}".format(out_prefix, env_val), cbar_format)

def _check_only_ok(Bratio, Tratio, only_B, only_T):
    if only_B != None and only_T != None:
        eps = 1e-9
        b_ok = abs(Bratio - only_B) < eps
        T_ok = abs(Tratio - only_T) < eps
        return b_ok and T_ok
    else:
        return True

def _main():
    parser = ArgumentParser(description="Construct phase diagram")
    parser.add_argument('--base_env_path', type=str, help="Base environment file path",
            default="twodof_phase_diagram_env.json")
    parser.add_argument('--read_prefix', type=str, default=None,
            help="If specified, read solved envs from here instead of solving new systems")
    parser.add_argument('--out_prefix', type=str, help="Output file path prefix",
            default="twodof_out_phase_diagram")
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
    parser.add_argument('--only_B', type=float, default=None,
            help="If specified with only_T and plot_spectrum or plot_dos, make the spectrum or dos plots only for the given (B, T) point")
    parser.add_argument('--only_T', type=float, default=None,
            help="If specified with only_B and plot_spectrum or plot_dos, make the spectrum or dos plots only for the given (B, T) point")
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

    M_plot_args = [[min_envs, args.out_prefix, "M01", "$|m_{0,1}|$", "abs", "%.2f"],
        [min_envs, args.out_prefix, "M11", "$|m_{1,1}|$", "abs", "%.2f"],
        [min_envs, args.out_prefix, "M02", "$|m_{0,2}|$", "abs", "%.2f"],
        [min_envs, args.out_prefix, "M12", "$|m_{1,2}|$", "abs", "%.2f"]]

    M_phase_plot_args = [[min_envs, args.out_prefix + "_phase", "M01", "$|m_{0,1}| > 0$", "phase", "%.2f"],
        [min_envs, args.out_prefix + "_phase", "M11", "$|m_{1,1}| > 0$", "phase", "%.2f"],
        [min_envs, args.out_prefix + "_phase", "M02", "$|m_{0,2}| > 0$", "phase", "%.2f"],
        [min_envs, args.out_prefix + "_phase", "M12", "$|m_{1,2}| > 0$", "phase", "%.2f"]]

    W_plot_args = [[min_envs, args.out_prefix, "W01", "$w_{0,1}$", None, "%.2f"],
        [min_envs, args.out_prefix, "W11", "$w_{1,1}$", None, "%.2f"],
        [min_envs, args.out_prefix, "W02", "$w_{0,2}$", None, "%.2f"],
        [min_envs, args.out_prefix, "W12", "$w_{1,2}$", None, "%.2f"]]

    plot_args = []
    plot_args.extend(M_plot_args)
    plot_args.extend(M_phase_plot_args)
    plot_args.extend(W_plot_args)
    plot_args.append([min_envs, args.out_prefix, "FreeEnergy", "$F$", None, None])

    if not args.ions:
        plot_args.append([min_envs, args.out_prefix, "Mu", "$\mu$", None, None])

    with Pool() as pool:
        pool.starmap(_make_BT_plot, plot_args)

if __name__ == "__main__":
    _main()
