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

def phase_sample(base_env, num_Bs, num_Ts, bounds=None):
    '''Generate a set of envs over which the phase diagram will be sampled,
    with constant properties taken from base_env.
    '''
    # Generate set of (B, T) values.
    base_Jbe = Jbe(base_env)

    #Bratio_start, Bratio_stop = 0.01, 0.7
    #Tratio_start, Tratio_stop = 0.01, 0.8
    Bratio_start, Bratio_stop = 0.01, 0.6
    Tratio_start, Tratio_stop = 0.01, 0.8
    if bounds != None:
        Bratio_start = bounds[0][0]
        Bratio_stop = bounds[0][1]
        Tratio_start = bounds[1][0]
        Tratio_stop = bounds[1][1]

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

def min_envs_from_base(base_path, ions, num_Bs, num_Ts, npar, bounds=None):
    base_env = read_env_file(base_path)
    sample = phase_sample(base_env, num_Bs, num_Ts, bounds)
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

def _make_val_diagram(Bs, Ts, vals, val_label, out_prefix, cbar_format=None, clim_vals=None, cbar_tick_labels=None):
    plt.xlabel("$b_{x}/4J_{b}$", fontsize='x-large')
    plt.ylabel("$T/4J_{b}$", fontsize='x-large')
    plt.title(val_label, fontsize='x-large')

    plt.xlim(0.0, max(Bs))
    plt.ylim(0.0, max(Ts))

    #plt.scatter(Bs, Ts, c=vals, cmap='gnuplot', s=100, edgecolors="none") # 10x10
    #plt.scatter(Bs, Ts, c=vals, cmap='gnuplot', s=15, edgecolors="none") # 100x100
    #plt.scatter(Bs, Ts, c=vals, cmap='Set1', s=15, edgecolors="none") # 100x100
    #plt.scatter(Bs, Ts, c=vals, cmap='Paired', s=15, edgecolors="none") # 100x100
    plt.scatter(Bs, Ts, c=vals, cmap='viridis', s=15, edgecolors="none") # 100x100

    cbar_ticks = None
    if clim_vals != None and cbar_tick_labels == None:
        cbar_ticks = np.linspace(clim_vals[0], clim_vals[1], 11, endpoint=True)
    if cbar_tick_labels != None:
        cbar_ticks = cbar_tick_labels[0]

    cbar = None
    if cbar_format != None:
        if clim_vals != None:
            cbar = plt.colorbar(ticks=cbar_ticks, format=cbar_format)
        else:
            cbar = plt.colorbar(format=cbar_format)
    else:
        if clim_vals != None:
            cbar = plt.colorbar(ticks=cbar_ticks)
        else:
            cbar = plt.colorbar()

    if cbar_tick_labels != None:
        cbar.ax.set_yticklabels(cbar_tick_labels[1])

    if out_prefix == None:
        plt.show()
    else:
        plt.savefig(out_prefix + '.png', bbox_inches='tight', dpi=500)
        plt.savefig(out_prefix + '.eps', bbox_inches='tight', dpi=500)

    plt.clf()

def _make_BT_plot(min_envs, out_prefix, env_val, env_val_label, func=None, cbar_format=None, clim_vals=None):
    Bs, Ts, Ms = _collect_BT_env_val(min_envs, env_val, func)
    _make_val_diagram(Bs, Ts, Ms, env_val_label, "{}_{}".format(out_prefix, env_val), cbar_format, clim_vals)

def _multival_phase_plot(min_envs, out_prefix, env_val_1, env_val_2, val_label, func=None, cbar_format=None, clim_vals=None, cbar_tick_labels=None):
    Bs, Ts, val1s = _collect_BT_env_val(min_envs, env_val_1, func=None)
    Bs, Ts, val2s = _collect_BT_env_val(min_envs, env_val_2, func=None)
    combined_vals = []

    for val1, val2 in zip(val1s, val2s):
        if func == "phase_incl_m2":
            def phase_func(m1, m2):
                eps = 1e-6
                if abs(m1) > eps and abs(m2) > eps:
                    return 1.0
                elif (abs(m1) > eps and abs(m2) < eps) or (abs(m1) < eps and abs(m2) > eps):
                    return 0.5
                else:
                    return 0.0
            combined_vals.append(phase_func(val1, val2))
        else:
            raise ValueError("func unsupported")

    _make_val_diagram(Bs, Ts, combined_vals, val_label, "{}_{}".format(out_prefix, val_label), cbar_format, clim_vals, cbar_tick_labels)

def _near_M_b_cutoff_plot(min_envs, out_prefix, env_val_1, env_val_2, env_val_label_1, env_val_label_2, delta_B, aspect=None):
    # Find B cutoff.
    max_val_B = None
    for this_env in min_envs:
        if this_env == None:
            continue
        this_Jbe = Jbe(this_env)
        Bratio = this_env["Bxy0"] / this_Jbe

        eps = 1e-6
        val_1 = this_env[env_val_1]
        val_2 = this_env[env_val_2]

        if abs(val_1) > eps or abs(val_2) > eps:
            if max_val_B == None or Bratio > max_val_B:
                max_val_B = Bratio

    # Collect (T, val_1) and (T, val_2) for desired Bs.
    B_val1_data, B_val2_data = {}, {}
    for this_env in min_envs:
        if this_env == None:
            continue
        this_Jbe = Jbe(this_env)
        Bratio = this_env["Bxy0"] / this_Jbe
        Tratio = (1.0 / this_env["Beta"]) / this_Jbe

        if Bratio > max_val_B or max_val_B - Bratio > delta_B:
            continue

        if Bratio not in B_val1_data:
            B_val1_data[Bratio] = []

        if Bratio not in B_val2_data:
            B_val2_data[Bratio] = []

        val_1 = this_env[env_val_1]
        val_2 = this_env[env_val_2]

        B_val1_data[Bratio].append([Tratio, val_1])
        B_val2_data[Bratio].append([Tratio, val_2])

    # Sort data and make plots.
    for B in B_val1_data.keys():
        Ts_1, vals_1, Ts_2, vals_2 = [], [], [], []

        val1_data = B_val1_data[B]
        val2_data = B_val2_data[B]

        sorted_val1 = sorted(val1_data, key=lambda x: x[0])
        sorted_val2 = sorted(val2_data, key=lambda x: x[0])

        for T, val in val1_data:
            Ts_1.append(T)
            vals_1.append(val)

        for T, val in val2_data:
            Ts_2.append(T)
            vals_2.append(val)

        if aspect != None:
            w, h = 10.0 * aspect[0], 10.0 * aspect[1]
            plt.figure(figsize=(w, h))

        plt.xlabel("$T/4J_{b}$", fontsize='x-large')

        plt.xlim(min(min(Ts_1), min(Ts_2)), max(max(Ts_1), max(Ts_2)))
        plt.ylim(min(min(vals_1), min(vals_2)), max(max(vals_1), max(vals_2)))

        plt.plot(Ts_1, vals_1, 'r-', label=env_val_label_1, linewidth=4)
        plt.plot(Ts_2, vals_2, 'k--', label=env_val_label_2, linewidth=4)

        plt.legend(loc=0, fontsize='x-large', title="$b_x/4J_b =$ {:.4f}".format(B))

        plt.savefig(out_prefix + '_Bxy_{:.4f}.png'.format(B), bbox_inches='tight', dpi=500)
        plt.savefig(out_prefix + '_Bxy_{:.4f}.eps'.format(B), bbox_inches='tight', dpi=500)
        plt.clf()

        if aspect != None:
            plt.close('all')

def _check_only_ok(Bratio, Tratio, only_B, only_T):
    if only_B != None and only_T != None:
        eps = 1e-9
        b_ok = abs(Bratio - only_B) < eps
        T_ok = abs(Tratio - only_T) < eps
        return b_ok and T_ok
    else:
        return True

def _set_M_avgs(min_envs):
    for env in min_envs:
        env["M1_avg"] = 0.5 * (env["M01"] + env["M11"])
        env["M2_avg"] = 0.5 * (env["M02"] + env["M12"])

def _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds=None):
    min_envs = None
    if read_prefix == None:
        min_envs, all_fenvs = min_envs_from_base(base_env_path, ions, num_Bs, num_Ts, npar, bounds)
        _save_min_envs(min_envs, out_prefix)
        _save_all_envs(all_fenvs, out_prefix)
    else:
        min_envs = _read_min_envs(read_prefix)

    return min_envs

def _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None):
    # m_{p \alpha} vs b, T plots
    M_plot_args = [[min_envs, out_prefix, "M01", "$|m_{0,1}|$", "abs", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "M11", "$|m_{1,1}|$", "abs", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "M02", "$|m_{0,2}|$", "abs", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "M12", "$|m_{1,2}|$", "abs", "%.2f", (0.0, 1.0)]]

    # m_{0 \alpha}, m_{1 \alpha} average vs b, T plots
    for env in min_envs:
        mode1_avg = (env["M01"] + env["M11"]) / 2.0
        mode2_avg = (env["M02"] + env["M12"]) / 2.0
        env["mode1_avg"] = abs(mode1_avg)
        env["mode2_avg"] = abs(mode2_avg)

    avg_plot_args = [[min_envs, out_prefix, "mode1_avg", "$|m_{1}|$", None, "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "mode2_avg", "$|m_{2}|$", None, "%.2f", (0.0, 1.0)]]

    # m phase diagram vs b, T plots
    M_phase_plot_args = [[min_envs, out_prefix + "_phase", "M01", "$|m_{0,1}| > 0$", "phase", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix + "_phase", "M11", "$|m_{1,1}| > 0$", "phase", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix + "_phase", "M02", "$|m_{0,2}| > 0$", "phase", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix + "_phase", "M12", "$|m_{1,2}| > 0$", "phase", "%.2f", (0.0, 1.0)]]

    # w vs b, T plots
    W_plot_args = [[min_envs, out_prefix, "W01", "$w_{0,1}$", None, "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "W11", "$w_{1,1}$", None, "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "W02", "$w_{0,2}$", None, "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "W12", "$w_{1,2}$", None, "%.2f", (0.0, 1.0)]]

    plot_args = []
    plot_args.extend(M_plot_args)
    plot_args.extend(avg_plot_args)
    plot_args.extend(M_phase_plot_args)
    plot_args.extend(W_plot_args)
    plot_args.append([min_envs, out_prefix, "FreeEnergy", "$F$", None, None])

    # electron-specific plots
    if not ions:
        plot_args.append([min_envs, out_prefix, "Mu", "$\mu$", None, None])

    for this_plot_args in plot_args:
        _make_BT_plot(*this_plot_args)

    # phase diagram in b, T: M1 or M2 or R
    phase_tick_labels = [[0.0, 0.5, 1.0], ["R", "M2", "M1"]]
    _multival_phase_plot(min_envs, out_prefix + "_phase_combine_M01_M02", "M01", "M02", "", func="phase_incl_m2", cbar_format="%.2f", clim_vals=(0.0, 1.0), cbar_tick_labels=phase_tick_labels)
    _multival_phase_plot(min_envs, out_prefix + "_phase_combine_M11_M12", "M11", "M12", "", func="phase_incl_m2", cbar_format="%.2f", clim_vals=(0.0, 1.0), cbar_tick_labels=phase_tick_labels)

    # m vs T for fixed b
    delta_B = 0.1
    _near_M_b_cutoff_plot(min_envs, out_prefix + "_M_T_p0", "M01", "M02", "$m_{0,1}$", "$m_{0,2}$", delta_B, aspect)
    _near_M_b_cutoff_plot(min_envs, out_prefix + "_M_T_p1", "M11", "M12", "$m_{1,1}$", "$m_{1,2}$", delta_B, aspect)

    _set_M_avgs(min_envs)
    _near_M_b_cutoff_plot(min_envs, out_prefix + "_M_T_avg", "M1_avg", "M2_avg", "$\\frac{1}{2}(m_{01}+m_{11})$", "$\\frac{1}{2}(m_{02}+m_{12})$", delta_B, aspect)

def _march16_plots():
    base_env_path = "march16_quad_plot.json"
    read_prefix = None
    out_prefix = "march16_quad_bT"
    ions = True
    num_Bs = 100
    num_Ts = 100
    npar = None
    plot_spectrum = False
    plot_dos = False
    only_B = None
    only_T = None
    bounds = [[0.01, 0.6], [0.01, 0.8]]

    # Jxz = 0, m vs b, T
    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Jxz = 0, m vs T at fixed b
    out_prefix = "march16_quad_fixed_b"
    aspect = [0.4, 1.0]
    bounds = [[0.01, 0.6], [0.01, 0.5]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect)

    # Jxz = 0, Jc = 0.2, m vs b, T
    base_env_path = "march16_quad_Jc_plot.json"
    out_prefix = "march16_quad_Jc_bT"
    bounds = [[0.01, 0.6], [0.01, 0.8]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Jxz = 0, Jc = 0.2, m vs T at fixed b
    out_prefix = "march16_quad_Jc_fixed_b"
    aspect = [0.4, 1.0]
    bounds = [[0.01, 0.6], [0.01, 0.5]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect)

    # Jxz = 0.1, m vs b, T
    base_env_path = "march16_quart_plot.json"
    out_prefix = "march16_quart_bT"
    bounds = [[0.01, 0.6], [0.01, 0.8]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Jxz = 0.1, m vs T at fixed b
    out_prefix = "march16_quart_fixed_b"
    bounds = [[0.01, 0.6], [0.01, 0.5]]
    aspect = [0.4, 1.0]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect)

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
    parser.add_argument('--march16', action='store_true', help="Make March16 plots")
    args = parser.parse_args()

    # TODO -- add only_B, only_T options to fix (B, T) point for generating
    # a few spectra/dos

    # TODO - don't assume run in own directory

    if not args.march16:
        min_envs = _get_min_envs(args.base_env_path, args.read_prefix, args.out_prefix, args.ions, args.num_Bs, args.num_Ts, args.npar)
        _make_plots(min_envs, args.out_prefix, args.ions, args.plot_spectrum, args.plot_dos, args.only_B, args.only_T)
        return
    else:
        _march16_plots()

if __name__ == "__main__":
    _main()
