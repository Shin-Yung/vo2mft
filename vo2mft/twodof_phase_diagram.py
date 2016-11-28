from argparse import ArgumentParser
from copy import deepcopy
from multiprocessing import Pool
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.font_manager import FontProperties
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

def min_envs_from_sample(sample_envs, ions, npar, initial_vals=None):
    eps = 1e-8
    sample_env_eps = []
    for env in sample_envs:
        twodof, twodof_body_indep = True, True
        sample_env_eps.append((env, eps, ions, twodof, twodof_body_indep, initial_vals))

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

def min_envs_from_base(base_path, ions, num_Bs, num_Ts, npar, bounds=None, initial_vals=None):
    base_env = read_env_file(base_path)
    sample = phase_sample(base_env, num_Bs, num_Ts, bounds)
    min_envs, all_fenvs = min_envs_from_sample(sample, ions, npar, initial_vals)
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
        #plt.savefig(out_prefix + '.eps', bbox_inches='tight', dpi=500)

    plt.clf()

def _make_BT_plot(min_envs, out_prefix, env_val, env_val_label, func=None, cbar_format=None, clim_vals=None):
    Bs, Ts, Ms = _collect_BT_env_val(min_envs, env_val, func)
    _make_val_diagram(Bs, Ts, Ms, env_val_label, "{}_{}".format(out_prefix, env_val), cbar_format, clim_vals)

def _multival_phase_plot(min_envs, out_prefix, env_val_1, env_val_2, val_label, func=None, cbar_format=None, clim_vals=None, cbar_tick_labels=None, env_val_1p=None, env_val_2p=None):
    Bs, Ts, val1s = _collect_BT_env_val(min_envs, env_val_1, func=None)
    Bs, Ts, val2s = _collect_BT_env_val(min_envs, env_val_2, func=None)
    # If env_val_1p and env_val_2p are specified, average env_val_1 with env_val_1p and
    # env_val_2 with env_val_2p.
    if env_val_1p != None:
        Bs, Ts, val1ps = _collect_BT_env_val(min_envs, env_val_1p, func=None)
    if env_val_2p != None:
        Bs, Ts, val2ps = _collect_BT_env_val(min_envs, env_val_2p, func=None)
    combined_vals = []

    val_set = list(zip(val1s, val2s))
    if env_val_1p != None and env_val_2p != None:
        valp_set = list(zip(val1ps, val2ps))

    for i in range(len(val1s)):
        val1, val2 = val_set[i]
        use1, use2 = val1, val2
        if env_val_1p != None and env_val_2p != None:
            val1p, val2p = valp_set[i]
            use1, use2 = 0.5 * (val1 + val1p), 0.5 * (val2 + val2p)

        if func == "phase_incl_m2":
            def phase_func(m1, m2):
                eps = 1e-6
                if abs(m1) > eps and abs(m2) > eps and abs(abs(m1) - abs(m2)) < eps:
                    return 1.0
                elif abs(m1) > eps and abs(m2) > eps:
                    #return 1.0 - 0.5 * abs(abs(m1) - abs(m2))
                    return 0.75
                elif (abs(m1) > eps and abs(m2) < eps) or (abs(m1) < eps and abs(m2) > eps):
                    return 0.5
                else:
                    return 0.0
            combined_vals.append(phase_func(use1, use2))
        else:
            raise ValueError("func unsupported")

    _make_val_diagram(Bs, Ts, combined_vals, val_label, "{}_{}".format(out_prefix, val_label), cbar_format, clim_vals, cbar_tick_labels)

def _find_max_Bx(min_envs, env_val_A, env_val_B):
    # max_val_B will be set to the largest Bratio in which val_A or val_B is
    # nonzero for either min_envs_1 or min_envs_2.
    max_val_Bx = None
    for this_env in min_envs:
        if this_env == None:
            continue
        this_Jbe = Jbe(this_env)
        Bratio = this_env["Bxy0"] / this_Jbe

        eps = 1e-6
        val_A = this_env[env_val_A]
        val_B = this_env[env_val_B]

        if abs(val_A) > eps or abs(val_B) > eps:
            if max_val_Bx == None or Bratio > max_val_Bx:
                max_val_Bx = Bratio

    return max_val_Bx

def _collect_constB_data(min_envs, max_val_Bx, env_val_A, env_val_B, delta_Bx):
    B_valA_data, B_valB_data = {}, {}
    for this_env in min_envs:
        if this_env == None:
            continue
        this_Jbe = Jbe(this_env)
        Bratio = this_env["Bxy0"] / this_Jbe
        Tratio = (1.0 / this_env["Beta"]) / this_Jbe

        if Bratio > max_val_Bx or max_val_Bx - Bratio > delta_Bx:
            continue

        if Bratio not in B_valA_data:
            B_valA_data[Bratio] = []

        if Bratio not in B_valB_data:
            B_valB_data[Bratio] = []

        val_A = abs(this_env[env_val_A])
        val_B = abs(this_env[env_val_B])

        B_valA_data[Bratio].append([Tratio, val_A])
        B_valB_data[Bratio].append([Tratio, val_B])

    return B_valA_data, B_valB_data

def _sort_Bdata(B, B_valA_data, B_valB_data):
    Ts_A, vals_A, Ts_B, vals_B = [], [], [], []

    this_valA_data = B_valA_data[B]
    this_valB_data = B_valB_data[B]

    sorted_valA = sorted(this_valA_data, key=lambda x: x[0])
    sorted_valB = sorted(this_valB_data, key=lambda x: x[0])

    for T, val in sorted_valA:
        Ts_A.append(T)
        vals_A.append(val)

    for T, val in sorted_valB:
        Ts_B.append(T)
        vals_B.append(val)

    return Ts_A, vals_A, Ts_B, vals_B

def _all_min(xs_list):
    min_val = None
    for xs in xs_list:
        this_min = min(xs)
        if min_val == None or this_min < min_val:
            min_val = this_min

    return min_val

def _all_max(xs_list):
    max_val = None
    for xs in xs_list:
        this_max = max(xs)
        if max_val == None or this_max > max_val:
            max_val = this_max

    return max_val

def _near_M_b_cutoff_plot(min_envs, out_prefix, env_val_1, env_val_2, env_val_label_1, env_val_label_2, delta_B, aspect=None):
    # Find B cutoff.
    max_val_B = _find_max_Bx(min_envs, env_val_1, env_val_2)

    # Collect (T, val_1) and (T, val_2) for desired Bs.
    B_val1_data, B_val2_data = _collect_constB_data(min_envs, max_val_B,
            env_val_1, env_val_2, delta_B)

    # Sort data and make plots.
    for B in B_val1_data.keys():
        Ts_1, vals_1, Ts_2, vals_2 = _sort_Bdata(B, B_val1_data, B_val2_data)

        if aspect != None:
            w, h = 10.0 * aspect[0], 10.0 * aspect[1]
            plt.figure(figsize=(w, h))

        plt.xlabel("$T/4J_{b}$", fontsize='x-large')

        min_x, max_x = _all_min([Ts_1, Ts_2]), _all_max([Ts_1, Ts_2])
        min_y, max_y = _all_min([vals_1, vals_2]), _all_max([vals_1, vals_2])

        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)

        plt.plot(Ts_1, vals_1, 'r-', label=env_val_label_1, linewidth=4)
        plt.plot(Ts_2, vals_2, 'k--', label=env_val_label_2, linewidth=4)

        legend = plt.legend(loc=0, fontsize='x-large', title="$b_x/4J_b =$ {:.4f}".format(B))
        legend.get_title().set_fontsize('x-large')

        plt.savefig(out_prefix + '_Bxy_{:.4f}.png'.format(B), bbox_inches='tight', dpi=500)
        #plt.savefig(out_prefix + '_Bxy_{:.4f}.eps'.format(B), bbox_inches='tight', dpi=500)
        plt.clf()

        if aspect != None:
            plt.close('all')

def find_k_nearest(kt, d):
    '''Return the key of the dict d which is closest to kt.
    '''
    nearest = None
    for k in d.keys():
        if nearest is None or abs(kt - k) < abs(kt - nearest):
            nearest = k

    return nearest

def _near_M_b_cutoff_plot_multiple(min_envs_1, min_envs_2, out_prefix, env_val_A, env_val_B, env_val_labels_A, env_val_labels_B, delta_Bx, aspect=None, plot_all_Bs=False, fixed_xticks=None):
    # Find b cutoff.
    max_val_Bx = max(_find_max_Bx(min_envs_1, env_val_A, env_val_B),
            _find_max_Bx(min_envs_2, env_val_A, env_val_B))

    # Collect (T, val_A) and (T, val_B) for desired bs.
    B_valA_1_data, B_valB_1_data = _collect_constB_data(min_envs_1, max_val_Bx,
            env_val_A, env_val_B, delta_Bx)
    B_valA_2_data, B_valB_2_data = _collect_constB_data(min_envs_2, max_val_Bx,
        env_val_A, env_val_B, delta_Bx)

    # NOTE - Will assume each data set has the same set of B's
    # (that is, for each set at least 1 T exists at each b).
    if plot_all_Bs:
        Bs_to_plot = B_valA_1_data.keys()
        all_show_legend_labels = [True]*len(keys)
    else:
        near_keys = [0.4510, 0.4748, 0.4987]
        Bs_to_plot = [find_k_nearest(kt, B_valA_1_data) for kt in near_keys]
        all_show_legend_labels = [True, False, False]

    if aspect != None:
        w, h = 10.5 * aspect[0] * len(Bs_to_plot), 10.5 * aspect[1]
        fig = plt.figure(figsize=(w, h))
    else:
        fig = plt.figure()

    plt.subplots_adjust(wspace=0.001)

    x_limits, y_limits = None, None
    axes = []
    # Make plots vs T at fixed b.
    for B_index, (B, show_legend_labels) in enumerate(zip(Bs_to_plot, all_show_legend_labels)):
        Ts_A_1, vals_A_1, Ts_B_1, vals_B_1 = _sort_Bdata(B, B_valA_1_data, B_valB_1_data)
        Ts_A_2, vals_A_2, Ts_B_2, vals_B_2 = _sort_Bdata(B, B_valA_2_data, B_valB_2_data)

        ax = fig.add_subplot(1, len(Bs_to_plot), B_index+1)
        axes.append(ax)

        ax.set_xlabel("$T/4J_{b}$", fontsize='x-large')

        if B_index == 0:
            min_x = _all_min([Ts_A_1, Ts_B_1, Ts_A_2, Ts_B_2])
            max_x = _all_max([Ts_A_1, Ts_B_1, Ts_A_2, Ts_B_2])
            min_y = _all_min([vals_A_1, vals_B_1, vals_A_2, vals_B_2])
            max_y = _all_max([vals_A_1, vals_B_1, vals_A_2, vals_B_2])
            x_limits = (min_x, max_x)
            y_limits = (min_y, max_y)

            ax.set_xlim(min_x, max_x)
            ax.set_ylim(min_y, max_y)
        else:
            ax.set_xlim(x_limits[0], x_limits[1])
            ax.set_ylim(y_limits[0], y_limits[1])

        if fixed_xticks is not None:
            ax.set_xticks(fixed_xticks)

        label_A1, label_A2 = env_val_labels_A
        label_B1, label_B2 = env_val_labels_B

        if show_legend_labels:
            ax.plot(Ts_A_1, vals_A_1, 'k-', label=label_A1, linewidth=6)                            # solid black, m_A(1)
            ax.plot(Ts_B_1, vals_B_1, color='limegreen', linestyle='-', label=label_B1, linewidth=2)     # solid gray, m_B(1)
            ax.plot(Ts_A_2, vals_A_2, 'k--', label=label_A2, linewidth=6)                           # dashed black, m_A(2)
            ax.plot(Ts_B_2, vals_B_2, color='limegreen', linestyle='--', label=label_B2, linewidth=2)    # dashed gray, m_B(2)

            legend = ax.legend(loc=0, fontsize='x-large', title="$b_x/4J_b = {:.4f}$".format(B))
            legend.get_title().set_fontsize('x-large')
        else:
            ax.plot(Ts_A_1, vals_A_1, 'k-', linewidth=6)                            # solid black, m_A(1)
            ax.plot(Ts_B_1, vals_B_1, color='limegreen', linestyle='-', linewidth=2)     # solid gray, m_B(1)
            ax.plot(Ts_A_2, vals_A_2, 'k--', linewidth=6)                           # dashed black, m_A(2)
            ax.plot(Ts_B_2, vals_B_2, color='limegreen', linestyle='--', linewidth=2)    # dashed gray, m_B(2)

            abox = AnchoredText("$b_x/4J_b =$ {:.4f}".format(B), loc=3,
                    prop={"size": "x-large"}, frameon=True)
            ax.add_artist(abox)

    yticklabels = None
    for ax in axes[1:]:
        if yticklabels is None:
            yticklabels = ax.get_yticklabels()
        else:
            yticklabels = yticklabels + ax.get_yticklabels()

    plt.setp(yticklabels, visible=False)

    plt.savefig(out_prefix + '_Bxy_range.png', bbox_inches='tight', dpi=500)
    plt.savefig(out_prefix + '_Bxy_range.eps', bbox_inches='tight', dpi=500)
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
        env["M1_avg"] = 0.5 * (abs(env["M01"]) + abs(env["M11"]))
        env["M2_avg"] = 0.5 * (abs(env["M02"]) + abs(env["M12"]))

def _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds=None, initial_vals=None):
    min_envs = None
    if read_prefix == None:
        min_envs, all_fenvs = min_envs_from_base(base_env_path, ions, num_Bs, num_Ts, npar, bounds, initial_vals)
        _save_min_envs(min_envs, out_prefix)
        _save_all_envs(all_fenvs, out_prefix)
    else:
        min_envs = _read_min_envs(read_prefix)

    return min_envs

def _add_averages(min_envs):
    for env in min_envs:
        mode1_avg = (env["M01"] + env["M11"]) / 2.0
        mode2_avg = (env["M02"] + env["M12"]) / 2.0
        env["mode1_avg"] = abs(mode1_avg)
        env["mode2_avg"] = abs(mode2_avg)

def _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None):
    # m_{p \alpha} vs b, T plots
    M_plot_args = [[min_envs, out_prefix, "M01", "$|m_{0,1}|$", "abs", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "M11", "$|m_{1,1}|$", "abs", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "M02", "$|m_{0,2}|$", "abs", "%.2f", (0.0, 1.0)],
        [min_envs, out_prefix, "M12", "$|m_{1,2}|$", "abs", "%.2f", (0.0, 1.0)]]

    # m_{0 \alpha}, m_{1 \alpha} average vs b, T plots
    _add_averages(min_envs)

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
    phase_tick_labels = [[0.0, 0.5, 0.75, 1.0], ["R", "M2", "T-like", "M1"]]
    _multival_phase_plot(min_envs, out_prefix + "_phase_combine_M01_M02", "M01", "M02", "", func="phase_incl_m2", cbar_format="%.2f", clim_vals=(0.0, 1.0), cbar_tick_labels=phase_tick_labels)
    _multival_phase_plot(min_envs, out_prefix + "_phase_combine_M11_M12", "M11", "M12", "", func="phase_incl_m2", cbar_format="%.2f", clim_vals=(0.0, 1.0), cbar_tick_labels=phase_tick_labels)
    _multival_phase_plot(min_envs, out_prefix + "_phase_combine_avg", "M01", "M02", "", func="phase_incl_m2", cbar_format="%.2f", clim_vals=(0.0, 1.0), cbar_tick_labels=phase_tick_labels, env_val_1p="M11", env_val_2p="M12")

    # m vs T for fixed b
    delta_B = 0.1
    _near_M_b_cutoff_plot(min_envs, out_prefix + "_M_T_p0", "M01", "M02", "$|m_{0,1}|$", "$|m_{0,2}|$", delta_B, aspect)
    _near_M_b_cutoff_plot(min_envs, out_prefix + "_M_T_p1", "M11", "M12", "$|m_{1,1}|$", "$|m_{1,2}|$", delta_B, aspect)

    _set_M_avgs(min_envs)
    _near_M_b_cutoff_plot(min_envs, out_prefix + "_M_T_avg", "M1_avg", "M2_avg", "$\\frac{1}{2}(|m_{0,1}|+|m_{1,1}|)$", "$\\frac{1}{2}(|m_{0,2}|+|m_{1,2}|)$", delta_B, aspect)

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

    # Bxz = 0, m vs b, T
    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Bxz = 0, m vs T at fixed b
    out_prefix = "march16_quad_fixed_b"
    aspect = [0.4, 1.0]
    bounds = [[0.01, 0.6], [0.01, 0.5]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect)

    # Bxz = 0, Jc = 0, Fac_xy = 0.03, Poisson = 0, m vs b, T
    base_env_path = "march16_quad_Fxy_plot.json"
    out_prefix = "march16_quad_Fxy"
    bounds = [[0.01, 0.6], [0.01, 0.8]]
    initial_vals = "mode_symmetric"

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds, initial_vals)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Bxz = 0, Jc = 0, Fac_xy = 0.03, Poisson = 0, m vs T at fixed b
    base_env_path = "march16_quad_Fxy_plot.json"
    out_prefix = "march16_quad_Fxy_fixed_b"
    aspect = [0.5, 1.0]
    bounds = [[0.01, 0.6], [0.01, 0.5]]
    initial_vals = "mode_symmetric"

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds, initial_vals)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect)

    # Bxz = 0, Jc = 0, Fac_xy = 0.03, Poisson = 0.3, m vs b, T
    base_env_path = "march16_quad_Fxy_Poisson_plot.json"
    out_prefix = "march16_quad_Fxy_Poisson"
    bounds = [[0.01, 0.6], [0.01, 0.8]]
    initial_vals = "mode_symmetric"

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds, initial_vals)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Bxz = 0.1, Jc = 0, Fac_xy = 0.03, Poisson = 0, m vs b, T
    base_env_path = "march16_qq_Fxy_plot.json"
    out_prefix = "march16_qq_Fxy"
    bounds = [[0.01, 0.6], [0.01, 0.8]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds, initial_vals)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Bxz = 0, Jc = 0.2, m vs b, T
    base_env_path = "march16_quad_Jc_plot.json"
    out_prefix = "march16_quad_Jc_bT"
    bounds = [[0.01, 0.6], [0.01, 0.8]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Bxz = 0, Jc = 0.2, m vs T at fixed b
    out_prefix = "march16_quad_Jc_fixed_b"
    aspect = [0.4, 1.0]
    bounds = [[0.01, 0.6], [0.01, 0.5]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect)

    # Bxz = 0.1, m vs b, T
    base_env_path = "march16_quart_plot.json"
    out_prefix = "march16_quart_bT"
    bounds = [[0.01, 0.6], [0.01, 0.8]]

    min_envs = _get_min_envs(base_env_path, read_prefix, out_prefix, ions, num_Bs, num_Ts, npar, bounds)
    _make_plots(min_envs, out_prefix, ions, plot_spectrum, plot_dos, only_B, only_T, aspect=None)

    # Bxz = 0.1, m vs T at fixed b
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
    parser.add_argument('--multi_b_cutoff', action='store_true', help="Make b cutoff plots for multiple runs")
    args = parser.parse_args()

    # TODO -- add only_B, only_T options to fix (B, T) point for generating
    # a few spectra/dos

    # TODO - don't assume run in own directory

    if args.march16:
        _march16_plots()
    elif args.multi_b_cutoff:
        num_Bs, num_Ts = 100, 100
        bounds = [[0.01, 0.6], [0.01, 0.5]]
        ions = True
        npar = None

        # Check if runs have been performed.
        # If not, perform them.
        out_prefix = "run_multi_b_cutoff"
        out_prefix_bxz = "{}_bxz".format(out_prefix)

        if not os.path.exists(_min_env_filename(out_prefix_bxz)):
            read_prefix = None
            base_env_path_bxz = "multi_b_cutoff_bxz_env.json"
            min_envs_bxz = _get_min_envs(base_env_path_bxz, read_prefix, out_prefix_bxz,
                    ions, num_Bs, num_Ts, npar, bounds)
        else:
            min_envs_bxz = _read_min_envs(out_prefix_bxz)

        out_prefix_F = "{}_F".format(out_prefix)
        if not os.path.exists(_min_env_filename(out_prefix_F)):
            read_prefix = None
            base_env_path_F = "multi_b_cutoff_F_env.json"
            initial_vals = "mode_symmetric"

            min_envs_F = _get_min_envs(base_env_path_F, read_prefix, out_prefix_F,
                    ions, num_Bs, num_Ts, npar, bounds, initial_vals)
        else:
            min_envs_F = _read_min_envs(out_prefix_F)

        _add_averages(min_envs_bxz)
        _add_averages(min_envs_F)

        delta_Bx = 0.1
        aspect = [0.4, 1.0]
        # mode2-m_A and mode1-m_B correspondence reflects change between original formulation
        # and paper.
        labels_A = ["$m_A$ ($b_{xz}$)", "$m_A$ ($F_{xy}$)"]
        labels_B = ["$m_B$ ($b_{xz}$)", "$m_B$ ($F_{xy}$)"]
        _near_M_b_cutoff_plot_multiple(min_envs_bxz, min_envs_F, out_prefix,
                "mode2_avg", "mode1_avg", labels_A, labels_B, delta_Bx, aspect,
                fixed_xticks=[0.05, 0.15, 0.25, 0.35, 0.45])
    else:
        min_envs = _get_min_envs(args.base_env_path, args.read_prefix, args.out_prefix, args.ions, args.num_Bs, args.num_Ts, args.npar)
        _make_plots(min_envs, args.out_prefix, args.ions, args.plot_spectrum, args.plot_dos, args.only_B, args.only_T)

if __name__ == "__main__":
    _main()
