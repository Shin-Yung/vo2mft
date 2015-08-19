import numpy as np
import matplotlib.pyplot as plt
from vo2mft.elHamiltonian import ElHamiltonian

def _sc_kpath(alat):
    # k-path runs over simple cubic Brillouin zone:
    # G-X-M-G-R-X-M-R.
    labels = ["$\\Gamma$", "$X$", "$M$", "$\\Gamma$", "$R$", "$X$", "$M$", "$R$"]

    # High-symmetry points in reciprocal lattice coordinates.
    G = (0.0, 0.0, 0.0)
    M = (0.5, 0.5, 0.0)
    R = (0.5, 0.5, 0.5)
    X = (0.0, 0.5, 0.0)
    # Columns of Dlat = direct lattice vectors.
    Dlat = np.array([[alat, 0.0, 0.0], [0.0, alat, 0.0], [0.0, 0.0, alat]])
    # Rows of Rlat = reciprocal lattice vectors.
    Rlat = 2.0 * np.pi * np.linalg.inv(Dlat)

    # High-symmetry points in Cartesian coordinates.
    Gc = np.dot(G, Rlat)
    Mc = np.dot(M, Rlat)
    Rc = np.dot(R, Rlat)
    Xc = np.dot(X, Rlat)

    kpath = (Gc, Xc, Mc, Gc, Rc, Xc, Mc, Rc)
    return kpath, labels

def _interpolate_kpoints(kpath, kpoints_per_panel):
    interpolated = []
    for panel_index in range(len(kpath)):
        # Take panel from kpath[i-1] to kpath[i].
        if panel_index == 0:
            continue
        kstart = kpath[panel_index-1]
        kstop = kpath[panel_index]
        step = np.subtract(kstop, kstart) / (kpoints_per_panel - 1)
        k = kstart
        for k_index in range(kpoints_per_panel):
            # To avoid doubling k-points at panel ends, take panel start point
            # only on first panel.
            if k_index == 0:
                if panel_index == 1:
                    interpolated.append(k)
                continue
            # Past panel start point.
            k = np.add(k, step)
            interpolated.append(k)

    return interpolated

def _collect_ys(env, kpath):
    ys = None
    for k in kpath:
        H = ElHamiltonian(env, k)
        evals = sorted(np.linalg.eigvalsh(H))
        # Initialize ys to len(evals)-length list of lists.
        if ys == None:
            ys = []
            for i in range(len(evals)):
                ys.append([])
        # For each eval, put it in the corresponding list in ys based on its
        # sorted index.
        for i, this_eval in enumerate(evals):
            ys[i].append(this_eval)

    return ys

def plot_spectrum(env, plot_filename=None):
    alat = 1.0
    kpoints_per_panel = 50

    kpath, labels = _sc_kpath(alat)
    all_ks = _interpolate_kpoints(kpath, kpoints_per_panel)
    xs = range(len(all_ks)) # TODO - scale by distance between kpoints
    ys = _collect_ys(env, all_ks)

    # Set plot boundaries.
    plt.xlim(0, xs[-1])

    # Set symmetry point axis markers/lines.
    sym_xs = [0]
    for i in range(len(labels)):
        if i == 0:
            continue
        sym_xs.append(sym_xs[-1] + kpoints_per_panel - 1)

    for x in sym_xs:
        plt.axvline(x, color='k')
    plt.xticks(sym_xs, labels)

    # Plot data.
    for y_set in ys:
        plt.plot(xs, y_set, 'r')

    # Show or save plot.
    if plot_filename == None:
        plt.show()
    else:
        plt.savefig(plot_filename + '.png', bbox_inches='tight', dpi=500)
    plt.clf()
