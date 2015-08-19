import numpy as np

def _sc_kpath(alat):
    # k-path runs over simple cubic Brillouin zone:
    # G-X-M-G-R-X-M-R.

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

    return [Gc, Xc, Mc, Gc, Rc, Xc, Mc, Rc]

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

def plot_spectrum(env):
    alat = 1.0
    kpoints_per_panel = 3
    kpath = _interpolate_kpoints(_sc_kpath(alat), kpoints_per_panel)
    for k in kpath:
        print(k)
