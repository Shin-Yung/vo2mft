import numpy as np

def _cubic_R(a):
    '''Reciprocal lattice for simple cubic direct lattice of lattice
    constant a.
    '''
    D = np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    R = 2.0 * np.pi * np.linalg.inv(D)
    return R
