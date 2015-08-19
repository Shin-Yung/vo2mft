import numpy as np

def ElHamiltonian(env, k):
    '''Calculate 4x4 electronic Hamiltonian.
    Assumes that k is scaled such that kx*a, ky*a, kz*c --> kx, ky, kz.
    '''
    # KQ = k + Q
    KQ = [k[0] + np.pi, k[1] + np.pi, k[2] + np.pi]

    EpsAE = EpsilonAE(env, k)
    EpsBE = EpsilonBE(env, k)
    EpsBE_KQ = EpsilonBE(env, KQ)
    EpsAO = EpsilonAO(env, k)
    EpsBO = EpsilonBO(env, k)
    ikd = complex(0.0, k[0]/2.0+k[1]/2.0+k[2]/2.0)
    mu = complex(env["Mu"], 0.0)

    H = np.zeros((4, 4), dtype=np.complex128)
    H[0, 0] = EpsAE - mu
    H[1, 0] = 2.0 * EpsAO
    H[2, 0] = EpsBE * np.exp(-ikd)
    H[3, 0] = -cmplx.Conj(EpsBO) * np.exp(-ikd)

    H[0, 1] = -2.0 * EpsAO
    H[1, 1] = -EpsAE - mu
    H[2, 1] = cmplx.Conj(EpsBO) * np.exp(-ikd)
    H[3, 1] = complex(0.0, 1.0) * EpsBE_KQ * np.exp(-ikd)

    H[0, 2] = EpsBE * np.exp(ikd)
    H[1, 2] = EpsBO * np.exp(ikd)
    H[2, 2] = EpsAE - mu
    H[3, 2] = 2.0 * EpsAO

    H[0, 3] = -EpsBO * np.exp(ikd)
    H[1, 3] = complex(0.0, -1.0) * EpsBE_KQ * np.exp(ikd)
    H[2, 3] = -2.0 * EpsAO
    H[3, 3] = -EpsAE - mu

    return H

def EpsilonAE(env, k):
    '''Cubic axes, even symmetry (k, p; k, p)
    '''
    rp = -2.0 * (env["Tae"]*(np.cos(k[0])+np.cos(k[1])) + env["Tce"]*np.cos(k[2]))
    return complex(rp, 0.0)

def EpsilonBE(env, k):
    '''Body diagonal, even symmetry (k, p; k, pbar)
    '''
    rp = -8.0 * env["Tbe"] * np.cos(k[0]/2.0) * np.cos(k[1]/2.0) * np.cos(k[2]/2.0)
    return complex(rp, 0.0)

def EpsilonAO(env, k):
    '''Cubic axes, odd symmetry (k, p; k+Q, p)
    '''
    ip = -2.0 * env["M"] * (env["Tao"]*(np.sin(k[0])+np.sin(k[1])) + env["Tco"]*np.sin(k[2]))
    return complex(0.0, ip)

def EpsilonBO(env, k):
    '''Body diagonal, odd symmetry (k, p; k+Q, pbar)
    '''
    rp = -8.0 * env["M"] * env["Tbo"] * np.cos(k[0]/2.0) * np.cos(k[1]/2.0) * np.cos(k[2]/2.0)
    ip = 8.0 * env["M"] * env["Tbo"] * np.sin(k[0]/2.0) * np.sin(k[1]/2.0) * np.sin(k[2]/2.0)
    return complex(rp, ip)
