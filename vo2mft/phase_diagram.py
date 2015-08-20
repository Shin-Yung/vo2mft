from copy import deepcopy
import numpy as np
import vo2mft.environment

def phase_sample(base_env):
    '''Generate a set of envs over which the phase diagram will be sampled,
    with constant properties taken from base_env.
    '''
    # Generate set of (B, T) values.
    # TODO - generalize away from EpsilonM - EpsilonR = 0.
    QJ_ion = environment.QJ_ion(base_env)
    Bratio_start, Bratio_stop = 0.01, 0.75
    Tratio_start, Tratio_stop = 0.1, 10.0
    num_Bs, num_Ts = 20, 20

    Bratios = np.linspace(Bratio_start, Bratio_stop, num_Bs)
    Tratios = np.linspace(Tratio_start, Tratio_stop, num_Ts)

    B_T_vals = []
    for Br in Bratios:
        for Tr in Tratios:
            B_val = Br*QJ_ion
            T_val = Tr*QJ_ion
            B_T_vals.append([B_val, T_val])

    sample_envs = []
    for B, T in B_T_vals:
        this_env = deepcopy(base_env)
        this_env["B"] = B
        this_env["Beta"] = 1.0/T
        sample_envs.append(this_envs)

    return sample_envs
