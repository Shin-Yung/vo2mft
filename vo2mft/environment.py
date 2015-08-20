# Functions that operate on environments and final environments.
# Those that take 'env' operate on either;
# those that that 'fenv' operate only on final environments.
def QJ(fenv):
    return QJ_ion(fenv) + QJ_el(fenv)

def QJ_ion(env):
    return 4.0*env["Ja"] + 2.0*env["Jc"]

def QJ_el(fenv):
    return 4.0*fenv["Tao"]*fenv["Dao"] + 2.0*fenv["Tco"]*fenv["Dco"]

def DeltaS(env):
    # TODO - generalize away from n = 1?
    n = 1.0
    return env["B"] + (env["EpsilonM"] - env["EpsilonR"])*n
