import subprocess
import os
import json
from uuid import uuid4
from vo2mft.util import _solve_front_path, _twodof_solve_front_path, _twodof_body_fixed_solve_front_path

def solve(env, eps=1e-6, ions=False, flags=None, twodof=False, twodof_body_indep=False):
    '''Return the solved final env corresponding to the given env, solved to
    accuracy given by eps.
    '''
    solver_path = _solve_front_path()
    if twodof and twodof_body_indep:
        solver_path = _twodof_solve_front_path()
    elif twodof:
        solver_path = _twodof_body_fixed_solve_front_path()

    in_path, out_path = str(uuid4()), str(uuid4())

    write_env_file(env, in_path)

    # Run solver.
    solver_call = None
    if twodof:
        solver_call = [solver_path, "--eps", str(eps)]
        if ions:
            solver_call.append("--ions")
        if flags != None:
            solver_call.extend(flags)
        solver_call.extend([in_path, out_path])
    elif ions and not twodof:
        solver_call = [solver_path, "--eps", str(eps), "--ions"]
        if flags != None:
            solver_call.extend(flags)
        solver_call.extend([in_path, out_path])
    subprocess.call(solver_call)

    # Read solver output, if it exists.
    final_env_path = out_path + "_fenv.json"
    final_env = None
    try:
        final_env = read_env_file(final_env_path)
    except FileNotFoundError:
        pass

    # Clean up solver input/output.
    try:
        os.remove(in_path)
        os.remove(final_env_path)
    except FileNotFoundError:
        pass

    return final_env

def solve_set(envs, eps=1e-6, ions=False, flags=None, twodof=False, twodof_body_indep=False):
    '''Return a list of solved final envs corresponding to the given list of
    envs, solved to accuracy given by eps.

    The set of envs is solved serially (only one process is invoked).
    '''
    final_envs = []
    for i, initial_env in enumerate(envs):
        this_flags = None
        if flags != None:
            this_flags = flags[i]
        this_final_env = solve(initial_env, eps, ions, this_flags, twodof, twodof_body_indep)
        final_envs.append(this_final_env)
    return final_envs

def write_env_file(env, env_path):
    env_str = json.dumps(env)
    with open(env_path, 'w') as fp:
        fp.write(env_str)

def read_env_file(env_path):
    env = None
    with open(env_path, 'r') as fp:
        env_str = fp.read()
        env = json.loads(env_str)
    return env
