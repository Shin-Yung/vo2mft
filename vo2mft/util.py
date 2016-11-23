import os
import inspect

def _base_dir():
    frame = inspect.getfile(inspect.currentframe())
    this_dir = os.path.dirname(os.path.abspath(frame))
    base_dir = os.path.join(this_dir, "..")
    return os.path.normpath(base_dir)

def _solve_front_path():
    return os.path.join(_base_dir(), "vo2solve", "vo2solve_front", "vo2solve_front")

def _twodof_solve_front_path():
    return os.path.join(_base_dir(), "twodof", "vo2solve_front", "vo2solve_front")

def _twodof_body_fixed_solve_front_path():
    return os.path.join(_base_dir(), "twodofavg", "vo2solve_front", "vo2solve_front")

def _run_dos_path():
    return os.path.join(_base_dir(), "tetra_dos", "RunDosValues.out")
