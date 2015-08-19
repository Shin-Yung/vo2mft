import unittest
from solve import read_env_file
from min_free_energy import minimize_free_energy

class SolveEnvTest(unittest.TestCase):
    def test_system_test_env(self):
        # TODO - don't assume test run in its directory.
        env = read_env_file("../vo2solve/system_test_env.json")
        eps = 1e-6
        min_env = minimize_free_energy(env, eps)
        print(min_env)

if __name__ == "__main__":
    unittest.main()
