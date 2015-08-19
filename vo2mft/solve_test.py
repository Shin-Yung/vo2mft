import unittest
from solve import read_env_file, solve

class SolveEnvTest(unittest.TestCase):
    def test_system_test_env(self):
        # TODO - don't assume test run in its directory.
        env = read_env_file("../vo2solve/system_test_env.json")
        eps = 1e-6
        final_env = solve(env, eps)
        print(final_env)

if __name__ == "__main__":
    unittest.main()
