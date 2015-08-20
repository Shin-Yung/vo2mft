import unittest
from solve import read_env_file
from min_free_energy import minimize_free_energy

class SolveEnvTest(unittest.TestCase):
    def test_system_test_env(self):
        # TODO - don't assume test run in its directory.
        env = read_env_file("min_free_energy_test_env.json")
        eps = 1e-6
        env["EpsilonR"] = 0.05
        env["EpsilonM"] = 0.05

        M_expected_beta_0_1 = 0.0
        M_expected_beta_10 = 0.9999999999971557

        env["Beta"] = 0.1
        min_env_beta_0_1 = minimize_free_energy(env, eps)
        print("Minimum free energy env with beta = 0.1:", min_env_beta_0_1)
        self.assertTrue(abs(M_expected_beta_0_1 - min_env_beta_0_1["M"]) < 2.0*eps)

        env["Beta"] = 10.0
        min_env_beta_10 = minimize_free_energy(env, eps)
        print("Minimum free energy env with beta = 10.0:", min_env_beta_10)
        self.assertTrue(abs(M_expected_beta_10 - min_env_beta_10["M"]) < 2.0*eps)

if __name__ == "__main__":
    unittest.main()
