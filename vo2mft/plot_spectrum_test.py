import unittest
from vo2mft.solve import read_env_file, solve
from vo2mft.plot_spectrum import plot_spectrum

class PlotSpectrumTest(unittest.TestCase):
    def test_system_test_env(self):
        # TODO - don't assume test run in its directory.
        env = read_env_file("../vo2solve/system_test_env.json")
        eps = 1e-6
        env["EpsilonM"] = 0.05
        env["EpsilonR"] = 0.05
        #env["M"] = 0.0
        #env["W"] = 0.0

        mw_1_fenv = solve(env, eps)
        print("got env = {}".format(str(mw_1_fenv)))
        plot_spectrum(mw_1_fenv)


if __name__ == "__main__":
    unittest.main()
