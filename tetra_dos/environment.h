#ifndef TETRA_DOS_ENVIRONMENT
#define TETRA_DOS_ENVIRONMENT

#include <stdbool.h>

typedef struct {
	// Size of BZ on one edge (total number of BZ points is this cubed).
	int BZPointsPerDim;
	// Hopping parameters, even symmetry (a, c, diagonal axes).
	double Tae, Tce, Tbe;
	// Hopping parameters, odd symmetry (a, c, diagonal axes).
	double Tao, Tco, Tbo;
	// Order parameter <S>.
	double M;
	// Order parameter <S^2>.
	double W;
	// Chemical potential.
	double Mu;
	// Inverse temperature, 1 / (k_B * T).
	double Beta;
	// One-spin term for BEG model: coefficient for (S_i)^2.
	double B;
	// Exchange parameters for BEG model: coefficients to S_i dot S_j.
	// Jb is excluded since it does not contribute to results.
	double Ja, Jc;
	// Biquadratic exchange parameters for BEG model: coefficients to (S_i)^2 * (S_j)^2.
	double Ka, Kc, Kb;
	// On-site energies in M and R phases.
	double EpsilonM, EpsilonR;
	// Consider only ionic part of the problem:
	// only ions contribute to free energy; should solve
	// for (M, W).
	// If this is set to true, need to also set the following to 0:
	// Tae, Tce, Tbe, Tao, Tco, Tbo, EpsilonM, EpsilonR, Mu.
	// (maybe don't need to fix Mu = 0 -- large negative value could
	// be better).
	bool IonsOnly;
} Environment;

#endif // TETRA_DOS_ENVIRONMENT
