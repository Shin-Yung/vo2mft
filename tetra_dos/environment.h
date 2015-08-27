#ifndef TETRA_DOS_ENVIRONMENT
#define TETRA_DOS_ENVIRONMENT

#include <stdbool.h>

typedef struct {
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
	// On-site energies in M and R phases.
	double EpsilonM, EpsilonR;
} Environment;

#endif // TETRA_DOS_ENVIRONMENT
