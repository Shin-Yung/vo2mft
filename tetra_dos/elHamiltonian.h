#ifndef TETRA_DOS_ELHAMILTONIAN_H
#define TETRA_DOS_ELHAMILTONIAN_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "environment.h"

void ElHamiltonian_Recip(Environment *env, double k[3], gsl_matrix_complex *H);

void ElHamiltonian(Environment *env, double k[3], gsl_matrix_complex *H);

gsl_complex EpsilonAE(Environment *env, double k[3]);

gsl_complex EpsilonBE(Environment *env, double k[3]);

gsl_complex EpsilonAO(Environment *env, double k[3]);

gsl_complex EpsilonBO(Environment *env, double k[3]);

#endif // TETRA_DOS_ELHAMILTONIAN_H
