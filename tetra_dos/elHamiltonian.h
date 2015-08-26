#ifndef TETRA_DOS_ELHAMILTONIAN_H
#define TETRA_DOS_ELHAMILTONIAN_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "environment.h"

void ElHamiltonian_Recip(Environment *env, gsl_vector *k, gsl_matrix_complex *H);

void ElHamiltonian(Environment *env, gsl_vector *k, gsl_matrix_complex *H);

gsl_complex EpsilonAE(Environment *env, gsl_vector *k);

gsl_complex EpsilonBE(Environment *env, gsl_vector *k);

gsl_complex EpsilonAO(Environment *env, gsl_vector *k);

gsl_complex EpsilonBO(Environment *env, gsl_vector *k);

#endif // TETRA_DOS_ELHAMILTONIAN_H
