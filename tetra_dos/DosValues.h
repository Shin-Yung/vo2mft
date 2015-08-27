#ifndef TETRA_DOS_DOS_VALUES_H
#define TETRA_DOS_DOS_VALUES_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "ctetra/dos.h"
#include "elHamiltonian.h"
#include "environment.h"

gsl_matrix* cubicRecipLat(double a);

double* DosValues(Environment *env, int n, double *Es, double num_dos);

#endif // TETRA_DOS_DOS_VALUES_H
