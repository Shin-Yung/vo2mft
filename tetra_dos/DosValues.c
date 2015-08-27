#include "DosValues.h"

gsl_matrix* cubicRecipLat(double a) {
    gsl_matrix *R = gsl_matrix_calloc(3, 3);
    int i;
    for (i = 0; i < 3; i++) {
        gsl_matrix_set(R, i, i, 2.0*M_PI/a);
    }
    return R;
}

double* DosValues(Environment *env, int n, double *Es, double num_dos) {
    // Set up memory.
    int num_bands = 4;
    gsl_matrix_complex *Hk = gsl_matrix_complex_calloc(num_bands, num_bands);
    gsl_eigen_herm_workspace *work = gsl_eigen_herm_alloc(num_bands);

    // Reciprocal lattice: use simple cubic, a = 1.
    gsl_matrix *R = cubicRecipLat(1.0);

    // Use GCC nested function to make closure around H(env, k).
    // https://gcc.gnu.org/onlinedocs/gcc/Nested-Functions.html
    // Efn puts the eigenvalues of H(k) into energies.
    void Efn(double k[3], gsl_vector *energies) {
        // Set Hk = H(k).
        ElHamiltonian_Recip(env, k, Hk);
        // Calculate eigenvalues.
        gsl_eigen_herm(Hk, energies, work);
        gsl_sort_vector(energies);
    }

    // Calculate DOS values.
    double *dos_vals = Tetra_AllDosList(Efn, n, num_bands, R, Es, num_dos);

    // Clean up memory.
    gsl_matrix_free(R);
    gsl_eigen_herm_free(work);
    gsl_matrix_complex_free(Hk);

    return dos_vals;
}
