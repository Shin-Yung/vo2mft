#include "elHamiltonian.h"

// Set H equal to the 4x4 electronic Hamiltonian H(k).
// k is in the reciprocal lattice basis (i.e. k = (k_1, k_2, k_3) with
// corresponding Cartesian representation k_Cart = k_1 b_1 + k_2 b_2 + k_3 b_3).
void ElHamiltonian_Recip(Environment *env, gsl_vector *k, gsl_matrix_complex *H) {
    // k --> 2*pi*k
    double ki;
    int i;
    for (i = 0; i < 3; i++) {
        ki = gsl_vector_get(k, i);
        gsl_vector_set(k, i, 2.0*M_PI*ki);
    }
    ElHamiltonian(env, k, H);
}

// Set H equal to the 4x4 electronic Hamiltonian H(k).
// k is in the Cartesian basis, with each component scaled by the corresponding
// lattice constant; i.e. k = (a kx, a ky, c kz) and a kx, a ky, c kz range
// over [-pi, pi) and periodic copies of this interval.
void ElHamiltonian(Environment *env, gsl_vector *k, gsl_matrix_complex *H) {
    double kx = gsl_vector_get(k, 0);
    double ky = gsl_vector_get(k, 1);
    double kz = gsl_vector_get(k, 2);

    // Functions of k.
	gsl_complex EpsAE = EpsilonAE(env, k);
	gsl_complex EpsBE = EpsilonBE(env, k);
	gsl_complex EpsAO = EpsilonAO(env, k);
	gsl_complex EpsBO = EpsilonBO(env, k);

    gsl_complex ikd = gsl_complex_rect(0.0, kx/2.0+ky/2.0+kz/2.0);
    gsl_complex eikd = gsl_complex_exp(ikd);
    gsl_complex emikd = gsl_complex_exp(gsl_complex_mul_real(ikd, -1.0));

    // Functions of k + Q.
    // k --> k + Q
    double ki;
    int i;
    for (i = 0; i < 3; i++) {
        ki = gsl_vector_get(k, i);
        gsl_vector_set(k, i, ki + M_PI);
    }
	gsl_complex EpsBE_KQ = EpsilonBE(env, k);

    // Parts independent of k.
    double ident_re = (1.0-env->W)*env->EpsilonR+env->W*env->EpsilonM - env->Mu;
	gsl_complex ident_part = gsl_complex_rect(ident_re, 0.0);

    // Construct H(k).
	gsl_matrix_complex_set(H, 0, 0, gsl_complex_add(EpsAE, ident_part));
	gsl_matrix_complex_set(H, 1, 0, gsl_complex_mul_real(EpsAO, 2.0));
	gsl_matrix_complex_set(H, 2, 0, gsl_complex_mul(EpsBE, emikd));
	gsl_matrix_complex_set(H, 3, 0, gsl_complex_mul_real(gsl_complex_mul(gsl_complex_conjugate(EpsBO), emikd), -1.0));

	gsl_matrix_complex_set(H, 0, 1, gsl_complex_mul_real(EpsAO, -2.0));
	gsl_matrix_complex_set(H, 1, 1, gsl_complex_add(gsl_complex_mul_real(EpsAE, -1.0), ident_part));
	gsl_matrix_complex_set(H, 2, 1, gsl_complex_mul(gsl_complex_conjugate(EpsBO), emikd));
	gsl_matrix_complex_set(H, 3, 1, gsl_complex_mul(gsl_complex_rect(0.0, 1.0), gsl_complex_mul(EpsBE_KQ, emikd)));

	gsl_matrix_complex_set(H, 0, 2, gsl_complex_mul(EpsBE, eikd));
	gsl_matrix_complex_set(H, 1, 2, gsl_complex_mul(EpsBO, eikd));
	gsl_matrix_complex_set(H, 2, 2, gsl_complex_add(EpsAE, ident_part));
	gsl_matrix_complex_set(H, 3, 2, gsl_complex_mul_real(EpsAO, 2.0));

	gsl_matrix_complex_set(H, 0, 3, gsl_complex_mul_real(gsl_complex_mul(EpsBO, eikd), -1.0));
	gsl_matrix_complex_set(H, 1, 3, gsl_complex_mul(gsl_complex_rect(0.0, -1.0), gsl_complex_mul(EpsBE_KQ, eikd)));
	gsl_matrix_complex_set(H, 2, 3, gsl_complex_mul_real(EpsAO, -2.0));
	gsl_matrix_complex_set(H, 3, 3, gsl_complex_add(gsl_complex_mul_real(EpsAE, -1.0), ident_part));
}

// Cubic axes, even symmetry (k, p; k, p)
gsl_complex EpsilonAE(Environment *env, gsl_vector *k) {
    double kx = gsl_vector_get(k, 0);
    double ky = gsl_vector_get(k, 1);
    double kz = gsl_vector_get(k, 2);
    double rp = -2.0 * (env->Tae*(gsl_sf_cos(kx)+gsl_sf_cos(ky)) + env->Tce*gsl_sf_cos(kz));
    return gsl_complex_rect(rp, 0.0);
}

// Body diagonal, even symmetry (k, p; k, pbar)
gsl_complex EpsilonBE(Environment *env, gsl_vector *k) {
    double kx = gsl_vector_get(k, 0);
    double ky = gsl_vector_get(k, 1);
    double kz = gsl_vector_get(k, 2);
    double rp = -8.0 * env->Tbe * gsl_sf_cos(kx/2.0) * gsl_sf_cos(ky/2.0) * gsl_sf_cos(kz/2.0);
    return gsl_complex_rect(rp, 0.0);
}

// Cubic axes, odd symmetry (k, p; k+Q, p)
gsl_complex EpsilonAO(Environment *env, gsl_vector *k) {
    double kx = gsl_vector_get(k, 0);
    double ky = gsl_vector_get(k, 1);
    double kz = gsl_vector_get(k, 2);
    double ip = -2.0 * env->M * (env->Tao*(gsl_sf_sin(kx)+gsl_sf_sin(ky)) + env->Tco*gsl_sf_sin(kz));
    return gsl_complex_rect(0.0, ip);
}

// Body diagonal, odd symmetry (k, p; k+Q, pbar)
gsl_complex EpsilonBO(Environment *env, gsl_vector *k) {
    double kx = gsl_vector_get(k, 0);
    double ky = gsl_vector_get(k, 1);
    double kz = gsl_vector_get(k, 2);
	double rp = -8.0 * env->M * env->Tbo * gsl_sf_cos(kx/2.0) * gsl_sf_cos(ky/2.0) * gsl_sf_cos(kz/2.0);
	double ip = 8.0 * env->M * env->Tbo * gsl_sf_sin(kx/2.0) * gsl_sf_sin(ky/2.0) * gsl_sf_sin(kz/2.0);
    return gsl_complex_rect(rp, ip);
}
