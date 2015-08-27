#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "bstrlib/bstrlib.h"
#include "DosValues.h"

int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals);

int main(int argc, char *argv[]) {
    if (argc < 15) {
        printf("Usage: DosValues.out 'out_path' 'n' 'num_dos' 'Tae' 'Tce' 'Tbe' 'Tao' 'Tco' 'Tbo' 'EpsilonR' 'EpsilonM' 'M' 'W' 'Mu'\n");
        printf("Example: DosValues.out 'dos_vals' '8' '500' '0.2' '1.0' '0.4' '0.08' '0.4' '0.16' '0.05' '0.05' '1.0' '1.0' '-2.0'\n");
        return 2;
    }
    // Parse arguments.
    char *out_path = argv[1];
    int n = atoi(argv[2]);
    int num_dos = atoi(argv[3]);
    double Tae = atof(argv[4]);
    double Tce = atof(argv[5]);
    double Tbe = atof(argv[6]);
    double Tao = atof(argv[7]);
    double Tco = atof(argv[8]);
    double Tbo = atof(argv[9]);
    double EpsilonR = atof(argv[10]);
    double EpsilonM = atof(argv[11]);
    double M = atof(argv[12]);
    double W = atof(argv[13]);
    double Mu = atof(argv[14]);

    // Set up data.
    double *Es = malloc(num_dos * sizeof(double));
    Environment *env = malloc(sizeof(Environment));
    env->Tae = Tae;
    env->Tce = Tce;
    env->Tbe = Tbe;
    env->Tao = Tao;
    env->Tco = Tco;
    env->Tbo = Tbo;
    env->EpsilonM = EpsilonM;
    env->EpsilonR = EpsilonR;
    env->M = M;
    env->W = W;
    env->Mu = Mu;

    // Calculate DOS.
    double *dos_vals = DosValues(env, n, Es, num_dos);

    // Write out DOS.
    int write_err = write_dos_vals(out_path, num_dos, Es, dos_vals);
    if (write_err != 0) {
        printf("Error writing DOS values\n");
    }

    // Cleanup. 
    free(dos_vals);
    free(Es);
    return 0;
}

int write_dos_vals(char *outPath, int num_dos, double *Es, double *dos_vals) {
    FILE *fp = fopen(outPath, "w");
    if (fp == NULL) {
        return 1;
    }

    // Write header.
    fprintf(fp, "E\tDOS\n");
    // Write DOS data.
    int i;
    for (i = 0; i < num_dos; i++) {
        fprintf(fp, "%.10f\t%.10f\n", Es[i], dos_vals[i]);
    }

    fclose(fp);
    return 0;
}
