#ifndef MP_SOLVER_H
#define MP_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#define MAX_NOF_PATHS 3
#define MAX_NOF_PILOTS 32

// Define the structure for configuration parameters
typedef struct
{
    int nof_paths;  // P
    int nof_pilots; // K
    int m[MAX_NOF_PILOTS];
    int n[MAX_NOF_PILOTS];
    _Complex double y[MAX_NOF_PILOTS];
} mp_config_t; // input

typedef struct
{
    double tau[MAX_NOF_PATHS];
    double nu[MAX_NOF_PATHS];
    _Complex double h[MAX_NOF_PATHS];
} mp_profile_t; // output

int mp_solver(mp_config_t *mp_config, mp_profile_t *mp_profile);

#endif