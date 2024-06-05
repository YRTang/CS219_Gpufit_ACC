#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_multifit_nlinear.h>
#include <immintrin.h>
#include <chrono>
#include "mp_solver.h"
#include "csv_reader.h"
#include "/home/test/sarahTang/Gpufit/Gpufit/gpufit.h" // ../../Gpufit/gpufit.h

#define MAX_ITER 20
#define CSV_PATH "/home/test/sarahTang/Gpufit/examples/c++/myProject/sample_input.csv"

using namespace std;

static int
channel_eq_func(const gsl_vector *x, void *data,
                gsl_vector *f)
{
  mp_config_t *mp_config = (mp_config_t *)data;

  double tau[MAX_NOF_PATHS];
  double nu[MAX_NOF_PATHS];
  _Complex double h[MAX_NOF_PATHS];

  for (int p = 0; p < mp_config->nof_paths; ++p)
  {
    tau[p] = gsl_vector_get(x, p);
    nu[p] = gsl_vector_get(x, mp_config->nof_paths + p);
    h[p] = gsl_vector_get(x, 2 * mp_config->nof_paths + p) + gsl_vector_get(x, 3 * mp_config->nof_paths + p) * I;
  }

  for (int k = 0; k < mp_config->nof_pilots; ++k)
  {
    _Complex double Yk = 0;
    for (int p = 0; p < mp_config->nof_paths; ++p)
    {
      Yk += h[p] * cexp(-I * 2 * M_PI * (tau[p] * (mp_config->m[k] + 0.5) - nu[p] * (mp_config->n[k] + 0.5)));
    }
    gsl_vector_set(f, k, cabs(Yk - mp_config->y[k]));
  }

  return GSL_SUCCESS;
}


int mp_solver(mp_config_t *mp_config, mp_profile_t *mp_profile)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
      gsl_multifit_nlinear_default_parameters();

  double weights[MAX_NOF_PILOTS];
  double x_init[4 * MAX_NOF_PATHS];
  for (int p = 0; p < mp_config->nof_paths; ++p)
  {
    x_init[p] = mp_profile->tau[p];
    x_init[mp_config->nof_paths + p] = mp_profile->nu[p];
    x_init[2 * mp_config->nof_paths + p] = creal(mp_profile->h[p]);
    x_init[3 * mp_config->nof_paths + p] = cimag(mp_profile->h[p]);
  }

  gsl_vector_view x = gsl_vector_view_array(x_init, 4 * mp_config->nof_paths);
  gsl_vector_view wts = gsl_vector_view_array(weights, mp_config->nof_pilots);
  int status, info;

  const double xtol = 1e-6;
  const double gtol = 1e-6;
  const double ftol = 0.0;

  /* define the function to be minimized */
  fdf.f = channel_eq_func;
  fdf.df = NULL;  /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL; /* not using geodesic acceleration */
  fdf.n = mp_config->nof_pilots;
  fdf.p = 4 * mp_config->nof_paths;
  fdf.params = mp_config;

  /* this is the data to be fitted */
  for (int i = 0; i < mp_config->nof_pilots; ++i)
  {
    weights[i] = 1.0;
  };

  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc(T, &fdf_params, mp_config->nof_pilots, 4 * mp_config->nof_paths);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit(&x.vector, &wts.vector, &fdf, w);

  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(MAX_ITER, xtol, gtol, ftol,
                                       NULL, NULL, &info, w);

  /* read final solution */
  gsl_vector *res = gsl_multifit_nlinear_position(w);
  for (int p = 0; p < mp_config->nof_paths; ++p)
  {
    mp_profile->tau[p] = gsl_vector_get(res, p);
    mp_profile->nu[p] = gsl_vector_get(res, mp_config->nof_paths + p);
    double h_r = gsl_vector_get(res, 2 * mp_config->nof_paths + p);
    double h_i = gsl_vector_get(res, 3 * mp_config->nof_paths + p);
    mp_profile->h[p] = h_r + h_i * I;
  }

  gsl_multifit_nlinear_free(w);

  return 0;
}

int main()
{
  // read input
  mp_profile_t mp_profile;
  csvReader reader(CSV_PATH);
  reader.readData();
  mp_config_t mp_config = reader.getData();
  
  
  // mp_config.nof_pilots = 32;
  // mp_config.nof_paths = 3;


  auto t1 = chrono::high_resolution_clock::now();
  // for (int i = 0; i < 1000; ++i) { 
  //   mp_solver(&mp_config, &mp_profile);
  // } 
  auto t2 = chrono::high_resolution_clock::now();
  chrono::duration<double, std::milli> ms_double = t2 - t1;
  cout << "execution time:" << ms_double.count() << "ms" << endl;
  return 0;
}
