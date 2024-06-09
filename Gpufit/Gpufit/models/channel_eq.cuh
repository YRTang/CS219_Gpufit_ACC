#include <math.h>

#ifndef GPUFIT_CHANNEL_EQ_CUH_INCLUDED
#define GPUFIT_CHANNEL_EQ_CUH_INCLUDED

/***
 * parameters:
 *  tau1, nu1, h_real1, h_imaginary1, 
 *  tau2, nu2, h_real2, h_imaginary2,
 *  tau3, nu3, h_real3, h_imaginary3
 *  p[0], p[1], p[4], p[3]
 *
 * n_fits (number of fits)
 *  number of estimation tasks to be performed at 1 time (custom)
 *
 * n_points (number of data per fit)
 *  pilot number
 *
 * value
 *  real part of estimated y,
 *  imaginary part of estimated y
 *
 * user_info
 *  contains m, n
 *
 */
__device__ void calculate_channel_eq(
    REAL const *parameters,
    int const n_fits,
    int const n_points,
    float *value,
    float *derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char *user_info,
    std::size_t const user_info_size)
{
    ///////////////////////////// modify user_info //////////////////////////////
    REAL *user_info_float = (REAL *)user_info;
    
    REAL m, n = 0;
    REAL const *p = parameters;
    /**
     * path number set to 3 for convenience
     * To do: need to find a way to extract it from parameter size
     */
    int size = 3; 

    int const chunk_begin = chunk_index * n_fits * n_points * 2;
    int const fit_begin = fit_index * n_points;
    m = user_info_float[chunk_begin + fit_begin + point_index * 2];
    n = user_info_float[chunk_begin + fit_begin + point_index * 2 + 1];
    // printf("m=%.3f\n", m);
    // printf("n=%3f\n", n);
    // printf("\n");

    ///////////////////////////// value //////////////////////////////
    // for storing the intermediate periodic function calculations
    REAL * periodic = new REAL[size];
    // initialize real and imaginary part of estimated y to 0
    value[point_index * 2] = 0;
    value[point_index * 2 + 1] = 0;


    for(int i = 0; i < size; i++)
    {
        // cos and sin using tau_i+1, mu_i+1, and h_i+1
        periodic[2 * i] = cos(-2 * M_PI * (m * p[4 * i] - n * p[4 * i + 1]));
        periodic[2 * i + 1] = sin(-2 * M_PI * (m * p[4 * i] - n * p[4 * i + 1]));
        // real part of estimated y
        value[point_index * 2] +=  periodic[2 * i] * p[4 * i + 2] - periodic[2 * i + 1] * p[4* i + 3];     
        // imaginary part of estimated y
        value[point_index * 2 + 1] += periodic[2 * i] * p[4 * i + 3] + periodic[2 * i + 1] * p[4 * i + 2]; 
    }

    /////////////////////////// derivative ///////////////////////////
    REAL *current_derivative = derivative + point_index;

    /**
     * derivative of real part of y
     * dtau1_real, dtau1_img, dtau2_real, dtau2_img, dtau3_real, dtau3_img
     * dnu1_real, dnu1_img, dnu2_real, dnu2_img, dnu3_real, dnu3_img
     * dhreal1_real, dhreal1_img, dhreal2_real, dhreal2_img, dhreal3_real, dhreal3_img
     * dhimg1_real, dhimg1_img, dhimg2_real, dhimg2_img, dhimg3_real, dhimg3_img
     */

    // derivative wrt tau
    // derivative targeted index: 0, 
    int offset = 0;
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = -2 * M_PI * m * periodic[2 * i + 1] * p[4 * i + 2] 
                                                - 2 * M_PI * m * periodic[2 * i] * p[4 * i + 3];
    offset += 3;

    // derivative wrt nu
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = 2 * M_PI * n * periodic[2 * i + 1] * p[4 * i + 2] 
                                                    + 2 * M_PI * n * periodic[2 * i] * p[4 * i + 3];
    offset += 3;

    // derivative wrt h_real
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = periodic[2 * i];
    offset += 3;

    // derivative wrt h_img
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = -periodic[2 * i + 1];
    offset += 3;

    // derivative of imaginary part of y

    // derivative wrt tau
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = -2 * M_PI * m * periodic[2 * i + 1] * p[4 * i + 3] 
                                                    + 2 * M_PI * m * periodic[2 * i] * p[4 * i + 2];
    offset += 3;

    // derivative wrt nu
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = 2 * M_PI * n * periodic[2 * i + 1] * p[4 * i + 3] 
                                                    - 2 * M_PI * n * periodic[2 * i] * p[4 * i + 2];
    offset += 3;

    // derivative wrt h_real
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = periodic[2 * i];
    offset += 3;

    // derivative wrt h_img
    for(int i=0; i < 3; i++)
        current_derivative[(i+offset) * n_points] = periodic[2 * i + 1];
}

#endif