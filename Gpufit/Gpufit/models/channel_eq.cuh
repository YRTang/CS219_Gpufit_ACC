#include <math.h>

#ifndef GPUFIT_CHANNEL_EQ_CUH_INCLUDED
#define GPUFIT_CHANNEL_EQ_CUH_INCLUDED

/***
 * parameters:
 *  tau, nu, h_real, h_imaginary
 *  p[0], p[1], p[2], p[3]
 * 
 * n_fits (number of fits)
 *  number of estimation tasks to be performed at 1 time (custom)
 * 
 * n_points (number of data per fit)
 *  pilot number
 * 
 * value
 *  estimated y
 * 
 * user_info
 *  m, n, y_real, y_imaginary
 * 
 * 
*/
__device__ void channel_eq (
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    float * value,
    float * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
    ///////////////////////////// modify user_info //////////////////////////////
    REAL * user_info_float = (REAL*)user_info;
    REAL * m, n= 0;

    int const chunk_begin = chunk_index * n_fits * n_points * 2;
    int const fit_begin = fit_index * n_points;
    m = user_info_float[chunk_begin + fit_begin + point_index * 2];
    n = user_info_float[chunk_begin + fit_begin + point_index * 2 + 1];

    ///////////////////////////// value //////////////////////////////

    REAL const A = cos(2 * M_PI * (m * p[0] - n * p[1]));
    REAL const B = sin(2 * M_PI * (m * p[0] - n * p[1]));

    // calculate estimated y
    value[point_index * 2] = A * p[2] - B * p[3];              // real part of estimated y
    value[point_index * 2 + 1] = A * p[3] + B * p[2];          // imaginary part of estimated y

    /////////////////////////// derivative ///////////////////////////
    REAL * current_derivative = derivative + point_index;

    // derivative of real part of y

    // derivative wrt tau
    for (int i = 0; i < 3; i++)
        current_derivative[i * 2 * n_points] = -2 * M_PI * m * B * p[2] - 2 * M_PI * m * A * p[3]; 

    // derivative wrt nu
    for (int i = 0; i < 3; i++)
        current_derivative[(i + 3) * 2 * n_points ] = 2 * M_PI * n * B * p[2] + 2 * M_PI * n * A * p[3]; 
    
    // derivative wrt h_real
    for (int i = 0; i < 3; i++)
        current_derivative[(i + 6) * 2  * n_points] = A; 

    // derivative wrt h_img
        // derivative wrt h_real
    for (int i = 0; i < 3; i++)
        current_derivative[(i + 9) * 2  * n_points] = -B; 

    // derivative of imaginary part of y

        // derivative wrt tau
    for (int i = 0; i < 3; i++)
        current_derivative[(2 * i - 1) * n_points] = -2 * M_PI * m * B * p[3] + 2 * M_PI * m * A * p[2]; 

    // derivative wrt nu
    for (int i = 0; i < 3; i++)
        current_derivative[(2 * i + 5) * n_points] = 2 * M_PI * n * B * p[3] - 2 * M_PI * n * A * p[2]; 
    
    // derivative wrt h_real
    for (int i = 0; i < 3; i++)
        current_derivative[(2 * i + 11) * n_points] = A; 

    // derivative wrt h_img
    for (int i = 0; i < 3; i++)
        current_derivative[(2 * i + 17) * n_points] = B; 
}

#endif