#include <math.h>

#ifndef GPUFIT_CHANNEL_EQ_CUH_INCLUDED
#define GPUFIT_CHANNEL_EQ_CUH_INCLUDED

/***
 * parameters:
 *  tau, nu, h_real, h_imaginary
 *  p[0], p[1], p[2], p[3]
 *  p[4], p[5], p[6], p[7]
 *  p[8], p[9], p[10], p[11]
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
#define TAU(x) parameters[x * 4]
#define NU(x) parameters[x * 4 + 1]
#define H_REAL(x) parameters[x * 4 + 2]
#define H_IMAG(x) parameters[x * 4 + 3]

    bool real_part = point_index % 2 == 0;

    if (user_info_size / sizeof(REAL) == n_points)
    {
        if (real_part)
        {
            m = user_info_float[point_index];
            n = user_info_float[point_index + 1];
        }
        else
        {
            m = user_info_float[point_index - 1];
            n = user_info_float[point_index];
        }
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        if (real_part)
        {
            m = user_info_float[chunk_begin + fit_begin + point_index];
            n = user_info_float[chunk_begin + fit_begin + point_index + 1];
        }
        else
        {
            m = user_info_float[chunk_begin + fit_begin + point_index - 1];
            n = user_info_float[chunk_begin + fit_begin + point_index];
        }
    }

    ///////////////////////////// value //////////////////////////////

    /* REAL const A = cos(-2 * M_PI * (m * p[0] - n * p[1]));
    REAL const B = sin(-2 * M_PI * (m * p[0] - n * p[1])); */

    REAL cos_alpha[3] = {
        cos(-2 * M_PI * (m * TAU(0) - n * NU(0))),
        cos(-2 * M_PI * (m * TAU(1) - n * NU(1))),
        cos(-2 * M_PI * (m * TAU(2) - n * NU(2)))};

    REAL sin_alpha[3] = {
        sin(-2 * M_PI * (m * TAU(0) - n * NU(0))),
        sin(-2 * M_PI * (m * TAU(1) - n * NU(1))),
        sin(-2 * M_PI * (m * TAU(2) - n * NU(2)))};

    // calculate estimated y
    /* value[point_index * 2] = A * p[2] - B * p[3];     // real part of estimated y
    value[point_index * 2 + 1] = A * p[3] + B * p[2]; // imaginary part of estimated y */

    REAL tmp = 0;
    if (real_part)
    {
        for (int i = 0; i < 3; i++)
        {
            tmp += H_REAL(i) * cos_alpha[i] - H_IMAG(i) * sin_alpha[i];
        }
    }
    else
    {
        for (int i = 0; i < 3; i++)
        {
            tmp += H_REAL(i) * sin_alpha[i] + H_IMAG(i) + cos_alpha[i];
        }
    }

    value[point_index] = tmp;

    /////////////////////////// derivative ///////////////////////////
    REAL *current_derivative = derivative + point_index;

// derivative of real part of y
#define D_TAU_OFF(x) 4 * x *n_points
#define D_NU_OFF(x) (4 * x + 1) * n_points
#define D_H_REAL_OFF(x) (4 * x + 2) * n_points
#define D_H_IMAG_OFF(x) (4 * x + 3) * n_points

    if (real_part)
    {
        for (int i = 0; i < 3; i++)
        {
            current_derivative[D_TAU_OFF(i)] = 2 * M_PI *
                                               (H_REAL(i) * sin_alpha[i] + H_IMAG(i) * cos_alpha[i]) * m;

            current_derivative[D_NU_OFF(i)] = -2 * M_PI *
                                              (H_REAL(i) * sin_alpha[i] + H_IMAG(i) * cos_alpha[i]) * n;

            current_derivative[D_H_REAL_OFF(i)] = cos_alpha[i];

            current_derivative[D_H_IMAG_OFF(i)] = -sin_alpha[i];
        }
    }
    else
    {
        for (int i = 0; i < 3; i++)
        {
            current_derivative[D_TAU_OFF(i)] = -2 * M_PI *
                                               (H_REAL(i) * cos_alpha[i] - H_IMAG(i) * sin_alpha[i]) * m;

            current_derivative[D_NU_OFF(i)] = 2 * M_PI *
                                              (H_REAL(i) * cos_alpha[i] - H_IMAG(i) * sin_alpha[i]) * n;

            current_derivative[D_H_REAL_OFF(i)] = sin_alpha[i];

            current_derivative[D_H_IMAG_OFF(i)] = cos_alpha[i];
        }
    }
}

#endif