#include "../../../Gpufit/gpufit.h"
#include "mp_solver.h"
#include "csv_reader.h"
#include <vector>
#include <random>
#include <iostream>
#include <math.h>

#include <immintrin.h>
#include <chrono>
#include <complex.h>

#include <cstdlib>
#include <iostream>
#include <time.h>

#define CSV_PATH "../Gpufit/examples/c++/myProject/sample_input.csv"

void mp_solver_example(mp_profile_t *mp_profile, mp_config_t *mp_config)
{
	// number of fits, fit points and parameters
	std::size_t const n_fits = 1;
	std::size_t const n_points_per_fit = mp_config->nof_pilots; //=32
	std::size_t const n_model_parameters = 12;

	// custom x positions for the data points of every fit, stored in user info
	std::vector<REAL> user_info(n_points_per_fit * 2);

	for (std::size_t i = 0; i < n_points_per_fit; i++)
	{
		user_info[i * 2] = mp_config->m[i];
		user_info[i * 2 + 1] = mp_config->n[i];
	}
	// size of user info in bytes
	std::size_t const user_info_size = 2 * n_points_per_fit * sizeof(REAL);

	// initial parameters (randomized)
	std::vector<REAL> initial_parameters(n_fits * n_model_parameters);
	for (std::size_t i = 0; i != n_fits; i++)
	{
		for (int p = 0; p < mp_config->nof_paths; ++p)
		{
			// initial_parameters[i * n_model_parameters + p * 4 + 0] = mp_profile->tau[p]; // 4: t, v, h_real, h_imag
			// initial_parameters[i * n_model_parameters + p * 4 + 1] = mp_profile->nu[p];
			// initial_parameters[i * n_model_parameters + p * 4 + 2] = creal(mp_profile->h[p]);
			// initial_parameters[i * n_model_parameters + p * 4 + 3] = cimag(mp_profile->h[p]);

			initial_parameters[i * n_model_parameters + p * 4 + 0] = 0; // 4: t, v, h_real, h_imag
			initial_parameters[i * n_model_parameters + p * 4 + 1] = 0;
			initial_parameters[i * n_model_parameters + p * 4 + 2] = 1; // h[0]=h[1]=h[2]=1
			initial_parameters[i * n_model_parameters + p * 4 + 3] = 0;
		}
		// initial_parameters[i * n_model_parameters + 2] = 1; //h[0]=1
	}

	// generate data
	std::vector<REAL> data(n_points_per_fit * n_fits);
	for (std::size_t j = 0; j < n_fits; j++)
	{
		for (std::size_t i = 0; i < n_points_per_fit; i++)
		{
			/* data[j * n_fits + 2 * i] = creal(mp_config->y[i]);
			data[j * n_fits + 2 * i + 1] = cimag(mp_config->y[i]); */

			data[j * n_fits + i] = creal(mp_config->y[i]) + cimag(mp_config->y[i]);
		}
	}

	// tolerance
	REAL const tolerance = 0.001f;

	// maximum number of iterations
	int const max_number_iterations = 40;

	// estimator ID
	int const estimator_id = LSE;

	// model ID
	int const model_id = CHANNEL_EQ;

	// parameters to fit (all of them)
	std::vector<int> parameters_to_fit(n_model_parameters, 1);

	// output parameters
	std::vector<REAL> output_parameters(n_fits * n_model_parameters);
	std::vector<int> output_states(n_fits);
	std::vector<REAL> output_chi_square(n_fits);
	std::vector<int> output_number_iterations(n_fits);

	// call to gpufit (C interface)
	int const status = gpufit(
		n_fits,
		n_points_per_fit,
		data.data(),
		0,
		model_id,
		initial_parameters.data(),
		tolerance,
		max_number_iterations,
		parameters_to_fit.data(),
		estimator_id,
		user_info_size,
		reinterpret_cast<char *>(user_info.data()),
		output_parameters.data(),
		output_states.data(),
		output_chi_square.data(),
		output_number_iterations.data());

	// check status
	if (status != ReturnState::OK)
	{
		throw std::runtime_error(gpufit_get_last_error());
	}

	// examine output_parameters
	//  Matrix[3x4]
	cout << "size=" << output_parameters.size() << endl;
	for (int i = 0; i < 12; i++)
	{
		if (i % 4 == 0)
			cout << endl;
		cout << output_parameters.data()[i] << " ";
	}

	cout << endl
		 << "iteration number: " << endl;

	for (int i : output_number_iterations)
	{
		cout << i << endl;
	}

	cout << endl
		 << "chi-square" << endl;

	for (REAL i : output_chi_square)
	{
		cout << i << endl;
	}

	cout << endl
		 << "state" << endl;

	for (int i : output_states)
	{
		cout << i << endl;
	}

	cout << "=============== Convergence analysis ===================" << endl;
	cout << "In this part, we test for differnt iterations on different prarmeter initialization\n \
	 but the same data is used. To see whether the model really converge."
		 << endl
		 << endl;

	int iter_times_test[] = {1, 2, 3, 4, 5, 10, 15, 20, 30, 50};

	srand(time(0));
	for (int x = 0; x < 10; x++)
	{
		REAL tau, nu, h_real, h_imag;
		tau = rand();
		nu = rand();
		h_real = rand();
		h_imag = rand();

		for (int p = 0; p < mp_config->nof_paths; ++p)
		{
			// initial_parameters[i * n_model_parameters + p * 4 + 0] = mp_profile->tau[p]; // 4: t, v, h_real, h_imag
			// initial_parameters[i * n_model_parameters + p * 4 + 1] = mp_profile->nu[p];
			// initial_parameters[i * n_model_parameters + p * 4 + 2] = creal(mp_profile->h[p]);
			// initial_parameters[i * n_model_parameters + p * 4 + 3] = cimag(mp_profile->h[p]);

			initial_parameters[n_model_parameters + p * 4 + 0] = tau; // 4: t, v, h_real, h_imag
			initial_parameters[n_model_parameters + p * 4 + 1] = nu;
			initial_parameters[n_model_parameters + p * 4 + 2] = h_real; // h[0]=h[1]=h[2]=1
			initial_parameters[n_model_parameters + p * 4 + 3] = h_imag;
		}

		cout << "Initialized parameters as ("
			 << tau << ", " << nu << ", " << h_real << ", " << h_imag << ")" << endl;

		for (int i = 0; i < 10; i++)
		{
			int const status = gpufit(
				1,
				n_points_per_fit,
				data.data(),
				0,
				model_id,
				initial_parameters.data(),
				tolerance,
				iter_times_test[i],
				parameters_to_fit.data(),
				estimator_id,
				user_info_size,
				reinterpret_cast<char *>(user_info.data()),
				output_parameters.data(),
				output_states.data(),
				output_chi_square.data(),
				output_number_iterations.data());

			// check status
			if (status != ReturnState::OK)
			{
				throw std::runtime_error(gpufit_get_last_error());
			}

			cout << "The chi-square under max-iteration " << iter_times_test[i] << " is "
				 << output_chi_square[0] << endl;

			cout << "The output parameters are " << endl;

			// examine output_parameters
			//  Matrix[3x4]
			for (int k = 0; k < 12; k++)
			{
				if (k % 4 == 0 && k != 0)
					cout << endl;
				cout << output_parameters.data()[k] << " ";
			}

			cout << endl << endl;;
		}
	}
}

int main(int argc, char *argv[])
{
	// read input
	mp_profile_t mp_profile;
	csvReader reader(CSV_PATH);
	reader.readData();
	mp_config_t mp_config = reader.getData();

	auto t1 = chrono::high_resolution_clock::now();
	mp_solver_example(&mp_profile, &mp_config);
	auto t2 = chrono::high_resolution_clock::now();
	chrono::duration<double, std::milli> ms_double = t2 - t1;
	cout << "execution time:" << ms_double.count() << "ms" << endl;

	std::cout << std::endl
			  << "Example completed!" << std::endl;
	std::cout << "Press ENTER to exit" << std::endl;
	std::getchar();

	return 0;
}
