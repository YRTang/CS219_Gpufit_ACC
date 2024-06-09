#include "../../../Gpufit/gpufit.h"
#include "mp_solver.h"
#include "csv_reader.h"
#include <vector>
#include <random>
#include <iostream>
#include <math.h>

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

			initial_parameters[i * n_model_parameters + p * 4 + 0] = 0.5; // 4: t, v, h_real, h_imag
			initial_parameters[i * n_model_parameters + p * 4 + 1] = 0.5;
			initial_parameters[i * n_model_parameters + p * 4 + 2] = 0.5;
			initial_parameters[i * n_model_parameters + p * 4 + 3] = 0.5;
		}
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
	int const max_number_iterations = 20;

	// estimator ID
	int const estimator_id = LSE_COMPLEX;

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
}

int main(int argc, char *argv[])
{
	// read input
	mp_profile_t mp_profile;
	csvReader reader(CSV_PATH);
	reader.readData();
	mp_config_t mp_config = reader.getData();

	mp_solver_example(&mp_profile, &mp_config);

	std::cout << std::endl
			  << "Example completed!" << std::endl;
	std::cout << "Press ENTER to exit" << std::endl;
	std::getchar();

	return 0;
}
