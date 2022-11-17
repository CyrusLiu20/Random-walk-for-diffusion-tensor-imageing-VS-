#include <cstdio>
#include <chrono>
#include <iostream>

// external
#include <Eigen/Dense>
#include <matio.h>
#include <omp.h>

// Inputs
#include "montecarlofile.h"
#include "sequencefile.h"
#include "substratefile.h"

// MRI
#include "sequence.h"
#include "scansequence.h"
 #include "process_signal.h"

// Substrate
 #include "read_myocytes.h"
 #include "myocytes.h"
 #include "substrate.h"
 #include "transform.h"
 #include "transform_parameter.h"

// Monte Carlo
#include "walkers.h"
#include "simulation.h"
#include "particle_state.h"

// CPP files VScode debugging use
#include "saving_file.h"

int main(int argc, char* argv[]) {

	// Monitor the computational time
	auto begin_total = std::chrono::steady_clock::now();


	// Configurating sequence
	ScanSequence scanner;
	sequencefile sequence_raw;
	sequence sequence;
	sequence = scanner.create(sequence_raw);

	// Configurating substrate
	substratefile substrate_raw;
	read_myocytes myocyte_scanner(substrate_raw.filename);
	myocytes myocytes;
	if (myocyte_scanner.scanned()) {
		myocytes.load(myocyte_scanner);
		substrate_raw.substrate_type = "block";
	}
	else {
		// To do: load myocytes from substratefile directly
	}
	substrate substrate(substrate_raw);
	substrate.myocytes_scan(myocytes);
	substrate.buildcache();
	std::vector<polygon> myos = substrate.myos;

	// Configurating monte carlo simulator
	montecarlofile montecarlo_raw;
	walkers particles(montecarlo_raw);
	Eigen::VectorXd buffer = Eigen::VectorXd::Zero(6);
	Eigen::VectorXd seedbox;
	double simulation_time_total;

	int dim = substrate_raw.dimension.length();
	if (dim == 0) {
		printf("run::sim::No dimension has been selected");
	}

	if (montecarlo_raw.seedbox_type == "text") { // User used text input
		if (montecarlo_raw.seedbox_text.find("voxel") != std::string::npos) {
			seedbox = substrate_raw.voxel_range;
		}
		else if (montecarlo_raw.seedbox_text.find("origin")) {
			seedbox = Eigen::VectorXd::Zero(6);
		}
		else {
			printf("run_sim::Unknown seedbox option, such input is not supported");
		}

		if (montecarlo_raw.seedbox_text.find("buffer")) {
			simulation_time_total = sequence.get_total_time();
			for (int i = 0; i < dim; i++) {
				buffer({ i * 2,i * 2 + 1 }) << -1, 1;
			}
			buffer = buffer * std::sqrt(6 * substrate_raw.D_ecs * simulation_time_total);
			seedbox = seedbox + buffer;
		}
	}
	else {
		seedbox = montecarlo_raw.seedbox_vect;
	}

	simulation system(particles);
	Eigen::MatrixXd seedboxes = seedbox.transpose();
	system.seedParticlesInBox(seedboxes, -1);
	Eigen::MatrixXd pos0 = system.get_position();

	// Testing module
	// testing unit_test("unit_test.mat"); // Debugging use
	// system.load_test_module(unit_test); // Debugging use

	// Simulation of particles
	auto begin_simulation = std::chrono::steady_clock::now();
	std::cout << "Random Walk for cardiac Diffusion Tensor Imaging" << std::endl;
	std::cout << "Simulation of " << particles.get_N_p() << " particles starts" << std::endl;

	system.performScan(sequence, substrate);
	auto finish_simulation = std::chrono::steady_clock::now();

	std::vector<particle_state> states = system.get_states();
	process_signal post_processing(states, substrate, sequence, pos0);

	auto finish_total = std::chrono::steady_clock::now();
	double run_time_total = std::chrono::duration_cast<std::chrono::milliseconds>(finish_total - begin_total).count() / (double)1000;///1000000.0;
	double run_time_simulation = std::chrono::duration_cast<std::chrono::milliseconds>(finish_simulation - begin_simulation).count() / (double)1000;///1000000.0;

	// Save results
	saving_file result(states, pos0, post_processing);
	// saving_file result(states, unit_test.initial_positions, post_processing); // Debugging use

	printf("Simulation time (perform scan) : %lf (seconds)\n", run_time_simulation);
	printf("Total simulation time : %lf (secon1ds)\n", run_time_total);
}