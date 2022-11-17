// #include <eigen3/Eigen/Dense>
#include <cstdio>

// header file
#include "walkers.h"

walkers& walkers::operator =(const walkers& pSrc) {
    // Transfer of resources
    if (this != &pSrc) {
        N_p = pSrc.N_p;
        rng_seed = pSrc.rng_seed;
        stepType = pSrc.stepType;
        position = pSrc.position;
        phase = pSrc.phase;
        flag = pSrc.flag;
    }
    return *this;
}

walkers::walkers(montecarlofile parameters) {
    walkers::initialize(parameters.N_p, parameters.rngseed, parameters.stepType);
}

bool walkers::initialize(int N_p_input, int seed_input, std::string stepType_input) {
    // Check if Number of particles and the seed is specified
    N_p = N_p_input;
    rng_seed = seed_input;
    stepType = stepType_input;
    if (not(N_p == -1 and seed_input == -1)) {
        initialized = true;
    }
    else {
        std::printf("Initialization status: %s", initialized ? "true" : "false");
    }

    position = Eigen::MatrixXd::Zero(N_p, 3); // initially located at origin
    phase = Eigen::MatrixXd::Zero(N_p, 3); // initially no accuired phase
    flag = Eigen::VectorXd::Zero(N_p, 1); // initially unflagged

    return initialized ? true : false;
}

int walkers::get_N_p() {
    return N_p;
}

unsigned int walkers::get_rng_seed() {
    return rng_seed;
}

std::string walkers::get_steptype() {
    return stepType;
}
