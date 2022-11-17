#pragma once
#ifndef WALKERS_H
#define WALKERS_H

#include <Eigen/Dense>
#include <string>
#include "montecarlofile.h"

class walkers {

public:

    // Construtors
    walkers() = default;
    //walkers(const walkers& pSrc) = default;	// copying class
    //walkers(walkers&& pSrc) = default;	// moving class
    //~walkers() = default;

    // Constructors with parameters input
    walkers(montecarlofile parameters);

    // Operators overloading
    walkers& operator =(const walkers& pSrc);	// Transfer of resourses using equal sign

    // Initialize
    bool initialize(int N_p_input, int seed_input, std::string stepType);

    // Data retrieval
    int get_N_p();
    unsigned int get_rng_seed();
    std::string get_steptype();

    // [nP x DIM] array of [x, y, z]-positions of all particles
    // [nP x DIM] array of accumulated [x, y, z]-phase of all particles
    Eigen::MatrixXd position, phase;
    // {nP x 1} cell array of flags
    Eigen::VectorXd flag;

private:

    // Number of particles
    int N_p;

    // Random seed to all rng
    int rng_seed;

    // Constant time step
    std::string stepType = "constant";

    // Flag for initialization
    bool initialized = false;

};

#endif