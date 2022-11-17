#pragma once
#ifndef TESTING_H
#define TESTING_H

#include <matio.h>
#include "particle_state.h"
#include <fstream>
#include <Eigen/Dense>
#include <vector>

// Please note, this assumes only 1 counter per time step
class testing {

public:

    testing() = default;
    testing(std::string filename);
    bool initialize(std::string filename);
    bool generate_report(std::vector<particle_state>& states);

    int get_N_p();

    std::vector<Eigen::MatrixXd> directions;
    std::vector<Eigen::VectorXd> probabilities;

    Eigen::MatrixXd positions;
    Eigen::MatrixXd initial_positions;
    Eigen::MatrixXd phases;

    Eigen::VectorXd myoindex;

private:

    // Mat file pointer
    mat_t* test;

    // Reading mat file data
    bool read_fields();
    Eigen::MatrixXd matrix_conversion(matvar_t* field);
    Eigen::VectorXd vector_conversion(matvar_t* field);
    std::vector<std::string> comparison(std::vector<particle_state>& states);

    // Whether the testing module has been scanned
    bool flag = false;
    bool read_directions = false;
    bool read_probabilities = false;
    bool read_positions = false;
    bool read_phases = false;

    // The number of particles simulated in matlab
    int N_p_test;
    // The number of time steps simulated in matlab
    int N_t_test;

    // Tolerance between the two values
    double tolerance = 1e-6;
};

#endif