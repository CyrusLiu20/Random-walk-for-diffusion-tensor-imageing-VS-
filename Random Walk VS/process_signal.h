#pragma once
#ifndef PROCESS_SIGNAL_H
#define PROCESS_SIGNAL_H

#include <vector>
#include "substrate.h"
#include "particle_state.h"
#include "sequence.h"

class process_signal {

public:

    // Constructors
    process_signal() = default;
    process_signal(std::vector<particle_state>& states, substrate& substrate_input, sequence& sequence_input, Eigen::MatrixXd& pos0);

    Eigen::Matrix3d process_phase(Eigen::MatrixXd& phase, double& bvalue);
    void process_displacement(Eigen::MatrixXd& displacement, double& T);


    // Results
    double bvalue;
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero(3, 3);

    double MD, Dx, Dy, Dz;

    Eigen::MatrixXd phase_all;
    Eigen::MatrixXd phase; // Only in ICS and ECS
    Eigen::MatrixXd phase_ICS;
    Eigen::MatrixXd phase_ECS;

    Eigen::MatrixXd position_all;

    Eigen::MatrixXd displacement_all;
    Eigen::MatrixXd displacement;
    Eigen::MatrixXd displacement_ICS;
    Eigen::MatrixXd displacement_ECS;

    Eigen::Array<bool, Eigen::Dynamic, 1> valid;
    Eigen::Array<bool, Eigen::Dynamic, 1> insideECS;
    Eigen::Array<bool, Eigen::Dynamic, 1> insideVoxel;

    bool no_tensor = false;
private:

    // Finding array mask
    Eigen::MatrixXd find_element(Eigen::MatrixXd& target, Eigen::Array<bool, Eigen::Dynamic, 1> input);

    // Number of particles
    int N_p;

};
#endif