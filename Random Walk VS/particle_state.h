#pragma once
#ifndef PARTICLE_STATE_H
#define PARTICLE_STATE_H

#include <Eigen/Dense>

class particle_state {

public:

    // [nP x DIM] array of [x, y, z]-positions of all particles
    // [nP x DIM] array of accumulated [x, y, z]-phase of all particles
    Eigen::MatrixXd position_history, phase_history;
    // Current position
    Eigen::Vector3d position, phase;
    // {nP x 1} cell array of flags
    bool flag;
    // In which grid the particle is in
    int myoindex;
};

#endif