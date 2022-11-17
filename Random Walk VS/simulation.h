#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Dense>
#include <matio.h>
#include "walkers.h"
#include "substrate.h"
#include "sequence.h"
#include "particle_state.h"
#include "polygon.h"
#include "testing.h"
#include <random>

class simulation {

public:
    // Construtors
    simulation() = default;
    simulation(walkers obj_input);
    simulation(const simulation& pSrc) = default;	// copying class
    simulation(simulation&& pSrc) = default;	// moving class
    ~simulation() = default;

    bool seedParticlesInBox(Eigen::MatrixXd boundingboxes_input, int particlesPerBox_input);
    void performScan(sequence sequence_input, substrate substrate_input);
    particle_state onewalker(sequence& sequence_input, substrate& substrate_input, int i_particle, Eigen::Vector3d position, Eigen::Vector3d phase_input, bool flag);
    Eigen::VectorXd one_dt(Eigen::Vector3d& position, double& dt, substrate& substrate_input, int& myoindex, int& i_particle, int& timestep);
    Eigen::Vector3d getLimitedSteps(std::string& dim, double& maxStepLength);

    // A custom function for refilling missing particles
    Eigen::VectorXd refill(Eigen::VectorXd particlesPerBox, int missing_particles);

    // data retreival
    Eigen::MatrixXd get_position();
    std::vector<particle_state> get_states();

    // For debugging purposes
    void load_test_module(testing test);

private:

    // Creates a history file for particle position and phase
    void create_initial_states(long long unsigned int number_of_particles);
    // void save_history(Eigen::MatrixXd position, Eigen::MatrixXd phase, int i_particle);
    void save_history();

    // Compute the normal distance
    Eigen::Vector3d reflect(Eigen::Vector3d& oldstep, Eigen::MatrixXd& faceVertices);
    double computeNormalDistance(Eigen::MatrixXd& faceVertices, Eigen::Vector3d& step);

    // Final states of all the particles (position, phase, and flag)
    std::vector<particle_state> particle_states;

    // The history of each particle in each time step
    // std::vector<Eigen::MatrixXd> position_histories;
    // std::vector<Eigen::MatrixXd> phase_histories;

    // Random seed
    std::random_device rd;

    // The bounding box
    Eigen::MatrixXd boundingboxes;
    // Side length of the voxel (voxel)
    Eigen::MatrixXd sidelengths_raw, sidelengths;
    // Volumes of all boxes
    Eigen::VectorXd boxVolumes;
    // Particles per box both theoretical and preliminary
    Eigen::VectorXd particlesPerBox_theo, particlesPerBox_prel, particlesPerBox;

    // Initial position
    Eigen::MatrixXd pos0;

    // Initializing all the matrices and vectors in this class
    bool initialize(int N_b);
    // Randomly placing the particles in bounding boxes
    bool seeding();
    // Flag for any zero dimensions
    bool hasZeroDim;

    // System of particles
    walkers obj;
    // Number of columns in bounding box input
    int cols;
    // Dimension of bounding box
    int dimension = 0;
    // Number of bounding boxes
    int N_b;

    // Total box volumes
    double total_boxVolume;

    mat_t* intermediate;

    // For debugging purposes
    testing unit_test;
    bool counter_more_than_one = false;
    bool debugging = false;
};

#endif