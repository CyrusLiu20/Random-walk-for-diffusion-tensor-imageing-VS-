#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "myocytes.h"
#include "substratefile.h"
#include "polygon.h"
#include "intersection_info.h"
#include "transform.h"

class substrate {

public:

    substrate() = default;
    substrate(substratefile substratefile);
    void myocytes_scan(myocytes myocytes_input);

    // Build the bounding box for the block
    void buildcache();
    // Find the position of the particles in the myocytes
    int findMyocyte(Eigen::Vector3d position, unsigned int seed, std::string refFrame);
    int search_myocytes(Eigen::Vector3d position_local, unsigned int& seed);

    intersection_info intersectMyocytes(Eigen::Vector3d& position, Eigen::Vector3d& dxdydz, std::string refFrame);
    Eigen::Array<bool, Eigen::Dynamic, 1> needsChecking(Eigen::Vector3d& position_updated, Eigen::Vector3d& dxdydz, std::string refFrame);

    // Detailed version of myocytes with bounding box, volume etc
    std::vector<polygon> myos;

    // load other geometry parameters
    Eigen::Vector3d dxdydz;

    // Voxel range, volume and surface area
    boundingbox voxel;
    // substrate block range, volume and surface area
    polygon block_bb;

    Eigen::VectorXd myocytes_bbrange;

    // boundary type
    std::string boundary;

    // Transform parameters
    transform substrate_transform;

    // Diffusion properties
    std::string transit_model;
    double kappa;
    double D_e;
    double D_i;
    std::string dim;

private:

    // Stupid iteration method
    Eigen::VectorXd iteration(double a, double b, double step);
    // Very stupid finding index method
    Eigen::VectorXd find_index(Eigen::Array<bool, Eigen::Dynamic, 1>& input);
    // Custom absolute function
    double absolute(double input);

    // myocyte geometry
    myocytes myocytes_process;

    // Maximum and minimum of voxel length in y dimension
    Eigen::Vector2d y_extent = Eigen::Vector2d::Zero(2);
    // Realistic y slice
    Eigen::VectorXd y_minvals;
    Eigen::MatrixXd y_slice_minmax;

    // length in y dimension
    double Ly;
    // Number of myocytes
    int N_m;

    // substrate type
    std::string substrate_type;

    // // Transform parameters
    // transform substrate_transform;

    // To be deleted

    // double yextent_min;
    // double yextent_max;

    // int deg_rot_per_L_in_y; // degree of rotation per unit length in y
    // double z_amplitude; // amplitude in z for the sinusoidal transform
    // double x_frequency; // frquency in x for the sinusoidal transform
    // bool shift_block;  

    // bool isIdentity;

    // Eigen::Vector3d dxdydz_bb;
};

#endif