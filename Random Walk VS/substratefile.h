#ifndef SUBSTRATEFILE_H
#define SUBSTRATEFILE_H

#include <Eigen/Dense>
#include <string>
class substratefile {

public:

    substratefile() {
        voxel_range.resize(6);
        voxel_range << voxel_xmin, voxel_xmax, voxel_ymin, voxel_ymax, voxel_zmin, voxel_zmax;
    };

    // Geometry 
    // Transform
    // y extent
    double yextent_min = -2000;
    double yextent_max = 5000;

    double deg_rot_per_L_in_y = 1.5; // degree of rotation per unit length in y
    double z_amplitude = 0; // amplitude in z for the sinusoidal transform
    double x_frequency = 0; // frquency in x for the sinusoidal transform
    bool shift_block = false;

    std::string substrate_type = ""; // This will be inputted later when scanning myocytes

    // myocyte
    std::string filename = "geometry_1.mat";
    // size of the block
    double block_Lx = 495.3992; // x dimension
    double block_Ly = 392.3432; // y dimension
    double block_Lz = 126.5612; // z dimension


    // Domain
    // Imaging voxel
    Eigen::VectorXd voxel_range;
    double voxel_xmin = 0;
    double voxel_xmax = 494.899;
    double voxel_ymin = 0;
    double voxel_ymax = 391.843;
    double voxel_zmin = 58.0310;
    double voxel_zmax = 68.0310;


    // Transit model and permeability for the membranes
    std::string transit_model = "HybridModel"; // Type of model
    double permeability = 0.9;
    // Compartment specific diffusivity
    // Diffusitivity
    double D_ics = 1.5; // Intra-cellular diffusivity
    double D_ecs = 2.5; // Extra-cellular diffusivity


    //Dimension
    std::string dimension = "xyz";
};

#endif
