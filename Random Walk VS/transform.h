#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Eigen/Dense>
#include "transform_parameter.h"

class transform {

public:

    // constructor
    transform() = default;

    transform_parameter global2local(Eigen::Vector3d& position);
    Eigen::Vector3d local2global(Eigen::Vector3d& pos_local, int& iX, int& iY, int& iZ);
    // properties
    bool isIdentity = true;
    Eigen::Vector3d dxdydz_bb = Eigen::Vector3d::Zero(3);
    Eigen::MatrixXd y_slice_minmax;
    double deg_rot_per_L_in_y = 0;
    double z_amplitude = 0;
    double x_frequency = 1;
    bool shift_block = true;

    // Will be made private
    Eigen::Vector3d rotate_y(Eigen::Vector3d& position, double& theta);
    double sin(double& dx, double& dz, double& x);
    double custom_mod(double a, double b);
    double find_yslice(double& y_position); // Apply sin(x) displacement in z' 

private:

    // stupid absolute function
    double absolute(double& input);
    // custom stupid vector mod function
    // Eigen::Vector3d vect_mod(Eigen::Vector3d &original, Eigen::Vector3d &vector_mod)
};

#endif
