#ifndef TRANSFORM_PARAMETER_H
#define TRANSFORM_PARAMETER_H

#include <Eigen/Dense>
class transform_parameter {

public:

    transform_parameter() = default;

    Eigen::Vector3d position_local;
    double angle;
    double angle_reverse; // The rotation angle in rotate_y function

    int iX = 0;
    int iY = 0;
    int iZ = 0;

    bool transform_inverse = false;

};

#endif
