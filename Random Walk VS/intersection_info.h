#ifndef INTERSECTION_INFO_H
#define INTERSECTION_INFO_H

#include <Eigen/Dense>

class intersection_info {

public:

    intersection_info() = default;
    intersection_info operator =(const intersection_info& pSrc) {
        // Transfer of resources
        if (this != &pSrc) {
            t = pSrc.t;
            myoindex = pSrc.myoindex;
            vertices = pSrc.vertices;
            empty = pSrc.empty;

        }
        return *this;
    }

    int myoindex;
    double t;
    Eigen::MatrixXd vertices;
    bool empty = true;

};

#endif
