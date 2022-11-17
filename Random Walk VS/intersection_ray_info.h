#ifndef INTERSECTION_RAY_INFO_H
#define INTERSECTION_RAY_INFO_H

#include <Eigen/Dense>

class intersection_ray_info {

public:

    intersection_ray_info() = default;
    intersection_ray_info& operator =(const intersection_ray_info& pSrc) {
        // Transfer of resources
        if (this != &pSrc) {
            intersect = pSrc.intersect;
            t = pSrc.t;
            u = pSrc.u;
            v = pSrc.v;

        }
        return *this;
    };

    // Intersection info
    Eigen::VectorXd t, u, v;
    Eigen::Array<bool, Eigen::Dynamic, 1> intersect;
};

#endif
