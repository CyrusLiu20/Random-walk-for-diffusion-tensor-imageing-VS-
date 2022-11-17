#ifndef POLYGON_H
#define POLYGON_H

#include <Eigen/Dense>
#include "boundingbox.h"
#include "intersection_ray_info.h"
#include "intersection_info.h"
// #include <random>

// It's actually polyhedron
class polygon {

public:

    polygon() = default;
    polygon(Eigen::MatrixXd vertices_input, Eigen::MatrixXd faces_input);
    polygon(std::string box_type, Eigen::MatrixXd bounding_box);

    // Create substrate bounding box
    void initialize(std::string box_type, Eigen::MatrixXd bounding_box);

    // Inputs the particle position
    bool containsPoint(Eigen::Vector3d& point, unsigned int& seed);
    // Whether particle intersects boundary
    intersection_info intersection(Eigen::Vector3d& orig, Eigen::Vector3d& dir);

    // Vertices of one myocyte
    Eigen::MatrixXd Vertices;
    // Faces of one myocyte
    Eigen::MatrixXd Faces;
    // Transposed vertices of one myocyte
    Eigen::MatrixXd Vertices_t;
    // Transposed faces of one myocyte
    Eigen::MatrixXd Faces_t;

    boundingbox boundingbox;


    // Number of faces
    int n_faces;
    // Number of vertices
    int n_vertices;
    // Myocyte volume
    double volume;
    // Myocyte surface area
    double surface_area;
    // Number of bytes
    int bytes;

    // Data retreival
    Eigen::Vector3d get_minXYZ();

protected:
    Eigen::MatrixXd crossMat(Eigen::MatrixXd& a, Eigen::MatrixXd& b);

private:

    Eigen::MatrixXd get_vertices(int column);
    intersection_ray_info TriangleRayIntersection(Eigen::Vector3d& point, Eigen::Vector3d& dir);
    bool intersection_is_certain(intersection_ray_info& ray, bool test_end, double eps);

    // Compute the mesh property
    void volume_compute(Eigen::MatrixXd& vertices_input, Eigen::MatrixXd& faces_input);
    void surface_compute(Eigen::MatrixXd& vertices_input, Eigen::MatrixXd& faces_input);

    // Stupid finding element method
    Eigen::MatrixXd find_element(Eigen::MatrixXd& target, Eigen::Array<bool, Eigen::Dynamic, 1>& input);
    Eigen::VectorXd find_element(Eigen::VectorXd& target, Eigen::Array<bool, Eigen::Dynamic, 1>& input);
    Eigen::VectorXd find_index(Eigen::Array<bool, Eigen::Dynamic, 1>& input);
    int find_first_index(Eigen::VectorXd& input, double& target);
    bool is_unique(Eigen::VectorXd& input);

    // Random seed
    // std::random_device rd; 

    // The minimum values of the vertices
    Eigen::Vector3d minXYZ = Eigen::Vector3d::Zero(3);
    // The maximum values of the vertices
    Eigen::Vector3d maxXYZ = Eigen::Vector3d::Zero(3);

    // Bounding box range
    Eigen::VectorXd bb_range = Eigen::VectorXd::Zero(6);

    // Mean of the vertices
    Eigen::Vector3d vertices_mean = Eigen::Vector3d::Zero(3);
    // A shifted vertices
    Eigen::MatrixXd vertices_shifted;
    // Volumes of each tetrahedron
    Eigen::VectorXd volume_tetra;
    // Areas of each tetrahedron
    Eigen::VectorXd area_surface;

    // Calculate contains point
    Eigen::MatrixXd V1;
    Eigen::MatrixXd V2;
    Eigen::MatrixXd V3;

    // Whether the vertices has been extracted through faces
    bool vertices_extracted = false;

};

#endif
