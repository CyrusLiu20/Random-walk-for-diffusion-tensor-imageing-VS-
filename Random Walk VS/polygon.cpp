#include "polygon.h"
#include <iostream>
#include <random>

polygon::polygon(Eigen::MatrixXd vertices_input, Eigen::MatrixXd faces_input) {

    // Assigning vertices and faces
    Vertices = vertices_input;
    Vertices_t = vertices_input.transpose();
    Faces = faces_input;
    Faces_t = faces_input.transpose();

    // Number of faces and vertices
    n_vertices = vertices_input.rows();
    n_faces = faces_input.rows();

    polygon::volume_compute(vertices_input, faces_input);
    polygon::surface_compute(vertices_input, faces_input);

    // For intersection functions
    V1 = polygon::get_vertices(0);
    V2 = polygon::get_vertices(1);
    V3 = polygon::get_vertices(2);
    vertices_extracted = true;

    // Computing the bounding box range
    bb_range.resize(6);
    minXYZ << vertices_input(Eigen::indexing::all, 0).minCoeff(), vertices_input(Eigen::indexing::all, 1).minCoeff(), vertices_input(Eigen::indexing::all, 2).minCoeff();
    maxXYZ << vertices_input(Eigen::indexing::all, 0).maxCoeff(), vertices_input(Eigen::indexing::all, 1).maxCoeff(), vertices_input(Eigen::indexing::all, 2).maxCoeff();
    bb_range << minXYZ(0), maxXYZ(0), minXYZ(1), maxXYZ(1), minXYZ(2), maxXYZ(2);

    boundingbox.initialize(bb_range);
}

// Substrate block bounding box
void polygon::initialize(std::string box_type, Eigen::MatrixXd bounding_box) {
    Eigen::VectorXd rangespec(6);
    Eigen::MatrixXd vertices(8, 3);
    Eigen::MatrixXd faces(12, 3);
    Eigen::Vector3d origin, range;
    if (box_type == "cuboid") {
        rangespec << bounding_box(0, 0), bounding_box(1, 0), bounding_box(0, 1), bounding_box(1, 1), bounding_box(0, 2), bounding_box(1, 2);
        faces << 1, 3, 4, 1, 4, 2, 5, 6, 8, 5, 8, 7, 2, 4, 8, 2, 8, 6,
            1, 5, 7, 1, 7, 3, 1, 2, 6, 1, 6, 5, 3, 7, 8, 3, 8, 4;
        vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1;
        vertices = vertices.array() - 0.5;
        origin = (bounding_box(1, Eigen::indexing::all) + bounding_box(0, Eigen::indexing::all)) * 0.5;
        range = bounding_box(1, Eigen::indexing::all) - bounding_box(0, Eigen::indexing::all);

        vertices = vertices.array().rowwise() * range.transpose().array();
        vertices = vertices.rowwise() + origin.transpose();

        // Assigning vertices and faces
        Vertices = vertices;
        Vertices_t = vertices.transpose();
        Faces = faces;
        Faces_t = faces.transpose();

        // Number of faces and vertices
        n_vertices = Vertices.rows();
        n_faces = Faces.rows();

        // Compute mesh property
        polygon::volume_compute(vertices, faces);
        polygon::surface_compute(vertices, faces);

        // Computing the bounding box range
        bb_range.resize(6);
        minXYZ << vertices(Eigen::indexing::all, 0).minCoeff(), vertices(Eigen::indexing::all, 1).minCoeff(), vertices(Eigen::indexing::all, 2).minCoeff();
        maxXYZ << vertices(Eigen::indexing::all, 0).maxCoeff(), vertices(Eigen::indexing::all, 1).maxCoeff(), vertices(Eigen::indexing::all, 2).maxCoeff();
        bb_range << minXYZ(0), maxXYZ(0), minXYZ(1), maxXYZ(1), minXYZ(2), maxXYZ(2);

        boundingbox.initialize(bb_range);

    }
}

void polygon::volume_compute(Eigen::MatrixXd& vertices_input, Eigen::MatrixXd& faces_input) {

    // Compute the volume and surface areas of the tetrahedron
    vertices_mean << vertices_input(Eigen::indexing::all, 0).mean(), vertices_input(Eigen::indexing::all, 1).mean(), vertices_input(Eigen::indexing::all, 2).mean();
    vertices_shifted = vertices_input.rowwise() - vertices_mean.transpose();
    // Since each vertex is associated with different face
    volume_tetra.resize(n_faces);
    Eigen::Matrix3d tetra; // To compute the volume of each tetrahedron
    for (int i = 0; i < n_faces; i++) {
        tetra = vertices_shifted({ faces_input(i, 0) - 1, faces_input(i, 1) - 1, faces_input(i, 2) - 1 }, Eigen::indexing::all);
        volume_tetra(i) = tetra.determinant() / 6;
    }

    volume = volume > 0 ? volume : -volume;
}

bool polygon::containsPoint(Eigen::Vector3d& point, unsigned int& seed) {
    // A polyhedron contains a point if and only if a ray eminating from that point
    // intersects the faces of the polyhedron an odd number of times.

    if (not((point.array() >= minXYZ.array()) and (point.array() <= maxXYZ.array())).all()) {
        return false;
    }

    // Ensure vertices has been extracted
    if (not(vertices_extracted)) {
        // For intersection functions
        V1 = polygon::get_vertices(0);
        V2 = polygon::get_vertices(1);
        V3 = polygon::get_vertices(2);
        vertices_extracted = true;
    }

    Eigen::VectorXd dists = (Vertices.rowwise() - point.transpose()).rowwise().norm();
    double maxdist = dists.maxCoeff();

    bool odd;
    bool certain = false;
    int counter = 0;
    // Not sure if the monte carlo seed is needed
    std::mt19937 gen(seed);

    // std::cout << dists.transpose() << std::endl;

    while (not(certain)) {
        counter++;
        if (counter > 50) {
            printf("Polyhedron:containsPoint:counter', 'Too many attempts\n");
            break;
        }

        // determine ray
        std::uniform_real_distribution<> direction(0, 1);
        Eigen::Vector3d dir;
        dir << direction(gen), direction(gen), direction(gen); // three random variables from
        dir = dir.array() - 0.5; // pick random direction
        dir = dir.array() / dir.norm(); // unit vector
        dir = dir * maxdist * 10; // long enough to fully pass through polyhedron

        // compute intersections and test
        intersection_ray_info ray = polygon::TriangleRayIntersection(point, dir);
        int nIntersects = ray.intersect.cast<double>().sum(); // Number of intersections
        odd = nIntersects % 2 > 0; // inside if odd

        // make sure ray stays away fron surface triangle edges
        certain = polygon::intersection_is_certain(ray, false, 1e-6);

    }

    return odd;
}

intersection_info polygon::intersection(Eigen::Vector3d& orig, Eigen::Vector3d& dir) {
    intersection_info output;

    // Ensure vertices has been extracted
    if (not(vertices_extracted)) {
        // For intersection functions
        V1 = polygon::get_vertices(0);
        V2 = polygon::get_vertices(1);
        V3 = polygon::get_vertices(2);
        vertices_extracted = true;
    }

    intersection_ray_info ray = polygon::TriangleRayIntersection(orig, dir);
    if (not(ray.intersect.any())) {
        return output;
    }

    // make sure ray stays away from face edges (includes vertices)
    if (not(polygon::intersection_is_certain(ray, true, 1e-6))) {
        throw std::logic_error("Polyhedron:intersection:uncertain', 'Too close to edge/vertex/face\n");
        // std::cout << "Polyhedron:intersection:uncertain', 'Too close to edge/vertex/face" << std::endl;
    }

    // find closest intersection and get info
    Eigen::VectorXd found_t = polygon::find_element(ray.t, ray.intersect);
    if (not(polygon::is_unique(found_t))) {
        std::cout << "Polyhedron:intersection:duplicate, Two equal t found" << std::endl;
    }

    if (found_t.size() == 0) {
        return output;
    }
    else {
        double min_t = found_t.minCoeff();
        int min_faceIDs = polygon::find_first_index(found_t, min_t);
        Eigen::VectorXd found_ID = polygon::find_index(ray.intersect);
        int faceID = found_ID(min_faceIDs);

        // Store
        output.empty = false;
        output.t = min_t;
        output.vertices.resize(3, 3);
        output.vertices(0, Eigen::indexing::all) = V1(faceID, Eigen::indexing::all);
        output.vertices(1, Eigen::indexing::all) = V2(faceID, Eigen::indexing::all);
        output.vertices(2, Eigen::indexing::all) = V3(faceID, Eigen::indexing::all);

    }

    return output;
}

void polygon::surface_compute(Eigen::MatrixXd& vertices_input, Eigen::MatrixXd& faces_input) {

    // Initialize area
    area_surface.resize(n_faces);
    Eigen::Matrix3d xyz, xy, yz, zx;
    Eigen::Vector3d ones = Eigen::VectorXd::Ones(3);
    Eigen::Vector3d x, y, z;

    xy = Eigen::Matrix3d::Ones(3, 3);
    yz = Eigen::Matrix3d::Ones(3, 3);
    zx = Eigen::Matrix3d::Ones(3, 3);


    for (int i = 0; i < n_faces; i++) {
        xyz = vertices_input({ faces_input(i, 0) - 1, faces_input(i, 1) - 1, faces_input(i, 2) - 1 }, { 0,1,2 });
        x = xyz(Eigen::indexing::all, 0);
        y = xyz(Eigen::indexing::all, 1);
        z = xyz(Eigen::indexing::all, 2);

        // x-y principal plane
        xy(Eigen::indexing::all, 0) = x;
        xy(Eigen::indexing::all, 1) = y;

        // y-z principal plane
        yz(Eigen::indexing::all, 0) = y;
        yz(Eigen::indexing::all, 1) = z;

        // z-x principal plane
        zx(Eigen::indexing::all, 0) = z;
        zx(Eigen::indexing::all, 1) = x;

        area_surface(i) = 0.5 * std::sqrt(std::pow(xy.determinant(), 2) + std::pow(yz.determinant(), 2) + std::pow(zx.determinant(), 2));

    }
    surface_area = area_surface.sum();
}

Eigen::MatrixXd polygon::get_vertices(int column) {
    return Vertices(Faces(Eigen::indexing::all, column).array() - 1, Eigen::indexing::all);
}

bool polygon::intersection_is_certain(intersection_ray_info& ray, bool test_end = false, double eps = 1e-6) {

    Eigen::MatrixXd bary(ray.intersect.rows(), 3);
    bary(Eigen::indexing::all, 0) = ray.u;
    bary(Eigen::indexing::all, 1) = ray.v;
    bary(Eigen::indexing::all, 2) = 1 - ray.u.array() - ray.v.array();

    // std::cout << polygon::find_element(bary, ray.intersect).array().abs() << std::endl;
    Eigen::VectorXd min_abs_bary = (polygon::find_element(bary, ray.intersect).array().abs()).rowwise().minCoeff();
    Eigen::VectorXd abs_t0_int = (polygon::find_element(ray.t, ray.intersect).array().abs());
    Eigen::VectorXd abs_t1_int = (1 - polygon::find_element(ray.t, ray.intersect).array()).abs();

    bool certain = (min_abs_bary.array() > eps).all() and (abs_t0_int.array() > eps).all();
    if (test_end) {
        certain = certain and (abs_t1_int.array() > eps).all();
    }

    return certain;
}

intersection_ray_info polygon::TriangleRayIntersection(Eigen::Vector3d& point, Eigen::Vector3d& direction) {
    // Modified version of TriangleRayIntersection.
    //
    // Original Author:
    //    Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
    //
    // License: BSD license (http://en.wikipedia.org/wiki/BSD_licenses)

    intersection_ray_info output;

    int Nverts = V1.rows();
    Eigen::MatrixXd orig = Eigen::MatrixXd::Ones(Nverts, 3).array().rowwise() * point.transpose().array();
    Eigen::MatrixXd dir = Eigen::MatrixXd::Ones(Nverts, 3).array().rowwise() * direction.transpose().array();

    // tolerances
    double eps = 1e-20;
    double zero = 0.0;

    output.intersect.resize(Nverts, 1);
    output.intersect.fill(false);
    // output.intersect = Eigen::VectorXd::Zero(Nverts);
    output.t = Eigen::VectorXd::Ones(Nverts) * -1; // infinite is replaced by -1 in matlab
    output.u = Eigen::VectorXd::Ones(Nverts) * -1;
    output.v = Eigen::VectorXd::Ones(Nverts) * -1;

    // some pre-calculations
    Eigen::MatrixXd edge1 = V2 - V1; // find vectors for two edges sharing V1
    Eigen::MatrixXd edge2 = V3 - V1;
    Eigen::MatrixXd tvec = orig - V1; // vector from V1 to ray origin

    Eigen::MatrixXd pvec = crossMat(dir, edge2);
    Eigen::MatrixXd qvec = crossMat(tvec, edge1);
    Eigen::VectorXd det = (edge1.array() * pvec.array()).rowwise().sum(); // determinant of the matrix M = dot(edge1, pvec)

    // find faces parallel to the ray
    // Eigen::Array<bool, Eigen::Dynamic, 1> angleOK = (det.array().abs2() > eps).array(); // if det ~ 0 then ray lies in the triangle plane   
    Eigen::Array<bool, Eigen::Dynamic, 1> angleOK = (det.array().abs() > eps).array(); // if det ~ 0 then ray lies in the triangle plane   

    if ((angleOK == false).all()) {
        printf("polygon::TriangleRayIntersection::angle not within tolerance\n");
        return output;
    }
    // To do : change to avoid division by zero

    // calculate all variables for all line/triangle pairs
    output.u = (tvec.array() * pvec.array()).rowwise().sum().array() / det.array();
    output.v = (dir.array() * qvec.array()).rowwise().sum().array() / det.array();
    output.t = (edge2.array() * qvec.array()).rowwise().sum().array() / det.array();

    // test if line/plane intersection is within the triangle
    Eigen::Array<bool, Eigen::Dynamic, 1> ok = (angleOK and (output.u.array() >= -zero).array() and (output.v.array() >= -zero).array() and ((output.u + output.v).array() <= 1.0 + zero).array());
    // std::cout << ok.cast<double>().sum() << std::endl;
    output.intersect = (ok and (output.t.array() >= -zero).array() and (output.t.array() <= 1.0 + zero).array());

    return output;
}

// Matrix cross product of every row of a and b
Eigen::MatrixXd polygon::crossMat(Eigen::MatrixXd& a, Eigen::MatrixXd& b) {
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(a.rows(), 3);
    for (int i = 0; i < a.rows(); i++) {
        output(i, { 0, 1, 2 }) = a(i, { 0, 1, 2 }).cross(b(i, { 0, 1, 2 }));
    }
    return output;
}

Eigen::Vector3d polygon::get_minXYZ() {
    return minXYZ;
}

// very Stupid finding element method
Eigen::MatrixXd polygon::find_element(Eigen::MatrixXd& target, Eigen::Array<bool, Eigen::Dynamic, 1>& input) {
    int rows = input.cast<double>().sum();
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, 3);
    if (rows != 0) {
        int counter = 0;
        for (int i = 0; i < input.rows(); i++) {
            if (input(i)) {

                output(counter, Eigen::indexing::all) = target(i, Eigen::indexing::all);
                counter++;

            }
        }
    }
    return output;
}

// very Stupid finding element method
Eigen::VectorXd polygon::find_element(Eigen::VectorXd& target, Eigen::Array<bool, Eigen::Dynamic, 1>& input) {
    int rows = input.cast<double>().sum();
    Eigen::VectorXd output = Eigen::VectorXd::Zero(rows);
    if (rows != 0) {
        int counter = 0;
        for (int i = 0; i < input.rows(); i++) {
            if (input(i)) {
                output(counter) = target(i);
                counter++;
            }
        }
    }
    return output;
}

// very Stupid finding element method
Eigen::VectorXd polygon::find_index(Eigen::Array<bool, Eigen::Dynamic, 1>& input) {
    int rows = input.cast<double>().sum();
    Eigen::VectorXd output = Eigen::VectorXd::Zero(rows);
    if (rows != 0) {
        int counter = 0;
        for (int i = 0; i < input.rows(); i++) {
            if (input(i)) {
                output(counter) = i;
                counter++;
            }
        }
    }
    return output;
}

int polygon::find_first_index(Eigen::VectorXd& input, double& target) {
    int output;
    for (int i = 0; i < input.rows(); i++) {
        if (input(i) == target) {
            output = i;
            break;
        }
    }

    return output;
}

// Stupid method to check if vector is unique
bool polygon::is_unique(Eigen::VectorXd& input) {
    bool flag = true;
    for (int i = 0; i < input.rows() - 1; i++) {
        for (int j = i + 1; j < input.rows(); j++) {
            if (input(i) == input(j)) {
                flag = false;
                break;
            }
        }
    }

    return flag;
}