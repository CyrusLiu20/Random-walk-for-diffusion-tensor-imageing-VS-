#include <algorithm>
#include <iostream>
#include "substrate.h"
#include <Eigen/Dense>

substrate::substrate(substratefile substrate) {

    // Voxel size
    y_extent << substrate.yextent_min, substrate.yextent_max;
    dxdydz << substrate.block_Lx, substrate.block_Ly, substrate.block_Lz;
    Ly = substrate.block_Ly;

    // y slicing
    Eigen::VectorXd low = substrate::iteration(0, substrate.yextent_min, Ly);
    Eigen::VectorXd high = substrate::iteration(0, substrate.yextent_max, Ly);
    y_minvals.resize(low.rows() + high.rows() - 1);
    y_minvals << low(Eigen::seq(0, Eigen::indexing::last - 1)), high;

    y_slice_minmax.resize(2, y_minvals.rows() - 1);
    y_slice_minmax(0, Eigen::indexing::all) = y_minvals(Eigen::seq(0, Eigen::indexing::last - 1));
    y_slice_minmax(1, Eigen::indexing::all) = y_minvals(Eigen::seq(1, Eigen::indexing::last));

    // Type of substrate partial or full
    substrate_type = substrate.substrate_type;

    // Initializing voxel
    voxel.initialize(substrate.voxel_range);
    transit_model = substrate.transit_model;
    kappa = substrate.permeability;
    D_i = substrate.D_ics;
    D_e = substrate.D_ecs;
    dim = substrate.dimension;

    // Copy transformation
    // yextent_min = substrate.yextent_min;
    // yextent_max = substrate.yextent_max;

    // deg_rot_per_L_in_y = substrate.deg_rot_per_L_in_y; // degree of rotation per unit length in y
    // z_amplitude = substrate.z_amplitude; // amplitude in z for the sinusoidal transform
    // x_frequency = substrate.x_frequency; // frquency in x for the sinusoidal transform
    // shift_block = substrate.shift_block;  


    substrate_transform.y_slice_minmax = y_slice_minmax;

    substrate_transform.deg_rot_per_L_in_y = substrate.deg_rot_per_L_in_y; // degree of rotation per unit length in y
    substrate_transform.z_amplitude = substrate.z_amplitude; // amplitude in z for the sinusoidal transform
    substrate_transform.x_frequency = substrate.x_frequency; // frquency in x for the sinusoidal transform
    substrate_transform.shift_block = substrate.shift_block;

    if (substrate_type == "block") {
        boundary = "periodic";
        substrate_transform.isIdentity = false;
        substrate_transform.dxdydz_bb = dxdydz;
    }

    // Creating substrate block bounding box

    // Parsing in Matlab (????????)
}

// DEBUG: try cross product (surface compute)? 
// Transform the raw matlab files into a standard vector of polyhedron
void substrate::myocytes_scan(myocytes myocytes_input) {
    N_m = myocytes_input.Vertices.size();
    myos.resize(N_m);

    for (int i = 0; i < N_m; i++) {
        polygon polygon_i(myocytes_input.Vertices[i], myocytes_input.Faces[i]);
        myos[i] = polygon_i;
        myos[i].bytes = sizeof(polygon_i);
    }

}

// Builds the bounding box for the substrate blok
void substrate::buildcache() {

    // Concatonating the boundary range of each bounding box into a huge vector
    Eigen::VectorXd tmp_bbs;
    tmp_bbs.resize(6);
    myocytes_bbrange.resize(N_m * 6);
    // myocytes_bbrange = Eigen::VectorXd::Zero(N_m*6);
    for (int i = 0; i < N_m; i++) {
        tmp_bbs = myos[i].boundingbox.bb_range;
        myocytes_bbrange({ i * 6, i * 6 + 1, i * 6 + 2, i * 6 + 3, i * 6 + 4, i * 6 + 5 }) = tmp_bbs;
    }

    // Creating the substrate block bounding box
    Eigen::MatrixXd cuboid(2, 3);
    cuboid << 0, 0, 0, dxdydz(0), dxdydz(1), dxdydz(2);
    block_bb.initialize("cuboid", cuboid);
    block_bb.bytes = sizeof(block_bb);
}

int substrate::findMyocyte(Eigen::Vector3d position, unsigned int seed, std::string refFrame = "global") {
    Eigen::Vector3d position_local;
    int myoindex = -1;
    if (substrate_type == "block" and refFrame == "global") { // Some reason local reference frame is not specified, we'll see later
        // Work in progress
        try {
            position_local = substrate_transform.global2local(position).position_local;
        }
        catch (const std::exception& ex) {
            throw;
        }
    }

    // Eigen::Vector3d temp; // For debugging use
    // temp << 89.0938  ,-19.6363   ,62.1991; // Please remember to remove this
    myoindex = substrate::search_myocytes(position_local, seed);
    // myoindex = substrate::search_myocytes(temp, seed);
    // std::cout << "myoindex : " << myoindex << std::endl;

    return myoindex;
}

int substrate::search_myocytes(Eigen::Vector3d position_local, unsigned int& seed) {
    int myoindex = -1;
    bool inside; // Whether a particle is in a myocyte

    // Eigen::Vector3d temp; // For debugging use
    // temp << 97.3409,  381.8844,   60.4011; // Please remember to remove this

    for (int i_myocytes = 0; i_myocytes < N_m; i_myocytes++) {

        inside = myos[i_myocytes].containsPoint(position_local, seed);
        // inside = myos[i_myocytes].containsPoint(temp, seed);
        try {
            if (inside) {
                if (not(myoindex == -1)) {
                    // printf("Substrate:search_myocytes:multiple, Point cannot be inside multiple myocytes\n");
                    throw std::logic_error("substrate:search_myocytes:multiple, Point cannot be inside multiple myocytes");
                }

                myoindex = i_myocytes;
                // std::cout << "Position : " << position_local.transpose() << ", found in index : " << myoindex << std::endl;
            }
        }
        catch (const std::exception& ex) {
            throw;
        }

    }

    return myoindex;
}

intersection_info substrate::intersectMyocytes(Eigen::Vector3d& position, Eigen::Vector3d& dxdydz, std::string refFrame = "global") {
    intersection_info intersect_info;

    Eigen::Vector3d position_updated;
    if (substrate_type == "block" and refFrame == "global") {
        position_updated = substrate_transform.global2local(position).position_local;
    }
    else {
        position_updated = position;
    }

    // find all that need proper checking - can do it in local since we checked above
    Eigen::Array<bool, Eigen::Dynamic, 1> needsChecks = substrate::needsChecking(position_updated, dxdydz, "local"); // ~(wontEnter && isOutside)
    Eigen::VectorXd checks_indices = substrate::find_index(needsChecks);

    // std::cout << checks_indices.transpose() << std::endl;
    for (int i = 0; i < checks_indices.rows(); i++) {
        intersection_info info = myos[checks_indices(i)].intersection(position_updated, dxdydz);

        if (not(info.empty)) {
            info.myoindex = checks_indices(i);
            if (intersect_info.empty) {
                intersect_info = info;
            }
            else {
                if (info.t < intersect_info.t) {
                    intersect_info = info;
                }
                else if (info.t == intersect_info.t) {
                    std::cout << "Substrate:intersectMyocytes:duplicate, Two intersections found" << std::endl;
                }
                else {
                    // leave it
                }
            }
        }
    }


    return intersect_info;
}

Eigen::Array<bool, Eigen::Dynamic, 1> substrate::needsChecking(Eigen::Vector3d& position_updated, Eigen::Vector3d& dxdydz, std::string refFrame) {
    Eigen::Array<bool, Eigen::Dynamic, 1> output;
    output.resize(N_m);
    // To do : transformation

    for (int i = 0; i < N_m; i++) {
        double absdist_step = 1.01 * dxdydz.norm();

        // we will use a bigger box, whose size depends on the step (one where being inside means a step could potentially hit the body)
        double bb_xmin = myocytes_bbrange[0 + i * 6] - absdist_step;
        double bb_xmax = myocytes_bbrange[1 + i * 6] + absdist_step;
        double bb_ymin = myocytes_bbrange[2 + i * 6] - absdist_step;
        double bb_ymax = myocytes_bbrange[3 + i * 6] + absdist_step;
        double bb_zmin = myocytes_bbrange[4 + i * 6] - absdist_step;
        double bb_zmax = myocytes_bbrange[5 + i * 6] + absdist_step;
        bool isInside = (position_updated[0] > bb_xmin && position_updated[0] < bb_xmax)  // inside in x, and
            && (position_updated[1] > bb_ymin && position_updated[1] < bb_ymax)  // inside in y, and
            && (position_updated[2] > bb_zmin && position_updated[2] < bb_zmax); // inside in z

        output(i) = isInside;
    }

    return output;
}

// very Stupid finding element method
Eigen::VectorXd substrate::find_index(Eigen::Array<bool, Eigen::Dynamic, 1>& input) {
    int rows = input.cast<double>().sum();
    Eigen::VectorXd output(rows);
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

// Stupid way to implement iteration (to calculate y_minval and y_slice_minmax)
Eigen::VectorXd substrate::iteration(double a, double b, double step) {
    Eigen::VectorXd output;
    int size = std::floor(substrate::absolute((a - b) / step)) + 1;
    output.resize(size);

    if (a > b) {
        for (int i = 0; i < size; i++) {
            output(i) = a - step * i;
        }
        return output.reverse();
    }
    else {
        for (int i = 0; i < size; i++) {
            output(i) = a + step * i;
        }
        return output;
    }
}

// custom stupid absolute function
double substrate::absolute(double input) {
    return input > 0 ? input : -input;
}