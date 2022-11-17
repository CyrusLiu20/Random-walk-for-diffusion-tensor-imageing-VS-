#include <iostream>
#include "process_signal.h"
#include <complex.h>
//#include <Eigen/MatrixFunctions>
#include <Eigen/LU>

// Post processing to create diffusion tensor and diffusion coefficients
process_signal::process_signal(std::vector<particle_state>& states, substrate& substrate_input, sequence& sequence_input, Eigen::MatrixXd& pos0) {

    // Initializing
    N_p = states.size();

    phase_all = Eigen::MatrixXd::Zero(N_p, 3);
    position_all = Eigen::MatrixXd::Zero(N_p, 3);
    displacement = Eigen::MatrixXd::Zero(N_p, 3);

    valid.resize(N_p, 1);
    valid.fill(false);
    insideVoxel.resize(N_p, 1);
    insideVoxel.fill(false);
    insideECS.resize(N_p, 1);
    insideECS.fill(false);

    // Check the status of the final position of each particle 
    for (int i = 0; i < N_p; i++) {
        phase_all(i, Eigen::indexing::all) = states[i].phase.transpose();
        position_all(i, Eigen::indexing::all) = states[i].position.transpose();

        // Check if there are any flags
        valid(i) = states[i].flag == 0;

        // read-out happens inside voxel only
        Eigen::Vector3d temp = position_all(i, Eigen::indexing::all);
        insideVoxel(i) = (substrate_input.voxel.containsPoint(temp));

        // Whether particle is ICS or ECS
        insideECS(i) = states[i].myoindex < 0;
    }

    // Sorting particles according to particle final position
    phase_ECS = process_signal::find_element(phase_all, (insideECS and insideVoxel and valid));
    phase_ICS = process_signal::find_element(phase_all, (not(insideECS) and insideVoxel and valid));
    displacement_all = position_all - pos0;
    displacement_ECS = process_signal::find_element(displacement_all, (insideECS and valid));
    displacement_ICS = process_signal::find_element(displacement_all, (not(insideECS) and valid));

    phase = process_signal::find_element(phase_all, (insideVoxel and valid));
    displacement = process_signal::find_element(displacement_all, (valid));


    // calculate diffusion tensor
    bvalue = sequence_input.bvalue;
    if (phase.rows() != 0) {
        tensor = process_signal::process_phase(phase, bvalue);
    }
    else {
        no_tensor = true;
        std::cout << "process_signal::there are no valid phase to compute diffusion tensor" << std::endl;
    }

    // calculate bulk diffusivity
    double T = sequence_input.get_dt_entire().sum(); // NOT equal to Delta, because we don't have that displ data
    if (displacement.rows() != 0) {
        process_signal::process_displacement(displacement, T);
    }
    else {
        std::cout << "process_signal::there are no valid displacement to compute diffusion tensor" << std::endl;
    }


}

Eigen::Matrix3d process_signal::process_phase(Eigen::MatrixXd& phase, double& bvalue) {
    // get diffusion tensor

    // gradient sampling directions
    Eigen::MatrixXd directions(6, 3);
    directions << 1, 1, 0,
        1, -1, 0,
        1, 0, 1,
        1, 0, -1,
        0, 1, 1,
        0, 1, -1;

    int ndirs = directions.rows();
    Eigen::MatrixXd b_matrix = Eigen::MatrixXd::Zero(ndirs, 6);
    Eigen::VectorXd signal_ratio = Eigen::VectorXd::Zero(ndirs, 1);
    for (int i = 0; i < ndirs; i++) {
        Eigen::Vector3d dir = directions(i, Eigen::indexing::all); // [Gx, Gy, Gz]
        // b-matrix (b_xx, b_yy, b_zz, b_xy, b_xz, b_yz)
        Eigen::VectorXd temporary(6);
        temporary << dir(0) * dir(0) * bvalue, dir(1)* dir(1)* bvalue, dir(2)* dir(2)* bvalue, 2 * dir(0) * dir(1) * bvalue, 2 * dir(0) * dir(2) * bvalue, 2 * dir(1) * dir(2) * bvalue;
        b_matrix(i, Eigen::indexing::all) = temporary.transpose();
        // attenuation vector
        Eigen::VectorXcd phi = (phase.array().rowwise() * dir.transpose().array()).rowwise().sum(); // combine components
        std::complex<double> imag(0, 1);
        signal_ratio(i) = std::abs((phi * -imag).array().exp().mean());

    }
    // least squares solution
    Eigen::MatrixXd A = b_matrix;
    Eigen::VectorXd b = signal_ratio.array().log();
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b) * -1; // solve

    Eigen::Matrix3d tensor_local;
    tensor_local << x(0), x(3), x(4),
        x(3), x(1), x(5),
        x(4), x(5), x(2);

    return tensor_local;
}

void process_signal::process_displacement(Eigen::MatrixXd& displacement, double& T) {
    // get diffusivity from RMS displacement
    int dim = displacement.cols();
    MD = (displacement.array().square().rowwise().sum()).mean() / (2 * dim * T);
    Dx = (displacement(Eigen::indexing::all, 0).array() * displacement(Eigen::indexing::all, 0).array()).mean() / (2 * T);
    Dy = (displacement(Eigen::indexing::all, 1).array() * displacement(Eigen::indexing::all, 1).array()).mean() / (2 * T);
    Dz = (displacement(Eigen::indexing::all, 2).array() * displacement(Eigen::indexing::all, 2).array()).mean() / (2 * T);

}

// very Stupid finding element method
Eigen::MatrixXd process_signal::find_element(Eigen::MatrixXd& target, Eigen::Array<bool, Eigen::Dynamic, 1> input) {
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