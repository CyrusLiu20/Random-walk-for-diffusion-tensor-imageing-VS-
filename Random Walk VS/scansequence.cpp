#include <cmath>
#include <iostream>
#include "sequencefile.h"
#include "sequence.h"
#include "scansequence.h"

sequence ScanSequence::Scan(Eigen::VectorXd dt_input, Eigen::MatrixXd gG_input) {
    sequence sequence_result;

    sequence_result.N = dt_input.rows();
    sequence_result.input_dt(dt_input);
    sequence_result.input_gG(gG_input);

    // Temporary containers to calculate b values
    Eigen::VectorXd t = ScanSequence::cumsum(dt_input);
    Eigen::VectorXd k = ScanSequence::cumtrapz(t, gG_input);
    Eigen::VectorXd k_squared = k.array().pow(2);
    sequence_result.bvalue = ScanSequence::trapz(t, k_squared);
    // sequence_result.bvalue = ScanSequence::trapz(t, k.array().pow(2));


    return sequence_result;
}

sequence ScanSequence::create(sequencefile seq) {

    // Matlab description will change later
    // % make a sequence object from the given parameters
    // %   NT : number of time steps (target)
    // %   dt_max : maximum time step
    // %       (1, 2) = (dt_free, dt_grad) - in case of real sequence
    // %       scalar = dt - in case of dummy sequence
    // %   SeqName : sequence name {PGSE, MCSE/M2SE, STEAM}
    // %   varargin : parameters specific to the given sequence with SeqName

    if (seq.type == "PGSE" or seq.type == "STEAM") {
        Eigen::VectorXd durations(9), ids(9);
        durations << seq.alpha90, seq.epsilon, seq.delta, seq.epsilon,
            seq.Delta - (2 * seq.epsilon + seq.delta),
            seq.epsilon, seq.delta, seq.epsilon, seq.alphaRO;
        ids << 0, 1, 2, 3, 0, -1, -2, -3, 0;
        ScanSequence::discretize(durations, ids, seq.N_t, seq.dt_max_free, seq.dt_max_grad, seq.Gmax);
    }
    else if (seq.type == "MCSE" or seq.type == "M2SE") {
        Eigen::VectorXd durations(15), ids(15);
        double del1 = seq.delta1 + 2 * seq.epsilon;
        double del2 = seq.delta2 + 2 * seq.epsilon;
        seq.Delta = (del2 * (-2 * del1 + seq.epsilon) + del1 * seq.epsilon) / (del1 - del2);
        durations << seq.alpha90, seq.epsilon, seq.delta1, seq.epsilon, seq.epsilon,
            seq.delta2, seq.epsilon, seq.Delta - (del1 + del2), seq.epsilon,
            seq.delta2, seq.epsilon, seq.epsilon, seq.delta1, seq.epsilon, seq.alphaRO;
        ids << 0, 1, 2, 3, -1, -2, -3, 0, 1, 2, 3, -1, -2, -3, 0;
        ScanSequence::discretize(durations, ids, seq.N_t, seq.dt_max_free, seq.dt_max_grad, seq.Gmax);
    }
    else {
        seq.gamma = 1;
        dt = Eigen::VectorXd::Ones(seq.N_t) * seq.dt_max_free;
        gG = Eigen::MatrixXd::Zero(seq.N_t, 1);
    }

    gG = gG * seq.gamma;

    sequence_output = ScanSequence::Scan(dt, gG);

    return sequence_output;
}

void ScanSequence::discretize(Eigen::VectorXd durations, Eigen::VectorXd ids, int NT, double dt_max_free, double dt_max_grad, double Gmax) {
    // Matlab description
    // % discretize the sequence
    // %   durations := length of intervals
    // %   ids := designtaion of interval
    // %       0 - flat (free, gradient OFF)
    // %       % gradients: +/- sign of id represents sign of Gmax at the end of the gradient
    // %       1 - (+/-) gradient ramp-up
    // %       2 - (+/-) gradient flat
    // %       3 - (+/-) gradient ramp-down
    Eigen::VectorXd Nt_intervals(ids.rows());

    // calculate the target time step dt
    double dt_aim, dt_free, dt_grad;
    dt_aim = durations.sum() / NT;
    dt_free = std::fmin(dt_aim, dt_max_free);
    dt_grad = std::fmin(dt_aim, dt_max_grad);

    for (int i = 0; i < ids.rows(); i++) {
        if (ids(i) == 0) {
            Nt_intervals(i) = std::ceil(durations(i) / dt_free);
        }
        else {
            Nt_intervals(i) = std::ceil(durations(i) / dt_grad);
        }
    }

    dt = Eigen::VectorXd::Zero(0);
    gG = Eigen::MatrixXd::Zero(0, 0);

    // Temporarily store data
    double gA, gB;
    for (int i = 0; i < durations.rows(); i++) {
        double Nt_i = Nt_intervals(i);
        double dt_i = durations(i) / Nt_i;
        // For repetition
        Eigen::VectorXd dt_temp = dt;
        Eigen::VectorXd dt_rep = Eigen::VectorXd::Ones(Nt_i) * dt_i;
        dt.resize(dt_temp.rows() + Nt_i);
        dt << dt_temp, dt_rep;

        // gradient
        double id_i = ids(i);
        switch (std::abs((int)id_i)) { // use absolute value for switch and use sign(id_i) inside cases
        case 0: // flat, gradient off
            gA = 0;
            gB = gA;
            break;
        case 1: // ramp-up gradient
            gA = 0;
            gB = ScanSequence::sign(id_i) * Gmax;
            break;
        case 2: // flat, gradient on
            gA = ScanSequence::sign(id_i) * Gmax;
            gB = gA;
            break;
        case 3: // ramp-down gradient
            gA = ScanSequence::sign(id_i) * Gmax;
            gB = 0;
            break;
        default:
            printf("ScanSequence::Ids input out of bound, please check the values of Ids");
        }
        Eigen::MatrixXd gG_temp = gG;
        Eigen::MatrixXd gvals = Eigen::VectorXd::LinSpaced(Nt_i, gA, gB - (gB - gA) / (Nt_i));
        // gvals.resize(Nt_i, 1); // Removing the last element
        gG.resize(gG_temp.rows() + (Nt_i), 1);
        if (i == 0) {
            gG << gvals;
        }
        else {
            gG << gG_temp, gvals;
        }
    }
}

// Custom built sign function
int ScanSequence::sign(double input) {
    if (input > 0) {
        return 1;
    }
    else if (input < 0) {
        return -1;
    }
    else {
        return 0;
    }
}

// Custom built cumulative sum function
Eigen::VectorXd ScanSequence::cumsum(Eigen::VectorXd input) {
    Eigen::VectorXd output(input.rows());
    output(0) = input(0);
    for (int i = 1; i < input.rows(); i++) {
        output(i) = input(i) + output(i - 1);
    }

    return output;
}

// Custom built cumulative sum function
Eigen::MatrixXd ScanSequence::cumtrapz(Eigen::VectorXd x, Eigen::MatrixXd y) {
    Eigen::VectorXd output(x.rows());
    output(0) = 0;
    if (y.cols() == 1) {
        if (y.rows() == x.rows()) {
            for (int i = 1; i < x.rows(); i++) {
                output(i) = (y(i) + y(i - 1)) * (x(i) - x(i - 1)) / 2 + output(i - 1);
            }
        }
        else {
            printf("ScanSequence::cumtrapz::input rows differ");
        }
    }
    else {
        printf("ScanSequence::cumtrapz::Matrix column number not equal to 1");
    }
    return output;
}

// Custom built trapezoid function
double ScanSequence::trapz(Eigen::VectorXd x, Eigen::MatrixXd y) {
    double sum = 0;
    if (y.cols() == 1) {
        if (y.rows() == x.rows()) {
            for (int i = 1; i < x.rows(); i++) {
                sum += (y(i) + y(i - 1)) * (x(i) - x(i - 1)) / 2;
            }
        }
        else {
            printf("ScanSequence::trapz::input rows differ");
        }
    }
    else {
        printf("ScanSequence::trap::Matrix column number not equal to 1");
    }
    return sum;
}



sequence ScanSequence::get_sequence() {
    return sequence_output;
}

Eigen::VectorXd ScanSequence::get_dt() {
    return dt;
}

Eigen::MatrixXd ScanSequence::get_gG() {
    return gG;
}
