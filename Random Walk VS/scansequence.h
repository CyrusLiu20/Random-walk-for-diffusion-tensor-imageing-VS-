#ifndef SCANSEQUENCE_H
#define SCANSEQUENCE_H

#include <Eigen/Dense>
#include "sequence.h"

class ScanSequence {

public:

    ScanSequence() = default;

    sequence Scan(Eigen::VectorXd dt_input, Eigen::MatrixXd gG_input);
    sequence create(sequencefile seq);
    void discretize(Eigen::VectorXd durations, Eigen::VectorXd ids, int NT, double dt_max_free, double dt_max_grad, double Gmax);

    // Data retrieval
    sequence get_sequence();
    Eigen::VectorXd get_dt();
    Eigen::MatrixXd get_gG();

private:

    // custom sign function
    int sign(double input);
    // Custom built cumulative sumation
    Eigen::VectorXd cumsum(Eigen::VectorXd input);
    Eigen::MatrixXd cumtrapz(Eigen::VectorXd x, Eigen::MatrixXd input);
    double trapz(Eigen::VectorXd x, Eigen::MatrixXd y);

    sequence sequence_output;

    Eigen::VectorXd dt;
    Eigen::MatrixXd gG;

};

#endif