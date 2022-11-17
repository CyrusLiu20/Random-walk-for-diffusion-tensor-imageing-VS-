#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <Eigen/Dense>

class sequence {

public:

    // Construtors
    sequence() = default;
    sequence(const sequence& pSrc) = default;	// copying class
    sequence(sequence&& pSrc) = default;	// moving class
    ~sequence() = default;

    // Operators overloading
    sequence& operator =(const sequence& pSrc) {
        // Transfer of resources
        if (this != &pSrc) {
            N = pSrc.N;
            bvalue = pSrc.bvalue;
            dt = pSrc.dt;
            gG = pSrc.gG;
        }
        return *this;
    };	// Transfer of resourses using equal sign

    // Number of dt
    int N;
    // What is this?
    double bvalue;

    // input data
    void input_dt(Eigen::VectorXd input) {
        dt = input;
    };

    void input_gG(Eigen::MatrixXd input) {
        gG = input;
    };

    // Data retrieval
    Eigen::VectorXd get_dt_entire() {
        return dt;
    };

    Eigen::MatrixXd get_gG_entire() {
        return gG;
    };

    double get_total_time() {
        return dt.sum();
    };

    double get_dt(int n) {
        return dt(n);
    };

    Eigen::VectorXd get_gG(int n) {
        return gG(n, Eigen::indexing::all);
    };

private:

    // Time step
    Eigen::VectorXd dt;
    // Something gradient
    Eigen::MatrixXd gG;

};

#endif
