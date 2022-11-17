#ifndef SEQUENCEFILE_H
#define SEQUENCEFILE_H

#include <cstdio>
#include <string>
class sequencefile {

public:

    // Construtor : initializing the type
    sequencefile() {
        if (type == "PGSE") {
            Gmax = 40.572;
            Delta = 20.547;
            epsilon = 0.676;
            delta = 9.921;
            alpha90 = 1;
            alphaRO = 1;
        }
        else if (type == "MCSE") {
            Gmax = 39.478;
            epsilon = 0.661;
            delta1 = 7.819;
            delta2 = 16.299;
            alpha90 = 1;
            alphaRO = 1;
        }
        else if (type == "STEAM") {
            Gmax = 35.690;
            Delta = 1000;
            epsilon = 0.595;
            delta = 0.977;
            alpha90 = 1;
            alphaRO = 1;
        }
        else {
            printf("sequencefile::type mismatched, please check your sequence type");
        }
    }

    std::string type = "MCSE";

    double Gmax;
    double Delta, delta1, delta2;
    double epsilon;
    double delta;
    double alpha90;
    double alphaRO;

    // rad/ms/mT (it contains 10^-6 because the gradient is in mT/m and positions in um)
    double gamma = 267.5e-6;
    // number of time steps
    int N_t = 100;
    // time step limit [dt_free, dt_grad]
    double dt_max_free = 0.01;
    double dt_max_grad = 0.01;

};

#endif