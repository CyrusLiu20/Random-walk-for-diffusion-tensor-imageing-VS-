#ifndef MYOCYTES_H
#define MYOCYTES_H

#include <vector>
#include <Eigen/Dense>
#include "read_myocytes.h"

class myocytes {

public:

    myocytes() = default;

    myocytes(read_myocytes scanner) {
        load(scanner);
    };

    void load(read_myocytes scanner) {
        if (scanner.scanned()) {
            Vertices = scanner.get_Vertices();
            Faces = scanner.get_Faces();
        }
        else {
            printf("myocytes::myocytes file has not been loaded\n");
        }
    };

    std::vector<Eigen::MatrixXd> Vertices;
    std::vector<Eigen::MatrixXd> Faces;
};

#endif
