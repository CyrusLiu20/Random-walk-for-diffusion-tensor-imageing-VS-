#ifndef READ_MYOCYTE_H
#define READ_MYOCYTE_H

#include <vector>
#include <Eigen/Dense>
#include <matio.h>

class read_myocytes {

public:

    read_myocytes(std::string filename);

    // Data retrieval
    std::vector<Eigen::MatrixXd> get_Vertices();
    std::vector<Eigen::MatrixXd> get_Faces();
    bool scanned();

    // substratefile get_myo

private:

    // Converting Matlab matrix to Eigen matrix
    Eigen::MatrixXd vertex_conversion(matvar_t* field);
    Eigen::MatrixXd face_conversion(matvar_t* field);

    // Mat file pointer
    mat_t* geometry;

    // Whether the myocytes mat file has been scanned;
    bool flag = false;
    // Number of myocytes
    int N_m;

    // Vertices
    std::vector<Eigen::MatrixXd> Vertices;
    // Faces
    std::vector<Eigen::MatrixXd> Faces;

};

#endif
