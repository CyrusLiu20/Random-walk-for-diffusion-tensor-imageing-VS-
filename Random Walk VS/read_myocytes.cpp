#include "read_myocytes.h"
#include "matio.h"

// type might be a problem
read_myocytes::read_myocytes(std::string filename) {

    // Compile the full path
    //std::string fullpath = "myocytes/" + filename;
    std::string fullpath = filename;
    const char* filepath = fullpath.c_str();

    // Read Mat file
    geometry = Mat_Open(filepath, MAT_ACC_RDONLY);
    matvar_t* matVar, * vertex, * face;

    matVar = Mat_VarReadInfo(geometry, "myocytes");
    if (matVar) { // If the variable exist

        // Number of myocytes
        N_m = matVar->dims[1];
        // Resizing Vertices and Faces to fit all myocytes
        Vertices.resize(N_m);
        Faces.resize(N_m);

        // Looping through all myocytes to retrieve vertices and faces
        for (int i = 0; i < N_m; i++) {
            vertex = Mat_VarGetStructFieldByName(matVar, "Vertices", i);
            face = Mat_VarGetStructFieldByName(matVar, "Faces", i);

            // Reading vertices
            bool vertex_read_error = Mat_VarReadDataAll(geometry, vertex);
            if (not(vertex_read_error)) {
                Vertices[i] = read_myocytes::vertex_conversion(vertex);
            }
            else {
                printf("read_myocytes::unable to read vertex");
                break;
            }

            // Reading Faces
            bool face_read_error = Mat_VarReadDataAll(geometry, face);
            if (not(face_read_error)) {
                Faces[i] = read_myocytes::face_conversion(face);
            }
            else {
                printf("read_myocytes::unable to read face");
                break;
            }

        }

        // myocytes file has been successfully loaded
        flag = true;

    }
    else {
        printf("read_myocytes::non existent variable, please check your .mat file\n");
    }

    Mat_Close(geometry);
}

Eigen::MatrixXd read_myocytes::vertex_conversion(matvar_t* field) {

    Eigen::MatrixXd output;
    int rows, cols;
    unsigned field_size;
    const double* Data = static_cast<const double*>(field->data);


    // Specify the matrix
    if (field->rank == 2) {
        rows = field->dims[0];
        cols = field->dims[1];
        field_size = field->nbytes / field->data_size;

        output.resize(rows, cols);
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) {
                output(i, j) = Data[rows * j + i];
            }
        }

    }
    else {
        printf("Matrix rank not equal to 2");
    }

    return output;
}

Eigen::MatrixXd read_myocytes::face_conversion(matvar_t* field) {

    Eigen::MatrixXd output;
    int rows, cols;
    unsigned field_size;
    const uint16_t* Data = static_cast<const uint16_t*>(field->data);


    // Specify the matrix
    if (field->rank == 2) {
        rows = field->dims[0];
        cols = field->dims[1];
        field_size = field->nbytes / field->data_size;

        output.resize(rows, cols);
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) {
                output(i, j) = Data[rows * j + i];
            }
        }

    }
    else {
        printf("Matrix rank not equal to 2");
    }

    return output;
}

bool read_myocytes::scanned() {
    return flag;
}

std::vector<Eigen::MatrixXd> read_myocytes::get_Vertices() {
    return Vertices;
}

std::vector<Eigen::MatrixXd> read_myocytes::get_Faces() {
    return Faces;
}
