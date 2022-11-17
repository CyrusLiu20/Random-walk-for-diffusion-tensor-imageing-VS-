#include "testing.h"
#include <iostream>
#include <cmath>

testing::testing(std::string filename) {
    testing::initialize(filename);
}

bool testing::initialize(std::string filename) {
    // Compile the full path
    std::string fullpath = "testing/" + filename;
    const char* filepath = fullpath.c_str();

    // Read Mat file
    test = Mat_Open(filepath, MAT_ACC_RDONLY);

    flag = testing::read_fields();

    Mat_Close(test);

    return flag;
}

bool testing::read_fields() {
    matvar_t* matVar_dirs, * matVar_probs, * matVar_states, * matVar_pos, * matVar_phis, * matVar_pos0;// struct variable
    matvar_t* dirs, * probs; // field variable
    matVar_dirs = Mat_VarReadInfo(test, "dxdydz_normaldistrib_raw_struct");
    matVar_probs = Mat_VarReadInfo(test, "probability_of_transition_struct");
    // matVar_states = Mat_VarReadInfo(test, "walker");
    matVar_pos = Mat_VarRead(test, "positions");
    matVar_phis = Mat_VarRead(test, "phases");
    matVar_pos0 = Mat_VarRead(test, "pos0");


    if (matVar_dirs and matVar_probs and matVar_pos and matVar_phis) { // If the variable exist

        try {
            // Number of particles
            if ((matVar_dirs->dims[1] != matVar_probs->dims[1]) and (matVar_probs->dims[1] != matVar_pos->dims[1]) and (matVar_pos->dims[1] != matVar_phis->dims[1])) {
                throw std::runtime_error("testing::Testing file unequal length (unequal number of particles>");
            }
            N_p_test = matVar_dirs->dims[1];
            directions.resize(N_p_test);
            probabilities.resize(N_p_test);
            positions = Eigen::MatrixXd::Zero(N_p_test, 4);
            phases = Eigen::MatrixXd::Zero(N_p_test, 3);
            initial_positions = Eigen::MatrixXd::Zero(N_p_test, 4);


            // Looping through all particles to retrieve directions, probability and states
            // Note they should have the same length
            for (int i = 0; i < N_p_test; i++) {
                dirs = Mat_VarGetStructFieldByName(matVar_dirs, "directions", i);
                probs = Mat_VarGetStructFieldByName(matVar_probs, "probability", i);

                // Reading vertices
                bool directions_read_error = Mat_VarReadDataAll(test, dirs);
                if (not(directions_read_error)) {
                    directions[i] = testing::matrix_conversion(dirs);
                }
                else {
                    throw std::runtime_error("testing::unable to read directions");
                    break;
                }

                // Reading vertices
                bool probabilities_read_error = Mat_VarReadDataAll(test, probs);
                if (not(probabilities_read_error)) {
                    probabilities[i] = testing::vector_conversion(probs);
                }
                else {
                    throw std::runtime_error("testing::unable to read probabilities");
                    break;
                }
            }

            read_probabilities = true;
            read_directions = true;

            positions = testing::matrix_conversion(matVar_pos);
            phases = testing::matrix_conversion(matVar_phis);
            initial_positions = testing::matrix_conversion(matVar_pos0);


            read_positions = true;
            read_phases = true;
        }
        catch (const std::exception& ex) {
            std::cout << ex.what() << std::endl;
        }
    }
    else {
        std::cout << "testing::non existent variable, please check your .mat file" << std::endl;
    }

    return read_directions and read_probabilities and read_phases and read_positions and read_probabilities;
}

Eigen::MatrixXd testing::matrix_conversion(matvar_t* field) {

    Eigen::MatrixXd output;
    int rows, cols;
    unsigned field_size;
    const double* Data = static_cast<const double*>(field->data);


    // Specify the matrix
    if (field->rank == 2) {
        rows = field->dims[0];
        cols = field->dims[1];
        // field_size = field->nbytes / field->data_size;

        output.resize(rows, cols);
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) {
                output(i, j) = Data[rows * j + i];
            }
        }

    }
    else {
        std::cout << "Matrix rank not equal to 2" << std::endl;
    }

    return output;
}

Eigen::VectorXd testing::vector_conversion(matvar_t* field) {

    Eigen::VectorXd output;
    int rows;
    unsigned field_size;
    const double* Data = static_cast<const double*>(field->data);

    rows = field->dims[0];
    field_size = field->nbytes / field->data_size;

    output.resize(rows);
    try {
        for (int i = 0; i < rows; i++) {
            // std::cout << i << std::endl;
            output(i) = Data[i];
        }
    }
    catch (const std::exception& ex) {
        std::cout << ex.what() << std::endl;
    }

    return output;
}

bool testing::generate_report(std::vector<particle_state>& states) {
    bool flag = false;
    std::ofstream report;
    report.open("testing/report.txt");
    if (report.is_open()) {
        std::vector<std::string> bugs = testing::comparison(states);
        for (int i = 0; i < bugs.size(); i++) {
            report << bugs[i] << std::endl;
        }

        flag = true;
    }
    else {
        std::cout << "Unable to generate report" << std::endl;
        flag = false;
    }
    report.close();

    return flag;
}

std::vector<std::string> testing::comparison(std::vector<particle_state>& states) {

    std::vector<std::string> results(7);
    std::vector<int> position_error;
    std::vector<int> phase_error;
    std::vector<int> myocytes_error;

    // Position 
    int counter_myocytes = 0;
    int counter_pos = 0;
    int counter_phase = 0;
    int counter_flag = 0;

    for (int i = 0; i < N_p_test; i++) {

        // Position imcompatible
        Eigen::Vector3d position_i = positions(i, { 0, 1, 2 });
        if (((states[i].position - position_i).array().abs() > tolerance).any()) {
            counter_pos++;
            position_error.push_back(i);
        }

        // Phase imcompatible
        Eigen::Vector3d phase_i = phases(i, { 0, 1, 2 });
        if (((states[i].phase - phase_i).array().abs() > tolerance).any()) {
            counter_phase++;
            phase_error.push_back(i);
        }

        // myocytes imcompatible
        if (states[i].myoindex == -1) {
            if (not(std::isnan(positions(i, 3)))) {
                counter_myocytes++;
                myocytes_error.push_back(i);
            }
        }
        else if (states[i].myoindex + 1 != positions(i, 3)) {
            counter_myocytes++;
            myocytes_error.push_back(i);
        }

        if (states[i].flag != 0) {
            counter_flag++;
        }
    }

    results[0] = "Random Walk for cardiac Diffusion Tensor Imaging";
    results[1] = "Number of particles simulated : " + std::to_string(states.size());


    results[2] = "Errors computing particle positions : " + std::to_string(counter_pos);
    results[3] = "Errors computing particle phase : " + std::to_string(counter_phase);
    results[4] = "Errors computing particle myocytes : " + std::to_string(counter_myocytes);
    results[5] = "Total number of flags : " + std::to_string(counter_flag);
    results[6] = " ";
    if (counter_flag + counter_myocytes + counter_phase + counter_pos == 0) {
        results[6] = "Run Successful : no bugs found";
    }
    else {
        results[6] = "Run Unsuccessful : underlying bugs exist";
    }

    return results;
}

int testing::get_N_p() {
    return N_p_test;
}