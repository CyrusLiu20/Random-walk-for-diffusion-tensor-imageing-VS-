#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <matio.h>
#include <omp.h>

// header file
#include "simulation.h"

simulation::simulation(walkers obj_input) {
    obj = obj_input;
}

// Initializes position, phase and flag and generates the position uniformly
bool simulation::seedParticlesInBox(Eigen::MatrixXd boundingboxes_input, int particlesPerBox_input = -1) {
    bool flag = true;

    // Initialize Mersenne Twister pseudo-random number generator
    // std::mt19937 gen(rd());
    std::mt19937 gen(obj.get_rng_seed());

    // Obtain the dimension
    cols = (int)boundingboxes_input.cols();
    // Obtain the number of bounding boxes
    N_b = (int)boundingboxes_input.rows();

    // initialization of matrices and vectors
    if (not(simulation::initialize(N_b))) {
        flag = false;
    }

    // Check the dimensions of the bounding boxes
    switch (cols) {
    case 2:
        boundingboxes(Eigen::indexing::all, Eigen::seqN(0, 2)) = boundingboxes_input;
        dimension = 1;
        break;
    case 4:
        boundingboxes(Eigen::indexing::all, Eigen::seqN(0, 4)) = boundingboxes_input;
        dimension = 2;
        break;
    case 6:
        boundingboxes = boundingboxes_input;
        dimension = 3;
        break;
    default:
        printf("simulation::seedParticlesInBox::inconsistent, BBs must be [M,N] with M=2,4,6 and N=N_p\n");
        flag = false;
    }

    // Check bounding box inconsistentcy
    for (int i = 0; i < N_b; i++) {
        bool inconsistent = (boundingboxes(i, 0) > boundingboxes(i, 1) or boundingboxes(i, 2) > boundingboxes(i, 3) or boundingboxes(i, 4) > boundingboxes(i, 5));
        if (inconsistent) {
            printf("ParticleWalker:seedParticlesInBox:inconsistent, BBs must be [[xmin;xmax;ymin;ymax;zmin;zmax]]\n");
            flag = false;
        }
    }

    if (particlesPerBox_input == -1) {
        // Calculate the side lengths of bounding boxes
        Eigen::MatrixXd max, min;
        max = boundingboxes(Eigen::indexing::all, { 1,3,5 });
        min = boundingboxes(Eigen::indexing::all, { 0,2,4 });
        sidelengths_raw = max - min;

        // Check if and dimension has zero length in specified dimensions (e.g. only x and y)
        hasZeroDim = (sidelengths_raw(Eigen::indexing::all, Eigen::seqN(0, dimension)).array() == 0.0).any();
        if (hasZeroDim) {
            printf("simulation::seedParticlesInBox::inconsistent, each dimension must be either zero or non-zero across all boxes");
            flag = false;
        }

        // Removing zero dimensions (?????)
        sidelengths = sidelengths_raw(Eigen::indexing::all, Eigen::seqN(0, dimension));

        // Calculating the box volume of all the boxes
        boxVolumes = sidelengths.rowwise().prod();
        total_boxVolume = boxVolumes.sum();

        // Calculating the particles per box
        particlesPerBox_theo = boxVolumes / total_boxVolume * obj.get_N_p();
        particlesPerBox_prel = Eigen::floor(particlesPerBox_theo.array());
        int missing_particles = obj.get_N_p() - particlesPerBox_prel.sum();

        // Refilling missing particles
        particlesPerBox = simulation::refill(particlesPerBox_prel, missing_particles);
    }
    else if (particlesPerBox_input > 0) {
        particlesPerBox.resize(N_b);
        particlesPerBox = Eigen::VectorXd::Ones(N_b, 1) * particlesPerBox;
    }

    // Check if there are still any missing particles
    // printf("Testing : difference in particle number : %d\n", particlesPerBox.sum() - obj.get_N_p());
    if (particlesPerBox.sum() != obj.get_N_p()) {
        printf("simulation::seedParticlesInBox::missing_particles, please check your preliminary particles number");
        flag = false;
    }

    // Seeding particles in bounding boxes
    if (not(simulation::seeding())) {
        printf("simulation::seedParticlesInBox::seeding");
    }

    return flag;
}

bool simulation::initialize(int N_b) {
    bool flag = true;
    // Initializing the bounding boxes, side lengths of bounding boxes, and box volumes
    boundingboxes = Eigen::MatrixXd::Zero(N_b, 6);
    sidelengths = Eigen::MatrixXd::Zero(N_b, 3);
    boxVolumes = Eigen::VectorXd::Zero(N_b, 1);
    // Initializing the preliminary and theoretical particles per box
    particlesPerBox_theo = Eigen::VectorXd::Zero(N_b, 1);
    particlesPerBox_prel = Eigen::VectorXd::Zero(N_b, 1);

    return flag;
}

// A custom function for refilling missing particles
Eigen::VectorXd simulation::refill(Eigen::VectorXd particlesPerBox, int missing_particles) {

    // Initialize Mersenne Twister pseudo-random number generator
    // std::mt19937 gen(rd());
    std::mt19937 gen(obj.get_rng_seed());


    int rand_index;
    int number = particlesPerBox.rows();

    std::uniform_int_distribution<> rand_index_gen(0, number - 1);
    for (int i = 0; i < missing_particles; i++) {
        rand_index = rand_index_gen(gen);
        particlesPerBox(rand_index) = particlesPerBox(rand_index) + 1;
    }

    return particlesPerBox;
}

bool simulation::seeding() {

    // Flag for seeding
    bool flag = true;

    // Initialize Mersenne Twister pseudo-random number generator
    std::mt19937 gen(obj.get_rng_seed());
    // std::mt19937 gen(rd());

    // Initialization
    Eigen::VectorXd boundingBox = Eigen::VectorXd::Zero(N_b, 1).transpose();
    int ip_SetLast = 0; // Index of the particles
    int index_first, index_last, nP_iB;

    for (int iB = 0; iB < N_b; iB++) {
        boundingBox = boundingboxes(iB, Eigen::indexing::all);
        nP_iB = particlesPerBox(iB);
        index_first = ip_SetLast;
        index_last = ip_SetLast + nP_iB;
        ip_SetLast += nP_iB;

        std::uniform_real_distribution<> x_interval(boundingBox(0), boundingBox(1));
        std::uniform_real_distribution<> y_interval(boundingBox(2), boundingBox(3));
        std::uniform_real_distribution<> z_interval(boundingBox(4), boundingBox(5));

        for (int j = index_first; j < index_last; j++) {

            obj.position(j, 0) = x_interval(gen);
            obj.position(j, 1) = y_interval(gen);
            obj.position(j, 2) = z_interval(gen);
        }
    }

    return flag;
}

void simulation::performScan(sequence sequence_input, substrate substrate_input) {

    // perform the scan using the provided sequence in the provided substrate
    int number_of_particles = obj.get_N_p();
    particle_states.resize(number_of_particles);
    // position_histories.resize(number_of_particles);
    // phase_histories.resize(number_of_particles);
    // simulation::create_initial_states(number_of_particles);
    // loop over all particles (all independent, thus parallel)
    Eigen::Vector3d position_i, phase_i;

    omp_set_num_threads(number_of_particles);
#pragma omp parallel private(phase_i, position_i, unit_test)
#pragma omp for
    for (int i_particle = 0; i_particle < number_of_particles; i_particle++) {
#pragma omp critical
        position_i = obj.position(i_particle, Eigen::indexing::all);
        phase_i = obj.phase(i_particle, Eigen::indexing::all);
        particle_states[i_particle] = simulation::onewalker(sequence_input, substrate_input, i_particle, position_i, phase_i, obj.flag(i_particle));
        // simulation::save_history();
    }

#pragma omp barrier

    std::cout << "Simulation ended" << std::endl;
    simulation::save_history();
}

// The simualtion of one particle, position, phase, and the flag is the state of that one particle
particle_state simulation::onewalker(sequence& sequence_input, substrate& substrate_input, int i_particle, Eigen::Vector3d position_input, Eigen::Vector3d phase_input, bool flag) {
    particle_state one_particle_state;

    // Initialize state history
    one_particle_state.position_history = Eigen::MatrixXd::Zero(sequence_input.N, 4);
    one_particle_state.phase_history = Eigen::MatrixXd::Zero(sequence_input.N, 3);

    int myoindex = -1;

    // Different psuedo number generator : please be careful
    std::mt19937 gen(obj.get_rng_seed());

    Eigen::Vector3d phase = phase_input;
    Eigen::Vector3d position = Eigen::Vector3d::Zero(3);
    if (not(debugging)) {
        position = position_input;
    }
    else {
        position = unit_test.initial_positions(i_particle, Eigen::indexing::all);
    }
    // Eigen::Vector3d position = position_debug;


    // To Do : implement try catch error control
    try {
        myoindex = substrate_input.findMyocyte(position, obj.get_rng_seed(), "global");

    }
    catch (const std::exception& ex) {
        std::string error_message = "Particle : " + std::to_string(i_particle);
        std::cout << error_message << " : " << ex.what() << std::endl;
        obj.flag[i_particle] = 1;
    }

    Eigen::Vector3d position_new;
    // To do: maybe a catch to break the for loop?
    for (int i = 0; i < sequence_input.N; i++) {

        // get sequence step values
        double dt = sequence_input.get_dt(i);
        Eigen::VectorXd gG = sequence_input.get_gG(i); // It may be 1D or 3D

        // phase
        if (gG.rows() == 1) {
            phase = phase + (position * gG(0) * dt);
            // one_particle_state.phase_history(i, Eigen::indexing::all) = phase.transpose();
            one_particle_state.phase_history(i, Eigen::indexing::all) << phase(0), phase(1), phase(2);
            // obj.phase(i_particle, Eigen::indexing::all) = (position_debug.array()*gG(0)*dt);
        }
        else if (gG.rows() == 3) {
            phase = phase.array() + (position.array() * gG.array() * dt);
            // obj.phase(i_particle, Eigen::indexing::all) = (position_debug.array()*gG.array()*dt);
        }
        else {
            std::cout << "Wrong matrix size" << std::endl;
        }


        // try step until success
        int counter = 0;
        bool step_success = false;
        while (not(step_success)) {
            try {

                counter++;
                if (counter > 50) {
                    std::cout << "One walker : Counter over 50" << std::endl;
                    obj.flag(i_particle) = 2;
                    break;
                }

                Eigen::VectorXd position_raw = simulation::one_dt(position, dt, substrate_input, myoindex, i_particle, i);
                position = position_raw({ 0,1,2 });
                myoindex = position_raw(3);

                // one_particle_state.position_history(i, Eigen::indexing::all) = position_raw.transpose();
                one_particle_state.position_history(i, Eigen::indexing::all) << position_raw(0), position_raw(1), position_raw(2), position_raw(3);
                step_success = true;

            }
            catch (const std::exception& ex) {
                std::string error_message = ex.what();
                if (error_message.find("Polyhedron:intersection:uncertain") != std::string::npos or error_message.find("ParticleWalker:one_dt:unfinished") != std::string::npos) {
                    // if (error_message == "Polyhedron:intersection:uncertain', 'Too close to edge/vertex/face\n" or error_message == "ParticleWalker:one_dt:unfinished, Step has not finished after 50 substeps"){
                    position = one_particle_state.position_history(i - 1, { 0, 1, 2 });
                    myoindex = one_particle_state.position_history(i - 1, 3);
                    step_success = false;
                }
                else {
                    // To do ： Sorting error exceptions
                    // std::cout << ex.what() << std::endl;
                    obj.flag(i_particle) = 1;
                    step_success = true; // Don't deal with this particle
                    break;
                }
            }
        }


    }

    // Storing the final state
    one_particle_state.phase = phase;
    one_particle_state.position = position;
    one_particle_state.flag = obj.flag(i_particle);
    one_particle_state.myoindex = myoindex;

    return one_particle_state;
}

// One time step
Eigen::VectorXd simulation::one_dt(Eigen::Vector3d& position, double& dt, substrate& substrate_input, int& myoindex, int& i_particle, int& timestep) {
    Eigen::VectorXd output = Eigen::VectorXd::Zero(4);;
    double maxStepLength = 5;
    std::mt19937 gen(rd());

    Eigen::Vector3d dxdydz_normaldistrib = Eigen::Vector3d::Zero(3);

    // Generates the direction of the particle
    if (not(debugging)) {
        dxdydz_normaldistrib = simulation::getLimitedSteps(substrate_input.dim, maxStepLength);
    }
    else {
        dxdydz_normaldistrib = unit_test.directions[i_particle](timestep, Eigen::indexing::all);
    }
    // Eigen::Vector3d dxdydz_debug{{-1, -1, 0}};
    double D_old;
    if (myoindex < 0) {
        D_old = substrate_input.D_e;
    }
    else {
        D_old = substrate_input.D_i;
    }
    double D_new = D_old; // D_new holds the new diffusivity for every sub-step
    Eigen::Vector3d dxdydz = dxdydz_normaldistrib * std::sqrt(2 * dt * D_old);
    // Eigen::Vector3d dxdydz = dxdydz_debug * std::sqrt(2*dt*D_old); // Please change

    // until no more step left
    double ZERO = 1e-12; // 1e-12[m] = 1e-6[um] (note: eps(1) == 2e-16)
    int counter = 0;

    double probability_of_transit, ds, term; // Distrance normal to boundary
    double D_low, l_low, p_fieremans, p_maruyama;
    transform_parameter transform_info;
    Eigen::Vector3d position_LOCAL = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d position_future = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d dxdydz_toIntersection = Eigen::Vector3d::Zero(3);
    try {
        while (dxdydz.norm() > ZERO) {
            D_old = D_new;
            counter++;
            // if (counter == 2){
            //     std::cout << "One dt:: counter = 2" << std::endl;
            // }
            if (counter > 50) {
                throw std::logic_error("ParticleWalker:one_dt:unfinished, Step has not finished after 50 substeps\n");
                break;
            }

            // substrate checks
            transform_info = substrate_input.substrate_transform.global2local(position);
            position_LOCAL = transform_info.position_local;
            dxdydz = substrate_input.substrate_transform.rotate_y(dxdydz, transform_info.angle_reverse);

            intersection_info intersectInfo = substrate_input.intersectMyocytes(position_LOCAL, dxdydz, "local");

            // intersectInfo now contains info about first encountered intersection
            double stepEps = 1e-8;
            if (intersectInfo.empty) {
                position_future = position_LOCAL + dxdydz;

                if (not(substrate_input.block_bb.boundingbox.containsPoint(position_future))) {
                    // may throw an error, just take it and flag particle
                    intersection_info  intersectInfoBB = substrate_input.block_bb.intersection(position_LOCAL, dxdydz);
                    if (intersectInfoBB.empty) {
                        throw std::logic_error("ParticleWalker:one_dt:bb_inconsistent, Empty intersection when there should be one");
                        // std::cout << "ParticleWalker:one_dt:bb_inconsistent, Empty intersection when there should be one" << std::endl;
                    }

                    dxdydz_toIntersection = dxdydz * intersectInfoBB.t;
                    dxdydz = dxdydz * (1 - intersectInfoBB.t);

                    // ENABLE if the particles cannot leave the bounding box
                    if (substrate_input.boundary == "reflect") {
                        dxdydz = simulation::reflect(dxdydz, intersectInfoBB.vertices);
                    }

                    position_LOCAL = position_LOCAL + dxdydz_toIntersection;
                    position_LOCAL = position_LOCAL + dxdydz * stepEps;
                    dxdydz = dxdydz * (1 - stepEps);
                }
                else {
                    position_LOCAL = position_LOCAL + dxdydz; // no need to worry about 'position_LOCAL'
                    dxdydz = Eigen::Vector3d::Zero(3);
                }
            }
            else { // If an intersection was encountered
                Eigen::Vector3d whole_step = dxdydz;
                dxdydz_toIntersection = dxdydz * intersectInfo.t;
                dxdydz = dxdydz * (1 - intersectInfo.t);

                if (substrate_input.transit_model == "constant") {
                    probability_of_transit = substrate_input.kappa;
                }
                else if (substrate_input.transit_model == "Fieremans2010") {
                    // Fieremans et al, 2010, NMR Biomed.
                    // Monte Carlo study of a two-compartment exchange model of diffusion
                    // DOI:10.1002/nbm.1577

                    ds = simulation::computeNormalDistance(intersectInfo.vertices, dxdydz_toIntersection);
                    term = (2 * ds * substrate_input.kappa) / D_old;
                    probability_of_transit = term / (1 + term);
                }
                else if (substrate_input.transit_model == "HybridModel") {
                    ds = simulation::computeNormalDistance(intersectInfo.vertices, dxdydz_toIntersection);
                    if (D_old == substrate_input.D_e) {
                        D_low = substrate_input.D_i;
                        l_low = ds * std::sqrt(D_low / D_old);
                        term = 2 * l_low * substrate_input.kappa / D_low;
                        p_fieremans = term / (term + 1);
                        p_maruyama = std::sqrt(substrate_input.D_i / D_old) > 1 ? 1 : std::sqrt(substrate_input.D_i / D_old); // equivalent to matlab min(1,sqrt(substrate.D_i/D_old))
                        probability_of_transit = p_fieremans * p_maruyama;
                    }
                    else {
                        term = 2 * ds * substrate_input.kappa / substrate_input.D_i;
                        p_fieremans = term / (term + 1);
                        // p_maruyama = min(1,sqrt(D_new/D_old)); (in matlab)
                        probability_of_transit = p_fieremans;
                    }
                }
                else {
                    std::cout << "Error:NotImplemented, Transit model not supported" << std::endl;
                }

                std::uniform_real_distribution<> transit(0, 1);
                double U;
                if (not(debugging)) {
                    U = transit(gen);
                }
                else {
                    // U = unit_test.probabilities[i_particle](timestep);
                    U = 0.2;
                }

                if (U < probability_of_transit) {
                    // %TODO: ensure boundary cases (rand==0 or rand==1) are handled correctly
                    // % --> go through

                    // % update index %TODO: ensure this is done correctly!!!

                    if (myoindex < 0) {
                        myoindex = intersectInfo.myoindex;
                    }
                    else {
                        myoindex = -1;
                    }

                    // change D
                    if (myoindex < 0) {
                        D_new = substrate_input.D_e;
                    }
                    else {
                        D_new = substrate_input.D_i;
                    }

                    dxdydz = dxdydz * std::sqrt(D_new / D_old);

                }
                else {
                    // --> reflect
                    dxdydz = simulation::reflect(dxdydz, intersectInfo.vertices);
                }
                // step a little bit to get off face
                // - this requires faces to be at least a certain distance away from each other (geometry check!)
                position_LOCAL = position_LOCAL + dxdydz_toIntersection; // initial sub-step to intersection side
                position_LOCAL = position_LOCAL + dxdydz * stepEps; // little Eps extra of new step to move away from face
                dxdydz = dxdydz * (1 - stepEps); // remove eps
            }

            // transform the position back from the local to the global frame
            position = substrate_input.substrate_transform.local2global(position_LOCAL, transform_info.iX, transform_info.iY, transform_info.iZ);
            dxdydz = substrate_input.substrate_transform.rotate_y(dxdydz, transform_info.angle); // rotate reverse
            output << position, myoindex;
            // std::cout << "Position : " << position.transpose() << std::endl;
        }
    }
    catch (const std::exception& ex) {
        throw;
    }

    return output;
}

Eigen::Vector3d simulation::getLimitedSteps(std::string& dim, double& maxStepLength) {

    double maxStep_squared = maxStepLength * maxStepLength;
    Eigen::Vector3d dxdydz = Eigen::Vector3d::Zero(3);
    bool needUpdate = true;
    int tries = 0;

    std::mt19937 gen(rd());
    std::uniform_real_distribution<> choiceVector(-1, 1);
    std::normal_distribution<double> normal(0, 1);

    while (needUpdate) {
        // check
        tries++;
        if (tries > 10) {
            std::cout << "ParticleWalker:getLimitedStep:tries', 'Cannot draw step with limited size" << std::endl;
            break;
        }

        // Dimensions must be in suquencial order x, y, z
        int dimension = dim.size();

        // Constant time step creates either 1 or -1 value 
        if (obj.get_steptype() == "constant") {
            for (int i = 0; i < dimension; i++) {
                dxdydz(i) = choiceVector(gen) > 0 ? 1 : -1;
            }
        }
        // A standard normal distribution
        else if (obj.get_steptype() == "normal") {
            for (int i = 0; i < dimension; i++) {
                dxdydz(i) = normal(gen);
            }
        }

        // needUpdate = dxdydz.norm() > maxStepLength;
        needUpdate = (dxdydz.array() * dxdydz.array()).sum() > maxStep_squared;
    }

    return dxdydz;
}

// When a particle hits a boundary and does not pass through
Eigen::Vector3d simulation::reflect(Eigen::Vector3d& oldstep, Eigen::MatrixXd& faceVertices) {
    Eigen::Vector3d V1 = faceVertices(0, Eigen::indexing::all);
    Eigen::Vector3d V2 = faceVertices(1, Eigen::indexing::all);
    Eigen::Vector3d V3 = faceVertices(2, Eigen::indexing::all);

    Eigen::Vector3d edge01 = V2 - V1;
    Eigen::Vector3d edge02 = V3 - V1;

    Eigen::Vector3d normal = edge01.cross(edge02);
    normal = normal / normal.norm();

    double step_magn = oldstep.norm();
    Eigen::Vector3d step_norm = oldstep / step_magn;
    Eigen::Vector3d step_norm_reflected = step_norm - normal * 2 * (step_norm.dot(normal));

    Eigen::Vector3d new_step = step_norm_reflected * step_magn;

    return new_step;
}

// Normal distance to the boundary
double simulation::computeNormalDistance(Eigen::MatrixXd& faceVertices, Eigen::Vector3d& step) {
    Eigen::Vector3d V1 = faceVertices(0, Eigen::indexing::all);
    Eigen::Vector3d V2 = faceVertices(1, Eigen::indexing::all);
    Eigen::Vector3d V3 = faceVertices(2, Eigen::indexing::all);

    Eigen::Vector3d edge01 = V2 - V1;
    Eigen::Vector3d edge02 = V3 - V1;

    Eigen::Vector3d normal = edge01.cross(edge02);

    // compute the dot product using normalised vectors
    normal = normal / normal.norm();
    Eigen::Vector3d step_normalized = step / step.norm();
    double dotproduct = step_normalized.dot(normal);

    // % determine the distance
    // %{
    // theta = acos(dotproduct);
    // if theta >= pi/2, theta = pi - theta; end
    // t = cos(theta);
    // %}

    double t = dotproduct > 0 ? dotproduct : -dotproduct; // this is the equivalent of the commented code above
    double distance = step.norm() * t;

    return distance;
}

void simulation::load_test_module(testing test) {
    unit_test = test;
    debugging = true;
    std::cout << "Please be alert : debugging mode activated" << std::endl;
}

// Save the particle trajectory of all particles
void simulation::save_history() {
    mat_t* intermediate = NULL;
    intermediate = Mat_CreateVer("intermediate_steps.mat", NULL, MAT_FT_MAT5);
    const char* structname = "histories";
    const char* fieldnames[2] = { "position", "phase" };
    size_t structdim[2] = { 1, (long long unsigned int)obj.get_N_p() }; // create 1 x n struct
    //main struct: Data with 2 fields
    matvar_t* matstruct = Mat_VarCreateStruct(structname, 2, structdim, fieldnames, 2);
    for (int i_particle = 0; i_particle < obj.get_N_p(); i_particle++) {
        // int i_particle = 9;
        size_t pos_dim[2] = { (long long unsigned int)particle_states[i_particle].position_history.rows(), (long long unsigned int)particle_states[i_particle].position_history.cols() }; //string dimension
        double* pos = particle_states[i_particle].position_history.data();
        matvar_t* pos_variable = Mat_VarCreate(fieldnames[0], MAT_C_DOUBLE, MAT_T_DOUBLE, 2, pos_dim, pos, 0);
        Mat_VarSetStructFieldByName(matstruct, fieldnames[0], i_particle, pos_variable); //insert Data(p).name (1 <= p <= n)

        size_t phase_dim[2] = { (long long unsigned int)particle_states[i_particle].phase_history.rows(), (long long unsigned int)particle_states[i_particle].phase_history.cols() }; //string dimension
        double* phase = particle_states[i_particle].phase_history.data();
        matvar_t* phase_variable = Mat_VarCreate(fieldnames[1], MAT_C_DOUBLE, MAT_T_DOUBLE, 2, phase_dim, phase, 0);
        Mat_VarSetStructFieldByName(matstruct, fieldnames[1], i_particle, phase_variable); //insert Data(p).name (1 <= p <= n)

    }
    Mat_VarWrite(intermediate, matstruct, MAT_COMPRESSION_NONE);
    Mat_Close(intermediate);
    Mat_VarFree(matstruct);
}

void simulation::create_initial_states(long long unsigned int number_of_particles) {
    const char* filename = "intermediate_steps.mat";
    mat_t* intermediate = NULL; //matfp contains pointer to MAT file or NULL on failure
    intermediate = Mat_CreateVer(filename, NULL, MAT_FT_MAT5); //or MAT_FT_MAT4 / MAT_FT_MAT73

    //Create a 1 x n struct 'Data' with fields: name, unit, value
    const char* structname = "histories";
    const char* fieldnames[2] = { "position","phase" };
    size_t structdim[2] = { 1, number_of_particles }; // create 1 x n struct

    //main struct: Data with 2 fields
    matvar_t* matstruct = Mat_VarCreateStruct(structname, 2, structdim, fieldnames, 2);
    Mat_VarWrite(intermediate, matstruct, MAT_COMPRESSION_NONE);
    Mat_VarFree(matstruct);

    Mat_Close(intermediate);
}

// Data retrieval

Eigen::MatrixXd simulation::get_position() {
    return obj.position;
}

std::vector<particle_state> simulation::get_states() {
    return particle_states;
}
// Notes
// dimension must be sequential