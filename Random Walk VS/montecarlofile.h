#ifndef MONTECARLOFILE_H
#define MONTECARLOFILE_H

#include <Eigen/Dense>
#include <string>
class montecarlofile {

public:

    // All valid constructors
    // Initiallisation
    montecarlofile() = default;
    montecarlofile(const montecarlofile& pSrc) = default;	// copying class
    montecarlofile(montecarlofile&& pSrc) = default;	// moving class
    ~montecarlofile() = default;

    // Initial seed for the random number generator
    unsigned int rngseed = 345676;
    // Number of walkers
    int N_p = 5;

    // stepping type used in the simulation
    std::string stepType = "normal";

    // Whether text or vector is used
    std::string seedbox_type = "text"; // either "text" or "numeric"
    //  Seeding location for walkers
    std::string seedbox_text = "voxel+box"; // 'voxel', 'voxel+buffer'
    // or a bounding box (which may be singular for some dimensions) [xmin,xmax,ymin,ymax,zmin,zmax]
    Eigen::VectorXd seedbox_vect{(0, 1, 0, 0, 0, 3)}; // = Eigen::VectorXd(6);
};

#endif
