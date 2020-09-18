#include <equations.h>
#include <iostream>
#include "headers/solvers.h"
#include "headers/searchers.h"


int main(){

//    Bounds spatialGrid{0., 92.1, 0.01};  // km
//    Bounds timeGrid{0., 79.1, 0.1};  // h

    Bounds spatialGrid{0., 100.01, 0.01};  // km
    Bounds timeGrid{0., 124.1, 0.1};  // h


//    Bounds spatialGrid{0., 3.14, 0.024};  // km
//    Bounds timeGrid{0., 2., 0.01};  // h





//    Bounds spatialGrid{37.6, 56.6, 0.01};  // km, for 0-5 put 0, 6 if the interval is 1
//    Bounds timeGrid{29.6, 74.6, 0.1};  // h, for 0-5 put 0, 6 if the interval is 1

//    Bounds advCoeffs{0.8, 1.11, 0.01}; // use minus later or sth
//    Bounds diffCoeffs{0.01, 0.11, 0.01};
    std::vector<double> rs{1,1,1,1,1,1}; // NEVER USED!

//    GridSearch gridSearch(spatialGrid, timeGrid);
//    gridSearch.search(advCoeffs, diffCoeffs, rs);

//    // TODO: grid 0:5 gives 0:4 etc
    Simulation simulation(spatialGrid, timeGrid);
    double adv_c = -1.;
    double diff_c = 0.03;
//    simulation.run_and_save(adv_c, diff_c, std::vector<double>{1, 1, 1, 1, 1, 1});
    simulation.run_and_save_measurepoints(adv_c, diff_c, std::vector<double>{1, 1, 1, 1, 1, 1});
}