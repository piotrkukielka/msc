#include <equations.h>
#include <iostream>
#include "headers/solvers.h"
#include "headers/searchers.h"


int main(){
//    Bounds spatialGrid{0., 10, 0.1};  // km
//    Bounds timeGrid{0., 15, 1};  // h

    Bounds spatialGrid{0., 92.1, 0.1};  // km
    Bounds timeGrid{0., 79.1, 1};  // h
    Bounds advCoeffs{0.1, 1, 0.1};
    Bounds diffCoeffs{0.1, 1, 0.1};
    std::vector<double> rs{1,1,1,1,1,1}; // NEVER USED!

    GridSearch gridSearch(spatialGrid, timeGrid);
    gridSearch.search(advCoeffs, diffCoeffs, rs);

//    Simulation simulation(spatialGrid, timeGrid);
//    simulation.run_and_save(0, 0, std::vector<double>{1, 1, 1, 1, 1, 1});
}