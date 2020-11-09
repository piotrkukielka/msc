#include <equations.h>
#include <iostream>
#include <fstream>
#include "headers/solvers.h"
#include "headers/searchers.h"


int main(){

    Bounds spatialGrid{0., 100.01, 0.01};  // km
    Bounds timeGrid{0., 124.1, 0.1};  // h
    std::vector<double> rs{1,1,1,1,1,1}; // NEVER USED!

//    // TODO: grid 0:5 gives 0:4 etc
    Simulation simulation(spatialGrid, timeGrid);
    std::string modeled_value = "O2";
    double adv_c = 1.14;
    double diff_c = 0.08;
    simulation.run_and_save(adv_c, diff_c, true, true, modeled_value);


//    std::ofstream outFile("../results/photosynthesis.dat");
//    ReactionTerm reactionTerm;
//    for (int i = 0; i < 1240; ++i) {
//        double dt = 0.1;
//        double t = i*dt;
//        outFile << t << " " << reactionTerm.photosynthesis(t, 0.87) << std::endl;
//    }
}