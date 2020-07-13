#include <equations.h>
#include <iostream>
#include "headers/solvers.h"
#include "headers/searchers.h"


int main(){
    Bounds advCoeffs{0.1, 2, 0.1};
    Bounds diffCoeffs{0.1, 2, 0.1};
    std::vector<double> rs{1,1,1,1,1,1};

    GridSearch gridSearch;
    gridSearch.search(advCoeffs, diffCoeffs, rs);
}