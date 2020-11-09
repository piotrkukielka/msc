//
// Created by Piotr on 02-Jul-20.
//

#ifndef MSC_0505_SEARCHERS_H
#define MSC_0505_SEARCHERS_H

#include "equations.h"
#include <string>

class Bounds {
public:
    double lower;
    double upper;

    int get_num_of_iters();

    double get_value(int i);

    Bounds(double lower, double upper, double interval);

    double interval;
};

class GridSearch {
    Bounds spatialGrid;
    Bounds timeGrid;

public:
    GridSearch(const Bounds &spatialGrid, const Bounds &timeGrid);

    void search(Bounds advCoeffs, Bounds diffCoeffs, std::string &modeled_value);

};

class Simulation {
    double dx, dt, starting_point_spatial, starting_point_time;
    int nx, nt;

    void save(std::vector<std::vector<double>> data, const std::string &path);

    std::vector<std::vector<double>>
    run(double advCoeff, double diffCoeff, bool is_transport, std::string &modeled_variable);

public:
    Simulation(Bounds spatialGrid, Bounds timeGrid);

    void run_and_save(double advCoeff, double diffCoeff, bool is_transport, bool is_only_measurepoints,
                      std::string &modeled_variable);

    std::string create_path(double adv_coeff, double diff_coeff);

    void save_only_measurepoints(std::vector<std::vector<double>> data, const std::string &path);
};


#endif //MSC_0505_SEARCHERS_H
