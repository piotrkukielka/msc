//
// Created by Piotr on 02-Jul-20.
//

#include <vector>
#include <fstream>
#include <iostream>
#include "searchers.h"
#include "equations.h"
#include "solvers.h"


void GridSearch::search() {
    Simulation simulation;
    std::vector<double> rs{1,1,1,1,1,1};
    Bounds advCoeffs{0.1, 2, 0.1};
    Bounds diffCoeffs{0.1, 2, 0.1};
    for (int i = 0; i < advCoeffs.num_of_iters; ++i) {
        for (int j = 0; j < diffCoeffs.num_of_iters; ++j) {
            simulation.run_and_save(advCoeffs.values[i], diffCoeffs.values[j], rs)
        }
    }
}

Bounds::Bounds(double lower, double upper, double interval) : lower(lower), upper(upper),
                                                              interval(interval),
                                                              num_of_iters(this->find_num_of_iters()) {}

int Bounds::find_num_of_iters() {
    this->num_of_iters = int((this->upper - this->lower) / this->interval);
    return this->num_of_iters;
}

std::vector<double> Bounds::find_values() {
    for (int i = 0; i < this->num_of_iters; ++i) {
        this->values[i] = this->lower + i * this->interval;
    }
    return this->values;
}

void Simulation::run_and_save(double adv_coeff_bounds, double diff_coeff_bounds, std::vector<double> rs) {
    std::vector<std::vector<double>> result = this->run()
    this->save(result)
}

void Simulation::save(std::vector<std::vector<double>> data) {
    std::ofstream outFile("../results.txt");
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[0].size(); ++j) {
            outFile << data[i][j] << " ";
        }
        outFile << std::endl;
    }
}

std::vector<std::vector<double>>
Simulation::run(double adv_coeff_bounds, double diff_coeff_bounds, std::vector<double> rs) {
    const double dt = 1.;
    const double dx = 0.1;
    const double xend = 40.;  // km?
    const double xstart = 0.;
    const double tend = 4.*24.;  // h
    const double tstart = 0.;
    int nx = int((xend - xstart) / dx);
    int nt = int((tend - tstart) / dt);
    ReactionTerm reactionTerm{};
    SolverRungeKutta6 solverRk6{dt, reactionTerm};

    double adv_coeff = -0.1/24.;  // km/h
    double diff_coeff = 0.1/24.;  // km^2/h
    SolverCrankNicolson solverCrankNicolson{dt, dx, adv_coeff, diff_coeff};

    double left_bound = 0.; // kg per m3

    std::vector<double> init_cond_adv_diff(nx);
    for (int i = 0; i < nx; ++i) {
        init_cond_adv_diff[i] = InitialCondition::evaluate(i * dx);
    }

    std::vector<std::vector<double>> result(nt, std::vector<double>(nx, 0.));
    result[0] = init_cond_adv_diff;

    for (int i = 1; i < nt; ++i) {
        result[i] = solverCrankNicolson.step_cn(nx, result[i - 1], left_bound);
        for (int j = 0; j < nx; ++j) {
            result[i][j] = solverRk6.step_rk6(result[i][j], i);
        }
    }
    return result;
}
