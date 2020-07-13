//
// Created by Piotr on 02-Jul-20.
//

#include <utility>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "searchers.h"
#include "equations.h"
#include "solvers.h"


void GridSearch::search(Bounds advCoeffs, Bounds diffCoeffs, std::vector<double> rs) {
    for (int i = 0; i < advCoeffs.get_num_of_iters(); ++i) {
        for (int j = 0; j < diffCoeffs.get_num_of_iters(); ++j) {
            std::cout<<advCoeffs.get_value(i)<<std::endl;
            Simulation simulation(advCoeffs.get_value(i), diffCoeffs.get_value(j), rs);
            simulation.run_and_save();
        }
    }
}

Bounds::Bounds(double lower, double upper, double interval) : lower(lower), upper(upper),
                                                              interval(interval) {}

int Bounds::get_num_of_iters() {
    int num_of_iters = int((this->upper - this->lower) / this->interval);
    return num_of_iters;
}

double Bounds::get_value(int i) {
    return this->lower + i * this->interval;
}

void Simulation::run_and_save() {
    std::vector<std::vector<double>> result = this->run();
    this->save(result);
}

void Simulation::save(std::vector<std::vector<double>> data) {
    std::string path = "../results/";
    std::string filename = std::to_string(this->diff_coeff) + "_" + std::to_string(this->adv_coeff);
    std::ofstream outFile(path + filename);
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[0].size(); ++j) {
            outFile << data[i][j] << " ";
        }
        outFile << std::endl;
    }
}

std::vector<std::vector<double>>
Simulation::run() {
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

    SolverCrankNicolson solverCrankNicolson{dt, dx, this->adv_coeff, this->diff_coeff};

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

Simulation::Simulation(double advCoeff, double diffCoeff, std::vector<double> rs) : adv_coeff(advCoeff),
                                                                                           diff_coeff(diffCoeff),
                                                                                           rs(std::move(rs)) {}
