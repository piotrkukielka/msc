//
// Created by Piotr on 02-Jul-20.
//

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <valarray>
#include "searchers.h"
#include "equations.h"
#include "solvers.h"


void GridSearch::search(Bounds advCoeffs, Bounds diffCoeffs, std::string &modeled_value) {
    Simulation simulation(this->spatialGrid, this->timeGrid);
    for (int i = 0; i < advCoeffs.get_num_of_iters(); ++i) {
        for (int j = 0; j < diffCoeffs.get_num_of_iters(); ++j) {
            std::cout << i * diffCoeffs.get_num_of_iters() + j + 1 << " out of "
                      << advCoeffs.get_num_of_iters() * diffCoeffs.get_num_of_iters() << std::endl;
            simulation.run_and_save(advCoeffs.get_value(i), diffCoeffs.get_value(j), true, true, modeled_value);
        }
    }
}

GridSearch::GridSearch(const Bounds &spatialGrid, const Bounds &timeGrid) : spatialGrid(spatialGrid),
                                                                            timeGrid(timeGrid) {}

Bounds::Bounds(double lower, double upper, double interval) : lower(lower), upper(upper),
                                                              interval(interval) {}

int Bounds::get_num_of_iters() {
    int num_of_iters = int(round((this->upper - this->lower) / this->interval));
    return num_of_iters;
}

double Bounds::get_value(int i) {
    return this->lower + i * this->interval;
}

void Simulation::run_and_save(double advCoeff, double diffCoeff, bool is_transport, bool is_only_measurepoints,
                              std::string &modeled_variable) {
    std::vector<std::vector<double>> result = this->run(advCoeff, diffCoeff, is_transport, modeled_variable);
    std::string path = this->create_path(advCoeff, diffCoeff);
    if (is_only_measurepoints) {
        this->save_only_measurepoints(result, path);
    } else {
        this->save(result, path);
    }
}

std::string Simulation::create_path(double adv_coeff, double diff_coeff) {
    // TODO: add dx and dt to name
    std::string path = "../results/";
    std::string filename = std::to_string(adv_coeff) + "_" + std::to_string(diff_coeff);
    return path + filename;
}

void Simulation::save(std::vector<std::vector<double>> data, const std::string &path) {
    std::ofstream outFile(path);
    for (int j = 0; j < this->nx; ++j) {
        if (j % 10 == 0) { // saving every 10th time and spatial step
            for (int i = 0; i < this->nt; ++i) {
                if (i % 10 == 0) { // saving every 10th time and spatial step
                    outFile << i * this->dt << " " << j * this->dx << " " << data[i][j] << std::endl;
                }
            }
            outFile << std::endl;
        }
    }
}

void Simulation::save_only_measurepoints(std::vector<std::vector<double>> data, const std::string &path) {
    std::ofstream outFile(path + "measurepoints");
    std::valarray<double> spatial_measurepoints{
            0,
            0,
            9.2,
            0,
            9.2,
            9.2,
            0,
            15.7,
            37.6,
            37.6,
            37.6,
            43.9,
            43.9,
            50.9,
            56.5,
            56.5,
            56.5,
            65.9,
            56.5,
            73.4,
            56.5,
            56.5,
            56.5,
            90.2,
            90.2
    };
    spatial_measurepoints = (spatial_measurepoints - this->starting_point_spatial) / dx;
    std::valarray<double> time_measurepoints{
            0,
            1.91666666662786,
            3.66666666674428,
            5.66666666662786,
            7.33333333331393,
            8.25,
            9,
            13.5,
            29.5833333333721,
            30.8333333334303,
            31.8333333333721,
            35.9166666667443,
            37.1666666666279,
            38.25,
            41.8333333333139,
            46.4166666667443,
            51.9166666666861,
            53.4166666666861,
            55.9166666666279,
            60.5000000000582,
            62.3333333334304,
            66.5833333333139,
            74.4999999999418,
            76.9166666666279,
            79.0833333333721,
    };
    time_measurepoints = time_measurepoints + 24.; // TODO: jest prestart
    time_measurepoints = (time_measurepoints - this->starting_point_time) / dt;
    for (int i = 0; i < spatial_measurepoints.size(); ++i) {
        outFile << data[int(time_measurepoints[i])][int(spatial_measurepoints[i])] << std::endl;
    }
    outFile << std::endl;
}

std::vector<std::vector<double>>
Simulation::run(double advCoeff, double diffCoeff, bool is_transport, std::string &modeled_variable) {
    ReactionTerm reactionTerm{modeled_variable};
    SolverRungeKutta6 solverRk6{dt, reactionTerm};

    SolverCrankNicolson solverCrankNicolson{dt, dx, advCoeff, diffCoeff, modeled_variable};

    LeftBoundaryCondition leftBoundaryCondition{modeled_variable};
    std::vector<double> left_bound(nt); // g per m3
    for (int i = 0; i < nt; ++i) {
        left_bound[i] = leftBoundaryCondition.evaluate(i * dt);
    }

    RightBoundaryCondition rightBoundaryCondition{modeled_variable};
    std::vector<double> right_bound(nt); // g per m3
    for (int i = 0; i < nt; ++i) {
        right_bound[i] = rightBoundaryCondition.evaluate(i * dt);
    }

    InitialCondition initialCondition{modeled_variable};
    std::vector<double> init_cond_adv_diff(nx);
    for (int i = 0; i < nx; ++i) {
        init_cond_adv_diff[i] = initialCondition.evaluate(i * dx);
    }

    std::vector<std::vector<double>> result(nt, std::vector<double>(nx, 0.));
    result[0] = init_cond_adv_diff;

    for (int i = 1; i < nt; ++i) {
        result[i] = solverCrankNicolson.step_cn(nx, result[i - 1], left_bound[i], right_bound[i]);
        if (is_transport) {
            for (int j = 0; j < nx; ++j) {
                result[i][j] = solverRk6.step_rk6(result[i][j], i);
            }
        } else {
            result[i] = result[i - 1];
        }
    }
    return result;
}

Simulation::Simulation(Bounds spatialGrid, Bounds timeGrid) {
    this->dx = spatialGrid.interval;
    this->dt = timeGrid.interval;
    this->nx = spatialGrid.get_num_of_iters();
    this->nt = timeGrid.get_num_of_iters();
    this->starting_point_spatial = spatialGrid.lower;
    this->starting_point_time = timeGrid.lower;
}

