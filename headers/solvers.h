//
// Created by Piotr on 27-Jun-20.
//

#ifndef NAREW_SOLVERS_H
#define NAREW_SOLVERS_H

#include <vector>
#include "equations.h"
#include <string>

class SolverRungeKutta6 {
    double dt;
    ReactionTerm reactionTerm;

public:
    SolverRungeKutta6(double dt, ReactionTerm reactionTerm);

    double step_rk6(double previous_value, int iteration_number);

};

class SolverCrankNicolson {
    double dt;
    double dx;
    double adv_coeff;
    double diff_coeff;
    std::string modeled_variable;


    static std::vector<double> use_thomas(std::vector<double> a,
                                          std::vector<double> b,
                                          std::vector<double> c,
                                          std::vector<double> d,
                                          int n);

public:
    SolverCrankNicolson(double dt, double dx, double advCoeff, double diffCoeff, std::string &modeled_variable);

    std::vector<double> step_cn(int nx, std::vector<double> y_prev,
                                double left_boundary_cond, double right_boundary_cond);
};

#endif //NAREW_SOLVERS_H
