#include <equations.h>
#include <iostream>
#include "headers/solvers.h"
#include "headers/searchers.h"


int main(){
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

    std::vector<std::vector<double>> diss_oxygen(nt, std::vector<double>(nx, 0.));
    diss_oxygen[0] = init_cond_adv_diff;

    for (int i = 1; i < nt; ++i) {
        diss_oxygen[i] = solverCrankNicolson.step_cn(nx, diss_oxygen[i - 1], left_bound);
        for (int j = 0; j < nx; ++j) {
            diss_oxygen[i][j] = solverRk6.step_rk6(diss_oxygen[i][j], i);
        }
    }
}