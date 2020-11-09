//
// Created by Piotr on 27-Jun-20.
//

#include "../headers/solvers.h"

#include <utility>

SolverRungeKutta6::SolverRungeKutta6(double dt,
                                     ReactionTerm reactionTerm) : dt(dt),
                                                                  reactionTerm(std::move(reactionTerm)) {}

double SolverRungeKutta6::step_rk6(double previous_value, int iteration_number) {
    double t, next_value, k0, k1, k2, k3, k4, k5, k6, k7;
    t = dt * iteration_number;
    k0 = dt * reactionTerm.evaluate(t, previous_value);
    k1 = dt * reactionTerm.evaluate(t + 1. / 9. * dt, previous_value + 1. / 9. * k0);
    k2 = dt * reactionTerm.evaluate(t + 1. / 6. * dt, previous_value + 1. / 24. * (k0 + 3 * k1));
    k3 = dt * reactionTerm.evaluate(t + 1. / 3. * dt, previous_value + 1. / 6. * (k0 - 3 * k1 + 4 * k2));
    k4 = dt * reactionTerm.evaluate(t + 0.5 * dt, previous_value + 1. / 8. * (k0 + 3 * k3));
    k5 = dt * reactionTerm.evaluate(t + 2. / 3. * dt, previous_value + 1. / 9. * (17 * k0 - 63 * k1 + 51 * k2 + k4));
    k6 = dt *
         reactionTerm.evaluate(t + 5. / 6. * dt,
                               previous_value + 1. / 24. * (-22 * k0 + 33 * k1 + 30 * k2 - 58 * k3 + 34 * k4 + 3 * k5));
    k7 = dt * reactionTerm.evaluate(t + dt,
                                    previous_value +
                                    1. / 82. *
                                    (281 * k0 - 243 * k1 - 522 * k2 + 876 * k3 - 346 * k4 - 36 * k5 + 72 * k6));
    next_value = previous_value + 1. / 840. * (41 * (k0 + k7) + 216 * (k2 + k6) + 27 * (k3 + k5) + 272 * k4);
    return next_value;
}

std::vector<double>
SolverCrankNicolson::step_cn(int nx, std::vector<double> y_prev, double left_boundary_cond,
                             double right_boundary_cond) {
    std::vector<double> y(nx);
//    std::vector<double> a(nx, 0.25 * adv_coeff / dx - 0.5 * diff_coeff / dx / dx);
//    std::vector<double> b(nx, 1. / dt + diff_coeff / dx / dx);
//    std::vector<double> c(nx, -0.25 * adv_coeff / dx - 0.5 * diff_coeff / dx / dx);
    std::vector<double> a(nx, 0.25 * adv_coeff / dx - 0.5 * diff_coeff / dx / dx);
    std::vector<double> b(nx, 1. / dt + diff_coeff / dx / dx);
    std::vector<double> c(nx, -0.25 * adv_coeff / dx - 0.5 * diff_coeff / dx / dx);
    a[0] = 0;
    a[nx - 1] = 0;
    b[0] = 1;
    b[nx - 1] = 1;
    c[0] = 0;
    c[nx - 1] = 0;
    std::vector<double> d(nx);
    for (int j = 0; j < nx; ++j) {
        if (j == 0) {
            d[j] = left_boundary_cond;
        } else if (j == nx - 1) {
            // TODO: CO2
            if (this->modeled_variable == "DIC") {
                d[j] = y_prev[j - 1];
            }
            if (this->modeled_variable == "O2") {
                d[j] = right_boundary_cond;
            }
        } else {
            d[j] = (1. / dt - diff_coeff / dx / dx) * y_prev[j] +
                   (0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y_prev[j + 1] +
                   (-0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y_prev[j - 1];
        }
    }
    y = SolverCrankNicolson::use_thomas(a, b, c, d, nx);
    return y;
}

std::vector<double> SolverCrankNicolson::use_thomas(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                                                    std::vector<double> d, int n) {
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i] * d[i + 1];
    }
    return d;
}

SolverCrankNicolson::SolverCrankNicolson(double dt, double dx, double advCoeff, double diffCoeff,
                                         std::string &modeled_variable) : dt(dt), dx(dx),
                                                                          adv_coeff(-advCoeff),
                                                                          diff_coeff(diffCoeff),
                                                                          modeled_variable(modeled_variable) {}

