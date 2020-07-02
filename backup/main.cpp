#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

std::vector<double> use_thomas(std::vector<double> a,
                               std::vector<double> b,
                               std::vector<double> c,
                               std::vector<double> d,
                               int n);

std::vector<double> step_cn(double dt,
                            double dx,
                            int nx,
                            std::vector<double> y_prev,
                            double left_boundary_cond,
                            double adv_coeff,
                            double diff_coeff);

double init_cond_fun(double x);

double step_rk6(double step,
                double previous_value,
                int iteration_number,
                double(*fun)(double, double));

double reaction_term(double t, double c);

int main() {
    double dt = 1.;
    double dx = 0.1;
    double xend = 40.;  // km?
    double xstart = 0.;
    double tend = 4.*24.;  // h
    double tstart = 0.;
    int nx = int((xend - xstart) / dx);
    int nt = int((tend - tstart) / dt);
    std::vector<double> init_cond_adv_diff(nx);
    double left_bound = 0.; // kg per m3
    for (int i = 0; i < nx; ++i) {
        init_cond_adv_diff[i] = init_cond_fun(i * dx);
    }

//    double adv_coeff = -0.1/24.;  // km/h
    double adv_coeff = -0.5;  // km/h
//    double diff_coeff = 0.1/24.;  // km^2/h
    double diff_coeff = 0.5;  // km^2/h

    std::vector<std::vector<double>> diss_oxygen(nt, std::vector<double>(nx, 0.));
    diss_oxygen[0] = init_cond_adv_diff;

    for (int i = 1; i < nt; ++i) {
        diss_oxygen[i] = step_cn(dt, dx, nx, diss_oxygen[i - 1], left_bound, adv_coeff, diff_coeff);
        for (int j = 0; j < nx; ++j) {
            diss_oxygen[i][j] = step_rk6(dt, diss_oxygen[i][j], i, reaction_term);
        }
    }

    std::ofstream outFile("../results.txt");
//    for (int i = 0; i < diss_oxygen.size(); ++i) {
//        for (int j = 0; j < diss_oxygen[0].size(); ++j) {
//            outFile << diss_oxygen[i][j] << " ";
//        }
//        outFile << std::endl;
//    }
    for (int j = 0; j < nx; ++j) {
        for (int i = 0; i < nt; ++i) {
            outFile << i*dt << " " << j*dx << " " << diss_oxygen[i][j] << std::endl;
        }
        outFile << std::endl;
//        outFile << std::endl;
    }
    return 0;
}

double init_cond_fun(double x) {
//    return sin(x);
    return 2.; // g per m3
}

double temp(double t){
    // from 19 to 21, use wolfram system of eqs with 1.5, -0.5
    return 1.*(sin(2.*M_PI*t/24.)+0.5) - (-32./9.);
}

double photosynthesis(double t){
    double P_max = 1.75;  // g / m^3 / h
    return P_max*(sin(2.*M_PI*t/24.)+0.5);
//    return 1;  // g / m^3 / h
}

double reaction_term(double t, double c) {
    // TODO:
    // r1 with variable temp (of water) (variable in time AND space?)
    // r2 with variable depth (and temp?)
    // find what is weird about r2, is C missing?
    // check unit and get tkn for r3
    // c_sat (variable!) for r5
    // are we going to use the wikipedia formula for variable in r5?
    // where can I find R for respiration?
    // and P?
    // totally change r6



    double k20 = 0.23/24.;  // h^-1, or less
    double k_alpha = k20*pow(1.047, (temp(t) - 20));

    double const L = 4.5;  //  g/m^3

    double r1 = -k_alpha*L;

    ///////////////////////////////////////////////////////////////////
    double h = 0.75;  // m, probably should be variable MAYBE UNIT PROBLEM

    double r2 = 0.5/h;
    ///////////////////////////////////////////////////////////////////

    double kn = 0.25/24./24.;  //h^-2

    double tkn = 0.;  // to be found
    double ln = 4.5 + tkn;  //unit? g/m^3?

    double r3 = kn*ln;
    ///////////////////////////////////////////////////////////////////

    double r4 = 0.;  // supposed to stay like this

    ///////////////////////////////////////////////////////////////////
    double k_a = 0.4/24.;  // h^-1
    double c_sat = 8.5;  // g/m^3

    double r5 = - k_a * (c - c_sat);
    ///////////////////////////////////////////////////////////////////
    double const R =  1.;

    double r6 = photosynthesis(t) - R;


//    return r6;
    return r1 + r2 + r3 + r4 + r5 + r6;
}

std::vector<double> use_thomas(std::vector<double> a,
                               std::vector<double> b,
                               std::vector<double> c,
                               std::vector<double> d,
                               int n) {
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

std::vector<double> step_cn(double dt,
                            double dx,
                            int nx,
                            std::vector<double> y_prev,
                            double left_boundary_cond,
                            double adv_coeff,
                            double diff_coeff) {
    std::vector<double> y(nx);
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
            d[j] = y_prev[j-1];
//            d[j] = 0.;
        } else {
            d[j] = (1. / dt - diff_coeff / dx / dx) * y_prev[j] +
                   (0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y_prev[j + 1] +
                   (-0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y_prev[j - 1];
        }
    }
    y = use_thomas(a, b, c, d, nx);
    return y;
}

double step_rk6(double step,
                double previous_value,
                int iteration_number,
                double(*fun)(double, double)) {
    double t, next_value, k0, k1, k2, k3, k4, k5, k6, k7;
    t = step * iteration_number;
    k0 = step * fun(t, previous_value);
    k1 = step * fun(t + 1. / 9. * step, previous_value + 1. / 9. * k0);
    k2 = step * fun(t + 1. / 6. * step, previous_value + 1. / 24. * (k0 + 3 * k1));
    k3 = step * fun(t + 1. / 3. * step, previous_value + 1. / 6. * (k0 - 3 * k1 + 4 * k2));
    k4 = step * fun(t + 0.5 * step, previous_value + 1. / 8. * (k0 + 3 * k3));
    k5 = step * fun(t + 2. / 3. * step, previous_value + 1. / 9. * (17 * k0 - 63 * k1 + 51 * k2 + k4));
    k6 = step *
         fun(t + 5. / 6. * step,
             previous_value + 1. / 24. * (-22 * k0 + 33 * k1 + 30 * k2 - 58 * k3 + 34 * k4 + 3 * k5));
    k7 = step * fun(t + step,
                    previous_value +
                    1. / 82. * (281 * k0 - 243 * k1 - 522 * k2 + 876 * k3 - 346 * k4 - 36 * k5 + 72 * k6));
    next_value = previous_value + 1. / 840. * (41 * (k0 + k7) + 216 * (k2 + k6) + 27 * (k3 + k5) + 272 * k4);
    return next_value;
}
