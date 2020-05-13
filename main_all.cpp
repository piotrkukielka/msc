#include <boost/math/tools/roots.hpp>
#include <iostream>
#include <vector>
#include <fstream>


double _analitical_solution(double t);

double testfun(double t, double y);

struct implicit_fun_wrapper {
    implicit_fun_wrapper(double step, double t,
                         double y_curr, double (*fun)(double, double)) {
        _step = step;
        _t = t;
        _y_curr = y_curr;
        _fun = fun;
    }

    double operator()(double x) {
        double to_solve = _y_curr + _step * _fun(_t + _step, x) - x;
        return to_solve;
    }

private:
    double _step;
    double _t;
    double _y_curr;

    double (*_fun)(double, double);
};

std::vector<double> solve_with_forward_euler(double step,
                                             double start,
                                             double end,
                                             double initial_cond,
                                             double(*fun)(double, double)) {
    int n = int((end - start) / step);
    std::vector<double> y(n);
    y[0] = initial_cond;
    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = y[i] + step * fun(i * step, y[i]);  // or maybe i*step?
    }
    return y;
}

void print_diff(std::vector<double> vec1, std::vector<double> vec2) {
    for (int i = 0; i < vec1.size(); ++i) {
        std::cout << vec1[i] - vec2[i] << std::endl;
    }
}


std::vector<double> analitical_solution(double step,
                                        double start,
                                        double end) {
    int n = int((end - start) / step);
    std::vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = _analitical_solution(i * step);
    }
    return y;
}


std::vector<double> solve_with_midpoint(double step,
                                        double start,
                                        double end,
                                        double initial_cond,
                                        double(*fun)(double, double)) {
    int n = int((end - start) / step);
    std::vector<double> y(n);
    y[0] = initial_cond;
    for (int i = 0; i < n - 1; i++) {
        double first_call = fun(step * i, y[i]);
        double second_call = fun(step * i + 0.5 * step, y[i] + 0.5 * step * first_call);
        y[i + 1] = y[i] + step * second_call;
    }
    return y;
}

std::vector<double> solve_with_rk4(double step,
                                   double start,
                                   double end,
                                   double initial_cond,
                                   double(*fun)(double, double)) {
    int n = int((end - start) / step);
    std::vector<double> y(n);
    y[0] = initial_cond;
    double t;
    double k1, k2, k3, k4;
    for (int i = 0; i < n - 1; i++) {
        t = step * i;
        k1 = fun(t, y[i]);
        k2 = fun(t + 0.5 * step, y[i] + 0.5 * step * k1);
        k3 = fun(t + 0.5 * step, y[i] + 0.5 * step * k2);
        k4 = fun(t + step, y[i] + step * k3);
        y[i + 1] = y[i] + 1. / 6. * step * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return y;
}

std::vector<double> solve_with_rk6(double step,
                                   double start,
                                   double end,
                                   double initial_cond,
                                   double(*fun)(double, double)) {
    int n = int((end - start) / step);
    std::vector<double> y(n);
    y[0] = initial_cond;
    double t;
    double k0, k1, k2, k3, k4, k5, k6, k7;
    for (int i = 0; i < n - 1; i++) {
        t = step * i;
        k0 = step*fun(t, y[i]);
        k1 = step*fun(t + 1./9. * step, y[i] + 1./9. * k0);
        k2 = step*fun(t + 1./6. * step, y[i] + 1./24. * (k0 + 3*k1));
        k3 = step*fun(t + 1./3. * step, y[i] + 1./6. * (k0 - 3*k1 + 4*k2));
        k4 = step*fun(t + 0.5 * step, y[i] + 1./8.*(k0 + 3*k3));
        k5 = step*fun(t + 2./3.*step, y[i] + 1./9.*(17*k0 - 63*k1 + 51*k2 + k4));
        k6 = step*fun(t + 5./6.*step, y[i] + 1./24.*(-22*k0 + 33*k1 + 30*k2 - 58*k3 + 34*k4 + 3*k5));
        k7 = step*fun(t + step, y[i] + 1./82.*(281*k0 - 243*k1 - 522*k2 + 876*k3 - 346*k4 - 36*k5 + 72*k6));
        y[i + 1] = y[i] + 1./840.*(41*(k0 + k7) + 216*(k2+k6) + 27*(k3+k5) + 272*k4);
    }
    return y;
}

std::vector<double> solve_with_backward_euler(double step,
                                              double start,
                                              double end,
                                              double initial_cond,
                                              double(*fun)(double, double)) {
    int n = int((end - start) / step);
    std::vector<double> y(n);
    y[0] = initial_cond;
    for (int i = 0; i < n - 1; i++) {
        boost::math::tools::eps_tolerance<double> tol(32);
        implicit_fun_wrapper implicit_fun(step, i, y[i], fun);
        std::pair<double, double> bisect_result = boost::math::tools::bisect(implicit_fun, 0., 10., tol);
        y[i + 1] = 0.5 * (bisect_result.first + bisect_result.second);  // or maybe i*step?
    }
    return y;
}

double _analitical_solution(double t) {
//    return exp(t);
    return sin(t);
}

double testfun(double t, double y) {
//    return y;
    return sqrt(1 - y * y);
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

std::vector<std::vector<double>> solve_with_crank_nicolson(double dt,
                                                           double dx,
                                                           double xstart,
                                                           double xend,
                                                           double tstart,
                                                           double tend,
                                                           std::vector<double> initial_cond,
                                                           std::vector<double> left_boundary_cond,
                                                           std::vector<double> right_boundary_cond) {
    double adv_coeff = 0.9;
    double diff_coeff = 0.3;
    int nx = int((xend - xstart) / dx);
    int nt = int((tend - tstart) / dt);
    std::vector<std::vector<double>> y(nt, std::vector<double>(nx));
    y[0] = initial_cond;
    for (int i = 1; i < nt; i++) {
        y[i][0] = left_boundary_cond[i];
        y[i][nx - 1] = right_boundary_cond[i];
    }
    std::ofstream outFile("../bounds.txt");
    for(int i=0; i<y.size(); ++i){
        for(int j=0; j<y[0].size(); ++j){
            outFile << y[i][j] << ",";
        }
        outFile << std::endl;
    }
    std::vector<double> a(nx, 0.25 * adv_coeff / dx - 0.5 * diff_coeff / dx / dx);
    std::vector<double> b(nx, 1. / dt + diff_coeff / dx / dx);
    std::vector<double> c(nx, -0.25 * adv_coeff / dx - 0.5 * diff_coeff / dx / dx);
    std::vector<double> d(nx);
    for (int i = 1; i < nt; ++i) {
        for (int j = 0; j < nx; ++j) {
            if (j == 0) {
                d[j] = (1. / dt - diff_coeff / dx / dx) * y[i - 1][j] +
                       (0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y[i - 1][j + 1];
            } else if (j == nx) {
                d[j] = (1. / dt - diff_coeff / dx / dx) * y[i - 1][j] +
                       (-0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y[i - 1][j - 1];
            } else {
                d[j] = (1. / dt - diff_coeff / dx / dx) * y[i - 1][j] +
                       (0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y[i - 1][j + 1] +
                       (-0.25 * adv_coeff / dx + 0.5 * diff_coeff / dx / dx) * y[i - 1][j - 1];
            }
        }
        y[i] = use_thomas(a, b, c, d, nx);
        y[i][0] = left_boundary_cond[i];
        y[i][nx - 1] = right_boundary_cond[i];
    }
    return y;
}

double init_cond_fun(double x) {
    return sin(x);
}


int main() {
    double step = 0.0001;
    double start = 0.0;
    double end = M_PI / 2.;
    double initial_cond = 0.;
//    std::vector<double> result_analitical = analitical_solution(step, start, end);
//    std::vector<double> result_euler = solve_with_forward_euler(step, start, end, initial_cond, testfun);
//    std::vector<double> result_midpoint = solve_with_midpoint(step, start, end, initial_cond, testfun);
//    std::vector<double> result_rk4 = solve_with_rk4(step, start, end, initial_cond, testfun);
//    std::vector<double> result_implicit_euler = solve_with_backward_euler(step, start, end, initial_cond, testfun);
//    std::vector<double> result_rk6 = solve_with_rk6(step, start, end, initial_cond, testfun);
//    print_diff(result_analitical, result_euler);
//    print_diff(result_analitical, result_midpoint);
//    print_diff(result_analitical, result_rk4);
//    print_diff(result_analitical, result_implicit_euler);
//    print_diff(result_analitical, result_rk6);


//    std::ofstream outFile("../out.txt");
//    for (int i = 0; i < result_implicit_euler.size(); ++i) {
//        outFile << i << "," << result_analitical[i] << "," << result_implicit_euler[i] << "\n";
//    }



//    std::vector<double> a{ 0, -1, -1, -1 };
//    std::vector<double> b{ 4,  4,  4,  4 };
//    std::vector<double> c{-1, -1, -1,  0 };
//    std::vector<double> d{ 5,  5, 10, 23 };
//    std::vector<double> res = use_thomas(a, b, c, d, 4);



    double dt = 0.01;
    double dx = 0.01;
    double xend = M_PI;
    double xstart = 0.;
    double tend = 2.8;
    double tstart = 0.;
    int nx = int((xend - xstart) / dx);
    int nt = int((tend - tstart) / dt);
    std::vector<double> init_cond_adv_diff(nx);
    std::vector<double> left_bound_adv_diff(nt, 0.);
    std::vector<double> right_bound_adv_diff(nt, 0.);
    for (int i = 0; i < nx; ++i) {
        init_cond_adv_diff[i] = init_cond_fun(i * dx);
    }
    std::vector<std::vector<double>> res_adv_diff = solve_with_crank_nicolson(dt, dx, xstart, xend, tstart, tend,
                                                                              init_cond_adv_diff, left_bound_adv_diff,
                                                                              right_bound_adv_diff);
    std::ofstream outFile("../adv_diff_cn.txt");
    for(int i=0; i<res_adv_diff.size(); ++i){
        for(int j=0; j<res_adv_diff[0].size(); ++j){
            outFile << res_adv_diff[i][j] << ",";
        }
        outFile << std::endl;
    }
    // sprawdzic czy w ktoryms momencie nie sa zle war brzegowe
}

