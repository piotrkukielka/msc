//
// Created by Piotr on 02-Jul-20.
//
#include <cmath>
#include <iostream>
#include "equations.h"

const double O2_GRAMM3_TO_CO2_MOLELITER = -1. / 31.9988 / 1000.;


double ReactionTerm::evaluate(double t, double c) {
    double k20 = 0.23/24.;  // h^-1, or less
    double k_alpha = k20*pow(1.047, (temp(t) - 20.));
    double k_sed = 0.0558*108*pow(600., -2./3.);  // u* used in meters/hour, as 10percent of normal flow vel 0.3

    double const L = 4.5;  //  g/m^3
    double r1 = -k_alpha*L;
    ///////////////////////////////////////////////////////////////////
    double h = 0.93;  // m, probably should be variable MAYBE UNIT PROBLEM
    double r2 = -k_sed/h;
    ///////////////////////////////////////////////////////////////////
    double kn = 0.01;  //h^-2
    double tkn = 0.95;  // g/m^3
    double ln = 4.57 * tkn;  //unit? g/m^3?
    double r3 = -kn*ln;
    ///////////////////////////////////////////////////////////////////
    double r4 = 0.;  // supposed to stay like this
    ///////////////////////////////////////////////////////////////////
    double k_a = 0.04;  // h^-1
    double c_sat = 9.;  // g/m^3
    double r5 = k_a * (c_sat - c);
    ///////////////////////////////////////////////////////////////////
    double const R =  0.;
    double r6 = photosynthesis(t) - R;

    // only for CO2
//    r5 = -r5; //TODO
    return (r1 + r2 + r3 + r4 + r5 + r6);//*O2_GRAMM3_TO_CO2_MOLELITER;
}

double ReactionTerm::photosynthesis(double t) {
    // TODO: starting point according to first time point?
    double P_max = 0.43;  // g / m^3 / h
//    return P_max*(cos(2.*M_PI*t/24.)+0.5);  //  start 12AM, RIGHT?
    return P_max*(sin(2.*M_PI*(t)/24.)+0.5);  //  start 5AM
}

double ReactionTerm::temp(double t){
    // from 19 to 21, use wolfram system of eqs with 1.5, -0.5
//    return 1.*(cos(2.*M_PI*t/24.)+20);  // start 12AM, RIGHT?
    return 1.*(sin(2.*M_PI*(t)/24.)+20);  // start 5AM
}

double InitialCondition::evaluate(double x) {
    // o2, WACHNIEW, old
    double a = 1.968;
    double b = 0.004377;
    double c = 3.401;
    double d = 9.134;
    // co2, fucked
//    double a = 4.28078802e-05;
//    double b = 4.74534478e-03;
//    double c = 3.32547779e+00;
//    double d = 3.15280166e-03;
//    return a*sin(b*(31.8 * 60) + c) + d;  // TEN SINUS JEST NIEZALEZNY OD X, TAK MA BYC
// todo
//    return 3.1; // O2 bondary
//    return 0.002596297044703; // DIC bondary 2 pomiar
    // t is in hrs, convering to minutes
    // 31.8 is the time starting point
// to chyba wlasnie jest stale z tego sinusa
    return 3;
}

double LeftBoundaryCondition::evaluate(double t) {
    // o2, WACHNIEW, old
    double a = 1.968;
    double b = 0.004377;
    double c = 3.401;
    double d = 9.134;
    // co2, fucked
//    double a = 4.28078802e-05;
//    double b = 4.74534478e-03;
//    double c = 3.32547779e+00;
//    double d = 3.15280166e-03;
//    return a*sin(b*((t+31.8) * 60.) + c) + d;  // t is in hrs, converting to minutes
// todo
//    return 3.1; // O2 bondary
//    return 0.002596297044703; // DIC bondary 2 pomiar
    return 3;
}
