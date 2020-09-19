//
// Created by Piotr on 02-Jul-20.
//
#include <cmath>
#include <iostream>
#include "equations.h"

const double O2_GRAMM3_TO_CO2_MOLELITER = -1. / 31.9988 / 1000.;


double ReactionTerm::evaluate(double t, double c) {
    double temp_value = temp(t);
    double mean_velocity = 0.317; //m/s

    double k20 = 0.23/24.;  // h^-1, or less //OK  //max przedzialu
    double k_alpha = k20*pow(1.047, (temp(t) - 20.)); //OK
    double k_sed = 0.0558*(mean_velocity*100)*pow(480., -2./3.);  // u* used in meters/hour, as 10percent of normal flow vel 0.3m per s czyli powinno byc 30 // IDK
    // todo: bylo 600 dla co2

    double const L = 5.2;  //  g/m^3 //TODO bylo 4.5 // dwie rozne wartosci //5.2 to wios, 4.5 to grabinska
    double r1 = -k_alpha*L;
    ///////////////////////////////////////////////////////////////////
    double h = 0.87;  // m, probably should be variable MAYBE UNIT PROBLEM //OK
//    double r2 = -k_sed/h;
    double r2 = -k_sed/h * pow(1.065, temp_value-20);
    ///////////////////////////////////////////////////////////////////
    double kn = 0.01;  //h^-2 // srodek przedzialu
    double tkn = 0.95;  // g/m^3  // chyba OK
    double ln = 4.57 * tkn;  //unit? g/m^3?
    double r3 = -kn*ln;
    ///////////////////////////////////////////////////////////////////
    double r4 = 0.;  // supposed to stay like this
    ///////////////////////////////////////////////////////////////////
    double k_a = 0.0175;  // h^-1 //srodek przedzialu dla large river of low velocity
    double c_sat = 14.652 - 0.41022*temp_value + 0.007991*pow(temp_value, 2) - 7.7774E-5 * pow(temp_value, 3); // g/m^3
    double r5 = k_a * (c_sat - c);
    ///////////////////////////////////////////////////////////////////
    double const R =  0.; //TODO bylo 0.1 dla rownowagi
    double r6 = photosynthesis(t) - R;

    // only for CO2
//    r5 = -r5; //TODO
    return (r1 + r2 + r3 + r4 + r5 + r6);//*O2_GRAMM3_TO_CO2_MOLELITER;
//    return r2 + r4 + r5 + r6;
}

double ReactionTerm::photosynthesis(double t) {
    // TODO: starting point according to first time point?
    double P_max = 0.43;  // g / m^3 / h
//    return P_max*(cos(2.*M_PI*t/24.)+0.5);  //  start 12AM, RIGHT?
    return P_max*(sin(2.*M_PI*(t)/24.)+0.5);  //  start 5AM ??????
}

double ReactionTerm::temp(double t){
    // from 19 to 21, use wolfram system of eqs with 1.5, -0.5
//    return 1.*(cos(2.*M_PI*t/24.)+20);  // start 12AM, RIGHT?
//    return 1.*(sin(2.*M_PI*(t)/24.)+20);  // start 5AM ????
    double a = 1.51567690e+00;
    double b = 4.35966961e-03;
    double c = -1.63198211e+00;
    double d =  2.21646753e+01;
    return a*sin(b*(t * 60) + c) + d;
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
    return 3.1; // O2 bondary
//    return 0.002596297044703; // DIC bondary 2 pomiar
    // t is in hrs, convering to minutes
    // 31.8 is the time starting point
// to chyba wlasnie jest stale z tego sinusa
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
    return 3.1; // O2 bondary
//    return 0.002596297044703; // DIC bondary 2 pomiar
}

double RightBoundaryCondition::evaluate(double t) {
    double a = 1.96756534e+00;
    double b = 4.37692588e-03;
    double c = -7.83028109e+00;
    double d =  9.13417396e+00;
    return a*sin(b*(t * 60) + c) + d;
}
