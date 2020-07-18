//
// Created by Piotr on 02-Jul-20.
//
#include <cmath>
#include "equations.h"


double ReactionTerm::evaluate(double t, double c) {
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

    return r1 + r2 + r3 + r4 + r5 + r6;
}

double ReactionTerm::photosynthesis(double t) {
    double P_max = 1.75;  // g / m^3 / h
    return P_max*(sin(2.*M_PI*t/24.)+0.5);
}

double ReactionTerm::temp(double t){
    // from 19 to 21, use wolfram system of eqs with 1.5, -0.5
    return 1.*(sin(2.*M_PI*t/24.)+0.5) - (-32./9.);
}

double InitialCondition::evaluate(double x) {
//    return 2.; // g per m3
    return sin(x); // TODO
}
