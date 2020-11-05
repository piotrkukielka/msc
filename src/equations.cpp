//
// Created by Piotr on 02-Jul-20.
//
#include <cmath>
#include <iostream>
#include <fstream>
#include "equations.h"

const double O2_GRAMM3_TO_CO2_MOLELITER = 1. / 31.9988 / 1000.;


double ReactionTerm::evaluate(double t, double c) {
    double temp_value = temp(t);
    double mean_velocity = 0.317; //m/s

    double k20 = 0.23/24.;  // h^-1, or less //OK  //max przedzialu
    double k_alpha = k20*pow(1.047, (temp(t) - 20.)); //OK
//    double k_sed = 0.0558*(mean_velocity*10)*pow(480., -2./3.); //cm/s
//    k_sed = k_sed*36; //m/h
//    std::cout << k_sed << std::endl;
    // todo: bylo 600 dla co2
    double k_sed = 0.021;  //g/m^2/h

    double const L = 5.2;  //  g/m^3 //TODO bylo 4.5 // dwie rozne wartosci //5.2 to wios, 4.5 to grabinska
    double r1 = -k_alpha*L;
//    std::cout << r1 << std::endl;
    ///////////////////////////////////////////////////////////////////
    double h = 0.87;  // m, probably should be variable MAYBE UNIT PROBLEM //OK
    double r2 = -k_sed/h * pow(1.065, temp_value-20);
    ///////////////////////////////////////////////////////////////////
    double kn = 0.01;  //h^-1 // srodek przedzialu
    double tkn = 0.95;  // g/m^3  // chyba OK
    double ln = 4.57 * tkn;  //unit? g/m^3?
    double r3 = -kn*ln;
//    std::cout << r3 << std::endl;
    ///////////////////////////////////////////////////////////////////
//    double r4 = 0.000035;  // supposed to stay like this //
    double r4 = 0.;  // supposed to stay like this
    ///////////////////////////////////////////////////////////////////
    double k_a = 0.0175;  // h^-1 //srodek przedzialu dla large river of low velocity
    double c_sat = 14.652 - 0.41022*temp_value + 0.007991*pow(temp_value, 2) - 7.7774E-5 * pow(temp_value, 3); // g/m^3
//    double r5 = k_a * (c_sat - c);
//TODO: dla DIC
    double r5 = 0.02/h; //dwa milimole per h per mkw wiec dzielone przez glebokosc
    r5 = r5*0.001;  //per m^3 to per liter
    ///////////////////////////////////////////////////////////////////
    double r6 = photosynthesis(t, h);

//    return (r1 + r2 + r3 + r4 + r5 + r6); // DO
//    return (r2 + r4 + r5 + r6); // DO bez
//    return (r2 - r4 - r6)*O2_GRAMM3_TO_CO2_MOLELITER; // - r5; // DIC // nitryfikacja out, BOD out
    return (-r2 - r6)*O2_GRAMM3_TO_CO2_MOLELITER - r5 + r4; // - r5; // DIC // nitryfikacja out, BOD out
    // co z war brzeg, szczegolnie prawym?
}

double ReactionTerm::photosynthesis(double t, double h) {
    // old
//    // TODO: starting point according to first time point?
//    double P_max = 0.43;  // g / m^3 / h
////    return P_max*(cos(2.*M_PI*t/24.)+0.5);  //  start 12AM, RIGHT?
//    return P_max*(sin(2.*M_PI*(t)/24.)+0.5);  //  start 5AM ??????
//    // todo:adding t+
//    // todo: jak -0 to srodek 6h nocy jest w srodku faktycznej nocy
    //new
//    double R = 0.;
//    double P_max = 0.28/h + R;  // g / m^3  // should this be per h?
//    double P = P_max*sin(2.*M_PI*(t-1)/24.);  // -1 to move 6hrs back compared to the okres (half of okres)
//    // to trzeba poprawic, bo on nie jest rozicagniety, tylko ma 12h
//    if(int(2.*M_PI*(t-5)/24.) % 2 == 0){  // -5 because 5AM and it would give us plaskie w zlym miejscu, to tylko czas nie ustawianie sinusa (?)
//        return P-R;
//    }else{
//        return -R;
//    }
    double timeline_start = 5.;
    bool is_day = false;
    double time_of_day = fmod(t + timeline_start, 24.);
    if(time_of_day > 5 and time_of_day < 21){
        is_day = true;
    }

    double R = 0.09; //z lipca 0.09 g/ m^3
    double P_max = 0.32/h + R;  // g / m^3  // should this be per h? // bylo 0.32 z lipca
    double P = P_max*sin(2.*M_PI*(time_of_day-timeline_start)/32.);  // -1 to move 6hrs back compared to the okres (half of okres)
    // okres nie ma 24h, tylko 16*2=32 bo tyle trwa doba od 5 do 21

    if(is_day){  // -5 because 5AM and it would give us plaskie w zlym miejscu, to tylko czas nie ustawianie sinusa (?)
        return P-R;
    }else{
        return -R;
    }
}

double ReactionTerm::temp(double t){
    double a = 1.51567690e+00;
    double b = 4.35966961e-03;
    double c = -1.63198211e+00;
    double d =  2.21646753e+01;
    return a*sin(b*(t * 60) + c) + d;
}

double InitialCondition::evaluate(double x) {
//    return 3.3525; // O2 bondary
    return 0.002595183535082; // DIC bondary mean
}

double LeftBoundaryCondition::evaluate(double t) {
//    return 3.3525; // O2 bondary mean
    return 0.002595183535082; // DIC bondary mean
}

double RightBoundaryCondition::evaluate(double t) {
    double a = 1.96756534e+00;
    double b = 4.37692588e-03;
    double c = -7.83028109e+00;
    double d =  9.13417396e+00;
//    return a*sin(b*(t * 60) + c) + d;
// commented cuz doing for co2
}
