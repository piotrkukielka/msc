//
// Created by Piotr on 02-Jul-20.
//
#include <cmath>
#include <utility>
#include "equations.h"

const double O2_GRAMM3_TO_CO2_MOLELITER = 1. / 31.9988 / 1000.;


double ReactionTerm::evaluate(double t, double c) {
    double r1 =0 , r2 = 0, r3 = 0, r4 = 0, r5 = 0, r6 = 0;
    double temp_value = temp(t);

    ///////////////////////////////////////////////////////////////////
    double k20 = 0.23/24.;  // h^-1, or less
    double k_alpha = k20*pow(1.047, (temp(t) - 20.)); //OK
    double const L = 5.2;  //  g/m^3, or less by Grabinska
    r1 = -k_alpha*L;
    ///////////////////////////////////////////////////////////////////
    double k_sed = 0.021;  //g/m^2/h
    double h = 0.87;  // m
    r2 = -k_sed/h * pow(1.065, temp_value-20);
    ///////////////////////////////////////////////////////////////////
    double kn = 0.01;  // h^-1
    double tkn = 0.95;  // g/m^3
    double ln = 4.57 * tkn;  // g/m^3
    r3 = -kn*ln;
    ///////////////////////////////////////////////////////////////////
    if (this->modeled_variable == "DIC") {
        r4 = 0.000035;
    }
    if (this->modeled_variable == "O2") {
        r4 = 0.;
    }
    ///////////////////////////////////////////////////////////////////
    double k_a = 0.0175;  // h^-1 large river of low velocity
    double c_sat = 14.652 - 0.41022*temp_value + 0.007991*pow(temp_value, 2) - 7.7774E-5 * pow(temp_value, 3); // g/m^3
    if (this->modeled_variable == "DIC") {
        r5 = 0.02/h;
        r5 = r5*0.001;  //per m^3 to per liter // DIC
    }
    if (this->modeled_variable == "O2") {
        r5 = k_a * (c_sat - c); //O2
    }
    ///////////////////////////////////////////////////////////////////
    r6 = photosynthesis(t, h);
    ///////////////////////////////////////////////////////////////////
    if (this->modeled_variable == "DIC") {
        return (-r2 - r6)*O2_GRAMM3_TO_CO2_MOLELITER - r5 + r4;
    } else if (this->modeled_variable == "O2") {
        return (r2 + r4 + r5 + r6);
    } else {
        return 0;
    }
}

double ReactionTerm::photosynthesis(double t, double h) {
    double timeline_start = 5.;
    bool is_day = false;
    double time_of_day = fmod(t + timeline_start, 24.);
    if(time_of_day > 5 and time_of_day < 21){
        is_day = true;
    }

    double R = 0.09; // g/ m^3, from July
    double P_max = 0.32/h + R;  // g / m^3
    double P = P_max*sin(2.*M_PI*(time_of_day-timeline_start)/32.);  // -1 to move 6hrs back compared to the half of the period
    // the sines period is not 24, but 32 because half of it should be the length of the day

    if(is_day){
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

ReactionTerm::ReactionTerm(std::string modeledVariable) : modeled_variable(std::move(modeledVariable)) {}

double InitialCondition::evaluate(double x) {
    if (this->modeled_variable == "DIC") {
        return 0.002595183535082; // DIC bondary mean
    } else if (this->modeled_variable == "O2") {
        return 3.3525; // O2 bondary
    } else {
        return 0;
    }
}

InitialCondition::InitialCondition(std::string modeledVariable) : modeled_variable(std::move(modeledVariable)) {}

double LeftBoundaryCondition::evaluate(double t) {
    if (this->modeled_variable == "DIC") {
        return 0.002595183535082; // DIC bondary mean
    } else if (this->modeled_variable == "O2") {
        return 3.3525; // O2 bondary mean
    } else {
        return 0;
    }
}

LeftBoundaryCondition::LeftBoundaryCondition(std::string modeledVariable) : modeled_variable(std::move(modeledVariable)) {}

double RightBoundaryCondition::evaluate(double t) {
    double a = 1.96756534e+00;
    double b = 4.37692588e-03;
    double c = -7.83028109e+00;
    double d =  9.13417396e+00;
    return a*sin(b*(t * 60) + c) + d;
}

RightBoundaryCondition::RightBoundaryCondition(std::string modeledVariable) : modeled_variable(std::move(modeledVariable)) {}
