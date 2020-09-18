//
// Created by Piotr on 02-Jul-20.
//

#ifndef MSC_0505_EQUATIONS_H
#define MSC_0505_EQUATIONS_H

#include <vector>

class Equation{
public:
    virtual double evaluate(double t, double c) = 0;
};

class ReactionTerm : public Equation{
    static double photosynthesis(double t);
    static double temp(double t);

public:
    double evaluate(double t, double c) override;
};

class InitialCondition{
public:
    static double evaluate(double x);
};

class LeftBoundaryCondition{
public:
    static double evaluate(double t);
};

class RightBoundaryCondition{
public:
    static double evaluate(double t);
};

#endif //MSC_0505_EQUATIONS_H
