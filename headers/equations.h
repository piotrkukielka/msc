//
// Created by Piotr on 02-Jul-20.
//

#ifndef MSC_0505_EQUATIONS_H
#define MSC_0505_EQUATIONS_H

#include <vector>
#include <string>

class Equation{
public:
    virtual double evaluate(double t, double c) = 0;
};

class ReactionTerm : public Equation{
    static double temp(double t);
    std::string modeled_variable;

public:
    explicit ReactionTerm(std::string modeledVariable);

    double evaluate(double t, double c) override;

    double photosynthesis(double t, double h);
};

class InitialCondition{
    std::string modeled_variable;
public:
    explicit InitialCondition(std::string modeledVariable);

    double evaluate(double x);
};

class LeftBoundaryCondition{
    std::string modeled_variable;
public:
    explicit LeftBoundaryCondition(std::string modeledVariable);

    double evaluate(double t);
};

class RightBoundaryCondition{
    std::string modeled_variable;
public:
    explicit RightBoundaryCondition(std::string modeledVariable);

    static double evaluate(double t);
};

#endif //MSC_0505_EQUATIONS_H
