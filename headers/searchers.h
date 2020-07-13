//
// Created by Piotr on 02-Jul-20.
//

#ifndef MSC_0505_SEARCHERS_H
#define MSC_0505_SEARCHERS_H

class Bounds{
    double lower;
    double upper;
    double interval;

public:
    int get_num_of_iters();
    double get_value(int i);
    Bounds(double lower, double upper, double interval);

};

class GridSearch {

public:
    void search(Bounds advCoeffs, Bounds diffCoeffs, std::vector<double> rs);

};

class Simulation{
    double adv_coeff;
    double diff_coeff;
    std::vector<double> rs;

    void save(std::vector<std::vector<double>> data);
    std::vector<std::vector<double>> run();
public:
    Simulation(double advCoeff, double diffCoeff, std::vector<double> rs);

    void run_and_save();

};



#endif //MSC_0505_SEARCHERS_H
