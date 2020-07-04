//
// Created by Piotr on 02-Jul-20.
//

#ifndef MSC_0505_SEARCHERS_H
#define MSC_0505_SEARCHERS_H

class Bounds{
    double lower;
    double upper;
    double interval;

    int find_num_of_iters();
    std::vector<double> find_values();

public:
    int num_of_iters;
    std::vector<double> values;
    Bounds(double lower, double upper, double interval);

};

class GridSearch {

public:
    void search();

};

class Simulation{
    void save(std::vector<std::vector<double>> data);
    std::vector<std::vector<double>> run(double adv_coeff_bounds, double diff_coeff_bounds, std::vector<double> rs);
public:
    void run_and_save(double adv_coeff_bounds, double diff_coeff_bounds, std::vector<double> rs);

};



#endif //MSC_0505_SEARCHERS_H
