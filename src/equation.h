#ifndef EQUATION_H
#define EQUATION_H

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cassert>

///@todo TODO add Expects and Ensures
namespace Equation {

    inline double linear(const double a, const double b, const int age)
    {
        return (a * static_cast<double>(age) + b);
    }

    inline double exponential(const double a, const double b, const int age)
    {
        return std::exp(a * static_cast<double>(age) + b);
    }

    inline double logistic_growth(const double a, const double b, const double c, const int age)
    {
        return a / (1.0 + std::exp(c * (b - static_cast<double>(age))));
    }

    double mean(const std::vector<int>& v);
    double mean(const std::vector<long long int>& v);
    double mean(const std::vector<double>& v);
    double sd(const std::vector<int>& v);
    double sd(const std::vector<long long int>& v);
    double sd(const std::vector<double>& v);

    void mean_and_sd(const std::vector<int>& v, double& mean, double& sd);
}

#endif // EQUATION_H
