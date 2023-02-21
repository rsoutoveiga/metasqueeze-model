#include "equation.h"

//-----------------------------------------------------------------------------

double Equation::mean(const std::vector<int>& v)
{
    double result {0.0};

    if (!std::empty(v))
    {
        result = std::accumulate(
                    std::begin(v)
                    , std::end(v)
                    , 0.0
                    , std::plus<double>()
                    )
                / std::size(v);
    }

    // Ensures
    assert(result >= 0.0);
    return result;
}

//-----------------------------------------------------------------------------

double Equation::mean(const std::vector<long long int>& v)
{
    double result {0.0};

    if (!std::empty(v))
    {
        result = std::accumulate(
                    std::begin(v)
                    , std::end(v)
                    , 0.0
                    , std::plus<double>()
                    )
                / std::size(v);
    }

    // Ensures
    assert(result >= 0.0);
    return result;
}

//-----------------------------------------------------------------------------

double Equation::mean(const std::vector<double>& v)
{
    double result {0.0};

    if (!std::empty(v))
    {        
        result = std::accumulate(
                    std::begin(v)
                    , std::end(v)
                    , 0.0
                    , std::plus<double>()
                    )
                / std::size(v);
    }

    // Ensures
    assert(result >= 0.0);
    return result;
}

//-----------------------------------------------------------------------------

double Equation::sd(const std::vector<int>& v)
{
    double stdev {0.0};

    if (!std::empty(v))
    {
        double mean {
            std::accumulate(
                        std::begin(v)
                        , std::end(v)
                        , 0.0
                        , std::plus<double>()
                        )
                    / std::size(v)
        };

        double accum {0.0};

        std::for_each(std::begin(v), std::end(v), [&](const double x) {
            accum += (x - mean) * (x - mean);
        });

        if (std::size(v) > 1)
        {
            stdev = std::sqrt(accum / (std::size(v) - 1));
        }
        else
        {
            stdev = 0.0;
        }
    }

    // Ensures
    assert(stdev >= 0.0);
    return stdev;
}

//-----------------------------------------------------------------------------

double Equation::sd(const std::vector<long long int>& v)
{
    double stdev {0.0};

    if (!std::empty(v))
    {
        double mean {
            std::accumulate(
                        std::begin(v)
                        , std::end(v)
                        , 0.0
                        , std::plus<double>()
                        )
                    / std::size(v)
        };

        double accum {0.0};

        std::for_each(std::begin(v), std::end(v), [&](const double x) {
            accum += (x - mean) * (x - mean);
        });

        if (std::size(v) > 1)
        {
            stdev = std::sqrt(accum / (std::size(v) - 1));
        }
        else
        {
            stdev = 0.0;
        }
    }

    // Ensures
    assert(stdev >= 0.0);
    return stdev;
}

//-----------------------------------------------------------------------------

double Equation::sd(const std::vector<double>& v)
{
    double stdev {0.0};

    if (!std::empty(v))
    {
        double mean {
            std::accumulate(
                        std::begin(v)
                        , std::end(v)
                        , 0.0
                        , std::plus<double>()
                        )
                    / std::size(v)
        };

        double accum {0.0};

        std::for_each(std::begin(v), std::end(v), [&](const double x) {
            accum += (x - mean) * (x - mean);
        });

        if (std::size(v) > 1)
        {
            stdev = std::sqrt(accum / (std::size(v) - 1));
        }
        else
        {
            stdev = 0.0;
        }
    }

    // Ensures
    assert(stdev >= 0.0);
    return stdev;
}

//-----------------------------------------------------------------------------

void Equation::mean_and_sd(const std::vector<int>& v, double& mean, double& sd)
{
    if (!std::empty(v))
    {
        mean = std::accumulate(
                    std::begin(v)
                    , std::end(v)
                    , 0.0
                    , std::plus<double>()
                    )
                / std::size(v);

        double accum {0.0};

        std::for_each(std::begin(v), std::end(v), [&](const double x) {
            accum += std::pow(x - mean, 2);
        });

        if (std::size(v) > 1)
        {
            sd = std::sqrt(accum / (std::size(v) - 1));
        }
        else
        {
            sd = 0.0;
        }
    }
}
