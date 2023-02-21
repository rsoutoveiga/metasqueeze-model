#include "fire_simulator.h"

int Fire_simulator::get_fire_interval_truncated_normal(
        const int mean, int lower_limit)
{
    int fire_interval {0};

    if (lower_limit >= mean)
    {
        std::cout << "WARNING: fire_lower_cut is greater than fire_mean\n"
                  << "fire_interval_lower_cut is truncated to fire_mean: \n"
                  << "fire_interval_lower_cut truncates from "
                  << lower_limit << " to " << mean;
      lower_limit = mean;
    }

    const int upper_limit { (mean - lower_limit) + mean };

    // range rule for the standard deviation
    double sd {static_cast<double>(upper_limit - lower_limit) / 4.0};

    if (sd < 0.0)
    {
      sd = 0.0;
    }

    int attempts {0};

    do
    {
        ++attempts;
        if (attempts > 1000000)
        {
            std::cerr << "too many attempts to find a suitable fire interval\n";
            std::exit(EXIT_FAILURE);
        }

        fire_interval = static_cast<int>(
                    Global::prng.gaussian(mean, sd) + 0.5);

    } while (fire_interval < lower_limit || fire_interval > upper_limit);

    assert(fire_interval >= 0);  // fire_interval cannot be zero or below zero

    return fire_interval;
}

//-----------------------------------------------------------------------------

int Fire_simulator::get_fire_interval_weibull(
        const double shape, const double scale, const int lower_limit)
{
    int fire_interval {
        static_cast<int>(Global::prng.weibull(shape, scale) + 0.5)
    };

    if (fire_interval < lower_limit)
    {
        fire_interval = lower_limit;
    }

    // Ensures(fire_interval >= 0)
    return fire_interval;
}

//-----------------------------------------------------------------------------

double Fire_simulator::get_fire_size_truncated_normal(
        const double mean, double lower_limit)
{
    double fire_size {0.0};

    if (lower_limit >= mean)
    {
        std::cout << "WARNING: fire_size_lower_limit is greater than mean\n"
                  << "fire_size_lower_limit is truncated to fire_size_mean: \n"
                  << "fire_size_lower_limit truncates from "
                  << lower_limit << " to " << mean;
      lower_limit = mean;
    }

    const double upper_limit { (mean - lower_limit) + mean };

    // range rule for the standard deviation
    double sd { (upper_limit - lower_limit) / 4.0 };

    if (sd < 0.0)
    {
      sd = 0.0;
    }

    int attempts {0};

    do
    {
        ++attempts;
        if (attempts > 1000000)
        {
            std::cerr << "too many attempts to find a suitable fire size\n";
            std::exit(EXIT_FAILURE);
        }

        fire_size = Global::prng.gaussian(mean, sd);

    } while (fire_size < lower_limit || fire_size > upper_limit);

    assert(fire_size >= 0.0);  // fire_interval cannot be zero or below zero

    return fire_size;
}
