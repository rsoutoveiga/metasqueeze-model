#ifndef FIRE_SIMULATOR_H
#define FIRE_SIMULATOR_H

#include <cassert>
#include <iostream>

#include "global.h"

namespace Fire_simulator
{
    int get_fire_interval_truncated_normal(const int mean, int lower_limit);
    int get_fire_interval_weibull(const double shape
                                  , const double scale
                                  , const int lower_limit);
    double get_fire_size_truncated_normal(const double mean
                                          , double lower_limit);
}

#endif // FIRE_SIMULATOR_H
