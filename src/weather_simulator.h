#ifndef WEATHER_SIMULATOR_H
#define WEATHER_SIMULATOR_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <memory>
#include <cctype>

#include "global.h"

class Weather_simulator
{
private:
    unsigned m_next_year {0};

    std::vector<std::vector<double>> m_climate_variables;

public:
    Weather_simulator()
        : m_next_year(0)
        , m_climate_variables()
    {}

    void read_climate(const std::string& file_name
                      , const bool header = true
                      , const char delimiter = ',');

    const std::vector<double>& get_random_year();
    const std::vector<double>& get_consecutive_year();

    void reset_weather_year() { m_next_year = 0; }
};

#endif // WEATHER_SIMULATOR_H
