/*
 *-----------------------------------------------------------------------------
 * This file is part of Interval Squeeze model.
 *
 * Interval Squeeze model is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * Interval Squeez model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Interval Squeeze model.  If not, see
 * <https://www.gnu.org/licenses/>.
 *-----------------------------------------------------------------------------
 */

#include "weather_simulator.h"

//-----------------------------------------------------------------------------

void Weather_simulator::read_climate(const std::string& file_name
                                     , const bool header
                                     , const char delimiter)
{
    m_climate_variables.clear();

    std::string file_path {Global::project_directory + "data/in/climate/" + file_name};
    std::ifstream ifile(file_path);
    std::string line {""};

    if (!ifile)
    {
      std::cerr << "ERROR! Could not open the climate input file:\n";
      std::cerr << file_path << '\n';
      std::cerr << '\a' << '\n';  // beep!
      std::exit(EXIT_FAILURE);
    }

    // skip header row
    if (header == true)
    {
        std::getline(ifile, line);
    }

    while (std::getline(ifile, line))
    {
        if (std::all_of(std::cbegin(line), std::cend(line), ::isspace))
        {
            break;
        }

        std::vector<double> weather_variables_set;

        std::istringstream iss{line};

        // read the tokens from current line separated by delimiter
        std::vector<std::string> tokens;  // excel row
        std::string token;                // excel cell

        while (std::getline(iss, token, delimiter))
        {
            tokens.push_back(token);
        }

        // Process the tokens
        ///@note NOTE cppcheck suggest to use std::transform, but this algorithm
        /// does not guarantee in-order application...
        for (const auto& elem : tokens)
        {
            weather_variables_set.push_back(std::stod(elem));
        }

        m_climate_variables.push_back(weather_variables_set);
    }
}

//-----------------------------------------------------------------------------

const std::vector<double>& Weather_simulator::get_random_year()
{
    assert(!std::empty(m_climate_variables));

    m_next_year = Global::prng.uniform_unsigned(
                      0
                      , std::size(m_climate_variables) - 1
                  );

    return m_climate_variables.at(m_next_year);
}

//-----------------------------------------------------------------------------

const std::vector<double>& Weather_simulator::get_consecutive_year()
{
    assert(!std::empty(m_climate_variables));

    // if last year, then select first year
    if (m_next_year > std::size(m_climate_variables) - 1)
    {
        m_next_year = 0;
    }

    const std::vector<double>& weather { m_climate_variables.at(m_next_year) };

    ++m_next_year;

    return weather;
}
