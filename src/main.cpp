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

///@todo TODO enable -Wconversion in CMakeLists to solve those problems

#include <iostream>
#include <algorithm>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <memory>
#include <cctype>
//#include <cstdio>     // std::exit()
//#include <cstdlib>    // std::getenv()
#include <map>
#include <string>
#include <iomanip>
#include <ctime>
#include <stdexcept>  // std::runtime_error()

#include "environment.h"
#include "parameter_reader.h"
#include "global.h"


int main(int argc, const char* argv[])
{
    // set parent directory
    {
        std::filesystem::path project_directory {argv[0]};
        project_directory.remove_filename();
        project_directory.append("../../");
        Global::project_directory = project_directory;
    }

    std::string name_sim_folder;
    std::string name_sim_file;
    std::string output_folder;

    // read arguments
    switch (argc)
    {
    case 1:    // no arguments, take default simulation
//        name_sim_folder = "default";
//        name_sim_file   = "sim_default.csv";
        name_sim_folder = "manuscript";
        name_sim_file = "test.csv";
        break;
    case 2:
        name_sim_file   = argv[1];
        break;
    case 3:
        name_sim_folder = argv[1];
        name_sim_file   = argv[2];
        break;
    case 4:
        name_sim_folder = argv[1];
        name_sim_file   = argv[2];
        output_folder   = argv[3];
        break;
    default:
        throw std::runtime_error(
                    "There are too many command line arguments:\n"
                    "  argument #1: simulation input file\n"
                    "  argumnet #2: output prefix (optional)\n"
                    "  argument #3: output_folder (optional)\n\n"
                    "NOTE: without arguments, it runs default simulation:\n"
                    "~/data/in/sim/sim_default.txt\n\a"
              );
    }

    const std::filesystem::path sim_file_path {
        Global::project_directory
                + "data/in/sim/"
                + name_sim_folder
                + '/'
                + name_sim_file
    };

    if (!std::filesystem::exists(sim_file_path))
    {
        std::cout << sim_file_path << '\n';
        throw std::runtime_error("The simulation input file does not exist!\n");
    }

    std::ifstream ifs(sim_file_path);

    if (!ifs)
    {
        std::cout << sim_file_path << '\n';
        throw std::runtime_error("The simulation file could not be opened!\n");
    }

    // take local time system to name the output folder and output files
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss_time;
    oss_time << std::put_time(&tm, "%Y_%m_%d_%H%M%S");

    if (std::empty(output_folder) == true)
    {
        output_folder = oss_time.str() + '_' + name_sim_folder;
    }

    std::string line;           // to read line from input simulation file
    std::getline(ifs, line);    // skip header

    while (std::getline(ifs, line))  // read one simulation line
    {
        // skip line if there are only spaces and tabs
        if (std::all_of(std::cbegin(line), std::cend(line), ::isspace))
        {
            continue;
        }

        // print in terminal the simulation number that is running
        std::string sim_number;
        std::istringstream iss(line);
        std::getline(iss, sim_number, ' ');

        std::cout << "Initialize! Simulation number: " << sim_number << '\n';

        // read simulation parameters
        auto parameter_reader = std::make_unique<Parameter_reader>();
        parameter_reader->read_simulation(line);

        // create global environment with the simulation parameters
        auto envi = Environment(
                        std::move(parameter_reader)
                        , output_folder
                        , name_sim_file
                        , oss_time.str()
                    );

        envi.setup_output(name_sim_folder);
        envi.read_climate();

        // run simulations: each study area replicate 'i' is repeated by the
        // number of scenario replicates, where 'j' is the replicate number
        for (unsigned i = 0; i < std::size(envi.get_study_replicates()); ++i)
        {
            envi.generate_study(i);

            for (int j = 0; j < envi.get_sim_repetitions(); ++j)
            {
                envi.one_run(j);
            }

            envi.evaluate_end_sim_repetitions();
        }
        std::cout << "Completed!  Simulation number: " << sim_number << '\n';
    }
    std::cout << "------------------------------\n";
    std::cout << "All simulations are completed!\n\n";
    return EXIT_SUCCESS;
}
