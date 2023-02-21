#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <map>
#include <string>
#include <string_view>
#include <cctype>   // for ::isspace

#include "population.h"
#include "random_generator.h"
#include "weather_simulator.h"
#include "fire_simulator.h"
#include "output.h"

#include "parameter_reader.h"
#include "parameters.h"

///
/// \brief The Environment class corresponds to the global evironment of the model.
///
/// \details The main flowchart of the model occurs in this class. This class hold
/// the other classes: Grid class for the spatially-explicit processes, Output class
/// for creating output files, Parameter class, and some private member variables,
/// such as simulation time (years), fire interval (time between fires), and time since last fire
///
class Environment
{
public:
    Environment(
            std::unique_ptr<Parameter_reader> parameter_reader
            , const std::string& output_folder
            , const std::string& name_sim_file
            , const std::string& time
    );

    void read_climate()
    {
        m_weather_simulator.read_climate(m_params.file_climate);
    }

    void generate_study(const unsigned study_num);

    void setup_output(const std::string& name_sim_folder)
    {
        m_output.setup_output(m_params, name_sim_folder);
    }

    void evaluate_end_sim_repetitions()
    {
        if (m_params.is_out_plants_per_pop)
        {
            m_output.print_plants_per_pop(m_params, m_study_num);
        }
    }

    void one_run(const int run_num);


    // getters
    const std::vector<unsigned>& get_study_replicates() const
    {
        return m_params.study_replicates;
    }

    auto get_sim_repetitions() const
    {
        return m_params.sim_repetitions;
    }

private:
    // private member functions
    void restart_run();

    bool is_metapopulation_extinct();

    void read_metapopulation_file();
    void restart_metapopulation();

    void set_fire_interval();
    void set_weather_fuzzy_values();

    void evaluate_fire_event();

    void run_first_year();

    void set_postfire_dispersal_first_year()
    {
        // Expects
        assert(m_year == 0);    // this method is only used in the first year
        m_grid.set_postfire_dispersal_first_year(m_params, m_metapopulation);
    }


    void set_postfire_dispersal(const Population& population)
    {
        // Expects
        assert(m_year > 0);    // this method cannot be used in the first year
        m_grid.set_postfire_dispersal(m_params, population);
    }

    void one_year();

    void fire_event_yes(Population& population);
    void fire_event_no(Population& population);


    void postfire_seed_dispersal();

    void print_one_year();
    void print_one_run();

    void locate_patchy_fire()
    {
        m_grid.locate_patchy_fire(m_params);
    }

    void evaluate_is_population_burned()
    {
        m_grid.evaluate_is_population_burned(m_params, m_metapopulation);
    }

    void evaluate_year();

    // private member variables
    Parameters m_params;
    Output m_output{"default", "", ""};
    Weather_simulator m_weather_simulator;
    std::map<int, Population> m_metapopulation;
    Grid m_grid;
    unsigned m_study_num {0};
    int m_run_num {0};
    std::vector<double> m_weather_conditions;

    int m_time_since_fire {0};
    int m_fire_interval   {0};
    int m_year            {0};

    // membership function values (i.e. probabilities) from fuzzy function
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> m_weather_fuzzy_values;

    std::vector<int> m_occupied_patches;

    std::vector<double> m_immigrants_percent_postfire;
    std::vector<double> m_residents_percent_postfire;
    std::vector<int> m_immigrants_postfire;
    std::vector<int> m_residents_postfire;

    ///@note NOTE these variables will be tested in future studies one gene flow
//    std::vector<double> m_immigrants_percent;
//    std::vector<double> m_residents_percent;
//    std::vector<int> m_immigrants;
//    std::vector<int> m_residents;
};

#endif // ENVIRONMENT_H
