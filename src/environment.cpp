#include "environment.h"

Environment::Environment(
        std::unique_ptr<Parameter_reader> parameter_reader
        , const std::string& output_folder
        , const std::string& name_sim_file
        , const std::string& time)

    : m_params{std::move(parameter_reader)}
    , m_output{output_folder, name_sim_file, time}
{}

//-----------------------------------------------------------------------------

void Environment::one_run(const int run_num)
{
    // Expects
    assert(run_num >= 0);

    m_run_num = run_num;

    restart_run();
    run_first_year();
    print_one_year();

    if (!is_metapopulation_extinct())
    {
        while (m_year < m_params.run_time_max)
        {
            ++m_year;
            one_year();

            if (is_metapopulation_extinct())
            {
                break;
            }

            print_one_year();
        }
    }

    if (m_params.is_out_recol_rescue || m_params.is_out_metapopulation_pops)
    {
        // increment rescue effects in the case there are remaining counts
        for ([[maybe_unused]] auto& [key, val] : m_metapopulation)
        {
            val.increment_demographic_rescue();
            val.increment_recolonization_rescue();
        }
    }

    print_one_run();

    // Ensures
    assert(m_year >= 0);
    assert(m_run_num >= 0);
}

//-----------------------------------------------------------------------------

void Environment::run_first_year()
{
    // Expects
    assert(m_year == 0);

    switch (m_params.weather_scenario)
    {
    case Weather_scenario::random:
        m_weather_conditions = m_weather_simulator.get_random_year();
        break;
    case Weather_scenario::consecutive:
        m_weather_conditions = m_weather_simulator.get_consecutive_year();
        break;
    }

    set_postfire_dispersal_first_year();

    for ([[maybe_unused]] auto& [key, val] : m_metapopulation)
    {
        val.fire_mortality();
    }

    postfire_seed_dispersal();
    evaluate_year();
}

//-----------------------------------------------------------------------------

void Environment::one_year()
{
    // Expects
    assert(m_year > 0);

    switch (m_params.weather_scenario)
    {
    case Weather_scenario::random:
        m_weather_conditions = m_weather_simulator.get_random_year();
        break;
    case Weather_scenario::consecutive:
        m_weather_conditions = m_weather_simulator.get_consecutive_year();
        break;
    }

    if (m_params.model_type == Model_type::plant_based)
    {
        set_weather_fuzzy_values();
    }

    evaluate_fire_event();

    for ([[maybe_unused]] auto& [key, val] : m_metapopulation)
    {
        if (val.get_is_fire_event())
        {
            fire_event_yes(val);
        }
        else
        {
            fire_event_no(val);
        }
    }

    postfire_seed_dispersal();

    evaluate_year();
}

//-----------------------------------------------------------------------------

void Environment::fire_event_yes(Population& population)
{
    // Expects
    assert(m_params.is_fire);

    population.evaluate_rescue_effects(m_params);
    set_postfire_dispersal(population);
    population.fire_mortality();
}

//-----------------------------------------------------------------------------

void Environment::fire_event_no(Population& population)
{
    population.plant_mortality(m_params, m_weather_conditions);
    population.aging_plants_cones(m_params);
    population.density_regulation(m_params);

    if (m_params.recruit_interfire > 0.0)
    {
        population.interfire_seed_dispersal(m_params);
    }

    switch (m_params.model_type)
    {
        case Model_type::cohort_based :
            population.cone_production_cohort_based(
                        m_params,
                        m_weather_conditions);
            break;

        case Model_type::plant_based :
            population.cone_production_plant_based(
                        m_params,
                        m_weather_fuzzy_values);
            break;
    }
}

//-----------------------------------------------------------------------------

void Environment::restart_run()
{
    m_time_since_fire = 0;
    m_fire_interval   = 0;
    m_year            = 0;
    m_weather_simulator.reset_weather_year();
    m_occupied_patches.clear();

    m_immigrants_percent_postfire.clear();
    m_residents_percent_postfire.clear();
    m_immigrants_postfire.clear();
    m_residents_postfire.clear();
//    m_immigrants_percent.clear();
//    m_residents_percent.clear();
//    m_immigrants.clear();
//    m_residents.clear();

    set_fire_interval();
    restart_metapopulation();
}

//-----------------------------------------------------------------------------

void Environment::evaluate_year()
{
    double total_residents_postfire  {0.0};
    double total_immigrants_postfire {0.0};

//    double total_residents {0.0};
//    double total_immigrants {0.0};

    for ([[maybe_unused]] auto& [key, val] : m_metapopulation)
    {
        ///@todo TODO print here connectivity_dynamics
        val.evaluate_year(m_params);

        if (m_params.is_out_metapopulation_full
            || m_params.is_out_metapopulation_pops)
        {
            // postfire
            {
                double residents_postfire  {0.0};
                double immigrants_postfire {0.0};

                val.genetics_residents_immigrants_postfire(residents_postfire,
                                                           immigrants_postfire);

                total_residents_postfire += residents_postfire;
                total_immigrants_postfire += immigrants_postfire;
            }

            ///@note NOTE this will be tested for future studies on gene flow
            // genetic
            {
//                double residents {0.0};
//                double immigrants {0.0};

//                val.genetics_residents_immigrants(residents, immigrants);
//                total_residents += residents;
//                total_immigrants += immigrants;
            }
        }
    }

    if (m_params.is_out_metapopulation_full
        || m_params.is_out_metapopulation_pops)
    {
        // postfire
        {
            double immigrants_perct_postfire {0.0};

            if (total_residents_postfire + total_immigrants_postfire  > 0.0)
            {
                immigrants_perct_postfire = total_immigrants_postfire /
                        (total_residents_postfire + total_immigrants_postfire) * 100.0;
            }

            assert(immigrants_perct_postfire >= 0.0 && immigrants_perct_postfire <= 100.0);

            const double residents_perct_postfire { 100.0 - immigrants_perct_postfire };
            m_immigrants_percent_postfire.push_back(immigrants_perct_postfire);
            m_residents_percent_postfire.push_back(residents_perct_postfire);

            m_immigrants_postfire.push_back(total_immigrants_postfire);
            m_residents_postfire.push_back(total_residents_postfire);
        }

        ///@note NOTE this will be tested for future studies on gene flow
        // genetic
//        {
//            double immigrants_perct {0.0};

//            if (total_immigrants + total_residents > 0.0)
//            {
//                immigrants_perct = total_immigrants /
//                        (total_immigrants + total_residents) * 100.0;
//            }

//            assert(immigrants_perct >= 0.0 && immigrants_perct <= 100.0);

//            const double residents_perct { 100.0 - immigrants_perct };
//            m_immigrants_percent.push_back(immigrants_perct);
//            m_residents_percent.push_back(residents_perct);

//            m_immigrants.push_back(total_immigrants);
//            m_residents.push_back(total_residents);
//        }
    }
}

//-----------------------------------------------------------------------------

void Environment::postfire_seed_dispersal()
{
    if (!m_grid.is_empty_dispersal_list())
    {
        if (m_params.birds_prop > 0.0)
        {
            m_grid.LDD_by_birds(m_params, m_metapopulation);
        }

        if (m_params.wind_prop > 0.0)
        {
            m_grid.LDD_by_wind(m_params, m_metapopulation);
        }

        m_grid.establish_SDD_cohorts(m_metapopulation);

        m_grid.remove_dispersal_list();
    }
}

//-----------------------------------------------------------------------------

void Environment::set_weather_fuzzy_values()
{
    m_weather_fuzzy_values.clear();

    for (const auto& [key, val] : m_params.fuzzy_set_list)
    {
        std::vector<double> tmp_values;
        for (const auto& mf : val)
        {
            // calculation of the fuzzy membership value computed using
            // triangle or trapezoid membership function
            const double value {
                mf->get_value(m_weather_conditions.at(mf->get_climate_index()))
            };

            tmp_values.push_back(value);
        }
        m_weather_fuzzy_values.emplace(key, tmp_values);
    }
}

//-----------------------------------------------------------------------------

bool Environment::is_metapopulation_extinct()
{
    unsigned total_extinct_patches {0};

    for ([[maybe_unused]] const auto& [key, val] : m_metapopulation)
    {
        total_extinct_patches += val.get_is_extinct();
    }

    if (m_params.is_out_metapopulation || m_params.is_out_metapopulation_full)
    {
        assert(total_extinct_patches <= std::size(m_metapopulation));

        m_occupied_patches.push_back(
                    std::size(m_metapopulation) - total_extinct_patches);
    }

    //terminate the run simulation if all populations are extinct
    if (total_extinct_patches == std::size(m_metapopulation))
    {
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------

void Environment::set_fire_interval()
{
    if (!m_params.is_fire)
    {
        return;
    }

    switch (m_params.fire_interval_scn)
    {
        case Fire_scenario::deterministic :
            m_fire_interval = m_params.fire_interval_mean;
            break;

        case Fire_scenario::truncated_normal :
            m_fire_interval = Fire_simulator::get_fire_interval_truncated_normal(
                        m_params.fire_interval_mean
                        , m_params.fire_interval_lower_cut);
            break;
        case Fire_scenario::weibull :
            m_fire_interval = Fire_simulator::get_fire_interval_weibull(
                        m_params.fire_a
                        , m_params.fire_b
                        , m_params.fire_interval_lower_cut);
            break;
    }
}

//-----------------------------------------------------------------------------

void Environment::evaluate_fire_event()
{
    ++m_time_since_fire;

    for ([[maybe_unused]] auto& [key, val] : m_metapopulation)
    {
        val.increment_time_since_fire();
    }

    if (m_params.is_fire)
    {
        if (m_time_since_fire >= m_fire_interval)
        {
            switch (m_params.fire_scale)
            {
                case Fire_scale::patchy :
                    locate_patchy_fire();
                    break;

                case Fire_scale::study_area :
                    // nothing happens
                    break;
            }
            evaluate_is_population_burned();
            m_time_since_fire = 0;
            set_fire_interval();

            if (m_params.is_out_fire)
            {
                m_output.print_fire_output(m_params
                                           , m_study_num
                                           , m_run_num
                                           , m_year
                                           , m_metapopulation);
            }
        }
    }
}

//-----------------------------------------------------------------------------

void Environment::read_metapopulation_file()
{
    m_metapopulation.clear();

    const std::string file_path {
        Global::project_directory
                + "data/in/metapopulation/"
                + m_params.file_metapopulation
    };

    std::ifstream ifs_meta(file_path);

    if (!ifs_meta)
    {
        std::cerr << "ERROR: The following file could not open.\n\a"
                  << file_path
                  << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    std::string data;
    std::getline(ifs_meta, data);  // skip header

    while (std::getline(ifs_meta, data))
    {
        // skip if line is empty
        if (std::all_of(std::cbegin(data), std::cend(data), ::isspace))
        {
            continue;
        }

        ///@note NOTE Unused vars: cohort_age, init_individuals, and init_seedbank
        int pop_id            {0};
        double patch_size_ha  {0.0};
        int age               {0};
        int init_individuals  {0};
        int init_seedbank     {0};
        std::string habitat_quality_input;

        std::stringstream ss(data);
        ss  >> pop_id
            >> patch_size_ha
            >> age
            >> init_individuals
            >> init_seedbank
            >> habitat_quality_input;

        assert(m_params.cell_size > 0);

        int num_cells {0};

        if (m_params.cell_size > 0)
        {
            num_cells = static_cast<int>(
                            (patch_size_ha * (100 / m_params.cell_size)) + 0.5);
        }

        if (num_cells <= 0)
        {
            std::cerr << "ERROR: population with id:" << pop_id
                      << " has zero cells (size too small)\n"
                      << file_path << "\n\a";
            std::exit(EXIT_FAILURE);
        }

        if (pop_id <= 0)
        {
            std::cerr << "ERROR: population id must be above zero\n"
                      << file_path << "\n\a";
            std::exit(EXIT_FAILURE);
        }

        m_metapopulation.emplace(
                    pop_id
                    , Population(
                        pop_id
                        , num_cells
                        , m_params.carrying_capacity * num_cells
                        , habitat_quality_input));
    }
}

//-----------------------------------------------------------------------------

void Environment::restart_metapopulation()
{
    assert(!std::empty(m_metapopulation));

    const std::string file_path {
        Global::project_directory
                + "data/in/metapopulation/"
                + m_params.file_metapopulation
    };

    constexpr bool is_postfire   {true};
    constexpr bool is_LDD_cohort {false};

    for (auto& [key, val]: m_metapopulation)
    {
        std::ifstream ifs(file_path);

        if (!ifs)
        {
            std::cerr << "ERROR: The following file could not open.\n\a"
                      << file_path << "\n\a";
            std::exit(EXIT_FAILURE);
        }

        std::string data;
        std::getline(ifs, data);  // skip header

        int pop_id            {0};
        double patch_size_ha  {0.0};
        int age               {0};
        int init_individuals  {0};
        int init_seedbank     {0};
        std::string habitat_quality_input;

        //------------------------
        bool is_population_found {false};

        ///@warning WARNING it is not evaluated if metapop files from read_metapopulation
        /// has the same number of populations (and ID number)
        while (std::getline(ifs, data))
        {
            std::stringstream ss(data);
            ss  >> pop_id
                >> patch_size_ha
                >> age
                >> init_individuals
                >> init_seedbank
                >> habitat_quality_input;

            if (pop_id == key)
            {
                is_population_found = true;
                break;
            }
        }

        if (!is_population_found)
        {
            std::cerr << "ERROR: Population id not found in Metapopulation file."
                      << "\n Population id not found: " << key
                      << "\n File: " << file_path << "\n\a";
            std::exit(EXIT_FAILURE);
        }

        val.remove_all_cohorts();

        ///@note NOTE init_seedbank is unused and age unimportant
        val.initialize_population(
                    habitat_quality_input
                    , Cohort(age
                             , init_individuals
                             , is_postfire
                             , Gene_flow(key)
                             , is_LDD_cohort, Dispersal_vector::SDD_wind));
        ifs.close();
    }
}

//-----------------------------------------------------------------------------

void Environment::generate_study(const unsigned study_num)
{
    m_study_num = study_num;
    m_grid.set_study_area(m_params);
    read_metapopulation_file();

    assert(!std::empty(m_metapopulation));

    if (m_params.is_out_plants_per_pop)
    {
        for ([[maybe_unused]] auto& [key, val] : m_metapopulation)
        {
            val.init_plants_per_population(m_metapopulation);
        }
        m_output.init_plants_per_population(m_metapopulation);
    }

    m_grid.locate_metapopulation(m_params.study_replicates.at(study_num), m_metapopulation);
    m_grid.set_distance_btw_populations(m_params, m_metapopulation);
    m_grid.set_cells_for_LDD_by_birds(m_params, m_metapopulation);

    if (m_params.is_out_study_area)
    {
        m_output.print_study_area(
                    m_grid.get_study_area()
                    , m_params
                    , m_study_num
        );
    }
}

//-----------------------------------------------------------------------------

void Environment::print_one_year()
{
    if (m_params.is_out_metapopulation_dynamics)
    {
        m_output.print_metapopulation_dynamics(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_year
                    , m_time_since_fire
                    , m_metapopulation
        );
    }

    if (m_params.is_out_seed_dynamics)
    {
        m_output.print_seed_dynamics(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_year
                    , m_time_since_fire
                    , m_metapopulation
        );
    }
}

//-----------------------------------------------------------------------------

void Environment::print_one_run()
{
    if (m_params.is_out_metapopulation)
    {
        m_output.print_metapopulation(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_year
                    , m_metapopulation
                    , m_occupied_patches);
    }

    if (m_params.is_out_metapopulation_full)
    {
        m_output.print_metapopulation_full(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_year
                    , m_metapopulation
                    , m_occupied_patches
                    , m_immigrants_percent_postfire
                    , m_residents_percent_postfire
                    , m_immigrants_postfire
                    , m_residents_postfire);
    }

    if (m_params.is_out_metapopulation_pops)
    {
        m_output.print_metapopulation_pops(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_year
                    , m_metapopulation);
    }

    if (m_params.is_out_recol_rescue)
    {
        m_output.print_rescue_effects(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_year
                    , m_metapopulation);
    }

    if (m_params.is_out_plants_per_pop)
    {
        m_output.add_plants_per_population(m_metapopulation);
        ///@note NOTE print is at the end of all sim repetitions in
        /// Environment::evaluate_end_sim_repetitions()
    }

    if (m_params.is_out_dispersal)
    {
        m_output.print_dispersal(
                    m_params
                    , m_study_num
                    , m_run_num
                    , m_metapopulation);
    }
}
