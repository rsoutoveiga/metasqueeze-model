#ifndef PARAMETER_READER_H
#define PARAMETER_READER_H

///@note NOTE read the F.7 about not passing smart pointers
/// https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#f7-for-general-use-take-t-or-t-arguments-rather-than-smart-pointers


#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <sstream>
#include <fstream>
#include <memory>
#include <map>
#include <tuple>
#include <cctype>   // for ::isspace
#include <cstdint>  // for int64_t
#include <cmath>    // std::fabs()
#include <limits>   // epsilon()

#include "model_options.h"
#include "plant_type.h"
#include "fuzzy_function.h"
#include "trapezoid.h"
#include "triangle.h"
#include "flower_distribution.h"
#include "habitat_quality.h"
#include "equation.h"
#include "global.h"

///
/// \brief The Parameter_reader class reads parameters from files, and checks
/// the validity of the input parameters values, and calculates some
/// derived paramters (e.g., available viable seeds per cone age)
///
class Parameter_reader
{
public:
    Parameter_reader() = default;

    void read_simulation(const std::string& line);

    // getters
    int get_sim_id() const { return m_sim_id; }
    bool get_is_study_replicates() const { return m_is_study_replicates; }
    const std::vector<unsigned>& get_study_replicates() const { return m_study_replicates; }
    int get_sim_repetitions() const { return m_sim_repetitions; }
    int get_run_time_max() const { return m_run_time_max; }
    bool get_is_out_study_area() const { return m_is_out_study_area; }
    bool get_is_out_fire() const { return m_is_out_fire; }
    bool get_is_out_seed_dynamics() const { return m_is_out_seed_dynamics; }
    bool get_is_out_metapopulation() const { return m_is_out_metapopulation; }
    bool get_is_out_metapopulation_full() const { return m_is_out_metapopulation_full; }
    bool get_is_out_metapopulation_dynamics() const { return m_is_out_metapopulation_dynamics; }
    bool get_is_out_metapopulation_pops() const { return m_is_out_metapopulation_pops; }
    bool get_is_out_recol_rescue() const { return m_is_out_recol_rescue; }
    bool get_is_out_plants_per_pop() const { return m_is_out_plants_per_pop; }
    bool get_is_out_dispersal() const { return m_is_out_dispersal; }

    Model_type get_model_type() const { return m_model_type; }
    int get_cell_size() const { return m_cell_size; }
    int get_study_size_x() const { return m_study_size_x; }
    int get_study_size_y() const { return m_study_size_y; }

    const std::string& get_file_metapop() const { return m_file_metapopulation; }
    int get_carrying_capacity() const { return m_carrying_capacity; }
    double get_init_cones_plant() const { return m_init_cones_plant;}
    double get_init_viable_seeds_plant() const { return m_init_viable_seeds_plant; }
    bool get_is_climate() const { return m_is_climate; }
    double get_long_term_rain() const { return m_long_term_rain; }

    Weather_scenario get_weather_scenario() const { return m_weather_scenario; }
    const std::string& get_file_climate() const { return m_file_climate; }

    bool get_is_fire() const { return m_is_fire; }
    Fire_scenario get_fire_interval_scn() const { return m_fire_interval_scn; }
    Fire_scenario get_fire_size_scn() const { return m_fire_size_scn; }
    Fire_scale get_fire_scale() const { return m_fire_scale; }

    double get_fire_size_x() const { return m_fire_size_x; }
    double get_fire_size_y() const { return m_fire_size_y; }

    int get_fire_interval_mean() const { return m_fire_interval_mean; }
    int get_fire_interval_lower_cut() const { return m_fire_interval_lower_cut; }
    int get_burned_lower_cut() const { return m_burned_lower_cut; }
    int get_burned_upper_cut() const { return m_burned_upper_cut; }
    double get_fire_a() const { return m_fire_a; }
    double get_fire_b() const { return m_fire_b; }

    // dispersal file
    double get_wind_prop() const { return m_wind_prop; }
    double get_wind_a() const { return m_wind_a; }
    double get_wind_b() const { return m_wind_b; }
    Wind_direction get_wind_direction() const { return m_wind_direction; }
    double get_birds_dist_max() const { return m_birds_dist_max; }
    double get_birds_prop() const { return m_birds_prop; }
    Distribution_type get_follicles_distr_type() const { return m_follicles_distr_type; }
    double get_follicles_distr_a() { return m_follicles_distr_a; }
    double get_follicles_distr_b() { return m_follicles_distr_b; }
    double get_postfire_follicles_open() const { return m_postfire_follicles_open; }

    // species file
    int get_plant_longevity() const { return m_plant_longevity; }
    int get_plant_age_young() const { return m_plant_age_young; }
    int get_plant_age_adult() const { return m_plant_age_adult;}

    Mortality_scenario get_mort_scenario() const { return m_mort_scenario; }

    double get_mort_recruit_postfire_min() const { return m_mort_recruit_postfire_min; }
    double get_mort_recruit_postfire_max() const { return m_mort_recruit_postfire_max; }
    double get_mort_recruit_postfire_mean() const { return m_mort_recruit_postfire_mean; }
    double get_recruit_interfire() const { return m_recruit_interfire; }
    double get_recruit_weather() const { return m_recruit_weather; }
    int get_senescence_age() const { return m_senescence_age; }
    double get_senescence_increase() const { return m_senescence_increase; }
    double get_mort_min() const { return m_mort_min; }
    double get_mort_a() const { return m_mort_a; }
    double get_mort_b() const { return m_mort_b; }
    double get_mort_c() const { return m_mort_c; }
    double get_mort_d() const { return m_mort_d; }
    double get_mort_e() const { return m_mort_e; }
    double get_mort_f() const { return m_mort_f; }

    unsigned get_cone_cycle() const { return m_cone_cycle; }
    unsigned get_seed_longevity() const { return m_seed_longevity; }
    Flower_scenario get_flower_scenario() const { return m_flower_scenario; }
    double get_flower_age_a() const { return m_flower_age_a;}
    double get_flower_age_b() const { return m_flower_age_b;}
    double get_flower_age_c() const { return m_flower_age_c;}
    double get_flower_weather_a() const { return m_flower_weather_a; }
    double get_flower_weather_b() const { return m_flower_weather_b;}
    double get_flower_weather_c() const { return m_flower_weather_c;}
    double get_flower_weather_d() const { return m_flower_weather_d;}
    double get_flower_weather_e() const { return m_flower_weather_e;}
    double get_pollination_success() const { return m_pollination_success; }
    double get_num_follicles() const { return m_num_follicles; }
    double get_num_seeds() const { return m_num_seeds; }
    double get_firm_seeds() const { return m_firm_seeds; }
    double get_viable_seeds() const { return m_viable_seeds; }

    const std::vector<double>& get_viable_seeds_cone_age() const {
        return m_viable_seeds_cone_age;
    }
    const std::vector<double>& get_viable_seeds_dispersed_cone_age() const {
        return m_viable_seeds_dispersed_cone_age;
    }

    const std::vector<double>& get_insect_cone_age() const {
        return m_insect_cone_age;
    }

    const std::vector<double>& get_decay_cone_age() const {
        return m_decay_cone_age;
    }

    const std::vector<double>& get_open_follicles_cone_age() const {
        return m_open_follicles_cone_age;
    }

    const std::map<unsigned, Plant_type>& get_plant_type_list() const {
        return m_plant_type_list;
    }

    const std::map<unsigned, double>& get_plant_type_share() const {
        return m_plant_type_share;
    }

    const std::map<std::pair<unsigned, unsigned>, std::vector<std::shared_ptr<Fuzzy_function>>>& get_fuzzy_set_list() const
    {
        return m_fuzzy_set_list;
    }

    const std::map<std::pair<unsigned, unsigned>, unsigned>& get_fuzzy_set_weather_var() const {
        return m_fuzzy_set_weather_var;
    }

    const std::map<std::tuple<unsigned, unsigned, unsigned>, Flower_distribution>& get_flower_distribution() const
    {
        return m_flower_dist;
    }

    const std::map<std::string, Habitat_quality>& get_habitat_quality() const
    {
        return m_habitat_quality;
    }



private:

    void set_study_replicates(const std::string& total_replicates);
    void read_study_replicates(const std::string& name_file);

    // this is the old derivated parameters function
    void read_scenario(const std::string& file_scenario);

    void read_study_size(const std::string& file_study_size);
    void read_fire(const std::string& file_fire);
    void read_dispersal(const std::string& file_dispersal);

    void read_species(const std::string& file_species);
    void calculate_canopy_params();
    void read_plant_type(const std::string& file_plant_type);
    void read_fuzzy_set_list(const std::string& file_fuzzy_set_list);
    void read_flower_distribution(const std::string& file_flower_dist);
    void read_habitat_quality(const std::string& file_habitat);


    // read individual parameters
    void read_changed_parameters(
            const std::vector<std::pair<std::string, std::string>>& parameters);

    // check if first line is empty, i.e. file is empty
    void check_file_is_empty(std::string_view folder, std::string_view file, std::string_view line)
    {
        if (std::all_of(std::cbegin(line), std::cend(line), ::isspace))
        {
            std::cerr << "ERROR! The following file is empty:\n"
                      << folder << file << "\n\a";
            std::exit(EXIT_FAILURE);
        }
    }

    // simulation file
    int m_sim_id                              {1};
    bool m_is_study_replicates                {true};
    std::vector<unsigned> m_study_replicates;  // seeds for random_generator to create study area
    int m_sim_repetitions                     {5};
    int m_run_time_max                        {500};  // [years]
    bool m_is_out_study_area                  {true};
    bool m_is_out_fire                        {true};
    bool m_is_out_seed_dynamics               {false};
    bool m_is_out_metapopulation              {true};
    bool m_is_out_metapopulation_full         {true};
    bool m_is_out_metapopulation_dynamics     {true};
    bool m_is_out_metapopulation_pops         {true};
    bool m_is_out_recol_rescue                {true};
    bool m_is_out_plants_per_pop              {true};
    bool m_is_out_dispersal                   {false};

//------------------------
    // scenario file
    Model_type m_model_type { Model_type::cohort_based };

    // study_size.txt
    std::string m_file_study_size;    // it has no getter
    int m_cell_size    {100}; // cell size length [m]
    int m_study_size_x {3000}; // [m] (He et al., 2010)
    int m_study_size_y {5000}; // [m] (He et al., 2010)

    std::string m_file_metapopulation {"metapopulation_default.txt"};

    int m_carrying_capacity             {2500};    // maximum number of adult plants in one hectare [#]
    double m_init_cones_plant           {80.0};    // initial fertile cones per plant individual [#]
    double m_init_viable_seeds_plant    {640.0};   // initial viable seeds per plant individual [#]
    bool m_is_climate                   {true};    // is climate activated?
    double m_long_term_rain             {453.848}; // long-term rainfall [mm]
    Weather_scenario m_weather_scenario {Weather_scenario::random};
    std::string m_file_climate          {"baseline.csv"};

    bool m_is_fire {true};  // true = fire occurs, false = fire does not occur

    std::string m_file_fire;              // it has no getter
    std::string m_file_dispersal;         // it has no getter
    std::string m_file_species;           // it has no getter
    std::string m_file_plant_type;        // it has no getter
    std::string m_file_fuzzy_set_list;    // it has no getter
    std::string m_file_flower_dist;       // it has no getter
    std::string m_file_habitat;           // it has no getter

    // fire file
    Fire_scenario m_fire_interval_scn {Fire_scenario::deterministic};
    Fire_scenario m_fire_size_scn     {Fire_scenario::deterministic};
    Fire_scale m_fire_scale           {Fire_scale::patchy};
    double m_fire_size_x              {2121.0};  // length of the semi x-axis of the fire ellipse [m].
    double m_fire_size_y              {2828.0};  // length of the semi y-axis of the fire ellipse [m].
    int m_fire_interval_mean          {17};      // mean of fire frequency [years]
    int m_fire_interval_lower_cut     {1};       // minimum time between fires [years]
    int m_burned_lower_cut            {3};       // minimum fire interval needed to be burned, i.e., minimum fuel load [years] (Groeneveld et al., 2008)
    int m_burned_upper_cut            {12};      // Fire interval in which the probability to be burned is 100% [years] (Groeneveld et al., 2008)
    double m_fire_a                   {11.7};    // NOT USED! scale param (Weibull distribution) for fire interval (Groeneveld 2008)
    double m_fire_b                   {2.1};     // NOT USED! shape param (Weibull distribution) for fire interval (Groeneveld 2008)

    // dispersal file
    double m_wind_prop {0.15};    // proportion of LDD of seeds by post-fire wind
    double m_wind_a    {6.787};   // logmean parameter of the log-normal probability density function, i.e., dispersal kernel by post-fire wind
    double m_wind_b    {0.676};   // log standard deviation parameter of the log-normal probability density function, i.e., dispersal kernel by post-fire wind

    Wind_direction m_wind_direction {Wind_direction::random};

    double m_birds_dist_max  {1250.0};  // maximum distance cone dispersed by birds
    double m_birds_prop {0.07};         // proportion of cones dispersed by birds
    Distribution_type m_follicles_distr_type {Distribution_type::nbinom_mu};  // probability density function type of number of follicles per fertile cone
    double m_follicles_distr_a {6.219};  // parameter 'a' for the probability density function of number of follicles per fertile cone.
    double m_follicles_distr_b {10.08};  // parameter 'b' for the probability density function of number of follicles per fertile cone.
    double m_postfire_follicles_open {0.5};  // proportion of follicles opened two hours after a fire

    // species file
    int m_plant_longevity {999}; // maximum plant age
    int m_plant_age_young {5};   // plant age where cones start to be fertile
    int m_plant_age_adult {15};  // plants reached maximum cone production and survival

    Mortality_scenario m_mort_scenario {Mortality_scenario::age_weather_relative};

    double m_mort_recruit_postfire_min  {0.95};   // minimum mortality postfire recruits
    double m_mort_recruit_postfire_max  {0.987};  // maximum mortality postfire recruits
    double m_mort_recruit_postfire_mean {0.969};  // survival of plants at the age of zero
    double m_recruit_interfire     {5};      // percentage of inter-fire recruitment
    double m_recruit_weather       {0.06};   // weather impacts on seedling recruitment

    int m_senescence_age         {25};  // age of increased mortality [years] (Enright et al., 1998)
    double m_senescence_increase {0.01};  // annual increase in mortality Probability due to senescence (Enright et al., 1998)

    double m_mort_min {0.02};  // minimum plant mortality probability (plant age ≥ 1 year).

    double m_mort_a {0.0553}; // parameter 'a' to calculate plant mortality (plant age ≥ 1 year)
    double m_mort_b {0.2645}; // parameter 'b' to calculate plant mortality (plant age ≥ 1 year)
    double m_mort_c {0.179};  // parameter 'c' to calculate plant mortality (plant age ≥ 1 year)

    double m_mort_d {-0.00061}; // parameter 'd' to calculate plant mortality (plant age ≥ 1 year)
    double m_mort_e {0.26873};  // parameter 'e' to calculate plant mortality (plant age ≥ 1 year)
    double m_mort_f {0.0};      // parameter 'f' to calculate plant mortality (plant age ≥ 1 year)

    unsigned m_cone_cycle     {1}; // time required for cones and seeds to reach maturity after successful pollination of the flowers [years].
    unsigned m_seed_longevity {12}; // maximum seed longevity [years] (Enright et al., 1996).

    Flower_scenario m_flower_scenario {Flower_scenario::baseline};

    double m_flower_age_a {9.179382}; // parameter 'a' to calculate the mean flower production per plant based on plant age
    double m_flower_age_b {8.710231}; // parameter 'b' to calculate the mean flower production per plant based on plant age
    double m_flower_age_c {0.6208263}; // parameter 'c' to calculate the mean flower production per plant based on plant age

    double m_flower_weather_a {0.053176};
    double m_flower_weather_b {-15.21237};
    double m_flower_weather_c {0.022191};
    double m_flower_weather_d {0.023055};
    double m_flower_weather_e {-33.344613};

    double m_pollination_success {90.0}; // [%]
    double m_num_follicles       {9.97};
    double m_num_seeds           {2.0};
    double m_firm_seeds          {0.83};
    double m_viable_seeds        {0.744};

    double m_insect_a {0.02};
    double m_insect_b {0.18};

    double m_decay_a {0.34};
    double m_decay_b {-5.95};

    double m_open_follicles_a {0.3896};
    double m_open_follicles_b {6.6873};
    double m_open_follicles_c {0.9588};

    // ----------------------------------------------------
    // derived parameters:
    // seeds per cone age after seed loss and dispersal
    // seeds dispersed per cone age
    std::vector<double> m_viable_seeds_cone_age;
    std::vector<double> m_viable_seeds_dispersed_cone_age;

    std::vector<double> m_insect_cone_age;
    std::vector<double> m_decay_cone_age;
    std::vector<double> m_open_follicles_cone_age;
    //-----------------------------------------------------

    // table of habitat quality values after reading habitat file
    // key = habitat label; value = Habitat_quality (label + habitat value)
    std::map<std::string, Habitat_quality> m_habitat_quality;

    //-----------------------------------------------------
    // parameters only for the plant-based model type

    // parameters after reading plant_type file
    // key = plant type; value =
    std::map<unsigned, Plant_type> m_plant_type_list;
    std::map<unsigned, double> m_plant_type_share;

    // parameters after reading fuzzy set file (model_type must be plant-based)
    std::map<std::pair<unsigned, unsigned>, std::vector<std::shared_ptr<Fuzzy_function>>> m_fuzzy_set_list;

    // pair <level, fuzzy set>, climate index column
    std::map<std::pair<unsigned, unsigned>, unsigned> m_fuzzy_set_weather_var;

    // parameter after reading flower distribution file
    std::map<std::tuple<unsigned, unsigned, unsigned>, Flower_distribution> m_flower_dist;
};


#endif // PARAMETER_READER_H
