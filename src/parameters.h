#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <memory>

#include "model_options.h"
#include "plant_type.h"
#include "fuzzy_function.h"
#include "flower_distribution.h"
#include "habitat_quality.h"

#include "parameter_reader.h"

///
/// \brief The Parameters struct holds all parameters used during the simulation
/// \note all parameters are const and the validity was checked in
/// Parameter_reader Class. Detail description of the parameters is in Parameter_reader
///
struct Parameters
{
    Parameters(std::unique_ptr<Parameter_reader> p)
        : sim_id{ p->get_sim_id() }
        , study_replicates{ p->get_study_replicates() }
        , sim_repetitions{ p->get_sim_repetitions() }
        , run_time_max{ p->get_run_time_max() }
        , is_out_study_area{ p->get_is_out_study_area() }
        , is_out_fire{ p->get_is_out_fire() }
        , is_out_seed_dynamics{ p->get_is_out_seed_dynamics() }
        , is_out_metapopulation{ p->get_is_out_metapopulation() }
        , is_out_metapopulation_full{ p->get_is_out_metapopulation_full() }
        , is_out_metapopulation_dynamics{ p->get_is_out_metapopulation_dynamics() }
        , is_out_metapopulation_pops{ p->get_is_out_metapopulation_pops() }
        , is_out_recol_rescue{ p->get_is_out_recol_rescue() }
        , is_out_plants_per_pop{ p->get_is_out_plants_per_pop() }
        , is_out_dispersal{ p->get_is_out_dispersal() }
        , model_type{ p->get_model_type() }
        , cell_size{ p->get_cell_size() }
        , study_size_x{ p->get_study_size_x() }
        , study_size_y{ p->get_study_size_y() }
        , file_metapopulation{ p->get_file_metapop() }
        , carrying_capacity{ p->get_carrying_capacity() }
        , init_cones_plant{ p->get_init_cones_plant() }
        , init_viable_seeds_plant{ p->get_init_viable_seeds_plant() }
        , is_climate{ p->get_is_climate() }
        , long_term_rain{ p->get_long_term_rain() }
        , weather_scenario{ p->get_weather_scenario() }
        , file_climate{ p->get_file_climate() }
        , is_fire{ p->get_is_fire() }
        , fire_interval_scn { p->get_fire_interval_scn() }
        , fire_size_scn { p->get_fire_size_scn() }
        , fire_scale{ p->get_fire_scale() }
        , fire_size_x{ p->get_fire_size_x() }
        , fire_size_y{ p->get_fire_size_y() }
        , fire_interval_mean{ p->get_fire_interval_mean() }
        , fire_interval_lower_cut{ p->get_fire_interval_lower_cut() }
        , burned_lower_cut{ p->get_burned_lower_cut() }
        , burned_upper_cut{ p->get_burned_upper_cut() }
        , fire_a{ p->get_fire_a() }
        , fire_b{ p->get_fire_b() }
        , wind_prop{ p->get_wind_prop() }
        , wind_a{ p->get_wind_a() }
        , wind_b{ p->get_wind_b() }
        , wind_direction{ p->get_wind_direction() }
        , birds_dist_max{ p->get_birds_dist_max() }
        , birds_prop{ p->get_birds_prop()}
        , follicles_distr_type{ p->get_follicles_distr_type() }
        , follicles_distr_a{ p->get_follicles_distr_a() }
        , follicles_distr_b{ p->get_follicles_distr_b() }
        , postfire_follicles_open{ p->get_postfire_follicles_open() }
        , plant_longevity{ p->get_plant_longevity() }
        , plant_age_young{ p->get_plant_age_young() }
        , plant_age_adult{ p->get_plant_age_adult() }
        , mort_scenario{ p->get_mort_scenario() }
        , mort_recruit_postfire_min{ p->get_mort_recruit_postfire_min() }
        , mort_recruit_postfire_max{ p->get_mort_recruit_postfire_max() }
        , mort_recruit_postfire_mean{ p->get_mort_recruit_postfire_mean() }
        , recruit_interfire{ p->get_recruit_interfire() }
        , recruit_weather{ p->get_recruit_weather() }
        , senescence_age{ p->get_senescence_age() }
        , senescence_increase{ p->get_senescence_increase() }
        , mort_min{ p->get_mort_min() }
        , mort_a{ p->get_mort_a() }
        , mort_b{ p->get_mort_b() }
        , mort_c{ p->get_mort_c() }
        , mort_d{ p->get_mort_d() }
        , mort_e{ p->get_mort_e() }
        , mort_f{ p->get_mort_f() }
        , cone_cycle{ p->get_cone_cycle() }
        , seed_longevity{ p->get_seed_longevity() }
        , flower_scenario{ p->get_flower_scenario() }
        , flower_age_a{ p->get_flower_age_a() }
        , flower_age_b{ p->get_flower_age_b() }
        , flower_age_c{ p->get_flower_age_c() }
        , flower_weather_a{ p->get_flower_weather_a() }
        , flower_weather_b{ p->get_flower_weather_b() }
        , flower_weather_c{ p->get_flower_weather_c() }
        , flower_weather_d{ p->get_flower_weather_d() }
        , flower_weather_e{ p->get_flower_weather_e() }
        , pollination_success{ p->get_pollination_success() }
        , num_follicles{ p->get_num_follicles() }
        , num_seeds{ p->get_num_seeds() }
        , firm_seeds{ p->get_firm_seeds() }
        , viable_seeds{ p->get_viable_seeds() }
        , viable_seeds_cone_age{ p->get_viable_seeds_cone_age() }
        , viable_seeds_dispersed_cone_age{ p->get_viable_seeds_dispersed_cone_age() }
        , insect_cone_age{ p->get_insect_cone_age() }
        , decay_cone_age{ p->get_decay_cone_age() }
        , open_follicles_cone_age{ p-> get_open_follicles_cone_age() }
        , habitat_quality{ p->get_habitat_quality() }
        , plant_type_list{ p->get_plant_type_list() }
        , plant_type_share{ p->get_plant_type_share() }
        , fuzzy_set_list{ p->get_fuzzy_set_list() }
        , fuzzy_set_weather_var{ p->get_fuzzy_set_weather_var() }
        , flower_distribution{ p->get_flower_distribution() }
    {}

    // parameter description is in Parameter_reader Class
    const int sim_id;
    const std::vector<unsigned> study_replicates;
    const int sim_repetitions;
    const int run_time_max;
    const bool is_out_study_area;
    const bool is_out_fire;
    const bool is_out_seed_dynamics;
    const bool is_out_metapopulation;
    const bool is_out_metapopulation_full;
    const bool is_out_metapopulation_dynamics;
    const bool is_out_metapopulation_pops;
    const bool is_out_recol_rescue;
    const bool is_out_plants_per_pop;
    const bool is_out_dispersal;
    const Model_type model_type;
    const int cell_size;
    const int study_size_x;
    const int study_size_y;
    const std::string file_metapopulation;
    const int carrying_capacity;
    const double init_cones_plant;
    const double init_viable_seeds_plant;
    const bool is_climate;
    const double long_term_rain;
    const Weather_scenario weather_scenario;
    const std::string file_climate;
    const bool is_fire;
    const Fire_scenario fire_interval_scn;
    const Fire_scenario fire_size_scn;
    const Fire_scale fire_scale;
    const double fire_size_x;
    const double fire_size_y;
    const int fire_interval_mean;
    const int fire_interval_lower_cut;
    const int burned_lower_cut;
    const int burned_upper_cut;
    const double fire_a;
    const double fire_b;
    const double wind_prop;
    const double wind_a;
    const double wind_b;
    const Wind_direction wind_direction;
    const double birds_dist_max;
    const double birds_prop;
    const Distribution_type follicles_distr_type;
    const double follicles_distr_a;
    const double follicles_distr_b;
    const double postfire_follicles_open;
    const int plant_longevity;
    const int plant_age_young;
    const int plant_age_adult;
    const Mortality_scenario mort_scenario;
    const double mort_recruit_postfire_min;
    const double mort_recruit_postfire_max;
    const double mort_recruit_postfire_mean;
    const double recruit_interfire;
    const double recruit_weather;
    const int senescence_age;
    const double senescence_increase;
    const double mort_min;
    const double mort_a;
    const double mort_b;
    const double mort_c;
    const double mort_d;
    const double mort_e;
    const double mort_f;
    const unsigned cone_cycle;
    const unsigned seed_longevity;
    const Flower_scenario flower_scenario;
    const double flower_age_a;
    const double flower_age_b;
    const double flower_age_c;
    const double flower_weather_a;
    const double flower_weather_b;
    const double flower_weather_c;
    const double flower_weather_d;
    const double flower_weather_e;
    const double pollination_success;
    const double num_follicles;
    const double num_seeds;
    const double firm_seeds;
    const double viable_seeds;
    const std::vector<double> viable_seeds_cone_age;
    const std::vector<double> viable_seeds_dispersed_cone_age;
    const std::vector<double> insect_cone_age;
    const std::vector<double> decay_cone_age;
    const std::vector<double> open_follicles_cone_age;
    const std::map<std::string, Habitat_quality> habitat_quality;
    const std::map<unsigned, Plant_type> plant_type_list;
    const std::map<unsigned, double> plant_type_share;
    const std::map<std::pair<unsigned, unsigned>, std::vector<std::shared_ptr<Fuzzy_function>>> fuzzy_set_list;
    const std::map<std::pair<unsigned, unsigned>, unsigned> fuzzy_set_weather_var;
    const std::map<std::tuple<unsigned, unsigned, unsigned>, Flower_distribution> flower_distribution;
};

#endif // PARAMETERS_H
