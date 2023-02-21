#ifndef COHORT_H
#define COHORT_H

#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <functional>
#include <memory>
#include <cassert>
#include <deque>
#include <map>


#include "cell.h"
#include "equation.h"
#include "parameters.h"
#include "gene_flow.h"
#include "plant.h"
#include "global.h"


///
/// \brief The Cohort class represents a group of plants with the same age.
/// Cohort class belongs to Population class.
///
class Cohort
{
public:

    Cohort(const int age
           , const int num_plants
           , const bool is_postfire
           , const Gene_flow& gene_flow
           , const bool is_LDD_cohort
           , const Dispersal_vector dispersal_vector
    );

    const std::vector<double>& getting_seedbank_vector() const {
        return m_seedbank_cones;
    }

    void init_mature_cohort(const Parameters& params);

    void plant_mortality(const Parameters& params
                         , const std::vector<double>& weather
                         , const int carrying_capacity);

    void aging_plants_cones(const Parameters& params);

    int interfire_seed_dispersal(const Parameters& params);

    // cone production cohort-based and individual-based functions
    void cone_production_cohort_based(const Parameters& params
                                      , const std::vector<double>& weather
                                      , const Habitat_quality& habitat_quality);

    void cone_production_plant_based(
            const Parameters& params
            , const std::map<std::pair<unsigned, unsigned>, std::vector<double>>& weather_fuzzy_values
            , const Habitat_quality& habitat_quality);

    void add_plant(const Parameters& params, const unsigned key)
    {
        m_plants_list.at(key).emplace_back(Plant(params));
    }

    // getters
    int get_age() const { return m_age; }
    int get_num_plants() const { return m_num_plants; }
    bool get_is_postfire() const { return m_is_postfire; }
    const Gene_flow& get_gene_flow() const { return m_gene_flow; }
    bool get_is_LDD_cohort() const { return m_is_LDD_cohort; }
    Dispersal_vector get_dispersal_vector() const { return m_dispersal_vector; }

    double get_seedbank_cones(const Parameters& params) const;
    int get_seedbank(const Parameters& params) const;
    bool is_there_seeds_in_canopy(const Parameters& params) const;

    // setters
    inline void set_num_plants(const int num_plants)
    {
        m_num_plants = num_plants;
    }

private:
    // # Sub-submodels
    // ## Seed production sub-submodels
    // Note: the sub-submodels firm seeds and viable seeds are just proportions
    // thus, they are placed directly in the SeedProduction() submodel
    double flower_production_age(const Parameters& params);
    double flower_age_impacts(const Parameters& params);
    double flower_weather_impacts(const Parameters& params,
                                  const std::vector<double>& weather);
//    void truncate_flowers(const Parameters& params, double& num_flowers);

    // plant mortality scenarios
    double plant_mort_age(const Parameters& params);
    double plant_mort_weather_relative(const Parameters& params
                                       , const std::vector<double>& weather);
    double plant_mort_weather_absolute(const Parameters& params
                                       , const std::vector<double>& weather);
    double plant_mort_spring_autumn(const Parameters& params
                                    , const std::vector<double>& weather);
    double plant_mort_SDD_autumn_LDD_spring(const Parameters& params
                                            , const std::vector<double>& weather);
    double plant_mort_resi_autumn_immi_spring(const Parameters& params
                                              , const std::vector<double>& weather);
    double plant_mort_lowest(const Parameters& params
                             , const std::vector<double>& weather);

    // upper and lower limit values of plant mortaliy
    double limit_mort(const Parameters& params, double mort_probability);

    // mortality cohort-based and individual-based functions
    void plant_mort_cohort_based(const double mort_probability
                                 , const double carrying_capacity);

    void plant_mort_plant_based(const Parameters& params
                                , const double mort_probability);

    // private setter
    void set_seedbank_to_zero(const Parameters& params);

    // private member variables
    int m_age        {0};
    int m_num_plants {0};
    std::vector<double> m_seedbank_cones;
    bool m_is_postfire    {true};
    Gene_flow m_gene_flow {1};
    bool m_is_LDD_cohort  {false};
    Dispersal_vector m_dispersal_vector {Dispersal_vector::SDD_wind};
    std::map<unsigned, std::deque<Plant>> m_plants_list;
};

#endif // COHORT_H
