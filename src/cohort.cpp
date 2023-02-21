#include "cohort.h"

Cohort::Cohort(const int age
               , const int num_plants
               , const bool is_postfire
               , const Gene_flow& gene_flow
               , const bool is_LDD_cohort
               , const Dispersal_vector dispersal_vector)
    : m_age{age}
    , m_num_plants{num_plants}
    , m_is_postfire{is_postfire}
    , m_gene_flow{gene_flow}
    , m_is_LDD_cohort{is_LDD_cohort}
    , m_dispersal_vector{dispersal_vector}
{}

//-----------------------------------------------------------------------------

void Cohort::plant_mortality(const Parameters& params
                             , const std::vector<double>& weather
                             , const int carrying_capacity)
{
    assert(m_age >= 0);

    if (m_age > params.plant_longevity)
    {
        m_num_plants = 0;

        if (params.model_type == Model_type::plant_based)
        {
            m_plants_list.clear();
        }
        std::cout << "Cohort reached maximum plant longevity!\n\a";
    }
    else
    {
        double mort_probability {0.0};

        switch (params.mort_scenario)
        {
        case Mortality_scenario::age_weather_relative :
            mort_probability = plant_mort_age(params);

            if (params.is_climate)
            {
                mort_probability *= plant_mort_weather_relative(params, weather);
            }
            break;
        case Mortality_scenario::age_weather_absolute :
            mort_probability = plant_mort_age(params);

            if (params.is_climate)
            {
                mort_probability += plant_mort_weather_absolute(params, weather);
            }
            break;
        case Mortality_scenario::spring_autumn_mean :
            mort_probability = plant_mort_spring_autumn(params, weather);
            break;
        case Mortality_scenario::SDD_autumn_LDD_spring :
            mort_probability = plant_mort_SDD_autumn_LDD_spring(params, weather);
            break;
        case Mortality_scenario::resi_autum_immi_spring :
            mort_probability = plant_mort_resi_autumn_immi_spring(params, weather);
            break;
        case Mortality_scenario::lowest_mortality :
            mort_probability = plant_mort_lowest(params, weather);
            break;
        }

        mort_probability = limit_mort(params, mort_probability);

        if (mort_probability > 1.0)
        {
            std::cerr << "mortality probability above 1.0\n"
                         "it could be only above 1.0 in a stress tests"
                         "and sensitivity analysis\n\a";

            std::exit(EXIT_FAILURE);
        }
        else if (mort_probability <= 0.0)
        {
            std::cerr << "mortality probability below 0.0\n"
                         "it could be only below 0.0 in a stress tests"
                         "and sensitivity analysis\n\a";

            std::exit(EXIT_FAILURE);
        }

        switch (params.model_type)
        {
            case Model_type::cohort_based :
                plant_mort_cohort_based(mort_probability, carrying_capacity);
                break;

            case Model_type::plant_based :
                if (m_age >= params.plant_age_young)
                {
                    plant_mort_plant_based(params, mort_probability);
                }
                else
                {
                    // at seedling stage plants have no seedbank and have the
                    // same mortality probability
                    plant_mort_cohort_based(mort_probability, carrying_capacity);
                }
                break;
        }
    }
}

//-----------------------------------------------------------------------------

void Cohort::plant_mort_cohort_based(const double mort_probability
                                     , const double carrying_capacity)
{
    // Expects
    assert(mort_probability >= 0.0 && mort_probability <= 1.0);
    assert(carrying_capacity >= 0.0);

    int dead_plants {0};

    if (m_num_plants > carrying_capacity)
    {
        dead_plants = Global::prng.binomial(m_num_plants
                                            , mort_probability);
    }
    else
    {
        for (int i = 0; i < m_num_plants; ++i)
        {
            if (Global::prng.uniform_01() <= mort_probability)
            {
                ++dead_plants;
            }
        }
    }

    if (m_num_plants > dead_plants)
    {
        m_num_plants -= dead_plants;
    }
    else
    {
        m_num_plants = 0;
    }

    // Ensures
    assert(m_num_plants >= 0);
}

//-----------------------------------------------------------------------------

void Cohort::plant_mort_plant_based(const Parameters& params
                                    , const double mort_probability)
{
    // Expects
    assert(params.model_type == Model_type::plant_based);

    for (auto& [key, val] : m_plants_list)
    {
        const double mort_plant_type {
            limit_mort(params
                       , mort_probability
                       + params.plant_type_list.at(key).m_mort_prob)
        };

        val.erase(
            std::remove_if(
                std::begin(val)
                , std::end(val)
                , [mort_plant_type]([[maybe_unused]] const Plant& plant)
                {
                    return Global::prng.uniform_01() <= mort_plant_type;
                }
            ),
            std::end(val)
        );
    }

    // update total number of plants in cohort
    m_num_plants = 0;

    for (const auto& [key, val] : m_plants_list)
    {
        m_num_plants += static_cast<int>(std::size(val));
    }
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_age(const Parameters& params)
{
    // Expects
    assert(m_age >= 0);

    double mort_probability {0.0};

    if (m_age == 0)
    {
      mort_probability = params.mort_recruit_postfire_mean;

    }
    else if (m_age > 0 && m_age <= params.plant_age_adult)
    {
      mort_probability = params.mort_a
                         + (params.mort_b / m_age)
                         + (params.mort_c / std::pow(m_age, 2));
    }
    else if (m_age > params.plant_age_adult)
    {
        mort_probability = params.mort_a
                           + (params.mort_b / params.plant_age_adult)
                           + (params.mort_c / std::pow(m_age, 2));

        if (m_age > params.senescence_age)
        {
            mort_probability = mort_probability
                               + (params.senescence_increase
                                  * (m_age - params.senescence_age));
        }

    }

    return mort_probability;
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_weather_absolute(
        const Parameters& params
        , const std::vector<double>& weather)
{
    // Expects
    assert(m_age >= 0);

    double weather_impacts {0.0};

    if (m_age == 0)
    {
        weather_impacts = ((params.long_term_rain - weather[0]) / 100.0)
                          * params.recruit_weather;
    }
    else
    {
        weather_impacts = params.mort_d * weather[0] + params.mort_e;
    }

    return weather_impacts;
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_weather_relative(
        const Parameters& params
        , const std::vector<double>& weather)
{
    // Expects
    assert(m_age >= 0);

    double weather_impacts {0.0};

    if (m_age == 0)
    {
        weather_impacts = 1 + (((params.long_term_rain - weather[0]) / 100.0)
                                * params.recruit_weather);
    }
    else
    {
        weather_impacts = 1 + (params.mort_d * weather[0] + params.mort_e);
    }

    return weather_impacts;
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_SDD_autumn_LDD_spring(
        const Parameters& params
        , const std::vector<double>& weather)
{
    // Expects
    assert(m_age >= 0);

    double mort_probability {0.0};

    if (m_age == 0)  // recruitment
    {
        if (m_is_LDD_cohort)
        {
            mort_probability = params.mort_recruit_postfire_min;
        }
        else
        {
            mort_probability = params.mort_recruit_postfire_mean;
        }

        const double weather_impacts {
            ((params.long_term_rain - weather[0]) / 100.0)
            * params.recruit_weather
        };
        mort_probability += weather_impacts;
    }
    else
    {
        if (m_is_LDD_cohort)
        {
            if (m_age <= params.plant_age_adult)
            {
                // srping site equation (Keith et al. 2014)
                const double lm_spring {
                    (params.mort_a / m_age) +
                            (params.mort_b * weather[0]) +
                            params.mort_c
                };

                mort_probability = lm_spring;
            }
            else
            {
                const double lm_spring {
                    (params.mort_a / params.plant_age_adult) +
                            (params.mort_b * weather[0]) +
                            params.mort_c
                };

                mort_probability = lm_spring;

                // mortality increases when plant age is greater than senescence
                if (m_age > params.senescence_age)
                {
                    mort_probability = mort_probability +
                                       (params.senescence_increase *
                                        (m_age - params.senescence_age));
                }
            }
        }
        else  // !is_LDD_cohort, i.e. cohort is from SDD
        {
            if (m_age <= params.plant_age_adult)
            {
                // autumn site equation (Keith et al. 2014)
                const double lm_autumn {
                    (params.mort_d / m_age) +
                            (params.mort_e * weather[0]) +
                            params.mort_f
                };

                mort_probability = lm_autumn;
            }
            else
            {
                const double lm_autumn {
                    (params.mort_d / params.plant_age_adult) +
                            (params.mort_e * weather[0]) +
                            params.mort_f
                };

                mort_probability = lm_autumn;

                if (m_age > params.senescence_age)
                {
                    mort_probability = mort_probability +
                                       (params.senescence_increase *
                                        (m_age - params.senescence_age));
                }
            }
        }
    }

    return mort_probability;
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_resi_autumn_immi_spring(
        const Parameters& params
        , const std::vector<double>& weather)
{
    // Expects
    assert(m_age >= 0);

    double mort_probability {0.0};

    if (m_age == 0)  // recruitment
    {
        if (m_gene_flow.get_previous() != m_gene_flow.get_current())
        {
            mort_probability = params.mort_recruit_postfire_min;
        }
        else
        {
            mort_probability = params.mort_recruit_postfire_mean;
        }

        const double weather_impacts {
            ((params.long_term_rain - weather[0]) / 100.0)
            * params.recruit_weather
        };
        mort_probability += weather_impacts;
    }
    else
    {
        if (m_gene_flow.get_previous() != m_gene_flow.get_current())
        {
            if (m_age <= params.plant_age_adult)
            {
                // spring site equation (Keith et al. 2014)
                const double lm_spring {
                    (params.mort_a / m_age) +
                            (params.mort_b * weather[0]) +
                            params.mort_c
                };

                mort_probability = lm_spring;
            }
            else
            {
                const double lm_spring {
                    (params.mort_a / params.plant_age_adult) +
                            (params.mort_b * weather[0]) +
                            params.mort_c
                };

                mort_probability = lm_spring;

                // mortality increases when plant age is greater than senescence
                if (m_age > params.senescence_age)
                {
                    mort_probability = mort_probability +
                                       (params.senescence_increase *
                                        (m_age - params.senescence_age));
                }
            }
        }
        else  // is resident cohort
        {
            if (m_age <= params.plant_age_adult)
            {
                // autumn site equation (Keith et al. 2014)
                const double lm_autumn {
                    (params.mort_d / m_age) +
                            (params.mort_e * weather[0]) +
                            params.mort_f
                };

                mort_probability = lm_autumn;
            }
            else
            {
                const double lm_autumn {
                    (params.mort_d / params.plant_age_adult) +
                            (params.mort_e * weather[0]) +
                            params.mort_f
                };

                mort_probability = lm_autumn;

                if (m_age > params.senescence_age)
                {
                    mort_probability = mort_probability +
                                       (params.senescence_increase *
                                        (m_age - params.senescence_age));
                }
            }
        }
    }

    return mort_probability;
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_spring_autumn(
        const Parameters& params
        , const std::vector<double>& weather)
{
    // Expects
    assert(m_age >= 0);

    double mort_probability {0.0};

    if (m_age == 0)  // recruitment
    {
        const double weather_impacts {
            ((params.long_term_rain - weather[0]) / 100.0)
            * params.recruit_weather
        };
        mort_probability = params.mort_recruit_postfire_mean;
        mort_probability += weather_impacts;
    }
    else
    {
        // arithmetic mean of spring and autumn (Keith et al. 2014)
        if (m_age <= params.plant_age_adult)
        {
            const double lm_spring {
                (params.mort_a / m_age)
                + (params.mort_b * weather[0])
                + params.mort_c
            };
            const double lm_autumn {
                (params.mort_d / m_age)
                + (params.mort_e * weather[0])
                + params.mort_f
            };

            mort_probability = (lm_spring + lm_autumn) / 2.0;
        }
        else
        {
            const double lm_spring {
                (params.mort_a / params.plant_age_adult)
                + (params.mort_b * weather[0])
                + params.mort_c
            };
            const double lm_autumn {
                (params.mort_d / params.plant_age_adult)
                + (params.mort_e * weather[0])
                + params.mort_f
            };

            mort_probability = (lm_spring + lm_autumn) / 2.0;

            if (m_age > params.senescence_age)
            {
                mort_probability = mort_probability
                                   + (params.senescence_increase
                                   * (m_age - params.senescence_age));
            }
        }
    }

    return mort_probability;
}

//-----------------------------------------------------------------------------

double Cohort::plant_mort_lowest(const Parameters& params
                                 , const std::vector<double>& weather)
{
    // Expects
    assert(m_age >= 0);

    double mort_probability {0.0};

    if (m_age == 0)
    {
        const double weather_impacts {
            ((params.long_term_rain - weather[0]) / 100.0)
            * params.recruit_weather
        };
        mort_probability = params.mort_recruit_postfire_min;
        mort_probability += weather_impacts;
    }
    else
    {
        if (m_age <= params.plant_age_adult)
        {
            const double lm_spring {
                (params.mort_a / m_age)
                + (params.mort_b * weather[0])
                + params.mort_c
            };
            const double lm_autumn {
                (params.mort_d / m_age)
                + (params.mort_e * weather[0])
                + params.mort_f
            };

            mort_probability = std::min(lm_spring, lm_autumn);
        }
        else
        {
            const double lm_spring {
                (params.mort_a / params.plant_age_adult)
                + (params.mort_b * weather[0])
                + params.mort_c
            };
            const double lm_autumn {
                (params.mort_d / params.plant_age_adult)
                        + (params.mort_e * weather[0])
                        + params.mort_f
            };

            mort_probability = std::min(lm_spring, lm_autumn);

            if (m_age > params.senescence_age)
            {
                mort_probability = mort_probability
                                    + (params.senescence_increase
                                    * (m_age - params.senescence_age));
            }
        }
    }

    return mort_probability;
}

//-----------------------------------------------------------------------------

double Cohort::limit_mort(const Parameters& params, double mort_probability)
{
    // Expects
    assert(m_age >= 0);

    if (m_age == 0)
    {
        if (mort_probability > params.mort_recruit_postfire_max)
        {
            mort_probability = params.mort_recruit_postfire_max;
        }
        else if (mort_probability < params.mort_recruit_postfire_min)
        {
            mort_probability = params.mort_recruit_postfire_min;
        }

        if (!m_is_postfire)
        {
            // convert mortality to survival
            double survival_probability {1.0 - mort_probability};

            if (survival_probability > 0.0)
            {
                survival_probability = survival_probability
                                       * (params.recruit_interfire / 100.0);

                // convert survival to mortality
                mort_probability = 1.0 - survival_probability;
            }
            else
            {
                mort_probability = 0.0;
            }
        }
    }
    else
    {
        if (mort_probability > params.mort_recruit_postfire_max)
        {
            mort_probability = params.mort_recruit_postfire_max;
        }
        else if (mort_probability < params.mort_min)
        {
            mort_probability = params.mort_min;
        }
    }

    // Ensures
    assert(mort_probability >= 0 && mort_probability <= 1.0);

    return mort_probability;
}

//-----------------------------------------------------------------------------

void Cohort::aging_plants_cones(const Parameters& params)
{
    // aging of the cohort, i.e., plants within the cohort are one year older
    ++m_age;

    if (m_age <= params.plant_age_young)
    {
        return;  ///>@note plants have no seedbank yet, i.e. vector is empty
                 ///>@note seedbank is initiated in density_regulation() in Population
    }

    // aging of cones
    switch (params.model_type)
    {
        case Model_type::cohort_based :
            assert(!std::empty(m_seedbank_cones));

            // aging of seeds: rotate the elements of vector to the right
            std::rotate(std::begin(m_seedbank_cones)
                        , std::end(m_seedbank_cones) - 1
                        , std::end(m_seedbank_cones));
            break;

        case Model_type::plant_based :
            for ([[maybe_unused]] auto& [key, val] : m_plants_list)
            {
                for (auto& plants : val)
                {
                    plants.aging_cones();
                }
            }
            break;
    }
}

//-----------------------------------------------------------------------------

int Cohort::interfire_seed_dispersal(const Parameters& params)
{
    int total_seeds {0};

    switch (params.model_type)
    {
        case Model_type::cohort_based :
        {
            double seeds_dispersed {0.0};

            for (unsigned i = params.cone_cycle; i < std::size(m_seedbank_cones); ++i)
            {
                seeds_dispersed += m_seedbank_cones.at(i) *
                                   params.viable_seeds_dispersed_cone_age.at(i);
            }

            total_seeds = static_cast<int>((m_num_plants * seeds_dispersed) + 0.5);
            break;
        }
        case Model_type::plant_based :
        {
            for ([[maybe_unused]] auto& [key, val] : m_plants_list)
            {
                for (auto& plants : val)
                {
                    total_seeds += plants.interfire_seed_dispersal(params);
                }
            }
            break;
        }
    }

    // Ensures
    assert(total_seeds >= 0);
    return total_seeds;
}

//-----------------------------------------------------------------------------

void Cohort::cone_production_cohort_based(const Parameters& params,
                                          const std::vector<double>& weather,
                                          const Habitat_quality& habitat_quality)
{
    // Expects
    assert(params.model_type == Model_type::cohort_based);
    assert(m_age >= params.plant_age_young);
    assert(!std::empty(m_seedbank_cones));

    double new_cones { flower_production_age(params) };

    if (params.is_climate == true)
    {
        new_cones *= flower_weather_impacts(params, weather);
    }

    new_cones *= habitat_quality.m_effect_flowers;

    new_cones *= params.pollination_success / 100.0;

    if (new_cones < 0.0) { new_cones = 0.0; }

    m_seedbank_cones.at(0) = new_cones;
}

//-----------------------------------------------------------------------------

void Cohort::cone_production_plant_based(
        const Parameters& params
        , const std::map<std::pair<unsigned, unsigned>
            , std::vector<double>>& weather_fuzzy_values
        , const Habitat_quality& habitat_quality)
{
    // Expects
    assert(params.model_type == Model_type::plant_based);
    assert(m_age >= params.plant_age_young);
    assert(std::empty(m_seedbank_cones));

    for (auto& [key, val] : m_plants_list)
    {
        for (auto& plant : val)
        {
            double total_cones {0.0};

            for (const auto& mf1 : params.fuzzy_set_list.at(std::make_pair(0, 0)))
            {
                const double proportion {
                    weather_fuzzy_values.at(
                                std::make_pair(0, 0)).at(mf1->get_mf_index())
                };

                assert(proportion >= 0.0 && proportion <= 1.0);

                if (proportion > 0.0)
                {
                    unsigned fuzzy_set_selected { mf1->get_mf_index() };

                    for (const auto& mf2 : params.fuzzy_set_list.at(
                             std::make_pair(1, fuzzy_set_selected)))
                    {
                        const double proportion_2 {
                            weather_fuzzy_values.at(
                                        std::make_pair(
                                            1
                                            , fuzzy_set_selected)).at(
                                        mf2->get_mf_index())
                        };

                        assert(proportion_2 >= 0.0 && proportion_2 <= 1.0);

                        if (proportion_2 > 0.0)
                        {
                            unsigned mf_selected { mf2->get_mf_index() };

                            // pick number of flowers from selected distribution at random
                            double cones {0.0};

                            const auto& flower_distr = params.flower_distribution.at(
                                        std::make_tuple(
                                            fuzzy_set_selected
                                            , mf_selected
                                            , key));

                            switch (flower_distr.m_dist_type)
                            {
                                case Distribution_type::null :
                                    cones = 0.0;
                                    break;
                                case Distribution_type::pois :
                                    cones = Global::prng.poisson(
                                                flower_distr.m_param_a);
                                    break;
                                case Distribution_type::nbinom :
                                    cones = Global::prng.negative_binomial(
                                                    flower_distr.m_param_a,
                                                    flower_distr.m_param_b);
                                    break;
                                case Distribution_type::nbinom_mu :
                                    cones = Global::prng.negative_binomial_mu(
                                                    flower_distr.m_param_a,
                                                    flower_distr.m_param_b);
                                    break;
                                case Distribution_type::geom :
                                    cones = Global::prng.geometric(
                                                flower_distr.m_param_a);
                                    break;
                            }

                            // apply proportion of the membership functions
                            total_cones += cones * proportion * proportion_2;
                        }
                    }
                }
            }

            total_cones *= flower_age_impacts(params);

            total_cones *= habitat_quality.m_effect_flowers;

            total_cones *= params.pollination_success / 100.0;

            if (total_cones < 0.0) { total_cones = 0.0; }

            plant.add_new_cones(static_cast<int>(total_cones + 0.5));
        }
    }
}

//-----------------------------------------------------------------------------

double Cohort::flower_production_age(const Parameters& params)
{
    return Equation::logistic_growth(params.flower_age_a
                                     , params.flower_age_b
                                     , params.flower_age_c
                                     , m_age);
}

//-----------------------------------------------------------------------------

double Cohort::flower_age_impacts(const Parameters& params)
{
    // Expects
    assert(params.model_type == Model_type::plant_based);

    double impact_prop {1.0};

    if (m_age < params.plant_age_adult)
    {
        const double flowers_current_age {
            Equation::logistic_growth(params.flower_age_a
                                      , params.flower_age_b
                                      , params.flower_age_c
                                      , m_age)
        };

        const double flowers_highest_production {
            Equation::logistic_growth(params.flower_age_a
                                      , params.flower_age_b
                                      , params.flower_age_c
                                      , params.plant_age_adult)
        };

        if (flowers_highest_production <= 0.0)
        {
            std::cerr << "ERROR! At least one parameter of flower age or the"
                      << "plant age at highest production is is wrong: \n"
                      << "flower_age_a = " << params.flower_age_a
                      << '\n'
                      << "flower_age_b = " << params.flower_age_b
                      << '\n'
                      << "flower_age_c = " << params.flower_age_c
                      << '\n'
                      << "plant_age_highest_prod = " << params.plant_age_adult
                      << "\n\a";
            std::exit(EXIT_FAILURE);
        }

        impact_prop = flowers_current_age / flowers_highest_production;
    }

    // Ensures
    assert(impact_prop >= 0.0 && impact_prop <= 1.0);
    return impact_prop;
}

//-----------------------------------------------------------------------------

double Cohort::flower_weather_impacts(const Parameters& params
                                      , const std::vector<double>& weather)
{
    double model_longterm {1.0};
    double model_scenario {1.0};

    model_longterm = params.flower_weather_a
                     * params.long_term_rain
                     + params.flower_weather_b;

    switch (params.flower_scenario)
    {
      case Flower_scenario::baseline :
        model_scenario = params.flower_weather_a * weather[0]
                         + params.flower_weather_b;
        break;

      case Flower_scenario::current :
        model_scenario = params.flower_weather_c * weather[1]
                         + params.flower_weather_d * weather[2]
                         + params.flower_weather_e;
        break;

      case Flower_scenario::current_nb :
        model_scenario = std::exp(
                params.flower_weather_c * weather[1]
                + params.flower_weather_d * weather[2]
                + params.flower_weather_e);
        break;
    }

    if (model_scenario < 0.0)
    {
        model_scenario = 0.0;
    }

    if (model_longterm <= 0.0)
    {
        throw std::runtime_error(
            "Division by zero, in flower production at longterm conditions\n"
        );
    }

    // Ensures
    assert(model_scenario >= 0.0);
    assert(model_longterm > 0.0);
    return model_scenario / model_longterm;
}


//-----------------------------------------------------------------------------

void Cohort::set_seedbank_to_zero(const Parameters& params)
{
    switch (params.model_type)
    {
        case Model_type::cohort_based :
            std::fill(std::begin(m_seedbank_cones), std::end(m_seedbank_cones), 0.0);
            break;

        case Model_type::plant_based :
            for ([[maybe_unused]] auto& [key, val] : m_plants_list)
            {
                for (auto& plant : val)
                {
                    plant.set_seedbank_to_zero();
                }
            }
            break;
    }
}

// ----------------------------------------------------------------------------

int Cohort::get_seedbank(const Parameters& params) const
{
    double total_viable_seeds {0.0};

    switch (params.model_type)
    {
    case Model_type::cohort_based :
    {
        if (!std::empty(m_seedbank_cones))
        {
            std::vector<double> seedbank(params.seed_longevity, 0.0);

            std::transform(std::begin(m_seedbank_cones)
                           + params.cone_cycle
                           , std::end(m_seedbank_cones)
                           , std::begin(params.viable_seeds_cone_age)
                           + params.cone_cycle
                           , std::begin(seedbank)
                           , std::multiplies<double>());

            total_viable_seeds = std::accumulate(std::cbegin(seedbank)
                                                 , std::cend(seedbank)
                                                 , 0.0
                                                 , std::plus<double>());

            total_viable_seeds *= m_num_plants;
        }
        break;
    }
    case Model_type::plant_based:
    {
        for ([[maybe_unused]] const auto& [key, val] : m_plants_list)
        {
            for (const auto& plant : val)
            {
                total_viable_seeds += plant.get_seedbank(params);
            }
        }
        break;
    }
    }

    // Ensures
    assert(total_viable_seeds >= 0.0);

    return static_cast<int>(total_viable_seeds + 0.5);
}

//-----------------------------------------------------------------------------

double Cohort::get_seedbank_cones(const Parameters& params) const
{
    double total_cones {0.0};

    switch (params.model_type)
    {
        case Model_type::cohort_based :
        {
            if (!std::empty(m_seedbank_cones))
            {
                total_cones = std::accumulate(
                            std::cbegin(m_seedbank_cones)
                            + params.cone_cycle
                            , std::cend(m_seedbank_cones)
                            , 0.0
                            , std::plus<double>());

                total_cones *= m_num_plants;
            }
            break;
        }
        case Model_type::plant_based :
        ///@todo TODO change to std::accumulate ??
        // https://stackoverflow.com/questions/31354947/adding-all-values-of-map-using-stdaccumulate

//        std::accumulate(std::begin(floor_plan)
//                      , std::end(floor_plan)
//                      , 0
//                      , [] (int value, const std::map<int, int>::value_type& p)
//                           { return value + p.second; }
//                       );
            for ([[maybe_unused]] const auto& [key, val] : m_plants_list)
            {
                for (const auto& plant : val)
                {
                    total_cones += plant.get_seedbank_cones(params);
                }
            }
            break;
    }

    assert(total_cones >= 0.0);

    return total_cones;
}

//-----------------------------------------------------------------------------

bool Cohort::is_there_seeds_in_canopy(const Parameters& params) const
{
    bool is_seed {false};

    switch (params.model_type)
    {
    case Model_type::cohort_based:
    {
        if (get_seedbank(params) >= 1.0) { is_seed = true; }
        break;
    }
    case Model_type::plant_based:
    {
        double total_viable_seeds {0.0};

        for ([[maybe_unused]] const auto& [key, val] : m_plants_list)
        {
            for (const auto& plants : val)
            {
                if (total_viable_seeds >= 1.0)
                {
                    is_seed = true;
                    break;
                }

                total_viable_seeds += plants.get_seedbank(params);
            }

            if (total_viable_seeds >= 1.0)
            {
                is_seed = true;
                break;
            }
        }
        break;
    }
    }
    return is_seed;
}

//-----------------------------------------------------------------------------

void Cohort::init_mature_cohort(const Parameters& params)
{
    if (m_age == params.plant_age_young)
    {
        switch (params.model_type)
        {
        case Model_type::cohort_based:
            m_seedbank_cones.assign(
                        params.cone_cycle
                        + params.seed_longevity
                        , 0.0);
            break;
        case Model_type::plant_based:
            for (const auto& [key, val] : params.plant_type_list)
            {
                m_plants_list.emplace(key, std::deque<Plant>());
            }
            break;
        }
    }
}

//-----------------------------------------------------------------------------

//int Cohort::postfire_viable_seeds(const Parameters& params) const
//{
//    /// @note NOTE in the future for other sp. the viable_seeds might be different
//    /// for each seed_age, i.e. viability decreases with seed age.

//    double total_viable_seeds {0.0};

//    switch (params.model_type)
//    {
//        case Model_type::cohort_based :
//            total_viable_seeds = get_seedbank(params);
//            break;

//        case Model_type::plant_based :
//            for (const auto& [key, val] : m_plants_list)
//            {
//                for (const auto& plants : val)
//                {
//                    total_viable_seeds += plants.get_seedbank(params);
//                }
//            }
//            break;
//    }
//    return static_cast<int>(total_viable_seeds + 0.5);
//}
