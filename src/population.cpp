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

#include "population.h"


//-----------------------------------------------------------------------------

Population::Population(
        const int id,
        const int num_cells,
        const int carrying_capacity,
        std::string_view habitat_quality)

    : m_id{id}
    , m_num_cells{num_cells}
    , m_carrying_capacity{carrying_capacity}
    , m_habitat_quality{habitat_quality}
{}

//-----------------------------------------------------------------------------

void Population::plant_mortality(const Parameters& params, const std::vector<double>& weather)
{
    for (auto& cohort : m_cohort_list)
    {
        cohort.plant_mortality(params, weather, m_carrying_capacity);
    }

    remove_empty_cohorts();
}

//-----------------------------------------------------------------------------

void Population::aging_plants_cones(const Parameters& params)
{
    for (auto& cohort : m_cohort_list)
    {
        cohort.aging_plants_cones(params);
    }
}

//-----------------------------------------------------------------------------
///@todo TODO I should disssentangle this function into several smaller functions??
void Population::density_regulation(const Parameters& params)
{
    int living_adults {0};
    int cohort_index  {0};
//    bool is_postfire  {false};

    // 1. get cohorts at maturity age
    std::vector<std::pair<int, int>> cohorts_at_maturity;

    for (auto& cohort : m_cohort_list)
    {
        if (cohort.get_age() > params.plant_age_young)
        {
            living_adults += cohort.get_num_plants();
        }

        if (cohort.get_age() == params.plant_age_young)
        {
            cohorts_at_maturity.emplace_back(
                        std::make_pair(
                            cohort_index,
                            cohort.get_num_plants() ));

            // evaluate if at least one cohort is at maturity age
            // for incrementing the number of generations (below)
//            if (cohort.get_is_postfire() == true)
//            {
//                is_postfire = true;
//            }
        }
        ++cohort_index;
    }

    // 2. count the total number of plants at maturity age
    int total_plants_at_maturity_age {0};

    for (const auto& mature_cohort : cohorts_at_maturity)
    {
        total_plants_at_maturity_age += mature_cohort.second;
    }

    int empty_places {0};

    ///@note NOTE I am applying density regulation respecting the carrying_capacity
    /// of the entire population. If I implement the cohorts spatially (i.e.
    /// track the cell locations of each cohort within population), then
    /// the density regulation should be regarding those cells.
    if (m_carrying_capacity > living_adults)
    {
        empty_places = m_carrying_capacity - living_adults;
    }

    ///@note NOTE now the generations are counted in evaluate_year() when
    /// the population/dune is burned instead of counting when the population
    /// reached to the mature time with living individuals. This makes
    /// comparable total of generations with total of rescue effects
    /// (i.e. immigrants source)
//    if (empty_places > 0 && is_postfire == true &&
//        total_plants_at_maturity_age > 0)
//    {
//        ++m_num_generations;
//    }

    // 3. actual density regulation of plants at maturity age
    if (total_plants_at_maturity_age > empty_places)
    {
        assert(!std::empty(cohorts_at_maturity));

        if (std::size(cohorts_at_maturity) == 1)
        {
            // there is no need to select plants at random for only one cohort
            for (auto& cohort : m_cohort_list)
            {
                if (cohort.get_age() == params.plant_age_young)
                {
                    cohort.set_num_plants(empty_places);
                    break;
                }
            }
        }
        else
        {
            std::deque<int> plants_at_maturity;

            // 1. include each plant individually to the container
            for (const auto& cohort : cohorts_at_maturity)
            {
                // '.first' refers to the cohort_index and '.second' to num_plants
                plants_at_maturity.insert(
                            plants_at_maturity.end(),
                            static_cast<unsigned>(cohort.second),
                            cohort.first);
            }

            // 2. shuffle all elements (i.e. plants) of the container
            std::shuffle(std::begin(plants_at_maturity),
                         std::end(plants_at_maturity),
                         Global::prng.get_engine());

            // 3. take the surviving plants that establish in the empty places
            std::deque<int> new_established_adults = {
                std::cbegin(plants_at_maturity),
                std::cbegin(plants_at_maturity) + empty_places
            };

            // set the new established plants in their respective cohort
            cohort_index = 0;

            for (auto& cohort : m_cohort_list)
            {
                if (cohort.get_age() == params.plant_age_young)
                {
                    const auto count_adults = std::count(
                                std::cbegin(new_established_adults),
                                std::cend(new_established_adults),
                                cohort_index);

                    cohort.set_num_plants(count_adults);
                }
                ++cohort_index;
            }
        }
    }

    remove_empty_cohorts();

    // create seedbank for cohort-based and plant list for plant-based
    for (auto& cohort : m_cohort_list)
    {
        if (cohort.get_age() == params.plant_age_young)
        {
            cohort.init_mature_cohort(params);
        }
    }

    // classify plants types in plant-based model
    if (params.model_type == Model_type::plant_based)
    {
        for (auto& cohort : m_cohort_list)
        {
            if (cohort.get_age() == params.plant_age_young)
            {
                const int num_plants { cohort.get_num_plants() };

                for (int i = 0; i < num_plants; ++i)
                {
                    const double random_num {
                        Global::prng.uniform_01()
                    };

                    for (const auto& [key, val] : params.plant_type_share)
                    {
                        if (random_num <= val)
                        {
                            cohort.add_plant(params, key);
                            break;
                        }
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
///@todo TODO Add Expects()
void Population::interfire_seed_dispersal(const Parameters& params)
{
    // Expects(params.recruit_interfire > 0.0 && params.recruit_interfire <= 1.0)

    constexpr int interfire_cohort_age {0};
    constexpr bool is_LDD_cohort       {false};
    constexpr bool is_postfire_cohort  {false};

    std::vector<Cohort> interfire_cohort_list;

    // create interfire cohorts from interfire seed dispersal
    for (auto& cohort : m_cohort_list)
    {
        if (cohort.get_age() <= params.plant_age_young)
        {
            continue;    // plants have no mature cones
        }

        const int viable_seeds_dispersed { cohort.interfire_seed_dispersal(params) };

        if (viable_seeds_dispersed > 0)
        {
            interfire_cohort_list.push_back(
                        Cohort(
                            interfire_cohort_age,
                            viable_seeds_dispersed,
                            is_postfire_cohort,
                            cohort.get_gene_flow(),
                            is_LDD_cohort,
                            cohort.get_dispersal_vector()));
        }
    }

    // add interfire cohorts in cohort list
    if (!std::empty(interfire_cohort_list))
    {
        m_cohort_list.reserve(std::size(m_cohort_list)
                              + std::size(interfire_cohort_list));
        m_cohort_list.insert(std::cend(m_cohort_list)
                             , std::cbegin(interfire_cohort_list)
                             , std::cend(interfire_cohort_list));
    }
}

//-----------------------------------------------------------------------------

void Population::cone_production_cohort_based(const Parameters& params, const std::vector<double>& weather)
{
    assert(params.model_type == Model_type::cohort_based);

    const auto& habitat_effect_flowers = params.habitat_quality.at(m_habitat_quality);

    for (auto& cohort : m_cohort_list)
    {
        if (cohort.get_age() >= params.plant_age_young)
        {
            cohort.cone_production_cohort_based(params, weather, habitat_effect_flowers);
        }
    }
}

//-----------------------------------------------------------------------------

void Population::cone_production_plant_based(
        const Parameters& params,
        const std::map<std::pair<unsigned, unsigned>, std::vector<double>>& weather_fuzzy_values)
{
    assert(params.model_type == Model_type::plant_based);
    ///@note NOTE should I save this variable as member variable?
    const auto& habitat_effect_flowers = params.habitat_quality.at(m_habitat_quality);

    for (auto& cohort : m_cohort_list)
    {
        if (cohort.get_age() >= params.plant_age_young)
        {
            cohort.cone_production_plant_based(params, weather_fuzzy_values, habitat_effect_flowers);
        }
    }
}

//-----------------------------------------------------------------------------

void Population::fire_mortality()
{
    m_time_since_fire = 0;
    m_is_fire_event   = false;
    m_cohort_list.clear();
}

//-----------------------------------------------------------------------------

void Population::evaluate_year(const Parameters& params)
{
    if (!std::empty(m_cohort_list))
    {
        ++m_persistence_time;

        if (m_is_extinct)
        {
            m_is_extinct     = false;
            m_is_recolonized = true;

            reset_recolonizations_tracker();
        }
        else
        {
            if (m_time_since_fire == 0)
            {
                ++m_total_generations;
            }
        }
    }
    else
    {
        if (!m_is_extinct)
        {
            ++m_extinctions;
        }
        m_is_extinct = true;
    }

    evaluate_rescue_effects(params);
}

//-----------------------------------------------------------------------------

void Population::add_new_cohort(const int viable_seeds,
                                const Gene_flow& gene_flow,
                                const bool is_postfire_cohort,
                                const bool is_LDD_cohort,
                                const Dispersal_vector dispersal_vector)
{
    m_cohort_list.push_back(
                Cohort(
                    0,
                    viable_seeds,
                    is_postfire_cohort,
                    gene_flow,
                    is_LDD_cohort,
                    dispersal_vector));
}

//-----------------------------------------------------------------------------

void Population::initialize_population(std::string_view habitat_quality
                                       , const Cohort& init_cohort)
{
    m_cohort_list.clear();
    m_cohort_list.emplace_back(init_cohort);

    m_habitat_quality = habitat_quality;  ///@todo TODO why I save habitat quality in cohort and populations??

    m_time_since_fire  = 0;
    m_total_generations = 0;
    m_persistence_time = 0;
    m_extinctions      = 0;

    if (init_cohort.get_num_plants() > 0)
    {
        m_is_extinct = false;
    }
    else
    {
        m_is_extinct = true;
    }

    m_is_recolonized = false;
    m_is_fire_event  = false;
    m_total_fires    = 0;

    m_total_recolonizations.germination = 0;
    m_total_recolonizations.seedling    = 0;
    m_total_recolonizations.adult       = 0;
    m_total_recolonizations.source      = 0;

    m_total_LDD_residents.germination = 0;
    m_total_LDD_residents.seedling    = 0;
    m_total_LDD_residents.adult       = 0;
    m_total_LDD_residents.source      = 0;

    m_total_immigrants.germination = 0;
    m_total_immigrants.seedling    = 0;
    m_total_immigrants.adult       = 0;
    m_total_immigrants.source      = 0;

    reset_recolonizations_tracker();
    reset_LDD_residents_tracker();
    reset_immigrants_tracker();

    m_total_LDD_seeds_by_wind.outside    = 0;
    m_total_LDD_seeds_by_wind.unsuitable = 0;
    m_total_LDD_seeds_by_wind.immigrant  = 0;
    m_total_LDD_seeds_by_wind.resident   = 0;

    m_total_LDD_seeds_by_birds.immigrant = 0;
    m_total_LDD_seeds_by_birds.resident  = 0;

    m_genetic_pops_postfire.clear();
    m_genetic_pops.clear();

    m_genetic_pops_postfire_occupied.clear();
    m_genetic_pops_occupied.clear();

    m_immigrants_postfire_perc.clear();
    m_immigrants_perc.clear();

    for ([[maybe_unused]] auto& [key, val] : m_plants_per_population)
    {
        val.germination = 0;
        val.seedling    = 0;
        val.adult       = 0;
        val.source      = 0;
    }
}

//-----------------------------------------------------------------------------

void Population::remove_empty_cohorts()
{
  m_cohort_list.erase(
      std::remove_if(
          std::begin(m_cohort_list),
          std::end(m_cohort_list),
          [](const Cohort& cohort) {
              return cohort.get_num_plants() <= 0;
          }
      ),
      std::end(m_cohort_list));
}

//-----------------------------------------------------------------------------

void Population::evaluate_rescue_effects(const Parameters& params)
{
    ///@note NOTE demographic rescue must be before recolonization
    if (params.is_out_recol_rescue || params.is_out_metapopulation_pops)
    {
        evaluate_demographic_rescue(params);
        evaluate_recolonization_rescue(params);
    }

    if (params.is_out_plants_per_pop)
    {
        evaluate_plants_per_population(params);
    }
}

//-----------------------------------------------------------------------------

void Population::evaluate_demographic_rescue(const Parameters& params)
{
    assert(params.is_out_recol_rescue || params.is_out_metapopulation_pops);

    if (!m_is_recolonized && !m_is_extinct)
    {
        for (const auto& cohort : m_cohort_list)
        {
            if (cohort.get_is_LDD_cohort())
            {
                const int cohort_age { cohort.get_age() };

                if (cohort.get_gene_flow().get_previous() == m_id)
                {
                    if (cohort_age > params.plant_age_young)
                    {
                        // nothing
                    }
                    else if (cohort_age == params.plant_age_young)
                    {
                        m_LDD_residents_tracker.adult = true;
                    }
                    else if (cohort_age == 1)
                    {
                        m_LDD_residents_tracker.seedling = true;
                    }
                    else if (cohort_age == 0)
                    {
                        m_LDD_residents_tracker.germination = true;
                    }
                }
                else
                {
                    if (cohort_age > params.plant_age_young)
                    {
                        // nothing
                    }
                    else if (cohort_age == params.plant_age_young)
                    {
                        m_immigrants_tracker.adult = true;
                    }
                    else if (cohort_age == 1)
                    {
                        m_immigrants_tracker.seedling = true;
                    }
                    else if (cohort_age == 0)
                    {
                        m_immigrants_tracker.germination = true;
                    }
                }
            }
        }

        if (m_is_fire_event)
        {
            bool is_source_immigrants    {false};
            bool is_source_LDD_residents {false};

            for (auto& coh : m_cohort_list)
            {
                if (coh.get_is_LDD_cohort())
                {
                    if (coh.get_gene_flow().get_previous() == m_id)
                    {
                        if (!is_source_LDD_residents)
                        {
                            is_source_LDD_residents = coh.is_there_seeds_in_canopy(params);
                        }
                    }
                    else
                    {
                        if (!is_source_immigrants)
                        {
                            is_source_immigrants = coh.is_there_seeds_in_canopy(params);
                        }
                    }
                }

                if (is_source_LDD_residents && is_source_immigrants)
                {
                    break;
                }
            }

            m_LDD_residents_tracker.source = is_source_LDD_residents;
            m_immigrants_tracker.source    = is_source_immigrants;
        }
    }

    if (!m_is_recolonized && (m_is_fire_event || m_is_extinct))
    {
        increment_demographic_rescue();
    }
}

//-----------------------------------------------------------------------------
///@todo TODO ??change the total_recolonization variables (counter and tracker) to
/// an array?? traker --> std::array<4, bool>
void Population::evaluate_recolonization_rescue(const Parameters& params)
{
    assert(params.is_out_recol_rescue || params.is_out_metapopulation_pops);

    if (m_is_recolonized && !m_is_extinct)
    {
        for (const auto& cohort : m_cohort_list)
        {
            const int cohort_age { cohort.get_age() };

            if (cohort_age > params.plant_age_young)
            {
                // nothing
            }
            else if (cohort_age == params.plant_age_young)
            {
                m_recolonizations_tracker.adult = true;
            }
            else if (cohort_age == 1)
            {
                m_recolonizations_tracker.seedling = true;
            }
            else if (cohort_age == 0)
            {
                m_recolonizations_tracker.germination = true;
            }
        }

        if (m_is_fire_event)
        {
            bool is_source { false };

            for (auto& cohort : m_cohort_list)
            {
                if (!is_source)
                {
                    is_source = cohort.is_there_seeds_in_canopy(params);
                }

                if (is_source)
                {
                    break;
                }
            }
            m_recolonizations_tracker.source = is_source;
        }
    }

    if (m_is_recolonized && (m_is_fire_event || m_is_extinct))
    {
        m_is_recolonized = false;
        increment_recolonization_rescue();
    }
}

//-----------------------------------------------------------------------------

void Population::evaluate_plants_per_population(const Parameters& params)
{
    assert(params.is_out_plants_per_pop);

    if (!m_is_recolonized && !m_is_extinct)
    {
        for (const auto& cohort : m_cohort_list)
        {
            if (cohort.get_age() > params.plant_age_young)
            {
                // skip
            }
            else if (cohort.get_age() == params.plant_age_young)
            {
                m_plants_per_population.at(
                            cohort.get_gene_flow().get_previous()).adult += cohort.get_num_plants();
            }
            else if (cohort.get_age() == 1)
            {
                m_plants_per_population.at(
                            cohort.get_gene_flow().get_previous()).seedling += cohort.get_num_plants();
            }
            else if (cohort.get_age() == 0)
            {
                m_plants_per_population.at(
                            cohort.get_gene_flow().get_previous()).germination += cohort.get_num_plants();
            }

            if (m_is_fire_event)
            {
                if (cohort.is_there_seeds_in_canopy(params))
                {
                    m_plants_per_population.at(
                                cohort.get_gene_flow().get_previous()).source += cohort.get_num_plants();
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

void Population::increment_demographic_rescue()
{    
    m_total_LDD_residents.germination += m_LDD_residents_tracker.germination;
    m_total_LDD_residents.seedling    += m_LDD_residents_tracker.seedling;
    m_total_LDD_residents.adult       += m_LDD_residents_tracker.adult;
    m_total_LDD_residents.source      += m_LDD_residents_tracker.source;

    m_total_immigrants.germination += m_immigrants_tracker.germination;
    m_total_immigrants.seedling    += m_immigrants_tracker.seedling;
    m_total_immigrants.adult       += m_immigrants_tracker.adult;
    m_total_immigrants.source      += m_immigrants_tracker.source;

    reset_LDD_residents_tracker();
    reset_immigrants_tracker();
}

//-----------------------------------------------------------------------------

void Population::increment_recolonization_rescue()
{
    m_total_recolonizations.germination += m_recolonizations_tracker.germination;
    m_total_recolonizations.seedling    += m_recolonizations_tracker.seedling;
    m_total_recolonizations.adult       += m_recolonizations_tracker.adult;
    m_total_recolonizations.source      += m_recolonizations_tracker.source;

    reset_recolonizations_tracker();
}

//-----------------------------------------------------------------------------

void Population::reset_recolonizations_tracker()
{
    m_recolonizations_tracker.germination = false;
    m_recolonizations_tracker.seedling    = false;
    m_recolonizations_tracker.adult       = false;
    m_recolonizations_tracker.source      = false;
}

//-----------------------------------------------------------------------------

void Population::reset_LDD_residents_tracker()
{
    m_LDD_residents_tracker.germination = false;
    m_LDD_residents_tracker.seedling    = false;
    m_LDD_residents_tracker.adult       = false;
    m_LDD_residents_tracker.source      = false;
}

//-----------------------------------------------------------------------------

void Population::reset_immigrants_tracker()
{
    m_immigrants_tracker.germination = false;
    m_immigrants_tracker.seedling    = false;
    m_immigrants_tracker.adult       = false;
    m_immigrants_tracker.source      = false;
}

//-----------------------------------------------------------------------------

void Population::genetics_residents_immigrants_postfire(
        double& residents,
        double& immigrants)
{
    assert(residents == 0.0);
    assert(immigrants == 0.0);

    std::unordered_set<int> genetic_pops;

    for (const auto& cohort : m_cohort_list)
    {
        if (cohort.get_gene_flow().get_previous() == m_id)
        {
            residents += cohort.get_num_plants();
        }
        else
        {
            immigrants += cohort.get_num_plants();
        }

        genetic_pops.emplace(cohort.get_gene_flow().get_previous());
    }
    m_genetic_pops_postfire.push_back(
                static_cast<int>(std::size(genetic_pops)));

    if (!std::empty(genetic_pops))
    {
        m_genetic_pops_postfire_occupied.push_back(
                    static_cast<int>(std::size(genetic_pops)));


        assert(residents + immigrants > 0);

        const double immi_perc {
            (immigrants / (residents + immigrants)) * 100.0
        };

        assert(immi_perc >= 0.0 && immi_perc <= 100.00);

        m_immigrants_postfire_perc.push_back(immi_perc);
    }
}


//-----------------------------------------------------------------------------

///@brief store the total number of unique genetic populations each year
void Population::genetics_residents_immigrants(
        double& residents,
        double& immigrants)
{
    assert(residents == 0.0);   // Expects residents equal to zero
    assert(immigrants == 0.0);  // Expects immigrants equal to zero

    std::unordered_set<int> genetic_pops;

    for (const auto& cohort : m_cohort_list)
    {
        if (cohort.get_gene_flow().get_genetic() == m_id)
        {
            residents += cohort.get_num_plants();
        }
        else
        {
            immigrants += cohort.get_num_plants();
        }

        genetic_pops.emplace(cohort.get_gene_flow().get_genetic());
    }

    // store number of populations within patch only when there is at least one
    if (!std::empty(genetic_pops))
    {
        m_genetic_pops.push_back(static_cast<int>(std::size(genetic_pops)));
    }


    if (!std::empty(genetic_pops))
    {
        m_genetic_pops_occupied.push_back(
                    static_cast<int>(std::size(genetic_pops)));

        assert(residents + immigrants > 0);

        const double immi_perc {
            (immigrants / (residents + immigrants)) * 100.0
        };

        assert(immi_perc >= 0.0 && immi_perc <= 100.00);

        m_immigrants_perc.push_back(immi_perc);
    }
}
