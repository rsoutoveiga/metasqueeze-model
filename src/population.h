#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <unordered_set>

#include <string>
#include <string_view>

#include "point.h"
#include "cell.h"
#include "cohort.h"
#include "gene_flow.h"
#include "global.h"
#include "parameters.h"


struct LDD_seeds_by_wind
{
    long long int outside    {0};
    long long int unsuitable {0};
    long long int immigrant  {0};
    long long int resident   {0};
};

struct LDD_seeds_by_birds
{
    long long int resident  {0};
    long long int immigrant {0};
};

template<typename T>
struct Plant_stages
{
    T germination {0};
    T seedling    {0};
    T adult       {0};
    T source      {0};
};

typedef Plant_stages<long long int> Total_plants_stage;
typedef Plant_stages<int> Plant_stages_counter;
typedef Plant_stages<bool> Plant_stages_tracker;

//-----------------------------------------------------------------------------


class Population
{
public:

    Population(const int id
               , const int num_cells
               , const int carrying_capacity
               , std::string_view habitat_quality);


    // getters
    int get_num_cells() const { return m_num_cells; }

    // setters
    void set_cell_locations(const std::vector<Cell>& cell_locations)
    {
        m_cell_locations.clear();

        for (const auto& cell : cell_locations)
        {
            m_cell_locations.emplace_back(cell);
        }
    }

    void init_plants_per_population(const std::map<int, Population> metapopulation)
    {
        assert(std::empty(m_plants_per_population));

        constexpr Total_plants_stage plants {0, 0, 0, 0};

        for (const auto& [key, val] : metapopulation)
        {
            m_plants_per_population.emplace(key, plants);
        }
    }

    const std::map<int, Total_plants_stage>& get_plants_per_pop() const
    {
        return m_plants_per_population;
    }

    void initialize_population(std::string_view  habitat_quality
                               , const Cohort& init_cohort);

    void plant_mortality(const Parameters& params, const std::vector<double>& weather);
    void aging_plants_cones(const Parameters& params);
    void density_regulation(const Parameters& params);
    void interfire_seed_dispersal(const Parameters& params);
    void cone_production_cohort_based(
            const Parameters& params
            , const std::vector<double>& weather);
    void cone_production_plant_based(
            const Parameters& params
            , const std::map<std::pair<unsigned, unsigned>, std::vector<double>>& weather_fuzzy_values);

    void fire_mortality();

    void evaluate_year(const Parameters& params);

    void remove_all_cohorts() { m_cohort_list.clear(); }
    void remove_empty_cohorts();

    void add_new_cohort(const int viable_seeds,
                        const Gene_flow& gene_flow,
                        const bool is_postfire_cohort,
                        const bool is_LDD_cohort,
                        const Dispersal_vector dispersal_vector);


    void remove_all_cell_locations() { m_cell_locations.clear(); }



    void increment_time_since_fire() { ++m_time_since_fire; }
    void increment_generations() { ++m_total_generations; }

    void evaluate_rescue_effects(const Parameters& params);

    void increment_demographic_rescue();
    void increment_recolonization_rescue();


    void increment_total_LDD_seeds_by_wind_outside() {
        ++m_total_LDD_seeds_by_wind.outside;
    }

    void increment_total_LDD_seeds_by_wind_unsuitable() {
        ++m_total_LDD_seeds_by_wind.unsuitable;
    }

    void increment_total_LDD_seeds_by_wind_immigrant() {
        ++m_total_LDD_seeds_by_wind.immigrant;
    }

    void increment_total_LDD_seeds_by_wind_resident() {
        ++m_total_LDD_seeds_by_wind.resident;
    }

    void increment_total_LDD_seeds_by_birds_resident(const int seeds) {
        m_total_LDD_seeds_by_birds.resident += seeds;
    }

    void increment_total_LDD_seeds_by_birds_immigrant(const int seeds) {
        m_total_LDD_seeds_by_birds.immigrant += seeds;
    }

    void add_cell_for_LDD_by_birds(const Cell& cell)
    {
        m_cell_locations_for_LDD_birds.emplace_back(cell);
    }

    void add_one_cell_location(const unsigned row, const unsigned col)
    {
        m_cell_locations.emplace_back(Cell(row, col));
    }

    const Cell& get_random_cell_location() const
    {
        assert(!std::empty(m_cell_locations));  // Expects not empty
        const auto cell_random = Global::prng.uniform_unsigned(
                    0, std::size(m_cell_locations) - 1);
        return m_cell_locations.at(cell_random);
    }

    const Cell& get_cell_random_LDD_birds() const
    {
        assert(std::size(m_cell_locations_for_LDD_birds) > 0);

        const auto cell_random = Global::prng.uniform_unsigned(
                                     0,
                                     std::size(m_cell_locations_for_LDD_birds) - 1);

        return m_cell_locations_for_LDD_birds.at(cell_random);
    }

    // setters
    void reset_time_since_fire() { m_time_since_fire = 0; }

    void reset_total_generations() { m_total_generations = 0; }

    void set_is_fire_event(const bool is_fire_event)
    {
        m_is_fire_event = is_fire_event;

        if (is_fire_event)
        {
            ++m_total_fires;
        }
    }

    // getters
    int get_id() const { return m_id; }

    bool get_is_fire_event() const { return m_is_fire_event; }

    int get_total_fires() const { return m_total_fires; }

    const std::vector<Cell>& get_cell_locations() const
    {
        assert(!std::empty(m_cell_locations));  // Expects not empty
        return m_cell_locations;
    }

    int get_carrying_capacity() const { return m_carrying_capacity; }
    int get_time_since_fire() const { return m_time_since_fire; }
    int get_extinctions() const { return m_extinctions; }
    const std::vector<Cohort>& get_cohort_list() const { return m_cohort_list; }
    int get_total_generations() const { return m_total_generations; }
    int get_persistence_time() const { return m_persistence_time; }
    bool get_is_extinct() const { return m_is_extinct; }

    int get_total_recolonizations_germination() const {
        return m_total_recolonizations.germination;
    }

    int get_total_recolonizations_seedling() const {
        return m_total_recolonizations.seedling;
    }

    int get_total_recolonizations_adult() const {
        return m_total_recolonizations.adult;
    }

    int get_total_recolonizations_source() const {
        return m_total_recolonizations.source;
    }

    int get_total_LDD_residents_germination() const {
        return m_total_LDD_residents.germination;
    }

    int get_total_LDD_residents_seedling() const {
        return m_total_LDD_residents.seedling;
    }

    int get_total_LDD_residents_adult() const {
        return m_total_LDD_residents.adult;
    }

    int get_total_LDD_residents_source() const {
        return m_total_LDD_residents.source;
    }

    int get_total_immigrants_germination() const {
        return m_total_immigrants.germination;
    }

    int get_total_immigrants_seedling() const {
        return m_total_immigrants.seedling;
    }

    int get_total_immigrants_adult() const {
        return m_total_immigrants.adult;
    }

    int get_total_immigrants_source() const {
        return m_total_immigrants.source;
    }

    bool get_is_recolonized() const { return m_is_recolonized; }


    long long int get_total_LDD_seeds_by_wind_outside() const {
        return m_total_LDD_seeds_by_wind.outside;
    }
    long long int get_total_LDD_seeds_by_wind_unsuitable() const {
        return m_total_LDD_seeds_by_wind.unsuitable;
    }
    long long int get_total_LDD_seeds_by_wind_immigrant() const {
        return m_total_LDD_seeds_by_wind.immigrant;
    }
    long long int get_total_LDD_seeds_by_wind_resident() const {
        return m_total_LDD_seeds_by_wind.resident;
    }
    long long int get_total_LDD_seeds_by_birds_resident() const {
        return m_total_LDD_seeds_by_birds.resident;
    }
    long long int get_total_LDD_seeds_by_birds_immigrant() const {
        return m_total_LDD_seeds_by_birds.immigrant;
    }

    //--------------------------------------------------------

    void add_dist_edge_to_edge(const int pop_id, const double distance)
    {
        m_dist_edge_to_edge.emplace(pop_id, distance);
    }

    void add_dist_midpoint_to_edge(const int pop_id, const double distance)
    {
        m_dist_midpoint_to_edge.emplace(pop_id, distance);
    }

    void add_dist_midpoint_to_midpoint(const int pop_id, const double distance)
    {
        m_dist_midpoint_to_midpoint.emplace(pop_id, distance);
    }

    double get_dist_edge_to_edge(const int genetic_id)
    {
        return m_dist_edge_to_edge.at(genetic_id);
    }

    double get_dist_midpoint_to_edge(const int genetic_id)
    {
        return m_dist_midpoint_to_edge.at(genetic_id);
    }

    double get_dist_midpoint_to_midpoint(const int genetic_id)
    {
        return m_dist_midpoint_to_midpoint.at(genetic_id);
    }

    void set_midpoint(const Point& midpoint)
    {
        m_midpoint = midpoint;
    }

    Point get_midpoint()
    {
        return m_midpoint;
    }


    //-----------------------------
    // new calculations for the new member variables for metapopulation output
    ///@brief store the total number of unique genetic populations each year
    void genetics_residents_immigrants_postfire(double& residents, double& immigrants);

    //-----------------------------
    // new calculations for the new member variables for metapopulation output
    ///@brief store the total number of unique genetic populations each year
    void genetics_residents_immigrants(double& residents, double& immigrants);

    const std::vector<int>& get_genetic_pops_occupied() const
    {
        return m_genetic_pops_occupied;
    }

    const std::vector<int>& get_genetic_pops_postfire_occupied() const
    {
        return m_genetic_pops_postfire_occupied;
    }

    const std::vector<double>& get_immigrants_postfire_perc() const
    {
        return m_immigrants_postfire_perc;
    }

    const std::vector<double>& get_immigrants_perc() const
    {
        return m_immigrants_perc;
    }

    const std::vector<int>& get_genetic_pops_postfire() const
    {
        return m_genetic_pops_postfire;
    }

    const std::vector<int>& get_genetic_pops() const
    {
        return m_genetic_pops;
    }


private:

    // private member functions
    void reset_recolonizations_tracker();
    void reset_LDD_residents_tracker();
    void reset_immigrants_tracker();

    void evaluate_recolonization_rescue(const Parameters& params);
    void evaluate_demographic_rescue(const Parameters& params);
    void evaluate_plants_per_population(const Parameters& params);

    // private member variables
    const int m_id {1};
    const int m_num_cells {1};
    const int m_carrying_capacity {2500};
    std::vector<Cell> m_cell_locations;
    int m_time_since_fire  {0};
    int m_total_generations {0};
    int m_persistence_time {0};
    int m_extinctions      {0};
    bool m_is_extinct      {false};
    bool m_is_recolonized  {false};
    bool m_is_fire_event   {false};
    int m_total_fires      {0};

    std::vector<Cohort> m_cohort_list;

    Plant_stages_counter m_total_recolonizations;
    Plant_stages_tracker m_recolonizations_tracker;
    Plant_stages_counter m_total_immigrants;
    Plant_stages_tracker m_immigrants_tracker;
    Plant_stages_counter m_total_LDD_residents;
    Plant_stages_tracker m_LDD_residents_tracker;

    LDD_seeds_by_wind m_total_LDD_seeds_by_wind;
    LDD_seeds_by_birds m_total_LDD_seeds_by_birds;

    std::vector<Cell> m_cell_locations_for_LDD_birds;

    Point m_midpoint{0.0, 0.0};

    std::map<const int, const double> m_dist_edge_to_edge;
    std::map<const int, const double> m_dist_midpoint_to_edge;
    std::map<const int, const double> m_dist_midpoint_to_midpoint;
    // midpoint is the mean point (x and y) of all center cells of the population.

    ///@note NOTE should I save this variable as member variable?
    std::string m_habitat_quality;

//    std::vector<int> m_residents;
//    std::vector<int> m_residents_postfire;
//    std::vector<int> m_immigrants;
//    std::vector<int> m_immigrants_postfire;
    std::vector<int> m_genetic_pops_postfire;
    std::vector<int> m_genetic_pops;

    std::vector<int> m_genetic_pops_postfire_occupied;
    std::vector<int> m_genetic_pops_occupied;

    std::vector<double> m_immigrants_postfire_perc;
    std::vector<double> m_immigrants_perc;

    std::map<int, Total_plants_stage> m_plants_per_population;

//    std::vector<int> m_genetic_pops_postfire;
};

#endif // POPULATION_H
