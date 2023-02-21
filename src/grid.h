#ifndef GRID_H
#define GRID_H

#include <cassert>
#include <vector>
#include <map>

#include "cell.h"
#include "parameters.h"
#include "population.h"
#include "random_generator.h"
#include "dispersal.h"
#include "fire_simulator.h"

/**
 * @headerfile grid.h "src/grid.h"
 * @brief The Grid class hols the spatially-explicit processes of the model.
 *
*/
class Grid
{
public:

    void set_study_area(const Parameters& params);

    void locate_metapopulation(const unsigned seed, std::map<int, Population>& metapopulation);

    const std::vector<std::vector<int>>& get_study_area() const
    {
        return m_study_area;
    }

    void set_distance_btw_populations(const Parameters& params,
                                      std::map<int, Population>& metapopulation);

    void set_cells_for_LDD_by_birds(const Parameters& params,
                                    std::map<int, Population>& metapopulation);


    void set_postfire_dispersal_first_year(const Parameters& params, const std::map<int, Population>& metapopulation);

    void set_postfire_dispersal(const Parameters& params, const Population& population);

    bool is_empty_dispersal_list()
    {
        return m_dispersal_list.empty();
    }

    void establish_SDD_cohorts(std::map<int, Population>& metapopulation);
    void LDD_by_birds(const Parameters& params, std::map<int, Population>& metapopulation);
    void LDD_by_wind(const Parameters& params, std::map<int, Population>& metapopulation);

    void remove_dispersal_list() { m_dispersal_list.clear(); }


    void locate_patchy_fire(const Parameters& params);

    bool is_inside_ellipse(const Cell& center_cell,
                           const Cell& target_cell,
                           const double x_radius,
                           const double y_radius,
                           const double radians);

    void evaluate_is_population_burned(const Parameters& params, std::map<int, Population>& metapopulation);

private:

    bool is_population_within_patchy_fire(const Population& population);

    double get_fire_probability(const Parameters& params, const int time_since_fire);

    void reset_dispersal_matrix()
    {
        for (auto& row : m_dispersal_matrix)
        {
            std::fill(std::begin(row), std::end(row), 0);
        }
    }

    void establish_LDD_cohorts(std::map<int, Population>& metapopulation,
                               const Gene_flow gene_flow,
                               const Dispersal_vector dispersal_vector);


    // these to prints are to check the behaviour visually in debug mode.
    void print_study_area_console();
    void print_console_patchy_fire_matrix();

    bool is_valid_cell(const std::vector<std::vector<int>>& study_area,
                       const int pop_id,
                       const int target_row,
                       const int target_col);

    double distance_edge_to_edge(const Population& pop_a, const Population& pop_b);

    double distance_midpoint_to_edge(const Point& midpoint_pop_a, const Population& pop_b);


//    double euclidean_distance(const int row_a, const int col_a, const int row_b, const int col_b);

//    double euclidean_distance(const double row_a,
//                              const double col_a,
//                              const double row_b,
//                              const double col_b);

//    double euclidean_distance(const Point& point_a, const Point& point_b);

    //-----------------------------------------------------------------------------

    inline double euclidean_distance(const int row_a, const int col_a, const int row_b, const int col_b)
    {
        return std::sqrt(std::pow(row_a - row_b, 2) + std::pow(col_a - col_b, 2));
    }

    //-----------------------------------------------------------------------------

    inline double euclidean_distance(
            const double row_a, const double col_a, const double row_b, const double col_b)
    {
        return std::sqrt(std::pow(row_a - row_b, 2) + std::pow(col_a - col_b, 2));
    }

    //-----------------------------------------------------------------------------

    inline double euclidean_distance(const Point& point_a, const Point& point_b)
    {
        return std::sqrt(std::pow(point_a.row - point_b.row, 2) +
                         std::pow(point_a.col - point_b.col, 2));
    }

    // private member variables
    std::vector<std::vector<int>> m_study_area;
    std::vector<std::vector<int>> m_dispersal_matrix;
    std::vector<Dispersal> m_dispersal_list;

    /**
     * @brief m_patchy_fire_matrix is a replicate of the study area to locate fires
     * @warning vector<bool> is not a standard std container. Some features does
     * not work. For example, it cannot work with references in a for-range loop.
     */
    std::vector<std::vector<bool>> m_patchy_fire_matrix;
};

#endif // GRID_H
