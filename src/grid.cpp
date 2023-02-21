#include "grid.h"


//-----------------------------------------------------------------------------

void Grid::set_study_area(const Parameters& params)
{
    // Expects
    ///@note NOTE these variables are evaluated in parameter_reader class
    assert(params.cell_size > 0);
    assert(params.study_size_x % params.cell_size == 0);
    assert(params.study_size_y % params.cell_size == 0);

    const unsigned rows {
        static_cast<unsigned>(params.study_size_y / params.cell_size)
    };
    const unsigned cols {
        static_cast<unsigned>(params.study_size_x / params.cell_size)
    };

    m_study_area.clear();
    m_study_area.assign(rows, std::vector<int>(cols, 0));
    m_dispersal_matrix.assign(rows, std::vector<int>(cols, 0));
    m_patchy_fire_matrix.assign(rows, std::vector<bool>(cols, false));

    // Ensures
    assert(std::size(m_study_area)    > 0);
    assert(std::size(m_study_area[0]) > 0);
    assert(std::size(m_study_area)    == std::size(m_dispersal_matrix));
    assert(std::size(m_study_area[0]) == std::size(m_dispersal_matrix[0]));
    assert(std::size(m_study_area)    == std::size(m_patchy_fire_matrix));
    assert(std::size(m_study_area[0]) == std::size(m_patchy_fire_matrix[0]));
//    print_study_area_console();
}

//-----------------------------------------------------------------------------

void Grid::print_study_area_console()
{
    for (unsigned row = 0; row < std::size(m_study_area); ++row)
    {
        for (unsigned col = 0; col < std::size(m_study_area[0]); ++col)
        {
            std::cout << m_study_area[row][col] << ' ';
        }
        std::cout << '\n';
    }

    std::cout << std::endl;
    std::cout << "----------------------------------------------------\n";
}

//-----------------------------------------------------------------------------

void Grid::print_console_patchy_fire_matrix()
{
    for (unsigned row = 0; row < std::size(m_patchy_fire_matrix); ++row)
    {
        for (unsigned col = 0; col < std::size(m_patchy_fire_matrix[0]); ++col)
        {
            std::cout << m_patchy_fire_matrix[row][col] << ' ';
        }
        std::cout << '\n';
    }

    std::cout << std::endl;
    std::cout << "----------------------------------------------------\n";
}

//-----------------------------------------------------------------------------

void Grid::locate_metapopulation(const unsigned seed, std::map<int, Population>& metapopulation)
{
    Random_generator prng;
    prng.set_seed(seed);

    int count_attemps {0};

    // to access the input populations file randomly
    std::vector<int> keys_metapop_random;

    for (const auto& [key, val] : metapopulation)
    {
        keys_metapop_random.push_back(key);
    }

    // mix randomly the population order of the population input file
    std::shuffle(std::begin(keys_metapop_random),
                 std::end(keys_metapop_random),
                 prng.get_engine());

    assert(!std::empty(m_study_area));
    assert(std::size(m_study_area[0]) > 0);

    const auto max_rows = static_cast<unsigned>(std::size(m_study_area) - 1);
    const auto max_cols = static_cast<unsigned>(std::size(m_study_area[0]) - 1);

    for (const auto& key_random : keys_metapop_random)
    {
        auto& population = metapopulation.at(key_random);

        bool is_patch_successful {false};

        while (!is_patch_successful)
        {
            unsigned row_random {0};
            unsigned col_random {0};

            do
            {
                ++count_attemps;

                ///@todo TODO I increased the number of attempts. I should try to
                /// find a better approach??
                if (count_attemps > 10000000)
                {
                    std::cout << "ERROR! There are not enough free spaces"
                                 " (cells) in the matrix to place the populations"
                                 ".\nPlease, reduce the number of populations or "
                                 "their sizes, or increase the size of the study area.\n\a";
                    std::exit(EXIT_FAILURE);
                }

                row_random = prng.uniform_unsigned(0, max_rows);
                col_random = prng.uniform_unsigned(0, max_cols);

            } while (m_study_area[row_random][col_random] != 0);

            std::vector<std::vector<int>> study_area_temp {m_study_area};
            study_area_temp[row_random][col_random] = key_random;

            int row_current {static_cast<int>(row_random)};
            int col_current {static_cast<int>(col_random)};
            int row_next    {0};
            int col_next    {0};

            if (is_valid_cell(study_area_temp, key_random, row_current, col_current) == false)
            {
                continue;
            }

            std::vector<Cell> cell_locations;
            cell_locations.emplace_back(Cell(row_random, col_random));

            int counter_step_attempts  {0};
//            bool is_cell_empty         {true};
            constexpr int max_attempts {500};

            is_patch_successful = true;    // this is necessary here to be able to
                                           // place populations with only one cell

            const unsigned population_num_cells {static_cast<unsigned>(population.get_num_cells())};

            // starts random walk
            // i starts with 1 (i = 1) to avoid entering with pops with size 1
            for (unsigned i = 1; i < population_num_cells; ++i)
            {
                is_patch_successful = false;
                bool is_cell_empty {false};

                if (counter_step_attempts >= max_attempts)
                {
                    is_patch_successful = false;
                    break;
                }

                counter_step_attempts = 0;

                while (!is_cell_empty && counter_step_attempts < max_attempts)
                {
                    ++counter_step_attempts;

                    //---------------------------------------
                    double x {static_cast<double>(col_current)};
                    double y {static_cast<double>(row_current)};

                    // 3. get random angle
//                    double angle_degrees {Random_generator::uniform_real(0.0, 360)};
//                    double angle_radians {angle_degrees * s_kPi / 180.0};
                    double angle_radians {prng.uniform_real(0.0, 2.0 * Global::pi)};

                    // 6. calculate x and y at the given angle and distance
                    x = x + std::sin(angle_radians) + 0.5;
                    y = y - std::cos(angle_radians) + 0.5;

                    ///@note NOTE I have to check in LDD by wind!!
                    /// why I evaluate if x and y are below -1 ???
                    // 7. check if x and y are within the study area
                    if (x < 0 || x >= static_cast<double>(std::size(m_study_area[0])) ||
                        y < 0 || y >= static_cast<double>(std::size(m_study_area)))
                    {
                        continue;
                    }

                    row_next = static_cast<int>(y);
                    col_next = static_cast<int>(x);

                    if (!is_valid_cell(study_area_temp, key_random, row_next, col_next))
                    {
                        continue;
                    }

                    unsigned row {static_cast<unsigned>(row_next)};
                    unsigned col {static_cast<unsigned>(col_next)};

                    if (study_area_temp.at(row).at(col) == key_random)
                    {
                        row_current = row_next;
                        col_current = col_next;
                    }
                    else if (study_area_temp.at(row).at(col) == 0)
                    {
                        row_current = row_next;
                        col_current = col_next;

                        is_cell_empty = true;

                        Cell cell {
                            static_cast<unsigned>(row_current),
                            static_cast<unsigned>(col_current)
                        };

                        cell_locations.emplace_back(cell);

                        study_area_temp.at(row).at(col) = key_random;
                    }
                    else
                    {
                        ///@note NOTE this should be assert??
                        std::cerr << "error in grid.cpp, random walk\n\a";
                        std::exit(EXIT_FAILURE);
                    }
                }
            }

            if (std::size(cell_locations) == population_num_cells)
            {
                is_patch_successful = true;
                m_study_area = study_area_temp;
                population.set_cell_locations(cell_locations);
            }
        }
    }
//    print_study_area_console();
}

//-----------------------------------------------------------------------------

bool Grid::is_valid_cell(const std::vector<std::vector<int>>& study_area,
                         const int pop_id,
                         const int target_row,
                         const int target_col)
{
    for (int i = -1; i < 2; ++i)
    {
        for (int j = -1; j < 2; ++j)
        {
//            // skip center cell
//            if (i == j)
//            {
//                continue;
//            }

            // skip rows out of range.
            if ((i + target_row) < 0 ||
                (i + target_row >= static_cast<int>(std::size(study_area))))
            {
                 continue;
            }

            // skip columns out of range.
            if ((j + target_col) < 0 ||
                (j + target_col >= static_cast<int>(std::size(study_area[0]))))
            {
                continue;
            }

            unsigned row {static_cast<unsigned>(i + target_row)};
            unsigned col {static_cast<unsigned>(j + target_col)};

            if (study_area.at(row).at(col) == 0 ||
                study_area.at(row).at(col) == pop_id)
            {
                continue;
            }
            else
            {
                return false;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------

void Grid::set_distance_btw_populations(const Parameters& params,
                                        std::map<int, Population>& metapopulation)
{
    assert(!std::empty(m_study_area));
    assert(!std::empty(metapopulation));

    // calculate midpoint for each population
    for (auto& [key, val] : metapopulation)
    {
        std::vector<double> row_values;
        std::vector<double> col_values;

        for (const auto& cell : val.get_cell_locations())
        {
            row_values.push_back(static_cast<double>(cell.row));
            col_values.push_back(static_cast<double>(cell.col));
        }

        const auto midpoint = Point(Equation::mean(row_values) + 0.5
                                    , Equation::mean(col_values) + 0.5);

        val.set_midpoint(midpoint);
    }

    // midpoint to edge
    for (auto& [key, val] : metapopulation)
    {
        for (auto& [key2, val2] : metapopulation)
        {
            if (key == key2) {
                continue;
            }

            double minimum_distance {
                distance_midpoint_to_edge(val.get_midpoint(), val2)
            };

            minimum_distance *= params.cell_size;
            val.add_dist_midpoint_to_edge(key2, minimum_distance);
        }
    }

    // distance edge to edge, and distance midpoint to midpoint
    std::vector<int> track_id;

    for (auto& [key, val] : metapopulation)
    {
        track_id.push_back(key);

        for (auto& [key2, val2] : metapopulation)
        {
            const auto result = std::find(std::cbegin(track_id)
                                          , std::cend(track_id)
                                          , key2);

            // skip population is already calculated because distance
            // edge to edge and midpoint to midpoint
            // is the same between populations
            if (result != std::cend(track_id)) {
                continue;
            }

            double minimum_distance { distance_edge_to_edge(val, val2) };
            minimum_distance *= params.cell_size;

            val.add_dist_edge_to_edge(key2, minimum_distance);
            val2.add_dist_edge_to_edge(key, minimum_distance);

            minimum_distance = euclidean_distance(val.get_midpoint(),
                                                  val2.get_midpoint());
            minimum_distance *= params.cell_size;

            val.add_dist_midpoint_to_midpoint(key2, minimum_distance);
            val2.add_dist_midpoint_to_midpoint(key, minimum_distance);
        }
    }
}

//-----------------------------------------------------------------------------

double Grid::distance_edge_to_edge(const Population& pop_a,
                                   const Population& pop_b)
{
    double minimum_distance {0.0};
    double current_distance {0.0};
    bool is_first_dist      {true};

    for (const auto& cell_a : pop_a.get_cell_locations())
    {
        const double row_a { static_cast<double>(cell_a.row) };
        const double col_a { static_cast<double>(cell_a.col) };

        for (const auto& cell_b : pop_b.get_cell_locations())
        {
            const double row_b { static_cast<double>(cell_b.row) };
            const double col_b { static_cast<double>(cell_b.col) };

            // check the four vertices of the cells (i.e. shift 0/1)
            for (int row_a_shift = 0; row_a_shift < 2; ++row_a_shift)
            {
                for (int col_a_shift = 0; col_a_shift < 2; ++col_a_shift)
                {
                    for (int row_b_shift = 0; row_b_shift < 2; ++row_b_shift)
                    {
                        for (int col_b_shift = 0; col_b_shift < 2; ++col_b_shift)
                        {
                            current_distance = euclidean_distance(
                                        row_a + row_a_shift
                                        , col_a + col_a_shift
                                        , row_b + row_b_shift
                                        , col_b + col_b_shift);

                            if (is_first_dist)
                            {
                                minimum_distance = current_distance;
                                is_first_dist = false;
                            }
                            else if (minimum_distance > current_distance)
                            {
                                minimum_distance = current_distance;
                            }
                        }
                    }
                }
            }
        }
    }

    // Ensures
    assert(minimum_distance >= 0);

    return minimum_distance;
}

//-----------------------------------------------------------------------------

double Grid::distance_midpoint_to_edge(const Point& midpoint_pop_a, const Population& pop_b)
{
    double minimum_distance {0.0};
    double current_distance {0.0};
    bool is_first_dist      {true};

    for (const auto& cell_b : pop_b.get_cell_locations())
    {
        const double row_b { static_cast<double>(cell_b.row) };
        const double col_b { static_cast<double>(cell_b.col) };

        for (int row_b_shift = 0; row_b_shift < 2; ++row_b_shift)
        {
            for (int col_b_shift = 0; col_b_shift < 2; ++col_b_shift)
            {
                current_distance = euclidean_distance(midpoint_pop_a.row
                                                      , midpoint_pop_a.col
                                                      , row_b + row_b_shift
                                                      , col_b + col_b_shift);

                if (is_first_dist)
                {
                    minimum_distance = current_distance;
                    is_first_dist = false;
                }
                else if (minimum_distance > current_distance)
                {
                    minimum_distance = current_distance;
                }
            }
        }
    }

    // Ensures
    assert(minimum_distance >= 0.0);

    return minimum_distance;
}

//-----------------------------------------------------------------------------

void Grid::set_cells_for_LDD_by_birds(const Parameters& params, std::map<int, Population>& metapopulation)
{
    assert(!std::empty(m_study_area));
    assert(!std::empty(metapopulation));

    for (unsigned row {0}; row < std::size(m_study_area); ++row)
    {
        for (unsigned col {0}; col < std::size(m_study_area[0]); ++col)
        {
            // each suitable cell is evaluated if it is within the maximum LDD
            // by birds for each population.
            if (m_study_area[row][col] > 0)    // only suitable cells
            {
                const double target_row { static_cast<double>(row) };
                const double target_col { static_cast<double>(col) };

                for ([[maybe_unused]] auto& [key, val] : metapopulation)
                {
                    double minimum_distance {0.0};
                    bool is_first_distance  {true};

                    for (const auto& cell : val.get_cell_locations())
                    {
                        const double cell_row { static_cast<double>(cell.row) };
                        const double cell_col { static_cast<double>(cell.col) };

                        for (int target_row_shift {0};
                             target_row_shift < 2;
                             ++target_row_shift)
                        {
                            for (int target_col_shift {0};
                                 target_col_shift < 2;
                                 ++target_col_shift)
                            {
                                for (int cell_row_shift {0};
                                     cell_row_shift < 2;
                                     ++cell_row_shift)
                                {
                                    for (int cell_col_shift {0};
                                         cell_col_shift < 2;
                                         ++cell_col_shift)
                                    {
                                        double distance {
                                            euclidean_distance(
                                                 target_row + target_row_shift,
                                                 target_col + target_col_shift,
                                                 cell_row + cell_col_shift,
                                                 cell_col + cell_col_shift)
                                        };

                                        if (is_first_distance)
                                        {
                                            minimum_distance = distance;
                                            is_first_distance = false;
                                        }
                                        else if (minimum_distance > distance)
                                        {
                                            minimum_distance = distance;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    minimum_distance *= params.cell_size;

                    if (minimum_distance <= params.birds_dist_max)
                    {
                        val.add_cell_for_LDD_by_birds(Cell(row, col));
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

void Grid::set_postfire_dispersal_first_year(const Parameters& params, const std::map<int, Population>& metapopulation)
{
    for ([[maybe_unused]] const auto& [key, val] : metapopulation)
    {
        for (const auto& cohort : val.get_cohort_list())
        {
            //-----------------------------------------------------
            // this is the only part that differs from set_postfire_dispersal(...)
            const int total_LDD_cones_by_birds {
                static_cast<int>(params.init_cones_plant
                                 * cohort.get_num_plants()
                                 * params.birds_prop
                                 + 0.5)
            };

            const int total_viable_seeds {
                static_cast<int>(params.init_viable_seeds_plant
                                 * cohort.get_num_plants()
                                 + 0.5)
            };
            //-----------------------------------------------------

            if (total_viable_seeds > 0)
            {
                m_dispersal_list.emplace_back(
                            Dispersal(cohort.get_gene_flow()
                                      , total_LDD_cones_by_birds
                                      , total_viable_seeds));
            }
        }
    }
}


//-----------------------------------------------------------------------------

void Grid::set_postfire_dispersal(const Parameters& params, const Population& population)
{
    for (const auto& cohort : population.get_cohort_list())
    {
        //-----------------------------------------------------------------
        // this is the only part that differs from set_postfire_dispersal_first_year(...)
        const int total_LDD_cones_by_birds {
            static_cast<int>(
                        (cohort.get_seedbank_cones(params) *
                        params.birds_prop) + 0.5)
        };

        const int total_viable_seeds { cohort.get_seedbank(params) };
        //-----------------------------------------------------------------

        if (total_viable_seeds > 0)
        {
            m_dispersal_list.emplace_back(
                        Dispersal(cohort.get_gene_flow()
                                  , total_LDD_cones_by_birds
                                  , total_viable_seeds));
        }
    }
}

//-----------------------------------------------------------------------------

void Grid::establish_SDD_cohorts(std::map<int, Population>& metapopulation)
{
    constexpr bool is_LDD_cohort {false};

    for (const auto& dispersal : m_dispersal_list)
    {
        if (dispersal.m_total_viable_seeds > 0)
        {
            auto& population = metapopulation.at(
                                    dispersal.m_gene_flow.get_current());

            auto gene_flow = dispersal.m_gene_flow;

            gene_flow.set_current(population.get_id());

            bool is_postfire_cohort {true};

            if (population.get_time_since_fire() > 0)
            {
                is_postfire_cohort = false;
            }

            population.add_new_cohort(dispersal.m_total_viable_seeds
                                      , gene_flow
                                      , is_postfire_cohort
                                      , is_LDD_cohort
                                      , Dispersal_vector::SDD_wind);
        }
    }
}

//-----------------------------------------------------------------------------

void Grid::LDD_by_birds(const Parameters& params, std::map<int, Population>& metapopulation)
{
    for (auto& dispersal : m_dispersal_list)
    {
        if (dispersal.m_total_viable_seeds <= 0)
        {
            continue;
        }

        reset_dispersal_matrix();

        auto& population = metapopulation.at(dispersal.m_gene_flow.get_current());

        for (int i {0}; i < dispersal.m_LDD_cones_by_birds; ++i)
        {
            if (dispersal.m_total_viable_seeds <= 0)
            {
                break;
            }

            int remaining_seeds_in_cone {0};

            {
                double seeds_in_cone {0.0};

                // 1. calculate number of follicles within cone
                switch (params.follicles_distr_type)
                {
                case Distribution_type::null:
                    throw std::runtime_error(
                                "The parameter 'follicles_distr_type' is not valid:\n"
                        "- value '1' for Poisson distribution\n"
                        "- value '3' for negative binomial with 'mu' parameter\n");
                    break;
                case Distribution_type::pois:
                    seeds_in_cone = Global::prng.poisson(params.follicles_distr_a);
                    break;
                case Distribution_type::nbinom:
                    throw std::runtime_error(
                                "The parameter 'follicles_distr_type' is not valid:\n"
                        "- value '1' for Poisson distribution\n"
                        "- value '3' for negative binomial with 'mu' parameter\n");
                    break;
                case Distribution_type::nbinom_mu:
                    seeds_in_cone = Global::prng.negative_binomial_mu(
                                params.follicles_distr_a
                                , params.follicles_distr_b);
                    break;
                case Distribution_type::geom:
                    throw std::runtime_error(
                                "The parameter 'follicles_distr_type' is not valid:\n"
                        "- value '1' for Poisson distribution\n"
                        "- value '3' for negative binomial with 'mu' parameter\n");
                    break;
                }

                assert(seeds_in_cone >= 0.0);

                // 2. potential seeds in cone
                seeds_in_cone *= params.num_seeds * params.firm_seeds;

                // 3. seeds lost in 1-year-old cone
                const double seeds_lost {
                    seeds_in_cone * (params.insect_cone_age.at(params.cone_cycle)
                                     + params.decay_cone_age.at(params.cone_cycle))
                };

                seeds_in_cone -= seeds_lost;

                if (seeds_in_cone < 0)
                {
                    seeds_in_cone = 0;
                }

                // 4. proportion of follicles open
                double follicles_open {
                    Global::prng.uniform_real(params.postfire_follicles_open, 1.0)
                };

                assert(follicles_open >= 0.0 && follicles_open <= 1.0);

                seeds_in_cone *= (1.0 - follicles_open);

                // 5. number of remaining viable seeds within 1-year-old cones
                remaining_seeds_in_cone = static_cast<int>(
                            seeds_in_cone * params.viable_seeds + 0.5);
            }

            // 6. disperse cone
            if (remaining_seeds_in_cone > 0)
            {
                if (dispersal.m_total_viable_seeds >= remaining_seeds_in_cone)
                {
                    dispersal.m_total_viable_seeds -= remaining_seeds_in_cone;
                }
                else if (dispersal.m_total_viable_seeds > 0)
                {
                    remaining_seeds_in_cone = dispersal.m_total_viable_seeds;
                    dispersal.m_total_viable_seeds = 0;
                }
                else
                {
                    remaining_seeds_in_cone = 0;
                    dispersal.m_total_viable_seeds = 0;
                }

                if (remaining_seeds_in_cone > 0)
                {
                    // effective dispersal
                    const auto cell = population.get_cell_random_LDD_birds();

                    m_dispersal_matrix[cell.row][cell.col] = remaining_seeds_in_cone;

                    if (params.is_out_dispersal)
                    {
                        if (population.get_id() == m_study_area.at(cell.row).at(cell.col))
                        {
                            population.increment_total_LDD_seeds_by_birds_resident(remaining_seeds_in_cone);
                        }
                        else
                        {
                            population.increment_total_LDD_seeds_by_birds_immigrant(remaining_seeds_in_cone);
                        }
                    }
                }
            }
        }

        establish_LDD_cohorts(metapopulation
                              , dispersal.m_gene_flow
                              , Dispersal_vector::LDD_bird);
    }
}

//-----------------------------------------------------------------------------

void Grid::LDD_by_wind(const Parameters& params, std::map<int, Population>& metapopulation)
{
    // Expects
    assert(params.wind_prop > 0.0 && params.wind_prop <= 1.0);
    assert(std::size(m_study_area) > 0);
    assert(std::size(m_study_area[0]) > 0);
    assert(params.cell_size > 0);

    for (auto& dispersal : m_dispersal_list)
    {
        if (dispersal.m_total_viable_seeds <= 0)
        {
            continue;
        }

        const int total_LDD_seeds_by_wind {
            static_cast<int>(dispersal.m_total_viable_seeds
                             * params.wind_prop
                             + 0.5)
        };

        if (dispersal.m_total_viable_seeds < total_LDD_seeds_by_wind)
        {
            std::cout << "ERROR: it might be that 'wind_prop' is greater than 1\n";
            std::exit(EXIT_FAILURE);
        }

        dispersal.m_total_viable_seeds -= total_LDD_seeds_by_wind;

        reset_dispersal_matrix();

        auto& population = metapopulation.at(dispersal.m_gene_flow.get_current());

        const auto pop_id = population.get_id();

        for (int i {0}; i < total_LDD_seeds_by_wind; ++i)
        {
            // 1. select a grid of population at random
            const auto cell_pop = population.get_random_cell_location();

            auto x = static_cast<double>(cell_pop.col);
            auto y = static_cast<double>(cell_pop.row);

            // 2. truncate the selected cell of the population
            x = x + Global::prng.uniform_01();
            y = y + Global::prng.uniform_01();

            // 3. get random angle
            const double angle_radians {
                Global::prng.uniform_real(0.0, 2.0 * Global::pi)
            };

            ///@note NOTE the commented code below is for the case it is added
            /// the wind direction parameter.
//            double angle_degrees {0.0};
//            double angle_radians {0.0};

//            switch (params.wind_direction)
//            {
//                case Wind_direction::random:
////                    angle_degrees = Global::prng.uniform_real(0.0, 360.0);
//                    angle_radians = Global::prng.uniform_real(
//                                        0.0, Global::pi * 2.0);
//                    break;
//                case Wind_direction::geraldton:
////                    angle_degrees = Wind_simulator::get_wind_direction();
//                    std::cout << "I do not want to use Geralton wind direction\n";
//                    std::exit(EXIT_FAILURE);
//                    break;
//            }

//            const double angle_radians { angle_degrees * Global::pi / 180.0 };

            // 4. get random distance
            double LDD_random {
                Global::prng.lognormal(params.wind_a, params.wind_b)
            };

            // 5. transform distance to unit cells
            LDD_random = LDD_random / params.cell_size;

            // 6. calculate x and y at the given angle and distance
            x = x + (LDD_random * std::sin(angle_radians));
            y = y - (LDD_random * std::cos(angle_radians));

            // 7. check if x and y are within the study area
            if (y < 0.0 || y >= static_cast<double>(std::size(m_study_area)) ||
                x < 0.0 || x >= static_cast<double>(std::size(m_study_area[0])))
            {
                if (params.is_out_dispersal)
                {
                    population.increment_total_LDD_seeds_by_wind_outside();
                }
            }
            else
            {                
                const unsigned col { static_cast<unsigned>(x) };
                const unsigned row { static_cast<unsigned>(y) };

                if (m_study_area.at(row).at(col) > 0)
                {
                    ++m_dispersal_matrix[row][col];
                }

                if (params.is_out_dispersal)
                {
                    if (m_study_area[row][col] == 0)
                    {
                        population.increment_total_LDD_seeds_by_wind_unsuitable();
                    }
                    else
                    {
                        if (m_study_area[row][col] == pop_id)
                        {
                            population.increment_total_LDD_seeds_by_wind_resident();
                        }
                        else
                        {
                            population.increment_total_LDD_seeds_by_wind_immigrant();
                        }
                    }
                }

            }
        }
        establish_LDD_cohorts(metapopulation
                              , dispersal.m_gene_flow
                              , Dispersal_vector::LDD_wind);
    }
}

//-----------------------------------------------------------------------------

void Grid::establish_LDD_cohorts(
        std::map<int, Population>& metapopulation
        , const Gene_flow gene_flow
        , const Dispersal_vector dispersal_vector)
{
    constexpr bool is_LDD_cohort {true};

    for (auto& [key, val] : metapopulation)
    {
        int LDD_seeds {0};

        for (const auto& cell : val.get_cell_locations())
        {
            LDD_seeds += m_dispersal_matrix[cell.row][cell.col];
        }

        if (LDD_seeds > 0)
        {
            bool is_postfire_cohort {true};

            if (val.get_time_since_fire() > 0)
            {
                is_postfire_cohort = false;
            }

            auto gene_flow_new = gene_flow;
            gene_flow_new.set_current(key);

            val.add_new_cohort(LDD_seeds
                               , gene_flow_new
                               , is_postfire_cohort
                               , is_LDD_cohort
                               , dispersal_vector);
        }
    }
}

//-----------------------------------------------------------------------------

void Grid::locate_patchy_fire(const Parameters& params)
{
    assert(params.fire_scale == Fire_scale::patchy);
    assert(std::empty(m_patchy_fire_matrix) == false);
    assert(std::empty(m_patchy_fire_matrix[0]) == false);

    for (auto& row : m_patchy_fire_matrix)
    {
        std::fill(std::begin(row), std::end(row), false);
    }

    double fire_size_x {0.0};
    double fire_size_y {0.0};

    switch (params.fire_size_scn)
    {
        case Fire_scenario::deterministic:
            fire_size_x = params.fire_size_x;
            fire_size_y = params.fire_size_y;
            break;
        case Fire_scenario::truncated_normal:
            fire_size_x = Fire_simulator::get_fire_size_truncated_normal(
                        params.fire_size_x, 0.0);
            fire_size_y = Fire_simulator::get_fire_size_truncated_normal(
                        params.fire_size_y, 0.0);
            break;
        case Fire_scenario::weibull:
            // there is no weibull distribution
            assert(false);
            break;
    }

    const Cell center_cell {
        Global::prng.uniform_unsigned(
                    0,
                    static_cast<unsigned>(std::size(m_patchy_fire_matrix)) - 1),
        Global::prng.uniform_unsigned(
                    0,
                    static_cast<unsigned>(std::size(m_patchy_fire_matrix[0])) - 1)
    };

    const double angle_radians {
        Global::prng.uniform_real(0.0, 2.0 * Global::pi)
    };

    for (unsigned row {0}; row < std::size(m_patchy_fire_matrix); ++row)
    {
        for (unsigned col {0}; col < std::size(m_patchy_fire_matrix[0]); ++col)
        {
            const Cell target_cell {row, col};

            m_patchy_fire_matrix[row][col] = is_inside_ellipse(center_cell
                                                               , target_cell
                                                               , fire_size_x
                                                               , fire_size_y
                                                               , angle_radians);
        }
    }
//    print_console_patchy_fire_matrix();
}

//-----------------------------------------------------------------------------

bool Grid::is_inside_ellipse(const Cell& center_cell,
                             const Cell& target_cell,
                             const double x_radius,
                             const double y_radius,
                             const double radians)
{
    const double k { static_cast<double>(center_cell.row) };
    const double h { static_cast<double>(center_cell.col) };
    const double y { static_cast<double>(target_cell.row) };
    const double x { static_cast<double>(target_cell.col) };

    return ((std::pow((x - h) * std::cos(radians) + (y - k) * std::sin(radians), 2)) /
            std::pow(x_radius, 2)) +
            ((std::pow((x - h) * std::sin(radians) - (y - k) * std::cos(radians), 2)) /
             std::pow(y_radius, 2)) <= 1.0;

}

//-----------------------------------------------------------------------------

void Grid::evaluate_is_population_burned(
        const Parameters& params
        , std::map<int, Population>& metapopulation)
{
    for ([[maybe_unused]] auto& [key, val] : metapopulation)
    {
        bool is_fire_event {false};

        switch (params.fire_scale)
        {
            case Fire_scale::patchy:
                is_fire_event = is_population_within_patchy_fire(val);
                break;
            case Fire_scale::study_area:
                is_fire_event = true;
                break;
        }

        if (is_fire_event)
        {
            const int time_since_fire {val.get_time_since_fire()};

            if (time_since_fire <= params.burned_lower_cut)
            {
                is_fire_event = false;
            }
            else if (time_since_fire < params.burned_upper_cut)
            {
                const double fire_probability {
                    get_fire_probability(params, time_since_fire)
                };

                if (Global::prng.uniform_01() > fire_probability)
                {
                    is_fire_event = false;
                }
            }
            else
            {
                is_fire_event = true;
            }
        }
        val.set_is_fire_event(is_fire_event);
    }
}

//-----------------------------------------------------------------------------

bool Grid::is_population_within_patchy_fire(const Population& population)
{
    for (const auto& cell : population.get_cell_locations())
    {
        if (m_patchy_fire_matrix.at(cell.row).at(cell.col))
        {
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------

double Grid::get_fire_probability(const Parameters& params, const int time_since_fire)
{
    // Expects
    assert(time_since_fire >= params.burned_lower_cut);

    const double tsf       { static_cast<double>(time_since_fire) };
    const double lower_cut { static_cast<double>(params.burned_lower_cut) };
    const double upper_cut { static_cast<double>(params.burned_upper_cut) };

    assert((upper_cut - lower_cut) > 0);

    double fire_prob {(tsf - lower_cut) / (upper_cut - lower_cut)};

    // Ensures
    assert(fire_prob >= 0);

    return fire_prob;
}
