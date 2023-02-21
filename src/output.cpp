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

#include "output.h"

//-----------------------------------------------------------------------------

Output::Output(
        const std::string& output_folder
        , const std::string& name_sim_file
        , const std::string& time)

    : m_output_folder{Global::project_directory + "data/out/" + output_folder}
    , m_name_sim_file{name_sim_file}
    , m_time{time}
{}

//-----------------------------------------------------------------------------

Output::~Output()
{
    cleanup();
}

//-----------------------------------------------------------------------------

void Output::setup_output(const Parameters& params
                          , const std::string& name_sim_folder)
{
    const std::string simulations_folder {"/in/sim/"};
    const std::string plots_folder {"/plots/"};

    if (!std::filesystem::exists(m_output_folder))
    {
        std::filesystem::create_directories(m_output_folder);
    }

    if (!std::filesystem::exists(m_output_folder + m_raw_folder))
    {
        std::filesystem::create_directories(m_output_folder + m_raw_folder);
    }

    if (!std::filesystem::exists(m_output_folder + plots_folder))
    {
        std::filesystem::create_directories(m_output_folder + plots_folder);
    }

    if (!std::filesystem::exists(m_output_folder + simulations_folder))
    {
        std::filesystem::create_directories(m_output_folder + simulations_folder);
    }

    if (!std::filesystem::exists(m_output_folder + simulations_folder + m_name_sim_file))
    {
        std::filesystem::copy_file(
                    Global::project_directory
                    + "data/in/sim/"
                    + name_sim_folder
                    + '/'
                    + m_name_sim_file
                    , m_output_folder + simulations_folder + m_name_sim_file
                    , std::filesystem::copy_options::skip_existing
        );
    }

    const std::filesystem::path name_sim_file_copy {m_name_sim_file};

    std::filesystem::path root_file_name {m_output_folder + m_raw_folder};
    root_file_name /= name_sim_file_copy.stem();
    root_file_name += '_';
    root_file_name += m_time;
    root_file_name += "_sim_";
    root_file_name += std::to_string(params.sim_id);

    if (params.is_out_fire)
    {
        assert(std::empty(m_fire_output_path));
        m_fire_output_path.assign(root_file_name);
        m_fire_output_path += "_fire.csv";

        check_file_exists(m_fire_output_path);
        m_fire_output_stream.open(m_fire_output_path, std::ios_base::app);

        check_file_is_open(m_fire_output_stream, m_fire_output_path);

        const std::vector<std::string> header_fire_output {
            "sim_id", "study_num", "run_num", "year",
            "burned_pops", "burned_pops_perct",
            "burned_area", "burned_area_perct"
        };

        print_header(header_fire_output, m_fire_output_stream);
    }

    if (params.is_out_seed_dynamics)
    {
        assert(std::empty(m_seed_dynamics_path));
        m_seed_dynamics_path.assign(root_file_name);
        m_seed_dynamics_path += "_seed_dynamics.csv";
        check_file_exists(m_metapop_dynamics_path);
        m_seed_dynamics_stream.open(m_seed_dynamics_path, std::ios_base::app);
        check_file_is_open(m_seed_dynamics_stream, m_seed_dynamics_path);

        const std::vector<std::string> header_seed_dynamics {
            "sim_id", "study_num", "run_num", "year",
            "pop", "tsf",
            "coh_age", "coh_LDD", "coh_postfire",
            "num_plants", "cones_plant", "seeds_plant",
            "total_cones", "total_seeds"
        };

        print_header(header_seed_dynamics, m_seed_dynamics_stream);
    }

    if (params.is_out_metapopulation)
    {
        assert(std::empty(m_metapop_path));
        m_metapop_path.assign(root_file_name);
        m_metapop_path += "_metapop.csv";

        check_file_exists(m_metapop_path);
        m_metapop_stream.open(m_metapop_path, std::ios_base::app);

        check_file_is_open(m_metapop_stream, m_metapop_path);

        const std::vector<std::string> header_metapop {
            "sim_id", "study_num", "run_num", "year_max",
            "persis_max", "persis_mean", "persis_sd",
            "occu_max", "occu_mean", "occu_sd",
            "gen_max", "gen_mean", "gen_sd"
        };

        print_header(header_metapop, m_metapop_stream);
    }

    if (params.is_out_metapopulation_full)
    {
        assert(std::empty(m_metapop_full_path));
        m_metapop_full_path.assign(root_file_name);
        m_metapop_full_path += "_metapop_full.csv";

        check_file_exists(m_metapop_full_path);
        m_metapop_full_stream.open(m_metapop_full_path, std::ios_base::app);

        check_file_is_open(m_metapop_full_stream, m_metapop_full_path);

        const std::vector<std::string> header_metapop_full {
            "sim_id", "study_num", "run_num", "year_max",
            "persis_max", "persis_mean", "persis_sd",
            "occu_max", "occu_mean", "occu_sd",
            "gen_max", "gen_mean", "gen_sd",
            "pop_mean", "pop_sd",
//            "pop_fire_mean", "pop_fire_sd",
            "immi_mean", "immi_sd", "resi_mean", "resi_sd",
            "immi_perc_mean", "immi_perc_sd", "resi_perc_mean", "resi_perc_sd",
//            "immi_fire_mean", "immi_fire_sd", "resi_fire_mean", "resi_fire_sd",
//            "immi_fire_perc_mean", "immi_fire_perc_sd",
//            "resi_fire_perc_mean", "resi_fire_perc_sd"
        };

        print_header(header_metapop_full, m_metapop_full_stream);
    }

    if (params.is_out_metapopulation_dynamics)
    {
        assert(std::empty(m_metapop_dynamics_path));
        m_metapop_dynamics_path.assign(root_file_name);
        m_metapop_dynamics_path += "_metapop_dynamics.csv";
        check_file_exists(m_metapop_dynamics_path);
        m_metapop_dynamics_stream.open(m_metapop_dynamics_path, std::ios_base::app);
        check_file_is_open(m_metapop_dynamics_stream, m_metapop_dynamics_path);

        const std::vector<std::string> metapopulation_dynamics_header {
            "sim_id", "study_num", "run_num", "year",
            "tsf", "burned_pop",
            "existing_pops", "total_plants",
            "immi", "immi_bird", "immi_wind",
            "immi_perc",
            "pops_mean", "pops_sd"
//            "mean_diversity", "sd_diversity"
        };

        print_header(metapopulation_dynamics_header, m_metapop_dynamics_stream);
    }

    if (params.is_out_metapopulation_pops)
    {
        assert(std::empty(m_metapop_pops_path));
        m_metapop_pops_path.assign(root_file_name);
        m_metapop_pops_path += "_metapop_pops.csv";

        check_file_exists(m_metapop_pops_path);
        m_metapop_pops_stream.open(m_metapop_pops_path, std::ios_base::app);

        check_file_is_open(m_metapop_pops_stream, m_metapop_pops_path);

//        const std::vector<std::string> header_metapop_pops {
//            "sim_id", "study_num", "run_num", "year_max",
//            "pop_id", "burned", "extinctions",
//            "persistence", "generations",
//            "pops_postfire_mean", "pops_postfire_sd",
//            "pops_mean", "pops_sd",
//            "immi_postfire_perc_mean", "immi_postfire_perc_sd",
//            "immi_perc_mean", "immi_perc_sd",
//            "resi_germi", "resi_seedling", "resi_adult", "resi_source",
//            "immi_germi", "immi_seedling", "immi_adult", "immi_source",
//            "recol_germi", "recol_seedling", "recol_adult", "recol_source"
//        };

        const std::vector<std::string> header_metapop_pops {
            "sim_id", "study_num", "run_num", "year_max",
            "pop_id", "burned", "extinctions",
            "persistence", "generations",
            "pops_mean", "pops_sd",
            "immi_mean", "immi_sd",
            "resi_germi", "resi_seedling", "resi_adult", "resi_source",
            "immi_germi", "immi_seedling", "immi_adult", "immi_source",
            "recol_germi", "recol_seedling", "recol_adult", "recol_source"
        };

        print_header(header_metapop_pops, m_metapop_pops_stream);
    }

    if (params.is_out_recol_rescue)
    {
        assert(std::empty(m_rescue_effects_path));
        m_rescue_effects_path.assign(root_file_name);
        m_rescue_effects_path += "_rescue.csv";
        check_file_exists(m_rescue_effects_path);
        m_rescue_effects_stream.open(m_rescue_effects_path, std::ios_base::app);
        check_file_is_open(m_rescue_effects_stream, m_rescue_effects_path);

        const std::vector<std::string> header_rescue {
            "sim_id", "study_num", "run_num", "year_max",
            "generations", "extinctions",
            "resi_germi", "resi_seedling", "resi_adult", "resi_source",
            "immi_germi", "immi_seedling", "immi_adult", "immi_source",
            "recol_germi", "recol_seedling", "recol_adult", "recol_source"
        };

        print_header(header_rescue, m_rescue_effects_stream);
    }

    if (params.is_out_plants_per_pop)
    {
        assert(std::empty(m_plants_per_pop_path));
        m_plants_per_pop_path.assign(root_file_name);
        m_plants_per_pop_path += "_plants_per_pop.csv";
        check_file_exists(m_plants_per_pop_path);
        m_plants_per_pop_stream.open(m_plants_per_pop_path, std::ios_base::app);
        check_file_is_open(m_plants_per_pop_stream, m_plants_per_pop_path);

        const std::vector<std::string> header_plants_per_pop {
            "sim_id", "study_num",
            "pop_id", "pop_source_id",
            "germination_mean",  "germination_sd",
            "seedling_mean", "seedling_sd",
            "adult_mean", "adult_sd",
            "source_mean", "source_sd"
        };

        print_header(header_plants_per_pop, m_plants_per_pop_stream);
    }

    if (params.is_out_dispersal)
    {
        assert(std::empty(m_dispersal_path));
        m_dispersal_path.assign(root_file_name);
        m_dispersal_path += "_dispersal.csv";
        check_file_exists(m_dispersal_path);
        m_dispersal_stream.open(m_dispersal_path, std::ios_base::app);
        check_file_is_open(m_dispersal_stream, m_dispersal_path);

        const std::vector<std::string> header_dispersal {
            "sim_id", "study_num", "run_num",
            "bird_immi", "bird_resi",
            "wind_immi", "wind_resi", "wind_out", "wind_unsui"
        };

        print_header(header_dispersal, m_dispersal_stream);
    }
}
//-----------------------------------------------------------------------------

void Output::check_file_exists(const std::string& file_name)
{
    std::fstream infile(file_name);
    ///@note NOTE I add this variable to close the fstream to check memory leaks.
    auto is_good = infile.good();
    infile.close();

    if (is_good)
    {
        std::cerr << "error: The following file already exists:\n"
                << file_name << "\n\a";
        std::exit(EXIT_FAILURE);
    }
}

//-----------------------------------------------------------------------------

void Output::check_file_exists(const std::filesystem::path& path)
{
    if (std::filesystem::is_regular_file(path))
    {
        std::cerr << "error: The following file already exists:\n"
                  << path.relative_path() << "\n\a";
        std::exit(EXIT_FAILURE);
    }
}

//-----------------------------------------------------------------------------

void Output::check_file_is_open(const std::ofstream& ofs, const std::filesystem::path& path)
{
    if (!ofs.is_open())
    {
        std::cerr << "error: The following file already exists:\n"
                  << path.relative_path() << "\n\a";
        std::exit(EXIT_FAILURE);
    }
}

//-----------------------------------------------------------------------------

void Output::check_file_is_empty(const std::filesystem::path& path)
{
    if (!std::filesystem::is_empty(path))
    {
        std::cerr << "error: The following file is not empty:\n"
                  << path.relative_path() << "\n\a";
        std::exit(EXIT_FAILURE);
    }
}

//-----------------------------------------------------------------------------

void Output::cleanup()
{
    if (m_seed_dynamics_stream.is_open())
    {
        m_seed_dynamics_stream.close();
        m_seed_dynamics_stream.clear();
    }

    if (m_fire_output_stream.is_open())
    {
        m_fire_output_stream.close();
        m_fire_output_stream.clear();
    }

    if (m_metapop_stream.is_open())
    {
        m_metapop_stream.close();
        m_metapop_stream.clear();
    }

    if (m_metapop_full_stream.is_open())
    {
        m_metapop_full_stream.close();
        m_metapop_full_stream.clear();
    }

    if (m_metapop_dynamics_stream.is_open())
    {
        m_metapop_dynamics_stream.close();
        m_metapop_dynamics_stream.clear();
    }

    if (m_rescue_effects_stream.is_open())
    {
        m_rescue_effects_stream.close();
        m_rescue_effects_stream.clear();
    }

    if (m_plants_per_pop_stream.is_open())
    {
        m_plants_per_pop_stream.close();
        m_plants_per_pop_stream.clear();
    }
}

//-----------------------------------------------------------------------------

void Output::print_row(const std::ostringstream& ss, std::ofstream& stream)
{
    assert(stream.good());

    stream << ss.str();
    stream.flush();
}

//-----------------------------------------------------------------------------

void Output::print_header(const std::vector<std::string>& header, std::ofstream& stream)
{
    assert(stream.good());

    std::ostringstream ss;

    std::copy(header.begin(),
              header.end() - 1,
              std::ostream_iterator<std::string>(ss, ","));

    ss << header.back();
    stream << ss.str() << '\n';
    stream.flush();
}

//-----------------------------------------------------------------------------

void Output::print_study_area(const std::vector<std::vector<int>>& study_area,
                              const Parameters& params,
                              const unsigned study_num)
{
    // Expects
    assert(std::size(study_area) > 0);
    assert(std::size(study_area[0]) > 0);

    if (!std::filesystem::exists(m_output_folder + m_raw_folder))
    {
        std::filesystem::create_directories(m_output_folder);
    }

    const auto study_seed = params.study_replicates.at(study_num);

    const std::filesystem::path name_sim_file_copy {m_name_sim_file};

    std::filesystem::path root_file_name {m_output_folder + m_raw_folder};
    root_file_name /= name_sim_file_copy.stem();
    root_file_name += '_';
    root_file_name += m_time;
    root_file_name += "_sim_";
    root_file_name += std::to_string(params.sim_id);
    root_file_name += "_num_";
    root_file_name += std::to_string(study_num);
    root_file_name += "_seed_";
    root_file_name += std::to_string(study_seed);
    root_file_name += "_study.csv";

    check_file_exists(root_file_name);

    std::ofstream ofs;
    ofs.open(root_file_name, std::ios_base::app);
    check_file_is_open(ofs, root_file_name);

    const unsigned col_max {
        static_cast<unsigned>(std::size(study_area[0]) - 1)
    };

    for (unsigned row {0}; row < std::size(study_area); ++row)
    {
        std::ostringstream ss;
        for (unsigned col {0}; col < std::size(study_area[0]); ++col)
        {
            ss << study_area[row][col];

            if (col < col_max)
            {
                ss << ',';
            }
            else
            {
                ss << '\n';
            }
        }
        print_row(ss, ofs);
    }
    ofs.close();
    ofs.clear();
}

//-----------------------------------------------------------------------------

void Output::print_fire_output(
        const Parameters& params,
        const unsigned study_num,
        const int run_num,
        const int year,
        const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(!std::empty(metapopulation));
    assert(params.cell_size > 0);
    assert(run_num >= 0);
    assert(year > 0);

    int populations_burned {0};
    double burned_area_meters {0.0};
    double suitable_area_meters {0.0};

    for (const auto& [key, val] : metapopulation)
    {
        suitable_area_meters += val.get_num_cells() * std::pow(params.cell_size, 2);
        if (val.get_is_fire_event())
        {
            ++populations_burned;
            assert(val.get_num_cells() > 0);
            burned_area_meters += val.get_num_cells() * std::pow(params.cell_size, 2);
        }
    }

    // convert square meters to hectares
    constexpr double m2_to_ha { 10000.0 };
    const double suitable_area_ha { suitable_area_meters / m2_to_ha };
    const double burned_area_ha { burned_area_meters / m2_to_ha };

    assert(burned_area_ha >= 0.0);
    assert(suitable_area_ha > 0.0 && suitable_area_ha >= burned_area_ha);

    const double burned_area_perct { burned_area_ha / suitable_area_ha * 100.0 };

    assert(burned_area_perct >= 0);

    const double total_populations {
        static_cast<double>(std::size(metapopulation))
    };

    assert(populations_burned <= total_populations);

    const double populations_burned_perct {
        populations_burned / total_populations * 100.0
    };

    assert(populations_burned_perct >= 0.0 && populations_burned_perct <= 100.0);

    std::ostringstream ss;

    ss << params.sim_id            << ','
       << study_num                << ','
       << run_num                  << ','
       << year                     << ','
       << populations_burned       << ','
       << populations_burned_perct << ','
       << burned_area_ha           << ','
       << burned_area_perct        << '\n';

    print_row(ss, m_fire_output_stream);
}

//-----------------------------------------------------------------------------

void Output::print_metapopulation(const Parameters& params,
                                  const unsigned study_num,
                                  const int run_num,
                                  const int year_max,
                                  const std::map<int, Population>& metapopulation,
                                  std::vector<int>& occupied_patches)
{
    // Expects
    assert(params.is_out_metapopulation);
    assert(params.run_time_max > 0);
    assert(!std::empty(metapopulation));


//    const unsigned run_time_max { static_cast<unsigned>(params.run_time_max) };

//    assert(std::size(occupied_patches) <= run_time_max);

    /**
     * @note NOTE I commmented the if statement below to calculate
     * the occupied patches until the metapopulation extinction
     */
    // add zeros in case metapopulation went extinct before run time max
//    if (std::size(occupied_patches) < run_time_max)
//    {
//        occupied_patches.resize(run_time_max);
//    }

    std::vector<int> meta_persistence;
    std::vector<int> meta_generations;

    for ([[maybe_unused]] const auto& [key, val] : metapopulation)
    {
        meta_persistence.push_back(val.get_persistence_time());
        meta_generations.push_back(val.get_total_generations());
    }

    const double num_patches {
        static_cast<double>(std::size(metapopulation))
    };

    double occupied_max {0.0};

    if (!std::empty(occupied_patches))
    {
        occupied_max = static_cast<double>(
                    *std::max_element(
                        std::cbegin(occupied_patches),
                        std::cend(occupied_patches))) / num_patches * 100.0;
    }

    const double occupied_mean {
        Equation::mean(occupied_patches) / num_patches * 100.0
    };

    const double occupied_sd {
        Equation::sd(occupied_patches) / num_patches * 100.0
    };

    assert(occupied_max >= 0.0 && occupied_max <= 100.0);
    assert(occupied_mean >= 0.0 && occupied_mean <= 100.0);
    assert(occupied_sd >= 0.0 && occupied_sd <= 100.0);

    int meta_persistence_max {0};

    assert(!std::empty(meta_persistence));

    if (!std::empty(meta_persistence))
    {
        meta_persistence_max = *std::max_element(
                    std::cbegin(meta_persistence),
                    std::cend(meta_persistence));
    }

    assert(meta_persistence_max >= 0);

    assert(!std::empty(meta_generations));

    int meta_generations_max {0};

    if (!std::empty(meta_generations))
    {
        meta_generations_max = *std::max_element(
                    std::cbegin(meta_generations), std::cend(meta_generations));
    }

    assert(meta_generations_max >= 0);

    std::ostringstream ss;
    ss << params.sim_id                    << ','
       << study_num                        << ','
       << run_num                          << ','
       << year_max                         << ','
       << meta_persistence_max             << ','
       << Equation::mean(meta_persistence) << ','
       << Equation::sd(meta_persistence)   << ','
       << occupied_max                     << ','
       << occupied_mean                    << ','
       << occupied_sd                      << ','
       << meta_generations_max             << ','
       << Equation::mean(meta_generations) << ','
       << Equation::sd(meta_generations)   << '\n';

    print_row(ss, m_metapop_stream);
}

//-----------------------------------------------------------------------------

void Output::print_metapopulation_full(
        const Parameters& params,
        const unsigned study_num,
        const int run_num,
        const int year_max,
        const std::map<int, Population>& metapopulation,
        const std::vector<int>& occupied_patches,
        const std::vector<double>& immigrants_percent_postfire,
        const std::vector<double>& residents_percent_postfire,
        const std::vector<int>& immigrants_postfire,
        const std::vector<int>& residents_postfire)
{
    // Expects
    assert(params.is_out_metapopulation_full);
    assert(params.run_time_max > 0);
    assert(!std::empty(metapopulation));

//    const unsigned run_time_max { static_cast<unsigned>(params.run_time_max) };

//    assert(std::size(occupied_patches) <= run_time_max);
//    assert(std::size(immigrants_percent_postfire) <= run_time_max);
//    assert(std::size(residents_percent_postfire) <= run_time_max);
//    assert(std::size(immigrants_postfire) <= run_time_max);
//    assert(std::size(residents_postfire) <= run_time_max);
//    assert(std::size(immigrants_percent) <= run_time_max);
//    assert(std::size(residents_percent) <= run_time_max);
//    assert(std::size(immigrants) <= run_time_max);
//    assert(std::size(residents) <= run_time_max);

    /**
     * @note NOTE I commmented the if statements below to calculate
     * the only those variables until the metapopulation extinction
     */
    // add zeros in case metapopulation went extinct before run time max
//    if (std::size(occupied_patches) < run_time_max)
//    {
//        occupied_patches.resize(run_time_max);
//    }

//    if (std::size(immigrants_percent_postfire) < run_time_max)
//    {
//        immigrants_percent_postfire.resize(run_time_max);
//        residents_percent_postfire.resize(run_time_max);
//        immigrants_postfire.resize(run_time_max);
//        residents_postfire.resize(run_time_max);

//        immigrants_percent.resize(run_time_max);
//        residents_percent.resize(run_time_max);
//        immigrants.resize(run_time_max);
//        residents.resize(run_time_max);
//    }

    std::vector<int> meta_persistence;
    std::vector<int> meta_generations;

    std::vector<double> meta_genetic_pops_postfire_mean;
    double meta_genetic_pops_postfire_max {0.0};

    std::vector<double> meta_genetic_pops_mean;
    double meta_genetic_pops_max {0.0};

    for (const auto& [key, val] : metapopulation)
    {
        meta_persistence.push_back(val.get_persistence_time());
        meta_generations.push_back(val.get_total_generations());

        // postfire
        {
//            std::vector<int> genetic_pops_postfire {
//                val.get_genetic_pops_postfire()
//            };
            std::vector<int> genetic_pops_postfire {
                val.get_genetic_pops_postfire_occupied()
            };


//            assert(std::size(genetic_pops_postfire) <= run_time_max);

            int genetic_pops_postfire_max {0};

            if (!std::empty(genetic_pops_postfire))
            {
                genetic_pops_postfire_max = *std::max_element(
                            std::cbegin(genetic_pops_postfire)
                            , std::cend(genetic_pops_postfire));
            }

            if (meta_genetic_pops_postfire_max < genetic_pops_postfire_max)
            {
                meta_genetic_pops_postfire_max = genetic_pops_postfire_max;
            }

            /**
             * @note NOTE I commmented the if statement below to calculate
             * the populations per dune until the metapopulation extinction
             */
            // add zeros if the vector genetic_pops is smaller than the max time
            // in order to calculate the mean in the entire simulation run
//            if (std::size(genetic_pops_postfire) < run_time_max)
//            {
//                genetic_pops_postfire.resize(run_time_max);
//            }

            meta_genetic_pops_postfire_mean.push_back(
                        Equation::mean(genetic_pops_postfire));
        }

        // genetic
        {
//            std::vector<int> genetic_pops { val.get_genetic_pops() };
            std::vector<int> genetic_pops { val.get_genetic_pops_occupied() };

//            assert(std::size(genetic_pops) <= run_time_max);

            int genetic_pops_max {0};

            if (!std::empty(genetic_pops))
            {
                genetic_pops_max = *std::max_element(
                            std::cbegin(genetic_pops)
                            , std::cend(genetic_pops));
            }

            if (meta_genetic_pops_max < genetic_pops_max)
            {
                meta_genetic_pops_max = genetic_pops_max;
            }

            // add zeros if the vector genetic_pops is smaller than the max time
            // in order to calculate the mean in the entire simulation run
//            if (std::size(genetic_pops) < run_time_max)
//            {
//                genetic_pops.resize(run_time_max);
//            }

            meta_genetic_pops_mean.push_back(Equation::mean(genetic_pops));
        }
    }

    const double num_patches {
        static_cast<double>(std::size(metapopulation))
    };

    double occupied_max {0.0};

    if (!std::empty(occupied_patches))
    {
        occupied_max = static_cast<double>(
                    *std::max_element(std::cbegin(occupied_patches),
                                      std::cend(occupied_patches))
                    ) / num_patches * 100.0;
    }

    const double occupied_mean {
        Equation::mean(occupied_patches) / num_patches * 100.0
    };

    const double occupied_sd {
        Equation::sd(occupied_patches) / num_patches * 100.0
    };

    assert(occupied_max >= 0.0 && occupied_max <= 100.0);
    assert(occupied_mean >= 0.0 && occupied_mean <= 100.0);
    assert(occupied_sd >= 0.0 && occupied_sd <= 100.0);

    int meta_persistence_max {0};

    assert(!std::empty(meta_persistence));

    if (!std::empty(meta_persistence))
    {
        meta_persistence_max = *std::max_element(
                    std::cbegin(meta_persistence),
                    std::cend(meta_persistence));
    }

    assert(meta_persistence_max >= 0);

    assert(!std::empty(meta_generations));

    int meta_generations_max {0};

    if (!std::empty(meta_generations))
    {
        meta_generations_max = *std::max_element(
                    std::cbegin(meta_generations), std::cend(meta_generations));
    }

    assert(meta_generations_max >= 0);

    std::ostringstream ss;
    ss << params.sim_id                                                 << ','
       << study_num                                                     << ','
       << run_num                                                       << ','
       << year_max                                                      << ','
       << meta_persistence_max                                          << ','
       << Equation::mean(meta_persistence)                              << ','
       << Equation::sd(meta_persistence)                                << ','
       << occupied_max                                                  << ','
       << occupied_mean                                                 << ','
       << occupied_sd                                                   << ','
       << meta_generations_max                                          << ','
       << Equation::mean(meta_generations)                              << ','
       << Equation::sd(meta_generations)                                << ','

//       << Equation::mean(meta_genetic_pops_mean)                        << ','
//       << Equation::sd(meta_genetic_pops_mean)                          << ','

       << Equation::mean(meta_genetic_pops_postfire_mean)               << ','
       << Equation::sd(meta_genetic_pops_postfire_mean)                 << ','

//       << Equation::mean(immigrants)                                    << ','
//       << Equation::sd(immigrants)                                      << ','
//       << Equation::mean(residents)                                     << ','
//       << Equation::sd(residents)                                       << ','

//       << Equation::mean(immigrants_percent)                            << ','
//       << Equation::sd(immigrants_percent)                              << ','
//       << Equation::mean(residents_percent)                             << ','
//       << Equation::sd(residents_percent)                               << ','

       << Equation::mean(immigrants_postfire)                           << ','
       << Equation::sd(immigrants_postfire)                             << ','
       << Equation::mean(residents_postfire)                            << ','
       << Equation::sd(residents_postfire)                              << ','

       << Equation::mean(immigrants_percent_postfire)                   << ','
       << Equation::sd(immigrants_percent_postfire)                     << ','
       << Equation::mean(residents_percent_postfire)                    << ','
       << Equation::sd(residents_percent_postfire)                      << '\n';

    print_row(ss, m_metapop_full_stream);
}

//-----------------------------------------------------------------------------

void Output::print_metapopulation_dynamics(
        const Parameters& params
        , const unsigned int study_num
        , const int run_num
        , const int year
        , const int time_since_fire
        , const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(params.is_out_metapopulation_dynamics);
    assert(params.run_time_max > 0);
    assert(!std::empty(metapopulation));


    long double immigrants_postfire {0.0L};
//    long double immigrants        {0.0L};

    long double immigrants_postfire_birds {0.0L};
    long double immigrants_postfire_wind  {0.0L};

//    long double immigrants_birds {0.0L};
//    long double immigrants_wind  {0.0L};

    long double total_plants {0.0L};
    int count_extinct_pop    {0};

    const double total_populations {
        static_cast<double>(std::size(metapopulation))
    };

    assert(total_populations >= 1.0);  // there must be at least one population

    std::vector<int> pops_per_patch;
//    std::vector<int> gene_diversity;

    int total_populations_burned {0};
    ///@todo TODO this commented variables will be placed in dispersal output

//    long double LDD_seeds_wind_resident   {0.0L};
//    long double LDD_seeds_wind_immigrant  {0.0L};
//    long double LDD_seeds_wind_unsuitable {0.0L};
//    long double LDD_seeds_wind_outside    {0.0L};
//    long double LDD_seeds_bird_resident   {0.0L};
//    long double LDD_seeds_bird_immigrant  {0.0L};

    for (const auto& [key, val] : metapopulation)
    {
        if (val.get_time_since_fire() == 0)
        {
            ++total_populations_burned;
//            LDD_seeds_wind_resident += val.get_total_LDD_seeds_by_wind_resident();
//            LDD_seeds_wind_immigrant += val.get_total_LDD_seeds_by_wind_immigrant();
//            LDD_seeds_wind_unsuitable += val.get_total_LDD_seeds_by_wind_unsuitable();
//            LDD_seeds_wind_outside += val.get_total_LDD_seeds_by_wind_outside();
//            LDD_seeds_bird_resident += val.get_total_LDD_seeds_by_birds_resident();
//            LDD_seeds_bird_immigrant += val.get_total_LDD_seeds_by_birds_immigrant();
        }

        count_extinct_pop += val.get_is_extinct();

//        std::unordered_set<int> pop_ids;
        std::unordered_set<int> pop_ids_generation;

        for (const auto& coh : val.get_cohort_list())
        {
//            pop_ids.emplace(coh.get_gene_flow().get_genetic());
            pop_ids_generation.emplace(coh.get_gene_flow().get_previous());

            if (coh.get_gene_flow().get_previous() != key)
            {
                immigrants_postfire += coh.get_num_plants();

                switch (coh.get_dispersal_vector())
                {
                    case Dispersal_vector::LDD_bird :
                        immigrants_postfire_birds += coh.get_num_plants();
                        break;

                    case Dispersal_vector::LDD_wind :
                        immigrants_postfire_wind += coh.get_num_plants();
                        break;

                    case Dispersal_vector::SDD_wind :
                        assert(false);    // SDD cohort should not be considered immigrant
                        break;            // e.g. An error was to assign interfire cohorts as SDD_wind
                }
            }

//            if (coh.get_gene_flow().get_genetic() != key)
//            {
//                immigrants += coh.get_num_plants();
//                switch (coh.get_dispersal_vector())
//                {
//                    case Dispersal_vector::LDD_bird :
//                        immigrants_birds += coh.get_num_plants();
//                        break;

//                    case Dispersal_vector::LDD_wind :
//                        immigrants_wind += coh.get_num_plants();
//                        break;
//                    case Dispersal_vector::SDD_wind :
//                        // nothing happens here
//                        break;
//                }
//            }
            total_plants += coh.get_num_plants();
        }

//        gene_diversity.push_back(static_cast<int>(std::size(pop_ids)));
        pops_per_patch.push_back(
                    static_cast<int>(std::size(pop_ids_generation)));
    }

    std::ostringstream ss;
    ss << params.sim_id                                           << ','
       << study_num                                               << ','
       << run_num                                                 << ','
       << year                                                    << ','
       << time_since_fire                                         << ','
       << total_populations_burned / total_populations * 100.0    << ','
       << (total_populations - count_extinct_pop) /
          total_populations * 100.0                               << ','
       << total_plants                                            << ','
       << immigrants_postfire                                       << ','
       << immigrants_postfire_birds                                 << ','
       << immigrants_postfire_wind                                  << ',';

    if (total_plants >= 1.0L)
    {
        ss << immigrants_postfire / total_plants * 100 << ',';
    }
    else
    {
        ss << 0 << ',';
    }
    ss << Equation::mean(pops_per_patch) << ','
       << Equation::sd(pops_per_patch)   << '\n';
//       << Equation::mean(gene_diversity)            << ','
//       << Equation::sd(gene_diversity)              << '\n';

    print_row(ss, m_metapop_dynamics_stream);
}

//-----------------------------------------------------------------------------

void Output::print_metapopulation_pops(
        const Parameters& params,
        const unsigned study_num,
        const int run_num,
        const int year_max,
        const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(params.is_out_metapopulation_pops);
    assert(params.run_time_max > 0);
    assert(!std::empty(metapopulation));

    for (const auto& [key, val] : metapopulation)
    {
          std::ostringstream ss;

          ss << params.sim_id                                            << ','
             << study_num                                                << ','
             << run_num                                                  << ','
             << year_max                                                 << ','
             << key                                                      << ','
             << val.get_total_fires()                                    << ','
             << val.get_extinctions()                                    << ','
             << val.get_persistence_time()                               << ','
             << val.get_total_generations()                              << ',';

          if (std::empty(val.get_genetic_pops_postfire_occupied()))
          {
              ss << '0' << ',';
              ss << '0' << ',';
          }
          else
          {
              ss << Equation::mean(val.get_genetic_pops_postfire_occupied()) << ','
                 << Equation::sd(val.get_genetic_pops_postfire_occupied())   << ',';

          }

          if (std::empty(val.get_immigrants_postfire_perc()))
          {
              ss << '0' << ',';
              ss << '0' << ',';
          }
          else
          {
              ss << Equation::mean(val.get_immigrants_postfire_perc()) << ','
                 << Equation::sd(val.get_immigrants_postfire_perc())   << ',';
          }

          ss << val.get_total_LDD_residents_germination()                << ','
             << val.get_total_LDD_residents_seedling()                   << ','
             << val.get_total_LDD_residents_adult()                      << ','
             << val.get_total_LDD_residents_source()                     << ','
             << val.get_total_immigrants_germination()                   << ','
             << val.get_total_immigrants_seedling()                      << ','
             << val.get_total_immigrants_adult()                         << ','
             << val.get_total_immigrants_source()                        << ','
             << val.get_total_recolonizations_germination()              << ','
             << val.get_total_recolonizations_seedling()                 << ','
             << val.get_total_recolonizations_adult()                    << ','
             << val.get_total_recolonizations_source()                   << '\n';

          print_row(ss, m_metapop_pops_stream);
    }
}

//-----------------------------------------------------------------------------

void Output::print_rescue_effects(
        const Parameters& params,
        const unsigned int study_num,
        const int run_num,
        const int year_max,
        const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(params.is_out_recol_rescue);
    assert(params.run_time_max > 0);
    assert(!std::empty(metapopulation));

    int resi_germi    {0};
    int resi_seedling {0};
    int resi_adult    {0};
    int resi_source   {0};

    int immi_germi    {0};
    int immi_seedling {0};
    int immi_adult    {0};
    int immi_source   {0};

    int recol_germi    {0};
    int recol_seedling {0};
    int recol_adult    {0};
    int recol_source   {0};

    int generations {0};
    int extinctions {0};

    for ([[maybe_unused]] const auto& [key, val] : metapopulation)
    {
        generations   += val.get_total_generations();
        extinctions   += val.get_extinctions();

        resi_germi    += val.get_total_LDD_residents_germination();
        resi_seedling += val.get_total_LDD_residents_seedling();
        resi_adult    += val.get_total_LDD_residents_adult();
        resi_source   += val.get_total_LDD_residents_source();

        immi_germi    += val.get_total_immigrants_germination();
        immi_seedling += val.get_total_immigrants_seedling();
        immi_adult    += val.get_total_immigrants_adult();
        immi_source   += val.get_total_immigrants_source();

        recol_germi    += val.get_total_recolonizations_germination();
        recol_seedling += val.get_total_recolonizations_seedling();
        recol_adult    += val.get_total_recolonizations_adult();
        recol_source   += val.get_total_recolonizations_source();
    }

    std::ostringstream ss;

    ss << params.sim_id     << ','
       << study_num         << ','
       << run_num           << ','
       << year_max          << ','
       << generations       << ','
       << extinctions       << ','
       << resi_germi        << ','
       << resi_seedling     << ','
       << resi_adult        << ','
       << resi_source       << ','
       << immi_germi        << ','
       << immi_seedling     << ','
       << immi_adult        << ','
       << immi_source       << ','
       << recol_germi       << ','
       << recol_seedling    << ','
       << recol_adult       << ','
       << recol_source      << '\n';

    print_row(ss, m_rescue_effects_stream);
}

//-----------------------------------------------------------------------------

void Output::print_plants_per_pop(const Parameters& params, const unsigned study_num)
{
    // Expects
    assert(params.is_out_plants_per_pop);

    for (const auto& [key, val] : m_plants_per_pop_out)
    {
        for (const auto& [key2, val2] : val)
        {
            std::vector<long long int> germination;
            std::vector<long long int> seedling;
            std::vector<long long int> adult;
            std::vector<long long int> source;

            for (const auto& v : val2)
            {
                germination.push_back(v.germination);
                seedling.push_back(v.seedling);
                adult.push_back(v.adult);
                source.push_back(v.source);
            }

            std::ostringstream ss;
            ss << params.sim_id               << ','
               << study_num                   << ','
               << key                         << ','
               << key2                        << ','
               << Equation::mean(germination) << ','
               << Equation::sd(germination)   << ','
               << Equation::mean(seedling)    << ','
               << Equation::sd(seedling)      << ','
               << Equation::mean(adult)       << ','
               << Equation::sd(adult)         << ','
               << Equation::mean(source)      << ','
               << Equation::sd(source)        << '\n';

            print_row(ss, m_plants_per_pop_stream);
        }
    }
}

//-----------------------------------------------------------------------------

void Output::add_plants_per_population(const std::map<int, Population>& metapopulation)
{
    assert(!std::empty(metapopulation));

    for (const auto& [key, val] : metapopulation)
    {
        assert(m_plants_per_pop_out.find(key) != m_plants_per_pop_out.end());  // not found

        for (const auto& [key2, val2] : val.get_plants_per_pop())
        {
            m_plants_per_pop_out.at(key).at(key2).push_back(val2);
        }
    }
}

//-----------------------------------------------------------------------------

void Output::init_plants_per_population(const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(!std::empty(metapopulation));

    m_plants_per_pop_out.clear();

    for (const auto& [key, val] : metapopulation)
    {
        std::map<int, std::vector<Total_plants_stage>> my_map;

        for (const auto& [key2, val2] : metapopulation)
        {
            my_map.emplace(key2, std::vector<Total_plants_stage>());
        }

        m_plants_per_pop_out.emplace(key, my_map);
    }

    assert(std::size(m_plants_per_pop_out) == std::size(metapopulation));
}

//-----------------------------------------------------------------------------

void Output::print_seed_dynamics(
        const Parameters& params,
        const unsigned study_num,
        const int run_num,
        const int year,
        const int time_since_fire,
        const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(params.is_out_seed_dynamics);
    assert(params.run_time_max > 0);
    assert(run_num >= 0);
    assert(year >= 0);
    assert(time_since_fire >= 0);
    assert(!std::empty(metapopulation));

    for (const auto& [key, val]: metapopulation)
    {
        for (const auto& cohort : val.get_cohort_list())
        {
            const auto num_plants = cohort.get_num_plants();
            const auto seedbank_cones = cohort.get_seedbank_cones(params);
            const auto seedbank = cohort.get_seedbank(params);
            assert(num_plants > 0);
            assert(seedbank_cones >= 0);
            assert(seedbank >= 0);

            std::ostringstream ss;

            ss << params.sim_id              << ','
               << study_num                  << ','
               << run_num                    << ','
               << year                       << ','
               << key                        << ','
               << time_since_fire            << ','
               << cohort.get_age()           << ','
               << cohort.get_is_LDD_cohort() << ','
               << cohort.get_is_postfire()   << ','
               << num_plants                 << ',';

            switch (params.model_type)
            {
            case Model_type::cohort_based:
                ss << seedbank_cones / num_plants << ','
                   << seedbank / num_plants       << ',';
                break;
            case Model_type::plant_based:
                // each plant has different number of cones and seeds
                ss << "NA" << ','
                   << "NA" << ',';
                break;
            }

            ss << seedbank_cones << ','
               << seedbank       << '\n';

            print_row(ss, m_seed_dynamics_stream);
        }
    }
}

// ----------------------------------------------------------------------------

void Output::print_dispersal(
        const Parameters& params
        , const unsigned study_num
        , const int run_num
        , const std::map<int, Population>& metapopulation)
{
    // Expects
    assert(params.is_out_dispersal);
    assert(run_num >= 0);
    assert(!std::empty(metapopulation));

    long long int birds_immigrant {0};
    long long int birds_resident {0};
    long long int wind_immigrant {0};
    long long int wind_resident {0};
    long long int wind_outside {0};
    long long int wind_unsuitable {0};

    for (const auto& [key, val]: metapopulation)
    {
          birds_immigrant += val.get_total_LDD_seeds_by_birds_immigrant();
          birds_resident += val.get_total_LDD_seeds_by_birds_resident();
          wind_immigrant += val.get_total_LDD_seeds_by_wind_immigrant();
          wind_resident += val.get_total_LDD_seeds_by_wind_resident();
          wind_outside += val.get_total_LDD_seeds_by_wind_outside();
          wind_unsuitable += val.get_total_LDD_seeds_by_wind_unsuitable();
    }

    std::ostringstream ss;

    ss << params.sim_id   << ','
       << study_num       << ','
       << run_num         << ','
       << birds_immigrant << ','
       << birds_resident  << ','
       << wind_immigrant  << ','
       << wind_resident   << ','
       << wind_outside    << ','
       << wind_unsuitable << '\n';

    print_row(ss, m_dispersal_stream);
}

//void Output::print_seed_dynamics(const std::shared_ptr<Population>& pop,
//                                 const int sim_replicate_num,
//                                 const int run_num,
//                                 const int year,
//                                 const int run_fire_interval_meam,
//                                 const double run_pollination_success_mean)
//{
////    assert(m_sim_output == Sim_output::seed_dynamics);

////    std::ostringstream ss;

////    if (pop->get_cohort_list().empty() == false)
////    {
////        for (const auto& coh : pop->get_cohort_list())
////        {
////            ss << pop->get_id()                                       << ','
////               << Params::Sim::id()                                   << ','
////               << sim_replicate_num                                   << ','
////               << Params::Scenario::id()                              << ','
////               << run_num                                             << ','
////               << year                                                << ','
////               << Params::Scenario::is_climate()                      << ','
////               << Params::Scenario::file_name_climate()               << ','
////               << Params::Scenario::is_fire()                         << ','
////               << static_cast<int>(Params::Scenario::fire_scenario()) << ','
////               << run_fire_interval_meam                              << ','
////               << run_pollination_success_mean                        << ','
////               << pop->get_time_since_fire()                          << ','
////               << pop->get_fire_interval()                            << ','
////               << coh.get_is_postfire()                              << ','
////               << coh.get_age()                                      << ','
////               << coh.get_num_plants()                               << ',';

////            switch (Params::Scenario::model_type())
////            {
/// WARNING: get_seedbank() and get_seedbank_cones() in cohort-based take all within cohort (before was per plant)
////                case Model_type::cohort_based :
////                    ss << coh.get_seedbank()                               << ','
////                       << coh.get_seedbank() * coh.get_num_plants()       << ','
////                       << coh.get_seedbank_cones()                         << ','
////                       << coh.get_seedbank_cones() * coh.get_num_plants() << '\n';
////                    break;

////                case Model_type::individual_based :
////                    ss << "NA"                      << ','
////                       << coh.get_seedbank()       << ','
////                       << "NA"                      << ','
////                       << coh.get_seedbank_cones() << '\n';
////                    break;
////            }
////        }
////    }
////    else
////    {
////        ss << pop->get_id()                                       << ','
////           << Params::Sim::id()                                   << ','
////           << sim_replicate_num                                   << ','
////           << Params::Scenario::id()                              << ','
////           << run_num                                             << ','
////           << year                                                << ','
////           << Params::Scenario::is_climate()                      << ','
////           << Params::Scenario::file_name_climate()               << ','
////           << Params::Scenario::is_fire()                         << ','
////           << static_cast<int>(Params::Scenario::fire_scenario()) << ','
////           << run_fire_interval_meam                              << ','
////           << run_pollination_success_mean                        << ','
////           << pop->get_time_since_fire()                          << ','
////           << pop->get_fire_interval()                            << ','
////           << 0                                                   << ','
////           << 0                                                   << ','
////           << 0                                                   << ','
////           << 0                                                   << ','
////           << 0                                                   << ','
////           << 0                                                   << ','
////           << 0                                                   << '\n';
////    }

////    print_row(ss, m_file_stream);
//}

////-----------------------------------------------------------------------------

//void Output::print_persistence(const std::shared_ptr<Population>& pop,
//                               const int sim_replicate_num,
//                               const int run_num,
//                               const int run_fire_interval_mean,
//                               const double run_pollination_success_mean)
//{
////    assert(m_sim_output == Sim_output::persistence);

////    std::ostringstream ss;

////    ss  << m_id                                                << ','
////        << Params::Sim::id()                                   << ','
////        << sim_replicate_num                                   << ','
////        << Params::Scenario::id()                              << ','
////        << run_num                                             << ','
////        << Params::Scenario::is_climate()                      << ','
////        << Params::Scenario::file_name_climate()               << ','
////        << Params::Scenario::is_fire()                         << ','
////        << static_cast<int>(Params::Scenario::fire_scenario()) << ','
////        << run_pollination_success_mean                        << ','
////        << run_fire_interval_mean                              << ','
////        << pop->get_is_extinct()                               << ','
////        << pop->get_persistence_time()                         << ','
////        << pop->get_num_generations()                          << std::endl;

////    print_row(ss, m_file_stream);
////}

//////-----------------------------------------------------------------------------

//void Output::print_metapopulation(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num)
//{
//    assert(m_sim_output == Sim_output::metapopulation);

//    const unsigned run_time_max {
//        static_cast<unsigned>(Params::Sim::run_time_max())
//    };

//    std::vector<int> meta_persistence;
//    std::vector<int> meta_generations;
//    std::vector<double> meta_genetic_pops;
//    int total_extinct {0};

//    for (const auto& population : metapopulation)
//    {

//        std::vector<int> genetic_pops {population->get_genetic_pops()};

//        assert(genetic_pops.size() <= run_time_max);

//        if (genetic_pops.size() < run_time_max)
//        {
//            genetic_pops.resize(run_time_max);
//        }

//        meta_persistence.push_back(population->get_persistence_time());
//        meta_generations.push_back(population->get_num_generations());
//        meta_genetic_pops.push_back(Equation::mean(genetic_pops));

//        total_extinct += population->get_is_extinct();
//    }

//    std::ostringstream ss;

//    ss << Params::Sim::id()                                 << ','
//       << sim_replicate_num                                 << ','
//       << Params::Scenario::id()                            << ','
//       << run_num                                           << ','
//       << Equation::mean(meta_persistence)                  << ','
//       << *std::max_element(std::begin(meta_persistence),
//                            std::end(meta_persistence))     << ','
//       << Equation::mean(meta_generations)                  << ','
//       << *std::max_element(std::begin(meta_generations),
//                            std::end(meta_generations))     << ','
//       << Equation::mean(meta_genetic_pops)                 << ','
//       << *std::max_element(std::begin(meta_genetic_pops),
//                            std::end(meta_genetic_pops))    << ','
//       << total_extinct                                     << '\n';

//    print_row(ss, m_file_stream);
//}

////-----------------------------------------------------------------------------

//void Output::print_metapopulation_pops(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num)
//{
////    assert(m_sim_output == Sim_output::metapopulation_pops);

////    const unsigned run_time_max {
////        static_cast<unsigned>(Params::Sim::run_time_max())
////    };

////    for (const auto& pop : metapopulation)
////    {
//////        std::vector<int> residents {pop->get_residents()};
//////        std::vector<int> immigrants {pop->get_immigrants()};
////        std::vector<int> genetic_pops {pop->get_genetic_pops()};

//////        assert(residents.size() <= run_time_max);
//////        assert(immigrants.size() <= run_time_max);
////        assert(genetic_pops.size() <= run_time_max);

//////        assert(residents.size() == immigrants.size());
//////        assert(residents.size() == genetic_pops.size());

////        if (genetic_pops.size() < run_time_max)
////        {
//////            residents.resize(run_time_max);
//////            immigrants.resize(run_time_max);
////            genetic_pops.resize(run_time_max);
////        }

////        std::ostringstream ss;

////        ss << Params::Sim::id()            << ','
////           << sim_replicate_num            << ','
////           << Params::Scenario::id()       << ','
////           << run_num                      << ','
////           << pop->get_id()                << ','
////           << Equation::mean(genetic_pops) << ','
////           << pop->get_persistence_time()  << ','
////           << pop->get_num_generations()   << ','
////           << pop->get_is_extinct()        << '\n';

////        print_row(ss, m_file_stream);
////    }
//}

////-----------------------------------------------------------------------------

////-----------------------------------------------------------------------------


//void Output::print_metapopulation_dynamics(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num,
//        const int run_year,
//        const int time_since_fire)
//{
////    assert(m_sim_output == Sim_output::metapopulation_dynamics);

////    long double immigrants_he2004 {0.0L};
////    long double immigrants        {0.0L};

////    long double immigrants_he2004_birds {0.0L};
////    long double immigrants_he2004_wind  {0.0L};

////    long double immigrants_birds {0.0L};
////    long double immigrants_wind  {0.0L};

////    long double total_plants            {0.0L};
////    int count_extinct_pop               {0};
////    const double total_populations {static_cast<double>(metapopulation.size())};

////    assert(total_populations >= 1.0);  // there must be at least one population

////    std::vector<int> gene_diversity_generation;
////    std::vector<int> gene_diversity;

////    int total_populations_burned {0};

////    long double LDD_seeds_wind_resident   {0.0L};
////    long double LDD_seeds_wind_immigrant  {0.0L};
////    long double LDD_seeds_wind_unsuitable {0.0L};
////    long double LDD_seeds_wind_outside    {0.0L};
////    long double LDD_seeds_bird_resident   {0.0L};
////    long double LDD_seeds_bird_immigrant  {0.0L};

////    for (const auto& pop : metapopulation)
////    {
////        if (pop->get_time_since_fire() == 0)
////        {
////            ++total_populations_burned;

////            LDD_seeds_wind_resident += pop->get_total_LDD_seeds_by_wind_resident();
////            LDD_seeds_wind_immigrant += pop->get_total_LDD_seeds_by_wind_immigrant();
////            LDD_seeds_wind_unsuitable += pop->get_total_LDD_seeds_by_wind_unsuitable();
////            LDD_seeds_wind_outside += pop->get_total_LDD_seeds_by_wind_outside();
////            LDD_seeds_bird_resident += pop->get_total_LDD_seeds_by_birds_resident();
////            LDD_seeds_bird_immigrant += pop->get_total_LDD_seeds_by_birds_immigrant();
////        }

////        count_extinct_pop += pop->get_is_extinct();

////        std::unordered_set<int> pop_ids;
////        std::unordered_set<int> pop_ids_generation;

////        for (const auto& coh : pop->get_cohort_list())
////        {
////            pop_ids.emplace(coh.get_gene_flow().get_genetic());
////            pop_ids_generation.emplace(
////                        coh.get_gene_flow().get_geographical_previous());

////            if (coh.get_gene_flow().get_geographical_previous() != pop->get_id())
////            {
////                immigrants_he2004 += coh.get_num_plants();
////                switch (coh.get_dispersal_vector())
////                {
////                    case Dispersal_vector::LDD_bird :
////                        immigrants_he2004_birds += coh.get_num_plants();
////                        break;

////                    case Dispersal_vector::LDD_wind :
////                        immigrants_he2004_wind += coh.get_num_plants();
////                        break;
////                    case Dispersal_vector::SDD_wind :
////                        std::cerr << "error: SDD_wind should not happen in immigrants he2004!\n\a";
////                        std::exit(EXIT_FAILURE);
////                }
////            }

////            if (coh.get_gene_flow().get_genetic() != pop->get_id())
////            {
////                immigrants += coh.get_num_plants();
////                switch (coh.get_dispersal_vector())
////                {
////                    case Dispersal_vector::LDD_bird :
////                        immigrants_birds += coh.get_num_plants();
////                        break;

////                    case Dispersal_vector::LDD_wind :
////                        immigrants_wind += coh.get_num_plants();
////                        break;
////                    case Dispersal_vector::SDD_wind :
////                        // nothing happens here
////                        break;
////                }
////            }
////            total_plants += coh.get_num_plants();
////        }

////        gene_diversity.push_back(
////                    static_cast<int>(pop_ids.size()));
////        gene_diversity_generation.push_back(
////                    static_cast<int>(pop_ids_generation.size()));
////    }

////    std::ostringstream ss;

////    ss << Params::Sim::id()                                       << ','
////       << sim_replicate_num                                       << ','
////       << Params::Scenario::id()                                  << ','
////       << run_num                                                 << ','
////       << run_year                                                << ','
////       << time_since_fire                                         << ','
////       << total_populations_burned                                << ','
////       << static_cast<int>(Params::Scenario::mort_scenario())     << ','
////       << Params::Scenario::LDD_wind_prop()                       << ','
////       << Params::Scenario::LDD_birds_prop()                      << ','
////       << Params::Scenario::LDD_birds_a()                         << ','
////       << LDD_seeds_wind_resident                                 << ','
////       << LDD_seeds_wind_immigrant                                << ','
////       << LDD_seeds_wind_unsuitable                               << ','
////       << LDD_seeds_wind_outside                                  << ','
////       << LDD_seeds_bird_resident                                 << ','
////       << LDD_seeds_bird_immigrant                                << ','
////       << total_plants                                            << ','
////       << immigrants_he2004                                       << ','
////       << immigrants_he2004_birds                                 << ','
////       << immigrants_he2004_wind                                  << ','
////       << immigrants                                              << ','
////       << immigrants_birds                                        << ','
////       << immigrants_wind                                         << ',';
////    if (total_plants >= 1.0L)
////    {
////        ss << immigrants_he2004 / total_plants * 100              << ','
////           << immigrants / total_plants * 100                     << ',';
////    }
////    else
////    {
////        ss << 0                                                   << ','
////           << 0                                                   << ',';
////    }
////    ss << Equation::mean(gene_diversity_generation)               << ','
////       << Equation::sd(gene_diversity_generation) << ','
////       << Equation::mean(gene_diversity)                          << ','
////       << Equation::sd(gene_diversity)            << ','
////       << (total_populations - count_extinct_pop) /
////          total_populations * 100                                 << '\n';

////    print_row(ss, m_file_stream);
//}

////-----------------------------------------------------------------------------

////void Output::print_connectivity(
////        const std::vector<std::shared_ptr<Population>>& metapopulation,
////        const int sim_replicate_num,
////        const int run_num)
////{
////    assert(m_sim_output == Sim_output::connectivity);

//////    for (const auto& pop : metapopulation)
//////    {
//////        for (const auto& cohort : pop->get_cohort_list())
//////        {
//////            auto v = cohort->get_genetic();
//////            remove_consecutive_duplicates(v);
//////            cohort->set_genetic(v);
//////        }

//////        int coh_a_index {0};
//////        int coh_b_index {0};

//////        // remove duplicate cohorts with the same genetics
//////        // after removing the consecutive duplicates
//////        for (const auto& coh_a : pop->get_cohort_list())
//////        {
//////            ++coh_a_index;
//////            for (const auto& coh_b : pop->get_cohort_list())
//////            {
//////                ++coh_b_index;
//////                if (coh_a_index < coh_b_index)
//////                {
//////                    if (coh_a->get_genetic() == coh_b->get_genetic())
//////                    {
//////                        coh_a->set_num_plants(0);
//////                    }
//////                }
//////            }
//////            coh_b_index = 0;
//////        }

//////        pop->remove_empty_cohorts();
//////    }

//////    int cohort_id    {0};
//////    int cohort_index {0};

//////    std::ostringstream ss;

//////    for (const auto& pop : metapopulation)
//////    {
//////        cohort_index = 0;
//////        for (const auto& coh : pop->get_cohort_list())
//////        {
//////            ++cohort_index;
//////            ++cohort_id;

//////            for (const auto& genetic : coh->get_genetic())
//////            {
//////                for (const auto& pop_2 : metapopulation)
//////                {
//////                    if (pop_2->get_id() == genetic)
//////                    {
//////                        Cell cell_rand = pop_2->get_random_cell_location();

//////                        ss << Params::Sim::id()                           << ','
//////                           << sim_replicate_num                           << ','
//////                           << Params::Scenario::id()                      << ','
//////                           << run_num                                     << ','
//////                           << pop->get_id()                               << ','
//////                           << cohort_index                                << ','
//////                           << cohort_id                                   << ','
//////                           << genetic                                     << ','
//////                           << pop_2->get_cell_locations().front().get_row() << ','
//////                           << pop_2->get_cell_locations().front().get_col() << ','
//////                           << cell_rand.get_row()                         << ','
//////                           << cell_rand.get_col()                         << '\n';

//////                        print_row(ss, m_file_stream);

//////                        // remove all the values of the row to read next row
//////                        ss.str("");
//////                        ss.clear();
//////                    }
//////                }
//////            }
//////        }
//////    }
////}

////-----------------------------------------------------------------------------

////void Output::print_connectivity_dynamics(
////        const std::shared_ptr<Population>& target_population,
////        const std::vector<std::shared_ptr<Population>>& metapopulation,
////        const int sim_replicate_num,
////        const int run_num,
////        const int run_year,
////        const std::string& rescue_effect,
////        const std::string& plant_stage)
////{
////    assert(m_sim_output == Sim_output::connectivity_dynamics);

////////    const std::vector<std::string> connectivity_dynamics_header
////////    {
////////         "sim_id", "sim_rep", "scn_id", "run_num", "run_year",
////////         "from_id", "to_id", "rescue", "stage",
////////         "from_row", "from_col", "to_row", "to_col"
////////     };

//////    for (const auto& coh : target_population->get_cohort_list())
//////    {
////////        for (const auto& genetic : coh->get_genetic())
////////        {
////////            if (coh->get_genetic().size() == 1)
////////            {
////////                continue;
////////            }

//////            if (coh->get_gene_flow().get_geographical_previous() != target_population->get_id())
//////            {
//////                std::ostringstream ss;
//////                for (const auto& population : metapopulation)
//////                {
////////                    if (population->get_id() == genetic)
////////                    {
////////                        ss << Params::Sim::id()                           << ','
////////                           << sim_replicate_num                           << ','
////////                           << Params::Scenario::id()                      << ','
////////                           << run_num                                     << ','
////////                           << run_year                                    << ',';

//////                        ///@todo TODO finish here!
//////                        /// Maybe try first to plot in R with fake data!!
//////                        /// frow row, from_col, etc.the place of each population
//////                        /// I could make an extra excel file with the row and col for
//////                        /// each population id, and then join the two tables to geom_path
//////                        /// also, I could make an extra landscape with only number in the first
//////                        /// cell of each population, OR, can I print the number
//////                        /// with the extra excel here mentioned???

//////                        print_row(ss, m_file_stream);

//////                        // remove all the values of the row to read next row
//////                        ss.str("");
//////                        ss.clear();
////////                    }
//////                }
//////            }
////////        }
//////    }
////}

////-----------------------------------------------------------------------------

//void Output::print_rescue_over_time(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num,
//        const int run_year,
//        const int time_since_fire)
//{
////    assert(m_sim_output == Sim_output::rescue_over_time);

////    int LDD_resident_populations       {0};
////    int LDD_immigrant_populations      {0};
////    int LDD_recolonization_populations {0};

////    long double LDD_resident_plants       {0.0L};
////    long double LDD_immigrant_plants      {0.0L};
////    long double LDD_recolonization_plants {0.0L};

////    long double total_plants {0.0L};

////    const int total_populations {
////        static_cast<int>(metapopulation.size())
////    };

////    int extinct_populations {0};

////    int total_populations_burned {0};

////    for (const auto& pop : metapopulation)
////    {
////        bool is_resident_pop       {false};
////        bool is_immigrant_pop      {false};
////        bool is_recolonization_pop {false};

////        if (pop->get_time_since_fire() == 0)
////        {
////            ++total_populations_burned;
////        }

////        if (pop->get_is_extinct() == true)
////        {
////            assert(pop->get_cohort_list().empty() == true);
////            ++extinct_populations;
////            continue;
////        }

////        for (const auto& coh : pop->get_cohort_list())
////        {
//////            assert(coh->get_genetic().empty() == false);    // it cannot be empty
////            total_plants += coh.get_num_plants();

////            if (coh.get_is_LDD_cohort() == true)
////            {
////                if (coh.get_gene_flow().get_geographical_previous() == pop->get_id())
////                {
////                    if (is_resident_pop == false)
////                    {
////                        ++LDD_resident_populations;
////                        is_resident_pop = true;
////                    }
////                    LDD_resident_plants += coh.get_num_plants();
////                }
////                else
////                {
////                    if (pop->get_is_recolonized() == true)
////                    {
////                        if (is_recolonization_pop == false)
////                        {
////                            ++LDD_recolonization_populations;
////                            is_recolonization_pop = true;
////                        }
////                        LDD_recolonization_plants += coh.get_num_plants();
////                    }
////                    else
////                    {
////                        if (is_immigrant_pop == false)
////                        {
////                            ++LDD_immigrant_populations;
////                            is_immigrant_pop = true;
////                        }
////                       LDD_immigrant_plants += coh.get_num_plants();
////                    }
////                }
////            }
////        }
////    }

////    std::ostringstream ss;

////    ss << Params::Sim::id()                                     << ','
////       << sim_replicate_num                                     << ','
////       << Params::Scenario::id()                                << ','
////       << run_num                                               << ','
////       << run_year                                              << ','
////       << time_since_fire                                       << ','
////       << total_populations_burned                              << ','
////       << static_cast<int>(Params::Scenario::mort_scenario())   << ','
////       << Params::Scenario::LDD_wind_prop()                     << ','
////       << Params::Scenario::LDD_birds_prop()                    << ','
////       << Params::Scenario::LDD_birds_a()                       << ','
////       << total_populations                                     << ','
////       << extinct_populations                                   << ','
////       << LDD_resident_populations                              << ','
////       << LDD_recolonization_populations                        << ','
////       << LDD_immigrant_populations                             << ','
////       << total_plants                                          << ','
////       << LDD_resident_plants                                   << ','
////       << LDD_recolonization_plants                             << ','
////       << LDD_immigrant_plants                                  << '\n';

////    print_row(ss, m_file_stream);
//}

////-----------------------------------------------------------------------------

//void Output::print_rescue_plant_stages(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num)
//{
////    assert(m_sim_output == Sim_output::rescue_plant_stages);

////    for (const auto& pop : metapopulation)
////    {
////        std::ostringstream ss;

////        ss << Params::Sim::id()                            << ','
////           << sim_replicate_num                            << ','
////           << Params::Scenario::id()                       << ','
////           << run_num                                      << ','
////           << static_cast<int>(Params::Scenario::mort_scenario()) << ','
////           << Params::Scenario::LDD_wind_prop()                 << ','
////           << Params::Scenario::LDD_birds_prop()                << ','
////           << Params::Scenario::LDD_birds_a()                   << ','
////           << pop->get_id()                                << ','
////           << pop->get_total_LDD_residents_germination()   << ','
////           << pop->get_total_LDD_residents_seedling()      << ','
////           << pop->get_total_LDD_residents_adult()         << ','
////           << pop->get_total_LDD_residents_source()        << ','
////           << pop->get_total_immigrants_germination()      << ','
////           << pop->get_total_immigrants_seedling()         << ','
////           << pop->get_total_immigrants_adult()            << ','
////           << pop->get_total_immigrants_source()           << ','
////           << pop->get_total_recolonizations_germination() << ','
////           << pop->get_total_recolonizations_seedling()    << ','
////           << pop->get_total_recolonizations_adult()       << ','
////           << pop->get_total_recolonizations_source()      << '\n';

////        print_row(ss, m_file_stream);
////    }
//}

////-----------------------------------------------------------------------------

//void Output::print_distance(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num)
//{
////    assert(m_sim_output == Sim_output::distance);

////    for (const auto& pop : metapopulation)
////    {
////        std::map<const int, int> joint_id_cohorts;

////        for (const auto& coh : pop->get_cohort_list())
////        {
////            if (coh.get_gene_flow().get_geographical_previous() != pop->get_id())
////            {
////                auto it = joint_id_cohorts.find(
////                              coh.get_gene_flow().get_genetic());

////                if (it != std::end(joint_id_cohorts))
////                {
////                    it->second += coh.get_num_plants();
////                }
////                else
////                {
////                    joint_id_cohorts.emplace(coh.get_gene_flow().get_genetic(),
////                                            coh.get_num_plants());
////                }
////            }
////        }

////        for (const auto& joint_coh : joint_id_cohorts)
////        {
////            std::ostringstream ss;

////            ss << Params::Sim::id()                                   << ','
////               << sim_replicate_num                                   << ','
////               << Params::Scenario::id()                              << ','
////               << run_num                                             << ','
////               << static_cast<int>(Params::Scenario::mort_scenario()) << ','
////               << Params::Scenario::LDD_wind_prop()                   << ','
////               << Params::Scenario::LDD_birds_prop()                  << ','
////               << Params::Scenario::LDD_birds_a()                     << ','
////               << pop->get_id()                                       << ','
////               << joint_coh.first                                     << ','
////               << joint_coh.second                                    << ','
////               << pop->get_dist_edge_to_edge(joint_coh.first)         << ','
////               << pop->get_dist_midpoint_to_edge(joint_coh.first)     << ','
////               << pop->get_dist_midpoint_to_midpoint(joint_coh.first) << '\n';

////            print_row(ss, m_file_stream);
////        }
////    }
//}

////-----------------------------------------------------------------------------

//void Output::print_distance_LDD_vector(
//        const std::vector<std::shared_ptr<Population>>& metapopulation,
//        const int sim_replicate_num,
//        const int run_num)
//{
////    assert(m_sim_output == Sim_output::distance_LDD_vector);

////    for (const auto& pop : metapopulation)
////    {
////        std::map<const int, int> joint_id_cohorts_birds;
////        std::map<const int, int> joint_id_cohorts_wind;

////        for (const auto& coh : pop->get_cohort_list())
////        {
////            if (coh.get_gene_flow().get_geographical_previous() != pop->get_id())
////            {
////                switch (coh.get_dispersal_vector())
////                {
////                    case Dispersal_vector::LDD_bird :
////                    {
////                        auto it = joint_id_cohorts_birds.find(
////                                        coh.get_gene_flow().get_genetic());

////                        if (it != std::end(joint_id_cohorts_birds))
////                        {
////                            it->second += coh.get_num_plants();
////                        }
////                        else
////                        {
////                            joint_id_cohorts_birds.emplace(
////                                        coh.get_gene_flow().get_genetic(),
////                                        coh.get_num_plants());
////                        }
////                        break;
////                    }
////                    case Dispersal_vector::LDD_wind :
////                    {
////                        auto it = joint_id_cohorts_wind.find(
////                                        coh.get_gene_flow().get_genetic());

////                        if (it != std::end(joint_id_cohorts_wind))
////                        {
////                            it->second += coh.get_num_plants();
////                        }
////                        else
////                        {
////                            joint_id_cohorts_wind.emplace(
////                                        coh.get_gene_flow().get_genetic(),
////                                        coh.get_num_plants());
////                        }
////                        break;
////                    }
////                    case Dispersal_vector::SDD_wind :
////                        // nothing
////                        break;
////                }
////            }
////        }

////        for (const auto& joint_coh : joint_id_cohorts_birds)
////        {
////            std::ostringstream ss;

////            ss << Params::Sim::id()                                   << ','
////               << sim_replicate_num                                   << ','
////               << Params::Scenario::id()                              << ','
////               << run_num                                             << ','
////               << static_cast<int>(Params::Scenario::mort_scenario()) << ','
////               << Params::Scenario::LDD_wind_prop()                   << ','
////               << Params::Scenario::LDD_birds_prop()                  << ','
////               << Params::Scenario::LDD_birds_a()                     << ','
////               << "birds"                                             << ','
////               << pop->get_id()                                       << ','
////               << joint_coh.first                                     << ','
////               << joint_coh.second                                    << ','
////               << pop->get_dist_edge_to_edge(joint_coh.first)         << ','
////               << pop->get_dist_midpoint_to_edge(joint_coh.first)     << ','
////               << pop->get_dist_midpoint_to_midpoint(joint_coh.first) << '\n';

////            print_row(ss, m_file_stream);
////        }

////        for (const auto& joint_coh : joint_id_cohorts_wind)
////        {
////            std::ostringstream ss;

////            ss << Params::Sim::id()                                   << ','
////               << sim_replicate_num                                   << ','
////               << Params::Scenario::id()                              << ','
////               << run_num                                             << ','
////               << static_cast<int>(Params::Scenario::mort_scenario()) << ','
////               << Params::Scenario::LDD_wind_prop()                   << ','
////               << Params::Scenario::LDD_birds_prop()                  << ','
////               << Params::Scenario::LDD_birds_a()                     << ','
////               << "wind"                                              << ','
////               << pop->get_id()                                       << ','
////               << joint_coh.first                                     << ','
////               << joint_coh.second                                    << ','
////               << pop->get_dist_edge_to_edge(joint_coh.first)         << ','
////               << pop->get_dist_midpoint_to_edge(joint_coh.first)     << ','
////               << pop->get_dist_midpoint_to_midpoint(joint_coh.first) << '\n';

////            print_row(ss, m_file_stream);
////        }
////    }
//}
//}

