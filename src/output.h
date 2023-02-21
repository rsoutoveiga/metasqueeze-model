#ifndef OUTPUT_H
#define OUTPUT_H


#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <cassert>
#include <unordered_set>
#include <map>
#include <cmath>
#include <stdexcept>  // std::runtime_error()

#include <filesystem>
#include <iomanip>
#include <ctime>

#include "parameters.h"
#include "population.h"
#include "grid.h"
#include "global.h"


///
/// \brief The Output class generate the output files and print results
///
class Output
{
public:
    Output(
        const std::string& output_folder
        , const std::string& name_sim_file
        , const std::string& time
    );

    ~Output();

    void setup_output(
            const Parameters& params
            , const std::string& name_sim_folder
         );

    void cleanup();

    void print_parameters(const Parameters& params);

    void print_study_area(
            const std::vector<std::vector<int>>& study_area
            , const Parameters & params
            , const unsigned study_num);

    void print_fire_output(
            const Parameters& params
            , const unsigned study_num
            , const int run_num
            , const int year
            , const std::map<int, Population>& metapopulation);

    void print_seed_dynamics(
            const Parameters& params
            , const unsigned study_num
            , const int run_num
            , const int year
            , const int time_since_fire
            , const std::map<int, Population>& metapopulation);

    void print_metapopulation(
            const Parameters& params
            , const unsigned study_num
            , const int run_num
            , const int year_max
            , const std::map<int, Population>& metapopulation
            , std::vector<int>& occupied_patches);

    void print_metapopulation_full(
            const Parameters& params
            , const unsigned study_num
            , const int run_num
            , const int year_max
            , const std::map<int, Population>& metapopulation
            , const std::vector<int>& occupied_patches
            , const std::vector<double>& immigrants_percent_postfire
            , const std::vector<double>& residents_percent_postfire
            , const std::vector<int>& immigrants_postfire
            , const std::vector<int>& residents_postfire);

    void print_metapopulation_dynamics(
            const Parameters& params,
            const unsigned study_num,
            const int run_num,
            const int year,
            const int time_since_fire,
            const std::map<int, Population>& metapopulation);

    void print_metapopulation_pops(
            const Parameters& params,
            const unsigned study_num,
            const int run_num,
            const int year_max,
            const std::map<int, Population>& metapopulation);

    void print_rescue_effects(
            const Parameters& params,
            const unsigned study_num,
            const int run_num,
            const int year_max,
            const std::map<int, Population>& metapopulation);

    void print_plants_per_pop(
            const Parameters& params
            , const unsigned study_num);

    void print_dispersal(
            const Parameters& params
            , const unsigned study_num
            , const int run_num
            , const std::map<int, Population>& metapopulation);

    void init_plants_per_population(
            const std::map<int, Population>& metapopulation);

    void add_plants_per_population(
                const std::map<int, Population>& metapopulation);


private:
    void check_file_exists(const std::string& file_name);
    void check_file_exists(const std::filesystem::path& path);
    void print_header(
            const std::vector<std::string>& header
            , std::ofstream& stream);
    void print_row(const std::ostringstream& ss, std::ofstream& stream);


    void check_file_is_open(
            const std::ofstream& ofs,
            const std::filesystem::path& path);
    void check_file_is_empty(const std::filesystem::path& path);


    const std::string m_output_folder;
    const std::string m_raw_folder {"/raw/"};
    const std::string m_name_sim_file;
    const std::string m_time;

//    std::map<int, std::ofstream> m_seed_dynamics_streams;
//    std::map<int, std::ofstream> m_persistence_streams;
    std::ofstream m_seed_dynamics_stream;
    std::ofstream m_fire_output_stream;
    std::ofstream m_metapop_stream;
    std::ofstream m_metapop_full_stream;
    std::ofstream m_metapop_dynamics_stream;
    std::ofstream m_metapop_pops_stream;
    std::ofstream m_rescue_effects_stream;
    std::ofstream m_plants_per_pop_stream;
    std::ofstream m_dispersal_stream;

//    std::map<int, std::filesystem::path> m_seed_dynamics_paths;
//    std::map<int, std::filesystem::path> m_persistence_paths;
    std::filesystem::path m_seed_dynamics_path;
    std::filesystem::path m_fire_output_path;
    std::filesystem::path m_metapop_path;
    std::filesystem::path m_metapop_full_path;
    std::filesystem::path m_metapop_dynamics_path;
    std::filesystem::path m_metapop_pops_path;
    std::filesystem::path m_rescue_effects_path;
    std::filesystem::path m_plants_per_pop_path;
    std::filesystem::path m_dispersal_path;


    //--------------------------------------------------------------------
//    std::map<int, std::vector<Total_plants_stage>> m_plants_per_pop_out;

    std::map<int, std::map<int, std::vector<Total_plants_stage>>> m_plants_per_pop_out;

//    bool is_file_exist(const char *file_name);

//    void print_parameters(const Parameters& params);


//    void remove_consecutive_duplicates(std::vector<int>& v)
//    {
//        auto last = std::unique(v.begin(), v.end());
//        v.erase(last, v.end());
//    }

//    void remove_all_duplicates(std::vector<int> &v)
//    {
//        std::unordered_set<int> s;
//        auto end = std::remove_if(v.begin(), v.end(),
//                                  [&s](const int &i) {
//            return !s.insert(i).second;
//        });

//        v.erase(end, v.end());
//    }

};

#endif // OUTPUT_H
