#include "parameter_reader.h"


void Parameter_reader::read_simulation(const std::string& line)
{
    int sim_id {0};
    bool is_study_replicates {false};
    std::string study_replicates;
    int sim_repetitions {0};
    int run_time_max    {0};
    bool is_out_study_area              {false};
    bool is_out_fire                    {false};
    bool is_out_seed_dynamics           {false};
    bool is_out_metapopulation          {false};
    bool is_out_metapopulation_full     {false};
    bool is_out_metapopulation_dynamics {false};
    bool is_out_metapopulation_pops     {false};
    bool is_out_recol_rescue            {false};
    bool is_out_plants_per_pop          {false};
    bool is_out_dispersal               {false};

    std::string file_scenario;

    std::stringstream ss_sim(line);

    ss_sim  >> sim_id
            >> is_study_replicates
            >> study_replicates
            >> sim_repetitions
            >> run_time_max
            >> is_out_study_area
            >> is_out_fire
            >> is_out_seed_dynamics
            >> is_out_metapopulation
            >> is_out_metapopulation_full
            >> is_out_metapopulation_dynamics
            >> is_out_metapopulation_pops
            >> is_out_recol_rescue
            >> is_out_plants_per_pop
            >> is_out_dispersal
            >> file_scenario;

    m_sim_id = sim_id;
    m_is_study_replicates = is_study_replicates;

    if (m_is_study_replicates)
    {
        read_study_replicates(study_replicates);
    }
    else
    {
        set_study_replicates(study_replicates);
    }

    m_sim_repetitions                = sim_repetitions;
    m_run_time_max                   = run_time_max;
    m_is_out_study_area              = is_out_study_area;
    m_is_out_fire                    = is_out_fire;
    m_is_out_seed_dynamics           = is_out_seed_dynamics;
    m_is_out_metapopulation          = is_out_metapopulation;
    m_is_out_metapopulation_full     = is_out_metapopulation_full;
    m_is_out_metapopulation_dynamics = is_out_metapopulation_dynamics;
    m_is_out_metapopulation_pops     = is_out_metapopulation_pops;
    m_is_out_recol_rescue            = is_out_recol_rescue;
    m_is_out_plants_per_pop          = is_out_plants_per_pop;
    m_is_out_dispersal               = is_out_dispersal;

    read_scenario(file_scenario);

    // read parameters to be changed
//    std::map<std::string, std::string> parameters;
    std::vector<std::pair<std::string, std::string>> parameters;

    {
        std::string buffer;  // temporal column value
        std::string param_name;
        bool is_param_value {false};

        while (ss_sim >> buffer)
        {
            if (is_param_value == true)
            {
                parameters.emplace_back(param_name, buffer);
                is_param_value = false;
                continue;
            }

            param_name = buffer;
            is_param_value = true;
        }
    }

    read_changed_parameters(parameters);   

    calculate_canopy_params();  // derived parameters

    ///@note NOTE fire size is in cell units!!!
    m_fire_size_x = m_fire_size_x / m_cell_size;
    m_fire_size_y = m_fire_size_y / m_cell_size;


    if (m_mort_recruit_postfire_min > m_mort_recruit_postfire_max)
    {
        m_mort_recruit_postfire_min = m_mort_recruit_postfire_max;
    }

    if (m_mort_recruit_postfire_max < m_mort_recruit_postfire_min)
    {
        m_mort_recruit_postfire_max = m_mort_recruit_postfire_min;
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::set_study_replicates(const std::string& total_replicates)
{
    // Expects
    assert(!m_is_study_replicates);
    assert(std::empty(m_study_replicates));

    const int total_study_replicates { std::stoi(total_replicates) };
    assert(total_study_replicates > 0);

    for (int i = 0; i < total_study_replicates; ++i)
    {
        unsigned study_seed {0};

        // find a seed that does not repeat itself
        do
        {
            study_seed = Global::prng.uniform_unsigned(
                            0
                            , std::numeric_limits<unsigned>::max()
                         );

        } while (std::find(
                    std::cbegin(m_study_replicates)
                    , std::cend(m_study_replicates)
                    , study_seed
                 ) != std::cend(m_study_replicates));

        m_study_replicates.push_back(study_seed);
    }

    // Ensures
    assert(std::size(m_study_replicates) == static_cast<unsigned>(total_study_replicates));
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_study_replicates(const std::string& name_file)
{
    // Expects
    assert(m_is_study_replicates);
    assert(std::empty(m_study_replicates));

    const std::string name_folder {Global::project_directory + "data/in/study_rep/"};
    std::ifstream ifs(name_folder + name_file);

    std::string line;
    std::getline(ifs, line);  // skip header

    // checking if the file has been opened correctly
    if (!ifs)
    {
        std::cerr << "ERROR! The following file could not be opened:\n"
                  << name_folder << name_file << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    while (std::getline(ifs, line))   // read simulation
    {
        bool is_random      {false};
        unsigned study_seed {0};

        std::stringstream ss(line);

        ss >> is_random
           >> study_seed;

        if (is_random == true)
        {
            // find a seed that does not repeat itself
            do
            {
                study_seed = Global::prng.uniform_unsigned(
                                0,
                                std::numeric_limits<unsigned>::max()
                             );

            } while (std::find(
                        std::cbegin(m_study_replicates)
                        , std::cend(m_study_replicates)
                        , study_seed
                     ) != std::cend(m_study_replicates));
        }

        m_study_replicates.push_back(study_seed);
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_scenario(const std::string& file_scenario)
{
    const std::string name_folder { Global::project_directory + "data/in/scenarios/" };
    std::ifstream ifs_scn(name_folder + file_scenario);

    // checking if the file has been opened correctly
    if (!ifs_scn)
    {
        std::cerr << "ERROR! The following scenario file could not be opened:\n"
                  << name_folder << file_scenario << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    std::string line_scn;
    std::getline(ifs_scn, line_scn);  // skip header
    std::getline(ifs_scn, line_scn);  // get line/data

    int model_type                 {0};
    std::string file_study_size;
    std::string file_metapop;
    int carrying_capacity          {0};
    double init_cones_plant        {0.0};
    double init_viable_seeds_plant {0.0};
    bool is_climate                {true};
    double long_term_rain          {0};
    int weather_scenario           {0};
    std::string file_climate;
    bool is_fire                   {true};
    std::string file_fire;
    std::string file_dispersal;
    std::string file_species;
    std::string file_habitat;
    std::string file_plant_type;
    std::string file_fuzzy_set_list;
    std::string file_flower_dist;

    std::stringstream ss_scenario(line_scn);

    ss_scenario >> model_type
                >> file_study_size
                >> file_metapop
                >> carrying_capacity
                >> init_cones_plant
                >> init_viable_seeds_plant
                >> is_climate
                >> long_term_rain
                >> weather_scenario
                >> file_climate
                >> is_fire
                >> file_fire
                >> file_dispersal
                >> file_species
                >> file_habitat
                >> file_plant_type
                >> file_fuzzy_set_list
                >> file_flower_dist;

    switch (model_type)
    {
        case 0 :
            m_model_type = Model_type::cohort_based;
            break;
        case 1 :
            m_model_type = Model_type::plant_based;
            break;
        default :
            std::cerr << "The parameter 'model_type' is not valid:\n"
                         "- value '0' for cohort-based modelling\n"
                         "- value '1' for plant-based modelling\n\a";
            std::exit(EXIT_FAILURE);
    }

    m_file_study_size = file_study_size;
    read_study_size(m_file_study_size);

    if (carrying_capacity <= 0)
    {
        std::cout << "WARNING: the parameter 'carrying_capacity' is zero\n\a";
        carrying_capacity = 0;
    }

    m_file_metapopulation     = file_metapop;
    m_carrying_capacity       = carrying_capacity;
    m_init_cones_plant        = init_cones_plant;
    m_init_viable_seeds_plant = init_viable_seeds_plant;
    m_is_climate              = is_climate;
    m_long_term_rain          = long_term_rain;

    switch (weather_scenario)
    {
    case 0:
        m_weather_scenario = Weather_scenario::random;
        break;
    case 1:
        m_weather_scenario = Weather_scenario::consecutive;
        break;
    default:
        throw std::runtime_error(
            "The parameter 'weather_scenario' is not valid:\n"
            "- value '0' for random scenario\n"
            "- value '1' for consecutive scenario \n\a"
        );
    }

    m_file_climate            = file_climate;
    m_is_fire                 = is_fire;
    m_file_fire               = file_fire;
    m_file_dispersal          = file_dispersal;
    m_file_species            = file_species;
    m_file_habitat            = file_habitat;
    m_file_plant_type         = file_plant_type;
    m_file_fuzzy_set_list     = file_fuzzy_set_list;
    m_file_flower_dist        = file_flower_dist;

    if (m_is_fire == true)
    {
        read_fire(m_file_fire);
    }

    read_dispersal(m_file_dispersal);
    read_species(m_file_species);

    ///@todo TODO add parameter whether the the user activate habitat quality or not?!
    read_habitat_quality(m_file_habitat);

    if (m_model_type == Model_type::plant_based)
    {
        read_plant_type(m_file_plant_type);
        read_fuzzy_set_list(m_file_fuzzy_set_list);
        read_flower_distribution(m_file_flower_dist);
    }
}


//-----------------------------------------------------------------------------

void Parameter_reader::read_study_size(const std::string& file_study_size)
{
    const std::string name_folder {Global::project_directory + "data/in/study_size/"};
    std::ifstream ifs_study(name_folder + file_study_size);

    // checking if the file has been opened correctly
    if (!ifs_study)
    {
        std::cerr << "ERROR! Could not open the fire file:\n"
                  << "Study size file: " << name_folder << file_study_size
                  << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    std::string line_study;
    std::getline(ifs_study, line_study);    // skip header
    std::getline(ifs_study, line_study);    // read data

    ///@note NOTE is this really checking the entire file?
    check_file_is_empty(name_folder, file_study_size, line_study);

    int cell_size    {0};
    int study_size_x {0};
    int study_size_y {0};

    std::stringstream ss_study(line_study);

    ss_study >> cell_size
             >> study_size_x
             >> study_size_y;

    if (study_size_x % cell_size != 0 || study_size_y % cell_size != 0)
    {
        std::cerr << "ERROR! Study size x and y must be "
                     "multiple of cell size\n"
                  << name_folder << file_study_size << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    m_cell_size = cell_size;
    m_study_size_x = study_size_x;
    m_study_size_y = study_size_y;
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_fire(const std::string& file_fire)
{
    const std::string name_folder {Global::project_directory + "data/in/fire/"};
    std::ifstream ifs_fire(name_folder + file_fire);

    // Checking if the file has been opened correctly
    if (!ifs_fire)
    {
        std::cerr << "ERROR! Could not open the fire file:\n"
                  << "Fire file: " << name_folder << file_fire << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    std::string line_fire;
    std::getline(ifs_fire, line_fire);    // skip header
    std::getline(ifs_fire, line_fire);    // read data

    ///@note NOTE is this really checking the entire file?
    check_file_is_empty(name_folder, file_fire, line_fire);

    int fire_interval_scn       {0};
    int fire_size_scn           {0};
    int fire_scale              {0};
    double fire_size_x          {0.0};
    double fire_size_y          {0.0};
    int fire_interval_mean      {0};
    int fire_interval_lower_cut {0};
    int burned_lower_cut        {0};
    int burned_upper_cut        {0};
    double fire_a               {0.0};
    double fire_b               {0.0};

    std::stringstream ss_fire(line_fire);

    ss_fire >> fire_interval_scn
            >> fire_size_scn
            >> fire_scale
            >> fire_size_x
            >> fire_size_y
            >> fire_interval_mean
            >> fire_interval_lower_cut
            >> burned_lower_cut
            >> burned_upper_cut
            >> fire_a
            >> fire_b;

    switch (fire_interval_scn)
    {
        case 0 :
            m_fire_interval_scn = Fire_scenario::deterministic;
            break;
        case 1 :
            m_fire_interval_scn = Fire_scenario::truncated_normal;
            break;
        case 2 :
            m_fire_interval_scn = Fire_scenario::weibull;
        break;
        default :
            std::cerr << "The parameter 'fire_scenario' is not valid:\n"
                         "- value '0' for deterministic fire events\n"
                         "- value '1' for truncated_normal fire events\n"
                         "- value '2' for Weibull distribution\n\a";
            std::exit(EXIT_FAILURE);
    }

    switch (fire_size_scn)
    {
        case 0 :
            m_fire_size_scn = Fire_scenario::deterministic;
            break;
        case 1 :
            m_fire_size_scn = Fire_scenario::truncated_normal;
            break;
        default :
            std::cerr << "The parameter 'fire_scenario' is not valid:\n"
                         "- value '0' for deterministic fire events\n"
                         "- value '1' for truncated_normal fire events\n\a";
            std::exit(EXIT_FAILURE);
    }

    switch (fire_scale)
    {
        case 0 :
            m_fire_scale = Fire_scale::patchy;
            break;
        case 1 :
            m_fire_scale = Fire_scale::study_area;
            break;
        default :
            std::cerr << "The parameter 'fire_scale' is not valid:\n"
                         "- value '0' for patchy fires\n"
                         "- value '1' for fire in the entire study area\n\a";
            std::exit(EXIT_FAILURE);
    }

    ///@note NOTE fire size is in cell units!!!
    m_fire_size_x = fire_size_x;
    m_fire_size_y = fire_size_y;

    m_fire_interval_mean = fire_interval_mean;
    m_fire_interval_lower_cut = fire_interval_lower_cut;

    if (burned_lower_cut < 0 || burned_upper_cut < 0)
    {
        throw std::runtime_error(
                "parameters 'burned_lower_cut' and 'burned_upper_cut' "
                "must be positive");
    }

    if (burned_lower_cut >= burned_upper_cut)
    {
        throw std::runtime_error(
                "parameters 'burned_lower_cut' must be lower than 'burned_upper_cut'");
    }

    m_burned_lower_cut = burned_lower_cut;
    m_burned_upper_cut = burned_upper_cut;

    m_fire_a = fire_a;
    m_fire_b = fire_b;
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_dispersal(const std::string& file_dispersal)
{
    const std::string name_folder {Global::project_directory + "data/in/dispersal/"};
    std::ifstream ifs_dispersal(name_folder + file_dispersal);

    // Checking if the file has been opened correctly
    if (!ifs_dispersal)
    {
        std::cerr << "ERROR! Could not open the fire file:\n"
                  << "Fire file: " << name_folder << file_dispersal << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    std::string line_dispersal;
    std::getline(ifs_dispersal, line_dispersal);    // skip header
    std::getline(ifs_dispersal, line_dispersal);    // read data

    ///@note NOTE is this really checking the entire file?
    check_file_is_empty(name_folder, file_dispersal, line_dispersal);

//    double wind_max          {0.0};
    double wind_prop         {0.0};
    double wind_a            {0.0};
    double wind_b            {0.0};
    int wind_direction       {0};
    double birds_dist_max    {0.0};
    double birds_prop        {0.02};         // proportion of cones dispersed by birds
    int follicles_distr_type {0};
    double follicles_distr_a {6.219};
    double follicles_distr_b {10.08};
//    int follicles_max        {22};
    double postfire_follicles_open {0.5};

    std::stringstream ss_dispersal(line_dispersal);

    ss_dispersal >> wind_prop
                 >> wind_a
                 >> wind_b
                 >> wind_direction
                 >> birds_dist_max
                 >> birds_prop
                 >> follicles_distr_type
                 >> follicles_distr_a
                 >> follicles_distr_b
//                 >> follicles_max
                 >> postfire_follicles_open;

//    m_wind_max = wind_max;
    m_wind_prop = wind_prop;
    m_wind_a = wind_a;
    m_wind_b = wind_b;

    switch (wind_direction)
    {
        case 0 :
            m_wind_direction = Wind_direction::random;
            break;
        case 1 :
            m_wind_direction = Wind_direction::geraldton;
            break;
        default :
            std::cerr << "The parameter 'wind_direction' is not valid:\n"
                         "- value '0' for random 360 degrees\n"
                         "- value '1' for wind roses from Geraldton station\n\a";
            std::exit(EXIT_FAILURE);
    }

    m_birds_dist_max = birds_dist_max;
    m_birds_prop = birds_prop;

    switch (follicles_distr_type)
    {
    case 1:
        m_follicles_distr_type = Distribution_type::pois;
        break;
    case 3:
        m_follicles_distr_type = Distribution_type::nbinom_mu;
        break;
    default:
        throw std::runtime_error(
                    "The parameter 'follicles_distr_type' is not valid:\n"
                    "- value '1' for Poisson distribution\n"
                    "- value '3' for negative binomial with 'mu' parameter\n");
    }

    m_follicles_distr_a = follicles_distr_a;
    m_follicles_distr_b = follicles_distr_b;
//    m_follicles_max = follicles_max;

    if (postfire_follicles_open < 0.0 || postfire_follicles_open > 1.0)
    {
        throw std::runtime_error(
                    "the parameter 'postfire_follicles_open'"
                    "must be double between 0 and 1");
    }
    m_postfire_follicles_open = postfire_follicles_open;
}
//-----------------------------------------------------------------------------

void Parameter_reader::read_species(const std::string& file_species)
{
    std::string name_folder {Global::project_directory + "data/in/species/"};
    std::ifstream ifs_species(name_folder + file_species);

    // checking if the file has been opened correctly
    if (!ifs_species)
    {
        std::cerr << "ERROR! Could not open the scenario file:\n"
                  << "Scenario file name: "
                  << name_folder << file_species << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    std::string line_species;
    std::getline(ifs_species, line_species); // skip header
    std::getline(ifs_species, line_species); // read data

    check_file_is_empty(name_folder, file_species, line_species);

    int plant_longevity {0};
    int plant_age_young {0};
    int plant_age_adult {0};

    int mort_scenario {0};

    double mort_recruit_postfire_min  {0.0};
    double mort_recruit_postfire_max  {0.0};
    double mort_recruit_postfire_mean {0.0};  // mort of plants at the age of zero
    double recruit_interfire {0.0};  // proportion of post-fire
    double recruit_weather   {0.0};  // weather impacts in seedling recruitment

    int senescence_age         {0};
    double senescence_increase {0.0};

    double mort_min {0.0};

//    double plant_mort_age_a {0.0};  // plant_mort_age_a
//    double plant_mort_age_b {0.0};  // plt_mort_age_b
//    double plant_mort_age_c {0.0};

//    double plant_mort_weather_a {0.0};
//    double plant_mort_weather_b {0.0};

    double mort_a {0.0553}; // old plant_mort_age_a
    double mort_b {0.2645}; // old plant_mort_age_b
    double mort_c {0.179};  // old plant_mort_age_c

    double mort_d {-0.00061}; // old plant_mort_weather_a
    double mort_e {0.26873};  // old plant_mort_weather_b
    double mort_f {0.0};      // new

    int cone_cycle     {0};
    int seed_longevity {0};

    int flower_scenario {0};
//    double flower_max   {0.0};

    double flower_age_a {0.0};
    double flower_age_b {0.0};
    double flower_age_c {0.0};

    double flower_weather_a {0.0};
    double flower_weather_b {0.0};
    double flower_weather_c {0.0};
    double flower_weather_d {0.0};
    double flower_weather_e {0.0};

    double pollination_success {0.0};
    double num_follicles       {0.0};
    double num_seeds           {0.0};
    double firm_seeds          {0.0};
    double viable_seeds        {0.0};

    double insect_a {0.0};
    double insect_b {0.0};

    double decay_a {0.0};
    double decay_b {0.0};

    double open_follicles_a {0.0};
    double open_follicles_b {0.0};
    double open_follicles_c {0.0};

    std::stringstream ss_species(line_species);

    ss_species  >> plant_longevity
                >> plant_age_young
                >> plant_age_adult
                >> mort_scenario
                >> mort_recruit_postfire_min
                >> mort_recruit_postfire_max
                >> mort_recruit_postfire_mean
                >> recruit_interfire
                >> recruit_weather
                >> senescence_age
                >> senescence_increase
                >> mort_min
                >> mort_a
                >> mort_b
                >> mort_c
                >> mort_d
                >> mort_e
                >> mort_f
                >> cone_cycle
                >> seed_longevity
                >> flower_scenario
//                >> flower_max
                >> flower_age_a
                >> flower_age_b
                >> flower_age_c
                >> flower_weather_a
                >> flower_weather_b
                >> flower_weather_c
                >> flower_weather_d
                >> flower_weather_e
                >> pollination_success
                >> num_follicles
                >> num_seeds
                >> firm_seeds
                >> viable_seeds
                >> insect_a
                >> insect_b
                >> decay_a
                >> decay_b
                >> open_follicles_a
                >> open_follicles_b
                >> open_follicles_c;

    ///@todo TODO evaluate input values!!

     m_plant_longevity            = plant_longevity;
     m_plant_age_young            = plant_age_young;
     m_plant_age_adult            = plant_age_adult;
     m_mort_recruit_postfire_min  = mort_recruit_postfire_min;
     m_mort_recruit_postfire_max  = mort_recruit_postfire_max;
     m_mort_recruit_postfire_mean = mort_recruit_postfire_mean;
     m_recruit_interfire          = recruit_interfire;
     m_recruit_weather            = recruit_weather;
     m_senescence_age             = senescence_age;
     m_senescence_increase        = senescence_increase;
     m_mort_min                   = mort_min;
     m_mort_a                     = mort_a;
     m_mort_b                     = mort_b;
     m_mort_c                     = mort_c;
     m_mort_d                     = mort_d;
     m_mort_e                     = mort_e;
     m_mort_f                     = mort_f;

     if (cone_cycle < 0)
     {
         std::cerr << "The parameter 'cone_cycle' cannot be negative:\n";
         std::exit(EXIT_FAILURE);
     }

     if (seed_longevity < 0)
     {
         std::cerr << "The parameter 'seed_longevity' cannot be negative:\n";
         std::exit(EXIT_FAILURE);
     }

     m_cone_cycle     = static_cast<unsigned>(cone_cycle);
     m_seed_longevity = static_cast<unsigned>(seed_longevity);

     switch (mort_scenario)
     {
     case 0 :
         m_mort_scenario = Mortality_scenario::age_weather_relative;
         break;
     case 1 :
         m_mort_scenario = Mortality_scenario::age_weather_absolute;
         break;
     case 2 :
         m_mort_scenario = Mortality_scenario::spring_autumn_mean;
         break;
     case 3 :
         m_mort_scenario = Mortality_scenario::SDD_autumn_LDD_spring;
         break;
     case 4 :
         m_mort_scenario = Mortality_scenario::resi_autum_immi_spring;
         break;
     case 5 :
         m_mort_scenario = Mortality_scenario::lowest_mortality;
         break;
     default :
         throw std::runtime_error(
                 "The parameter 'mort_scenario' is not valid:\n"
                 "- value '0': relative scenario from Keith et al. 2014\n"
                 "- value '1': absolute scenario from Keith et al. 2014\n"
                 "- value '2': arithmetic mean of spring and autumn recruit mortality\n"
                 "- value '3': SDD autumn and LDD spring\n"
                 "- value '4': residents autumn and immigrants spring\n"
                 "- value '5': lowest mortality rates\n\a");
     }

     switch (flower_scenario)
     {
         case 0 :
             m_flower_scenario = Flower_scenario::baseline;
             break;
         case 1 :
             m_flower_scenario = Flower_scenario::current;
             break;
         case 2 :
             m_flower_scenario = Flower_scenario::current_nb;
             break;
         default :
             std::cerr << "The parameter 'flower_scenario' is not valid:\n"
                          "- value '0' for baseline scenario\n"
                          "- value '1' for current scenario using lmer\n"
                          "- value '2' for current scenario using glmer.nb\n\a";
             std::exit(EXIT_FAILURE);
     }

//     m_flower_max          = flower_max;
     m_flower_age_a        = flower_age_a;
     m_flower_age_b        = flower_age_b;
     m_flower_age_c        = flower_age_c;
     m_flower_weather_a    = flower_weather_a;
     m_flower_weather_b    = flower_weather_b;
     m_flower_weather_c    = flower_weather_c;
     m_flower_weather_d    = flower_weather_d;
     m_flower_weather_e    = flower_weather_e;
     m_pollination_success = pollination_success;
     m_num_follicles       = num_follicles;
     m_num_seeds           = num_seeds;
     m_firm_seeds          = firm_seeds;
     m_viable_seeds        = viable_seeds;
     m_insect_a            = insect_a;
     m_insect_b            = insect_b;
     m_decay_a             = decay_a;
     m_decay_b             = decay_b;
     m_open_follicles_a    = open_follicles_a;
     m_open_follicles_b    = open_follicles_b;
     m_open_follicles_c    = open_follicles_c;
}

//-----------------------------------------------------------------------------

void Parameter_reader::calculate_canopy_params()
{
    // clear vectors to not accumulate elements each time reading scenario
    m_viable_seeds_cone_age.clear();
    m_viable_seeds_dispersed_cone_age.clear();
    // these vector will have the proportion of seed loss per cone age per fertile cone
    // These value do not change during the simulation run (i.e., based on cone age)
    // if this value change during the simulation (e.g. weather impacts as flower production and mortality)
    // then this vector should be included as in the old "interfireSeedLossandDispersal" in Cohort
    std::vector<double> insect_age_prop(m_cone_cycle + m_seed_longevity, 0.0);
    std::vector<double> decay_age_prop(m_cone_cycle + m_seed_longevity, 0.0);
    std::vector<double> open_age_prop(m_cone_cycle + m_seed_longevity, 0.0);

    int cone_age {0};
    // calculate the proportions per cone age for insect, decay, and open follicles sub-submodels
    // for (int i = 1; i <= Parameter_reader::seed_longevity; ++i)
    for (unsigned i {m_cone_cycle}; i < std::size(insect_age_prop); ++i)
    {
        ++cone_age;

        if (cone_age > 1)
        {
            insect_age_prop.at(i) =
                    Equation::linear(m_insect_a, m_insect_b, cone_age) -
                    Equation::linear(m_insect_a, m_insect_b, cone_age - 1);

            decay_age_prop.at(i) =
                    Equation::exponential(m_decay_a, m_decay_b, cone_age) -
                    Equation::exponential(m_decay_a, m_decay_b, cone_age - 1);

            open_age_prop.at(i) = Equation::logistic_growth(
                                      m_open_follicles_a,
                                      m_open_follicles_b,
                                      m_open_follicles_c,
                                      cone_age) -
                                  Equation::logistic_growth(
                                      m_open_follicles_a,
                                      m_open_follicles_b,
                                      m_open_follicles_c,
                                      cone_age - 1);
        }
        else if (cone_age == 1)
        {
            insect_age_prop.at(i) = Equation::linear(m_insect_a
                                                     , m_insect_b
                                                     , cone_age);

            decay_age_prop.at(i) = Equation::exponential(m_decay_a
                                                         , m_decay_b
                                                         , cone_age);

            open_age_prop.at(i) = Equation::logistic_growth(
                        m_open_follicles_a
                        , m_open_follicles_b
                        , m_open_follicles_c
                        , cone_age);
        }
    }

    m_insect_cone_age = insect_age_prop;
    m_decay_cone_age = decay_age_prop;
    m_open_follicles_cone_age = open_age_prop;

    // calculate the number of seeds per cone age per fertile cone, and
    // the number of seeds dispersed per cone age per fertile cone.
    m_viable_seeds_cone_age.clear();
    m_viable_seeds_dispersed_cone_age.clear();

    m_viable_seeds_cone_age.assign(m_cone_cycle + m_seed_longevity, 0.0);
    m_viable_seeds_dispersed_cone_age.assign(m_cone_cycle + m_seed_longevity, 0.0);

    const double potential_seeds_in_cone {
        m_num_follicles * m_num_seeds * m_firm_seeds
    };

    m_viable_seeds_cone_age.at(m_cone_cycle) = potential_seeds_in_cone;

//    double seeds_lost      {0.0};
//    double seeds_dispersed {0.0};

//    const unsigned vector_size {m_viable_seeds_cone_age.size()};

    for (unsigned i {m_cone_cycle}; i < std::size(m_viable_seeds_cone_age); ++i)
    {
        double seeds_lost {
            m_viable_seeds_cone_age.at(i) *
                    (insect_age_prop.at(i) + decay_age_prop.at(i))
        };

        if (seeds_lost > m_viable_seeds_cone_age.at(i))
        {
            m_viable_seeds_cone_age.at(i) = 0.0;
        }
        else
        {
            m_viable_seeds_cone_age.at(i) -= seeds_lost;

            double seeds_dispersed {
                m_viable_seeds_cone_age.at(i) * open_age_prop.at(i)
            };

            if (seeds_dispersed > m_viable_seeds_cone_age.at(i))
            {
              seeds_dispersed = m_viable_seeds_cone_age.at(i);
            }

            m_viable_seeds_cone_age.at(i) -= seeds_dispersed;

            // if not in the last element of the vector, copy the number of seeds to
            // the next element for next calculation (i.e. next cone age element)
            if (i < std::size(m_viable_seeds_cone_age) - 1)
            {
                m_viable_seeds_cone_age.at(i + 1) = m_viable_seeds_cone_age.at(i);
            }

            ///@note NOTE I add the next two lines in order to have only viable seeds
            m_viable_seeds_cone_age.at(i) *= m_viable_seeds;
            m_viable_seeds_dispersed_cone_age.at(i) = seeds_dispersed *
                                                          m_viable_seeds;
        }
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_plant_type(const std::string& file_plant_type)
{
    //    std::map<std::string, Plant_type> m_plant_type_list;
    //------------------------------
    std::string name_folder {Global::project_directory + "data/in/plant_types/"};
    std::ifstream ifs_plant_type(name_folder + file_plant_type);

    std::string header_plant_type;
    std::string data_plant_type;
    // Checking if the file has been opened correctly
    if (!ifs_plant_type)
    {
        std::cerr << "ERROR! Could not open the plant type file:\n";
        std::cerr << "data/in/plant_types/" << file_plant_type << "\n\a";
        std::exit(EXIT_FAILURE);
    }

    // skip the parameter names of the simulation file, i.e. header
    std::getline(ifs_plant_type, header_plant_type);

    m_plant_type_list.clear();

    while (std::getline(ifs_plant_type, data_plant_type))
    {
        if (std::all_of(std::cbegin(data_plant_type), std::cend(data_plant_type), ::isspace))
        {
            break;
        }

        std::stringstream ss_plant_type(data_plant_type);

        unsigned plant_type_name {0};
        double share_prob        {0.0};
        double mort_prob         {0.0};

        ss_plant_type >> plant_type_name
                      >> share_prob
                      >> mort_prob;

//        m_plant_type_list.push_back(Plant_type(plant_type_name,
//                                               share_prob,
//                                               mort_prob));
        m_plant_type_list.emplace(plant_type_name,
                                  Plant_type(plant_type_name,
                                             share_prob,
                                             mort_prob));

//        m_plant_type_list.push_back(Plant_type(plant_type_name,
//                                               share_prob,
//                                               mort_prob));

        //            m_plant_type_list.emplace(
        //                        plant_type_name,
        //                        Plant_type(plant_type_name, share_prob, mort_prob));
    }

    m_plant_type_share.clear();

    // accumulate the probability, i.e. add up probability to 1.0
    ///@todo TODO explain better this comment

    {
        double previous_value {0.0};

        for (const auto& [key, val] : m_plant_type_list)
        {
//            const double value {val.m_share_prob};

            if (std::empty(m_plant_type_share) == true)
            {
                m_plant_type_share.emplace(key, val.m_share_prob);
            }
            else
            {
                m_plant_type_share.emplace(key, val.m_share_prob + previous_value);
            }

            previous_value = val.m_share_prob;
        }

    }

//    for (const auto& plant : m_plant_type_list)
//    {
//        const double value {plant.m_share_prob};

//        if (m_plant_type_share.empty() == true)
//        {
//            m_plant_type_share.push_back(value);
//        }
//        else
//        {
//            m_plant_type_share.push_back(m_plant_type_share.back() + value);
//        }
//    }


    // accumulate the probability, i.e. add up probability to 1.0
    ///@todo TODO explain better this comment
//    m_plant_type_share.clear();

//    for (const auto& plant : m_plant_type_list)
//    {
//        const double value { plant.m_share_prob };

//        if (m_plant_type_share.empty() == true)
//        {
//            m_plant_type_share.push_back(value);
//        }
//        else
//        {
//            m_plant_type_share.push_back(value + m_plant_type_share.back());
//        }
//    }


//    auto it = std::max_element(std::begin(m_plant_type_share),
//                               std::end(m_plant_type_share));
//    if (*it != 1.0)
//    {
//        std::cerr << "ERROR! in the plant_type file:\n";
//        std::cerr << "'share_prob' does not add up to 1\n";
//        std::cerr << "data/in/plant_types/" << file_name_plant_type << "\n\a";
//        std::exit(EXIT_FAILURE);
//    }

    std::vector<double> check_share;

    for ([[maybe_unused]] const auto& [key, val] : m_plant_type_share)
    {
        check_share.push_back(val);
    }

    const auto it = std::max_element(std::cbegin(check_share), std::cend(check_share));

    if (*it != 1.0)
    {
        std::cerr << "ERROR! in the plant_type file:\n";
        std::cerr << "'share_prob' does not add up to 1\n";
        std::cerr << "data/in/plant_types/" << file_plant_type << "\n\a";
        std::exit(EXIT_FAILURE);
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_fuzzy_set_list(const std::string& file_fuzzy_set_list)
{
    std::ifstream ifs_fuzzy(Global::project_directory + "data/in/fuzzy_sets/" + file_fuzzy_set_list);

    // checking if the file has been opened correctly
    if (!ifs_fuzzy)
    {
        std::cerr << "ERROR! Could not open the fuzzy_set file:\n";
        std::cerr << "data/in/fuzzy_sets/" << file_fuzzy_set_list;
        std::exit(EXIT_FAILURE);
    }

    // Temporary strings
    std::string data_fuzzy;
    std::string header_fuzzy;

    // skip the parameter names of the simulation file, i.e. header
    std::getline(ifs_fuzzy, header_fuzzy);

    // remove the inner vector first because it has the smart pointers
    ///@note NOTE it might be unnecessary this step!!
    for ([[maybe_unused]] auto& [key, val] : m_fuzzy_set_list)
    {
        val.clear();
    }

    m_fuzzy_set_list.clear();
    m_fuzzy_set_weather_var.clear();

    while (std::getline(ifs_fuzzy, data_fuzzy))
    {
        if (std::all_of(std::cbegin(data_fuzzy), std::cend(data_fuzzy), ::isspace))
        {
            break;
        }

        std::stringstream ss_fuzzy(data_fuzzy);

        int id {0};
        std::string file_name_fuzzy_set;
        unsigned level {0};
        unsigned fuzzy_set_index;
        unsigned weather_variable {0};

        ss_fuzzy >> id
                 >> file_name_fuzzy_set
                 >> level
                 >> fuzzy_set_index
                 >> weather_variable;

//        m_fuzzy_set_weather_var.push_back(weather_variable);
//        m_fuzzy_set_weather_var.emplace(std::make_pair(fuzzy_set_index,
//                                                       weather_variable));

        m_fuzzy_set_weather_var.emplace(std::make_pair(level, fuzzy_set_index),
                                        weather_variable);

        std::ifstream ifs_fuzzy_set(Global::project_directory + "data/in/fuzzy_sets/" + file_name_fuzzy_set);

        // Checking if the file has been opened correctly
        if (!ifs_fuzzy_set)
        {
            std::cerr << "ERROR! Could not open the fuzzy_set file:\n";
            std::cerr << "data/in/fuzzy_sets/" << file_name_fuzzy_set;
            std::exit(EXIT_FAILURE);
        }

        // skip the parameter names of the simulation file, i.e. header
        std::getline(ifs_fuzzy_set, header_fuzzy);

        ///@todo TODO check if I can remove the shared_ptr and use normal object
        std::vector<std::shared_ptr<Fuzzy_function>> fuzzy_set;

        while (std::getline(ifs_fuzzy_set, data_fuzzy))
        {
            if (std::all_of(std::cbegin(data_fuzzy), std::cend(data_fuzzy), ::isspace))
            {
                break;
            }

            std::stringstream ss_fuzzy_set(data_fuzzy);

            int id_2 {0};
            std::string mf_shape;
            unsigned mf_index {0};
            double left {0.0};
            double middle_left {0.0};
            double middle_right {0.0};
            double right {0.0};

            ss_fuzzy_set >> id_2
                         >> mf_shape
                         >> mf_index
                         >> left
                         >> middle_left
                         >> middle_right
                         >> right;

            if (mf_shape == "trapezoid" || mf_shape == "1")
            {
                fuzzy_set.push_back(
                            std::make_shared<Trapezoid>(mf_index
                                                        , weather_variable
                                                        , left
                                                        , middle_left
                                                        , middle_right
                                                        , right));
            }
            else if (mf_shape == "triangle" || mf_shape == "2")
            {
                const auto a = middle_left;
                const auto b = middle_right + std::numeric_limits<double>::epsilon();

                if (std::fabs(a - b) > std::numeric_limits<double>::epsilon())
                {
                    std::cerr << "ERROR! 'middle_left' and 'middle_right'";
                    std::cerr << "must be equal in triangle mf shape:\n";
                    std::cerr << "File name: " << file_name_fuzzy_set << "\n";
                    std::cerr << "mf index: " << mf_index << "\n\a";
                    std::exit(EXIT_FAILURE);
                }

                fuzzy_set.push_back(
                            std::make_shared<Triangle>(mf_index
                                                       , weather_variable
                                                       , left
                                                       , middle_left
                                                       , right));
            }
            else
            {
                std::cerr << "ERROR! Wrong membership function shape:\n";
                std::cerr << "Available shapes are: \n";
                std::cerr << "- 'trapezoid' or '1'\n";
                std::cerr << "- 'triangle' or '2'\n";
                std::cerr << "File name: " << file_name_fuzzy_set << "\n\a";
                std::exit(EXIT_FAILURE);
            }
        }

        if (std::empty(fuzzy_set) == false)
        {
            m_fuzzy_set_list.emplace(
                        std::make_pair(level, fuzzy_set_index)
                        , fuzzy_set);
        }
        else
        {
            std::cerr << "ERROR! fuzzy set file is empty!:\n";
            std::cerr << "File name: " << file_name_fuzzy_set << "\n\a";
            std::exit(EXIT_FAILURE);
        }
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_flower_distribution(const std::string& file_flower_dist)
{
    const std::string folder_name { Global::project_directory + "data/in/flower_distributions/" };
    std::ifstream ifs_flower_dist(folder_name + file_flower_dist);

    std::string header_flower_dist;
    std::string data_flower_dist;

    // checking if the file has been opened correctly
    if (!ifs_flower_dist)
    {
        std::cerr << "ERROR! Could not open the flower distribution file:\n";
        std::cerr << folder_name << file_flower_dist;
        std::exit(EXIT_FAILURE);
    }

    // skip the parameter names of the simulation file, i.e. header
    std::getline(ifs_flower_dist, header_flower_dist);

    m_flower_dist.clear();

    while (std::getline(ifs_flower_dist, data_flower_dist))
    {
        if (std::all_of(std::cbegin(data_flower_dist),
                        std::cend(data_flower_dist),
                        ::isspace))
        {
            break;
        }

        std::stringstream ss_flower_dist(data_flower_dist);

        unsigned fuzzy_set;
        unsigned mf_name;
        unsigned flower_type;
        std::string dist_type;
        double param_a;
        double param_b;

        ss_flower_dist >> fuzzy_set
                       >> mf_name
                       >> flower_type
                       >> dist_type
                       >> param_a
                       >> param_b;

        Distribution_type distribution_type;

        if (dist_type == "null" || dist_type == "0")
        {
            distribution_type = Distribution_type::null;
        }
        else if (dist_type == "pois" || dist_type == "1")
        {
            distribution_type = Distribution_type::pois;
        }
        else if (dist_type == "nbinom" || dist_type == "2")
        {
            distribution_type = Distribution_type::nbinom;
        }
        else if (dist_type == "nbinom_mu" || dist_type == "3")
        {
            distribution_type = Distribution_type::nbinom_mu;
        }
        else if (dist_type == "geom" || dist_type == "4")
        {
            distribution_type = Distribution_type::geom;
        }
        else
        {
            std::cerr << "The parameter 'dist_type' is not valid:\n"
                         "- '0' or 'null' for no distribution (i.e. zero flowers)\n"
                         "- '1' or 'pois' for Poisson distribution\n"
                         "- '2' or 'nbinom' for Negative binomial (size and prob) distribution\n"
                         "- '3' or 'nbinom_mu' for Negative binomial (size and mu) distribution\n"
                         "- '4' or 'geom' for Geometric distribution\n\a";
            std::exit(EXIT_FAILURE);
        }

        m_flower_dist.emplace(
                    std::make_tuple(fuzzy_set
                                    , mf_name
                                    , flower_type)
                    , Flower_distribution(fuzzy_set
                                          , mf_name
                                          , flower_type
                                          , distribution_type
                                          , param_a
                                          , param_b));
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_habitat_quality(const std::string& file_habitat)
{
    const std::string folder_name_habitat { Global::project_directory + "data/in/habitat_quality/" };
    std::ifstream ifs_habitat(folder_name_habitat + file_habitat);

    std::string header_habitat;
    std::string data_habitat;

    // checking if the file has been opened correctly
    if (!ifs_habitat)
    {
        throw std::runtime_error(
            "ERROR! Could not open the flower distribution file:\n"
            + folder_name_habitat
            + file_habitat
            + "\n\a"
        );
    }

    // skip the parameter names of the simulation file, i.e. header
    std::getline(ifs_habitat, header_habitat);

    m_flower_dist.clear();

    while (std::getline(ifs_habitat, data_habitat))
    {
        if (std::all_of(std::cbegin(data_habitat), std::cend(data_habitat), ::isspace))
        {
            break;
        }

        std::stringstream ss_habitat(data_habitat);

        std::string habitat_label;
        double effect_flowers {0.0};

        ss_habitat >> habitat_label
                   >> effect_flowers;

        m_habitat_quality.emplace(habitat_label
                                  , Habitat_quality(habitat_label
                                                    , effect_flowers));
    }
}

//-----------------------------------------------------------------------------

void Parameter_reader::read_changed_parameters(
        const std::vector<std::pair<std::string, std::string>>& parameters)
{
    for (const auto& [key, val] : parameters)
    {
        if (key == "scenario_file" || key == "scn_file")
        {
            read_scenario(val);
        }
        // -- parameters from the scenario file ------------------
        else if (key == "model_type")
        {
            const int model_type { std::stoi(val) };

            switch (model_type)
            {
                case 0 :
                    m_model_type = Model_type::cohort_based;
                    break;
                case 1 :
                    m_model_type = Model_type::plant_based;
                    read_plant_type(m_file_plant_type);
                    read_fuzzy_set_list(m_file_fuzzy_set_list);
                    read_flower_distribution(m_file_flower_dist);
                    break;
                default :
                    throw std::runtime_error(
                        "The parameter 'model_type' is not valid:\n"
                        "- value '0' for cohort-based modelling\n"
                        "- value '1' for individual-based modelling\n\a"
                    );
            }
        }
        else if (key == "study_size")
        {
            m_file_study_size = val;
            read_study_size(val);
        }
        else if (key == "metapop_file")
        {
            m_file_metapopulation = val;
        }
        else if (key == "carrying_capacity")
        {
            m_carrying_capacity = std::stoi(val);
        }
        else if (key == "init_cones")
        {
            m_init_cones_plant = std::stod(val);
        }
        else if (key == "init_seeds")
        {
            m_init_viable_seeds_plant = std::stod(val);
        }
        else if (key == "is_climate")
        {
            if (val == "0" || val == "false")
            {
                m_is_climate = false;
            }
            else if (val == "1" || val == "true")
            {
                m_is_climate = true;
            }
            else
            {
                throw std::runtime_error(
                    "ERROR: " + val + " must be true (1) or false (0)\n"
                );
            }
        }
        else if (key == "long_term_rain")
        {
            m_long_term_rain = std::stod(val);
        }
        else if (key == "weather_scenario" || key == "weather_scn")
        {
            const int weather_scenario { std::stoi(val) };

            switch (weather_scenario)
            {
            case 0:
                m_weather_scenario = Weather_scenario::random;
                break;
            case 1:
                m_weather_scenario = Weather_scenario::consecutive;
                break;
            default:
                throw std::runtime_error(
                    "The parameter 'weather_scenario' is not valid:\n"
                    "- value '0' for random scenario\n"
                    "- value '1' for consecutive scenario \n\a"
                );
            }
        }
        else if (key == "climate_file")
        {
            m_file_climate = val;
        }
        else if (key == "is_fire")
        {
            if (val == "0" || val == "false")
            {
                m_is_fire = false;
            }
            else if (val == "1" || val == "true")
            {
                m_is_fire = true;
            }
            else
            {
                throw std::runtime_error(
                    "ERROR: " + val + " must be true (1) or false (0)\n"
                );
            }
        }
        else if (key == "fire_file")
        {
            m_file_fire = val;
            read_fire(val);
        }
        else if (key == "dispersal_file")
        {
            m_file_dispersal = val;
            read_dispersal(val);
        }
        else if (key == "species_file")
        {
            m_file_species = val;
            read_species(val);
        }
        else if (key == "habitat_file")
        {
            m_file_habitat = val;
            read_habitat_quality(val);
        }
        else if (key == "plant_type_file")
        {
            m_file_plant_type = val;
            read_plant_type(val);
        }
        else if (key == "fuzzy_set_file")
        {
            m_file_fuzzy_set_list = val;
            read_fuzzy_set_list(val);
        }
        else if (key == "flower_dist_file")
        {
            m_file_flower_dist = val;
            read_flower_distribution(val);
        }
        // -- parameters from the study_size file -------------------
        else if (key == "cell_size")
        {
            m_cell_size = std::stoi(val);
        }
        else if (key == "study_size_x")
        {
            m_study_size_x = std::stoi(val);
        }
        else if (key == "study_size_y")
        {
            m_study_size_y = std::stoi(val);
        }
        // -- parameters from fire file -------------------------------
        else if (key == "fire_interval_scn")
        {
            const int fire_interval_scn { std::stoi(val) };

            switch (fire_interval_scn)
            {
                case 0 :
                    m_fire_interval_scn = Fire_scenario::deterministic;
                    break;
                case 1 :
                    m_fire_interval_scn = Fire_scenario::truncated_normal;
                    break;
                case 2 :
                    m_fire_interval_scn = Fire_scenario::weibull;
                    break;
                default :
                    throw std::runtime_error(
                        "The parameter 'fire_interval_scn' in sim file is not valid:\n"
                        "- value '0' for deterministic fire events\n"
                        "- value '1' for truncated_normal fire events\n"
                        "- value '2' for Weibull distribution\n\a"
                    );
            }
        }
        // -- parameters from the fire file ----------------------
        else if (key == "fire_size_scn")
        {
            const int fire_size_scn { std::stoi(val) };

            switch (fire_size_scn)
            {
                case 0 :
                    m_fire_size_scn = Fire_scenario::deterministic;
                    break;
                case 1 :
                    m_fire_size_scn = Fire_scenario::truncated_normal;
                    break;
                default :
                    throw std::runtime_error(
                        "The parameter 'fire_scenario' in sim file is not valid:\n"
                        "- value '0' for deterministic fire events\n"
                        "- value '1' for truncated_normal fire events\n\a"
                    );
            }
        }
        else if (key == "fire_scale")
        {
            const int fire_scale { std::stoi(val) };

            switch (fire_scale)
            {
                case 0 :
                    m_fire_scale = Fire_scale::patchy;
                    break;
                case 1 :
                    m_fire_scale = Fire_scale::study_area;
                    break;
                default :
                    throw std::runtime_error(
                        "The parameter 'fire_scale' is not valid:\n"
                        "- value '0' for patchy fires\n"
                        "- value '1' for fire in the entire study area\n\a"
                    );
            }
        }
        else if (key == "fire_size_x")
        {
            m_fire_size_x = std::stod(val);
        }
        else if (key == "fire_size_y")
        {
            m_fire_size_y = std::stod(val);
        }
        else if (key == "fire_interval_mean")
        {
            m_fire_interval_mean = std::stoi(val);
        }
        else if (key == "fire_interval_lower_cut")
        {
            m_fire_interval_lower_cut = std::stoi(val);
        }
        else if (key == "burned_lower_cut")
        {
            m_burned_lower_cut = std::stoi(val);
        }
        else if (key == "burned_upper_cut")
        {
            m_burned_upper_cut = std::stoi(val);
        }
        else if (key == "fire_a")
        {
            m_fire_a = std::stod(val);
        }
        else if (key == "fire_b")
        {
            m_fire_b = std::stod(val);
        }
        // -- parameters from dispersal file ---------------------
//        else if (key == "wind_max")
//        {
//            m_wind_max = std::stod(val);
//        }
        else if (key == "wind_prop")
        {
            m_wind_prop = std::stod(val);
        }
        else if (key == "wind_a")
        {
            m_wind_a = std::stod(val);
        }
        else if (key == "wind_b")
        {
            m_wind_b = std::stod(val);
        }
        else if (key == "wind_direction")
        {
            const int wind_direction { std::stoi(val) };

            switch (wind_direction)
            {
                case 0 :
                    m_wind_direction = Wind_direction::random;
                    break;
                case 1 :
                    m_wind_direction = Wind_direction::geraldton;
                    break;
                default :
                    throw std::runtime_error(
                        "The parameter 'wind_direction' is not valid:\n"
                        "- value '0' for random 360 degrees\n"
                        "- value '1' for wind roses from Geraldton station\n\a"
                    );
            }
        }
        else if (key == "birds_dist_max")
        {
            m_birds_dist_max = std::stod(val);
        }
        else if (key == "birds_prop")
        {
            m_birds_prop = std::stod(val);
        }
        else if (key == "follicles_distr_type")
        {
            const int follicles_distr_type { std::stoi(val) };

            switch (follicles_distr_type)
            {
            case 1:
                m_follicles_distr_type = Distribution_type::pois;
                break;
            case 3:
                m_follicles_distr_type = Distribution_type::nbinom_mu;
                break;
            default:
                throw std::runtime_error(
                        "The parameter 'follicles_distr_type' is not valid:\n"
                        "- value '1' for Poisson distribution\n"
                        "- value '3' for negative binomial with 'mu' parameter\n");
            }
        }
        else if (key == "follicles_distr_a")
        {
            m_follicles_distr_a = std::stod(val);
        }
        else if (key == "follicles_distr_b")
        {
            m_follicles_distr_b = std::stod(val);
        }
//        else if (key == "follicles_max")
//        {
//            m_follicles_max = std::stoi(val);
//        }
        else if (key == "postfire_follicles_open")
        {
            m_postfire_follicles_open = std::stod(val);
        }
        // -- parameters from species file -------------------------
        else if (key == "longevity" || key == "plant_longevity")
        {
            m_plant_longevity = std::stoi(val);
        }
        else if (key == "young" || key == "plant_young" || key == "plant_age_young")
        {
            m_plant_age_young = std::stoi(val);
        }
        else if (key == "adult" || key == "plant_adult" || key == "plant_age_adult")
        {
            m_plant_age_adult = std::stoi(val);
        }
        else if (key == "mort_scn")
        {
            const int mort_scenario { std::stoi(val) };

            switch (mort_scenario)
            {
            case 0 :
                m_mort_scenario = Mortality_scenario::age_weather_relative;
                break;
            case 1 :
                m_mort_scenario = Mortality_scenario::age_weather_absolute;
                break;
            case 2 :
                m_mort_scenario = Mortality_scenario::spring_autumn_mean;
                break;
            case 3 :
                m_mort_scenario = Mortality_scenario::SDD_autumn_LDD_spring;
                break;
            case 4 :
                m_mort_scenario = Mortality_scenario::resi_autum_immi_spring;
                break;
            case 5 :
                m_mort_scenario = Mortality_scenario::lowest_mortality;
                break;
            default :
                throw std::runtime_error(
                        "The parameter 'mort_scenario' is not valid:\n"
                        "- value '0': relative scenario from Keith et al. 2014\n"
                        "- value '1': absolute scenario from Keith et al. 2014\n"
                        "- value '2': arithmetic mean of spring and autumn recruit mortality\n"
                        "- value '3': SDD autumn and LDD spring\n"
                        "- value '4': residents autumn and immigrants spring\n"
                        "- value '5': lowest mortality rates\n\a");
            }
        }
        else if (key == "mort_recruit_post_min" || key == "recruit_post_min")
        {
            m_mort_recruit_postfire_min = std::stod(val);

            if (m_mort_recruit_postfire_min < 0.0
                || m_mort_recruit_postfire_min > 1.0)
            {
                throw std::runtime_error(
                    "the parameter 'mort_recruit_postfire_min' is not "
                    "between 0 and 1");
            }
        }
        else if (key == "mort_recruit_post_max" || key == "recruit_post_max")
        {
            m_mort_recruit_postfire_max = std::stod(val);

            if (m_mort_recruit_postfire_max < 0.0
                || m_mort_recruit_postfire_max > 1.0)
            {
                throw std::runtime_error(
                    "the parameter 'mort_recruit_postfire_max' is not "
                    "between 0 and 1");
            }
        }
        else if (key == "mort_recruit_post_mean" || key == "recruit_post_mean")
        {
            m_mort_recruit_postfire_mean = std::stod(val);
        }
        else if (key == "recruit_interfire")
        {
            m_recruit_interfire = std::stod(val);
        }
        else if (key == "recruit_weather")
        {
            m_recruit_weather = std::stod(val);
        }
        else if (key == "senescence_age" || key == "plant_senescence_age")
        {
            m_senescence_age = std::stoi(val);
        }
        else if (key == "senescence_increase" || key == "plant_senescence_increase")
        {
            m_senescence_increase = std::stod(val);
        }
        else if (key == "mort_min" || key == "plant_mort_min")
        {
            m_mort_min = std::stod(val);

            if (m_mort_min < 0.0 || m_mort_min > 1.0)
            {
                throw std::runtime_error(
                    "the parameter 'mort_recruit_postfire_max' is not "
                    "between 0 and 1");
            }
        }
//        else if (key == "plant_mort_age_a")
//        {
//            m_plant_mort_age_a = std::stod(val);
//        }
//        else if (key == "plant_mort_age_b")
//        {
//            m_plant_mort_age_b = std::stod(val);
//        }
//        else if (key == "plant_mort_age_c")
//        {
//            m_plant_mort_age_c = std::stod(val);
//        }
//        else if (key == "plant_mort_weather_a")
//        {
//            m_plant_mort_weather_a = std::stod(val);
//        }
//        else if (key == "plant_mort_weather_b")
//        {
//            m_plant_mort_weather_b = std::stod(val);
//        }
        else if (key == "mort_a" || key == "plant_mort_a")
        {
            m_mort_a = std::stod(val);
        }
        else if (key == "mort_b" || key == "plant_mort_b")
        {
            m_mort_b = std::stod(val);
        }
        else if (key == "mort_c" || key == "plant_mort_c")
        {
            m_mort_c = std::stod(val);
        }
        else if (key == "mort_d" || key == "plant_mort_d")
        {
            m_mort_d = std::stod(val);
        }
        else if (key == "mort_e" || key == "plant_mort_e")
        {
            m_mort_e = std::stod(val);
        }
        else if (key == "mort_f" || key == "plant_mort_f")
        {
            m_mort_f = std::stod(val);
        }
        else if (key == "cone_cycle")
        {
            const int cone_cycle { std::stoi(val) };

            if (cone_cycle < 0)
            {
                throw std::runtime_error(
                    "ERROR: the parameter " + val + " must positive integer.\n"
                );
            }

            m_cone_cycle = static_cast<unsigned>(cone_cycle);
        }
        else if (key == "seed_longevity")
        {
            const int seed_longevity { std::stoi(val) };

            if (seed_longevity < 0)
            {
                throw std::runtime_error(
                    "ERROR: the parameter " + val + " must positive integer.\n"
                );
            }

            m_seed_longevity = static_cast<unsigned>(seed_longevity);
        }
        else if (key == "flower_scn")
        {
            const int flower_scenario { std::stoi(val) };

            switch (flower_scenario)
            {
                case 0 :
                    m_flower_scenario = Flower_scenario::baseline;
                    break;
                case 1 :
                    m_flower_scenario = Flower_scenario::current;
                    break;
                case 2 :
                    m_flower_scenario = Flower_scenario::current_nb;
                    break;
                default :
                    throw std::runtime_error(
                        "The parameter 'flower_scenario' is not valid:\n"
                        "- value '0' for baseline scenario\n"
                        "- value '1' for current scenario using lmer\n"
                        "- value '2' for current scenario using glmer.nb\n\a"
                    );
            }
        }
//        else if (key == "flower_max")
//        {
//            m_flower_max = std::stod(val);
//        }
        else if (key == "flower_age_a")
        {
            m_flower_age_a = std::stod(val);
        }
        else if (key == "flower_age_b")
        {
            m_flower_age_b = std::stod(val);
        }
        else if (key == "flower_age_c")
        {
            m_flower_age_c = std::stod(val);
        }
        else if (key == "flower_weather_a")
        {
            m_flower_weather_a = std::stod(val);
        }
        else if (key == "flower_weather_b")
        {
            m_flower_weather_b = std::stod(val);
        }
        else if (key == "flower_weather_c")
        {
            m_flower_weather_c = std::stod(val);
        }
        else if (key == "flower_weather_d")
        {
            m_flower_weather_d = std::stod(val);
        }
        else if (key == "flower_weather_e")
        {
            m_flower_weather_e = std::stod(val);
        }
        else if (key == "pollination")
        {
            m_pollination_success = std::stod(val);
        }
        else if (key == "follicles")
        {
            m_num_follicles = std::stod(val);
        }
        else if (key == "seeds")
        {
            m_num_seeds = std::stod(val);
        }
        else if (key == "firm_seeds")
        {
            m_firm_seeds = std::stod(val);
        }
        else if (key == "viable_seeds")
        {
            m_viable_seeds = std::stod(val);
        }
        else if (key == "insect_a")
        {
            m_insect_a = std::stod(val);
        }
        else if (key == "insect_b")
        {
            m_insect_b = std::stod(val);
        }
        else if (key == "decay_a")
        {
            m_decay_a = std::stod(val);
        }
        else if (key == "decay_b")
        {
            m_decay_b = std::stod(val);
        }
        else if (key == "open_a")
        {
            m_open_follicles_a = std::stod(val);
        }
        else if (key == "open_b")
        {
            m_open_follicles_b = std::stod(val);
        }
        else if (key == "open_c")
        {
            m_open_follicles_c = std::stod(val);
        }
        else
        {
            throw std::runtime_error(
                "ERROR: the following input parameter does not exist: " + key
            );
        }
    }
}
