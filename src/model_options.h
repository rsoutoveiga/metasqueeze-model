#ifndef MODEL_OPTIONS_H
#define MODEL_OPTIONS_H


///@brief all model options used in the simulation input files

//enum class Study_area_method { regular = 1, irregular};
enum class Weather_scenario  { random, consecutive };
enum class Flower_scenario   { baseline, current, current_nb };
enum class Fire_scenario     { deterministic, truncated_normal, weibull };
enum class Fire_scale        { patchy, study_area };

enum class Mortality_scenario { age_weather_relative,
                                age_weather_absolute,
                                spring_autumn_mean,
                                SDD_autumn_LDD_spring,
                                resi_autum_immi_spring,
                                lowest_mortality
                              };

enum class Dispersal_vector { SDD_wind, LDD_bird, LDD_wind };

enum class Wind_direction { random, geraldton };

enum class Flower_producer { low, high };

enum class Model_type { cohort_based, plant_based};

enum class Distribution_type { null, pois, nbinom, nbinom_mu, geom };


#endif // MODEL_OPTIONS_H
