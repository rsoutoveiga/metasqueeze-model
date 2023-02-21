#ifndef PLANT_H
#define PLANT_H

#include <vector>
#include <algorithm>

#include "parameters.h"

///
/// \brief The Plant class is to include intraspecific variability in
/// plant performance within patches (i.e., individual-based approach). Plants
/// belong to Cohort class.
///
class Plant
{
public:
    Plant(const Parameters& params);

    void aging_cones();

    inline void add_new_cones(const int new_cones)
    {
        m_seedbank_cones.at(0) = new_cones;
    }

    int interfire_seed_dispersal(const Parameters& params);


    // getters
    int get_seedbank_cones(const Parameters& params) const
    {
        return std::accumulate(
                    std::cbegin(m_seedbank_cones) + params.cone_cycle
                    , std::cend(m_seedbank_cones)
                    , 0
                    , std::plus<int>());
    }

    double get_seedbank(const Parameters& params) const
    {
        std::vector<double> seedbank(params.seed_longevity, 0.0);

        std::transform(std::begin(m_seedbank_cones) + params.cone_cycle
                       , std::end(m_seedbank_cones)
                       , std::begin(params.viable_seeds_cone_age)
                       + params.cone_cycle
                       , std::begin(seedbank)
                       , std::multiplies<double>());

        return std::accumulate(std::cbegin(seedbank)
                               , std::cend(seedbank)
                               , 0.0
                               , std::plus<double>());
    }

    // setters
    inline void set_seedbank_to_zero()
    {
        std::fill(std::begin(m_seedbank_cones)
                  , std::end(m_seedbank_cones)
                  , 0.0);
    }

private:
    std::vector<int> m_seedbank_cones;
};

#endif // PLANT_H
