#include "plant.h"

Plant::Plant(const Parameters& params)
    : m_seedbank_cones(params.cone_cycle + params.seed_longevity, 0)
{}

//-----------------------------------------------------------------------------

void Plant::aging_cones()
{
    // Expects
    assert(!std::empty(m_seedbank_cones));

    // aging of seeds: rotate the elements of vector to the right
    std::rotate(std::begin(m_seedbank_cones)
                , std::end(m_seedbank_cones) - 1
                , std::end(m_seedbank_cones));
}

//-----------------------------------------------------------------------------

int Plant::interfire_seed_dispersal(const Parameters& params)
{
    double seeds_dispersed {0.0};

    for (unsigned i {params.cone_cycle}; i < std::size(m_seedbank_cones); ++i)
    {
        seeds_dispersed += m_seedbank_cones.at(i)
                * params.viable_seeds_dispersed_cone_age.at(i);
    }

    // Ensures
    assert(seeds_dispersed >= 0.0);

    return static_cast<int>(seeds_dispersed + 0.5);
}

