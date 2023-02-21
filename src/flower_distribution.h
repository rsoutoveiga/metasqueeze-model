#ifndef FLOWER_DISTRIBUTION_H
#define FLOWER_DISTRIBUTION_H

#include "model_options.h"

///
/// \brief The Flower_distribution struct is used to store the
/// flower count probability density functions for poor and good producers.
///
struct Flower_distribution
{
    Flower_distribution(const unsigned fuzzy_set,
                        const unsigned mf_name,
                        const unsigned flower_type,
                        const Distribution_type dist_type,
                        const double param_a,
                        const double param_b)
        : m_fuzzy_set(fuzzy_set)
        , m_mf_name(mf_name)
        , m_flower_type(flower_type)
        , m_dist_type(dist_type)
        , m_param_a(param_a)
        , m_param_b(param_b)
    {}

    const unsigned m_fuzzy_set;
    const unsigned m_mf_name;
    const unsigned m_flower_type;
    const Distribution_type m_dist_type;
    const double m_param_a;
    const double m_param_b;
};


#endif // FLOWER_DISTRIBUTION_H
