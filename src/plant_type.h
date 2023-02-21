#ifndef PLANT_TYPE_H
#define PLANT_TYPE_H

///
/// \brief The Plant_type struct stores the plant trait characteristics
/// read in input file.
///
struct Plant_type
{
    Plant_type(const unsigned plant_type,
               const double share_prob,
               const double mort_prob)
        : m_plant_type(plant_type)
        , m_share_prob(share_prob)
        , m_mort_prob(mort_prob)
    {}

    const unsigned m_plant_type {0};
    const double m_share_prob {1.0};
    const double m_mort_prob {0.0};
};

#endif // PLANT_TYPE_H
