#ifndef HABITAT_QUALITY_H
#define HABITAT_QUALITY_H

#include <string>
#include <string_view>

///
/// \brief The Habitat_quality struct stores the habitat quality
/// coefficient of mean increase in flower production
///
struct Habitat_quality
{
    Habitat_quality(std::string_view habitat_label,
                    const double effect_flowers)

        : m_habitat_label{habitat_label}
        , m_effect_flowers{effect_flowers}
    {}

    const std::string m_habitat_label;
    const double m_effect_flowers;
};

#endif // HABITAT_QUALITY_H
