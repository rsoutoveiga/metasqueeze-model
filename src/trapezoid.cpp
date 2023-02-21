#include "trapezoid.h"

Trapezoid::Trapezoid(
        const unsigned mf_index,
        const unsigned climate_index,
        const double left,
        const double middle_left,
        const double middle_right,
        const double right)

    : m_mf_index{mf_index}
    , m_climate_index{climate_index}
    , m_left{left}
    , m_middle_left{middle_left}
    , m_middle_right{middle_right}
    , m_right{right}
{}
