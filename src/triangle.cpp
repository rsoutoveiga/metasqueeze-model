#include "triangle.h"

Triangle::Triangle(const unsigned mf_index,
                   const unsigned climate_index,
                   const double left,
                   const double middle,
                   const double right)
    : m_mf_index{mf_index}
    , m_climate_index{climate_index}
    , m_left{left}
    , m_middle{middle}
    , m_right{right}
{}
