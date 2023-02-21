#ifndef TRAPEZOID_H
#define TRAPEZOID_H


#include "fuzzy_function.h"

class Trapezoid : public Fuzzy_function
{
public:

    Trapezoid(const unsigned mf_index,
              const unsigned climate_index,
              const double left,
              const double middle_left,
              const double middle_right,
              const double right);


    double get_value(const double dot) const override
    {
        if (dot <= m_left)
        {
            return 0.0;
        }
        else if (dot < m_middle_left)
        {
            return (dot - m_left) / (m_middle_left - m_left);
        }
        else if (dot <= m_middle_right)
        {
            return 1.0;
        }
        else if (dot < m_right)
        {
            return (m_right - dot) / (m_right - m_middle_right);
        }
        else
        {
            return 0.0;
        }
    }

    unsigned get_mf_index() const override
    {
        return m_mf_index;
    }

    unsigned get_climate_index() const override
    {
        return m_climate_index;
    }

private:

    const unsigned m_mf_index {0};
    const unsigned m_climate_index {0};
    const double m_left {0.0};
    const double m_middle_left {0.0};
    const double m_middle_right {0.0};
    const double m_right {0.0};
};

#endif // TRAPEZOID_H
