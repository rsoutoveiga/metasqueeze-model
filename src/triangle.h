#ifndef TRIANGLE_H
#define TRIANGLE_H


#include "fuzzy_function.h"

class Triangle : public Fuzzy_function
{
public:

    Triangle(const unsigned mf_index,
             const unsigned climate_index,
             const double left,
             const double middle,
             const double right);

    double get_value(const double dot) const override
    {
        if (dot <= m_left)
        {
            return 0.0;
        }
        else if (dot < m_middle)
        {
            return (dot - m_left) / (m_middle - m_left);
        }
        else if (dot == m_middle)
        {
            return 1.0;
        }
        else if (dot < m_right)
        {
            return (m_right - dot) / (m_right - m_middle);
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
    const double m_middle {0.0};
    const double m_right {0.0};
};


#endif // TRIANGLE_H
