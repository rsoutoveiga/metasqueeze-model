#ifndef FUZZY_FUNCTION_H
#define FUZZY_FUNCTION_H


class Fuzzy_function
{

public:
    virtual ~Fuzzy_function();

    // calculation of the fuzzy membership value computed using
    // triangle or trapezoid membership function
    virtual double get_value(const double dot) const = 0;

    virtual unsigned get_mf_index() const = 0;
    virtual unsigned get_climate_index() const = 0;

};


#endif // FUZZY_FUNCTION_H
