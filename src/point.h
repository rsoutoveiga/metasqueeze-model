#ifndef POINT_H
#define POINT_H

struct Point
{
    Point(const double r,
          const double c)
        : row(r)
        , col(c)
    {}

    double row;
    double col;
};

#endif // POINT_H
