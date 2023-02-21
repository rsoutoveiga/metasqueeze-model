#ifndef CELL_H
#define CELL_H

struct Cell
{
    Cell(const unsigned r,
         const unsigned c)
        : row{r}
        , col{c}
    {}

    const unsigned row;
    const unsigned col;
};

#endif // CELL_H
