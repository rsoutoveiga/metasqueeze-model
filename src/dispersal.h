#ifndef DISPERSAL_H
#define DISPERSAL_H

#include <deque>

#include "gene_flow.h"

struct Dispersal
{
    Dispersal(const Gene_flow& gene_flow
              , const int LDD_cones_by_birds
              , const int total_viable_seeds)

        : m_gene_flow{gene_flow}
        , m_LDD_cones_by_birds{LDD_cones_by_birds}
        , m_total_viable_seeds{total_viable_seeds}
    {}

    const Gene_flow m_gene_flow;
    const int m_LDD_cones_by_birds;

    /** @brief m_total_viable_seeds seeds dispersed of a cohort.
     *  @note this variable is not assigned as 'const' because is updated
     *  when some number of seeds are disperse by a dispersal vector
     *  (e.g. LDD_by_birds, LDD_by_wind).
     *
     */
    int m_total_viable_seeds {0};
};
#endif // DISPERSAL_H
