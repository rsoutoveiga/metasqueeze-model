#ifndef GENE_FLOW_H
#define GENE_FLOW_H

///
/// \brief The Gene_flow class for each cohort
///
class Gene_flow
{
public:
    Gene_flow(const int population_id)
        : m_genetic{population_id}
        , m_previous{population_id}
        , m_current{population_id}
    {}

    // setters
    void set_current(const int population_id)
    {
        m_previous = m_current;
        m_current = population_id;
    }

    // getters
    int get_genetic() const { return m_genetic; }
    int get_previous() const { return m_previous; }
    int get_current() const { return m_current; }

private:
    int m_genetic  {1};  // original population id of the cohort
    int m_previous {1};  // population id before last fire
    int m_current  {1};  // current population id since last fire
};

#endif // GENE_FLOW_H
