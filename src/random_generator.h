#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H


#include <random>
#include <chrono>

class Random_generator
{
public:
    Random_generator()
    {
        std::random_device rd;
        std::seed_seq sd{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        m_engine.seed(sd);
    }

    // getter
    std::mt19937& get_engine() { return m_engine; }

    // setter
    void set_seed(const unsigned seed) { m_engine.seed(seed); }

    // distributions
    double uniform_01()
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        return distribution(m_engine);
    }

    int uniform_int(const int min, const int max)
    {
        std::uniform_int_distribution<int> distribution(min, max);

        return distribution(m_engine);
    }

    double uniform_real(const double min, const double max)
    {
        std::uniform_real_distribution<double> distribution(min, max);

        return distribution(m_engine);
    }

    double lognormal(const double log_mean, const double log_sd)
    {
        std::lognormal_distribution<double> distribution(log_mean, log_sd);

        return distribution(m_engine);
    }

    unsigned uniform_unsigned(const unsigned min, const unsigned max)
    {
        std::uniform_int_distribution<unsigned> distribution(min, max);

        return distribution(m_engine);
    }

//    long long unsigned uniform_unsigned(const long long unsigned min, const long long unsigned max)
//    {
//        std::uniform_int_distribution<long long unsigned> distribution(min, max);

//        return distribution(m_engine);
//    }

    int binomial(const int living_plants, const double prob_plant_mortality)
    {
        std::binomial_distribution<int>
                distribution(living_plants, prob_plant_mortality);

        return distribution(m_engine);
    }

    int poisson(const double mean)
    {
        std::poisson_distribution<int> distribution(mean);

        return distribution(m_engine);
    }

    double gaussian(const double mean, const double sd)
    {
        std::normal_distribution<double> distribution(mean, sd);

        return distribution(m_engine);
    }

    //------------------------------------------------------------------------------------------------
    ///@brief nbinom is taken from the R source (it is implemented in c)
    // https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/nmath/rnbinom.c
    //    double rnbinom(double size, double prob)
    //    {
    //        if(!R_FINITE(prob) || size <= 0 || prob <= 0 || prob > 1)
    //        /* prob = 1 is ok, PR#1218 */
    //        ML_ERR_return_NAN;
    //        if(!R_FINITE(size)) size = DBL_MAX / 2.; // '/2' to prevent rgamma() returning Inf
    //        return (prob == 1) ? 0 : rpois(rgamma(size, (1 - prob) / prob));
    //    }

    //    double rnbinom_mu(double size, double mu)
    //    {
    //        if(!R_FINITE(mu) || size <= 0 || mu < 0)
    //        ML_ERR_return_NAN;
    //        if(!R_FINITE(size)) size = DBL_MAX / 2.;
    //        return (mu == 0) ? 0 : rpois(rgamma(size, mu / size));
    //    }

    int negative_binomial(const double size, const double prob)
    {
        std::gamma_distribution<double> dist_gamma(size, (1 - prob) / prob);

        std::poisson_distribution<int> dist_pois(dist_gamma(m_engine));

        return dist_pois(m_engine);
    }

    int negative_binomial_mu(const double size, const double mu)
    {
        std::gamma_distribution<double> dist_gamma(size, mu / size);

        std::poisson_distribution<int> dist_pois(dist_gamma(m_engine));

        return dist_pois(m_engine);
    }

    int geometric(const double prob)
    {
        std::geometric_distribution<int> distribution(prob);

        return distribution(m_engine);
    }

    double weibull(const double shape, const double scale)
    {
        std::weibull_distribution<double> distribution(shape, scale);

        return distribution(m_engine);
    }

private:
    std::mt19937 m_engine;
};

#endif // RANDOM_GENERATOR_H
