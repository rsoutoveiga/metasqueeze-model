#include "wind_simulator.h"

double Wind_simulator::get_wind_direction()
{
    double random_num { Global::prng.uniform_01() };
    double angle_degrees {0.0};

    if (random_num < 0.5)  // wind roses, summer, 3 pm Summer at Geraldton station
    {
        random_num = Global::prng.uniform_01();
        assert(random_num >= 0.0 && random_num <= 1.0);

        if (random_num < 0.588)
        {
            random_num    = Global::prng.uniform_01();
            angle_degrees = Global::prng.uniform_real(0.0, 22.5);

            if (random_num < 0.5)
            {
                angle_degrees = 360.0 - angle_degrees;
            }
        }
        else if (random_num < 0.850)
        {
            angle_degrees = Global::prng.uniform_real(22.5, 67.5);
        }
        else if (random_num < 0.920)
        {
            angle_degrees = Global::prng.uniform_real(67.5, 112.5);
        }
        else if (random_num < 0.946)
        {
            angle_degrees = Global::prng.uniform_real(112.5, 157.5);
        }
        else if (random_num <= 1.0)
        {
            angle_degrees = Global::prng.uniform_real(247.5, 337.5);
        }
        else
        {
            assert(false && "not all cases were evaluated");
        }
    }
    else                   // wind roses, 9:00 a.m. Geraldton station
    {
        random_num = Global::prng.uniform_01();
        assert(random_num >= 0 && random_num <= 1.0);

        if (random_num < 0.316)
        {
            random_num    = Global::prng.uniform_01();
            angle_degrees = Global::prng.uniform_real(0.0, 22.5);

            if (random_num < 0.5)
            {
                angle_degrees = 360.0 - angle_degrees;
            }
        }
        else if (random_num < 0.548)
        {
            angle_degrees = Global::prng.uniform_real(292.5, 337.5);
        }
        else if (random_num < 0.680)
        {
            angle_degrees = Global::prng.uniform_real(247.5, 292.5);
        }
        else if (random_num < 0.800)
        {
            angle_degrees = Global::prng.uniform_real(202.5, 247.5);
        }
        else if (random_num < 0.827)
        {
            angle_degrees = Global::prng.uniform_real(157.5, 202.5);
        }
        else if (random_num < 0.862)
        {
            angle_degrees = Global::prng.uniform_real(112.5, 157.5);
        }
        else if (random_num < 0.898)
        {
            angle_degrees = Global::prng.uniform_real(67.5, 112.5);
        }
        else if (random_num < 0.965)
        {
            angle_degrees = Global::prng.uniform_real(22.5, 67.5);
        }
        else if (random_num <= 1.000)  // calm
        {
            angle_degrees = Global::prng.uniform_real(0.0, 360.0);
        }
        else
        {
            assert(false && "not all cases were evaluated");
        }
    }
    return angle_degrees;
}

